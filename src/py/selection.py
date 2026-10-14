#!/usr/bin/env python3
"""selection.py: the event selection, step by step on one event, with a figure per step.

    python3 src/py/selection.py RUN [--events 0-9] [--cm median|mean|interp|none]
            [--baseline end|pre] [--baseline-window 250] [--smooth 64]
            [--window 1000,4000] [--nsigma 5] [--from-step N] [--save NAME]

Interactive: Enter goes to the next step, n to the next event, q quits. See docs/selection_doc.tex.
"""
import argparse
import contextlib
import io
import sys
from pathlib import Path

import h5py
import numpy as np
import matplotlib.pyplot as plt

from wfview import DATA, PLOTS, find_run, parse_events, smooth

WRAP = 2 ** 24
CM_METHODS = ("median", "interp", "none")     # not the mean: it makes the sum identically zero
STEPS = ("read", "common-mode analysis", "locate the pulse (provisional)", "interpolate m across W",
         "subtract and sum", "baseline", "pulse search", "charge", "multiplicity", "verdict")


# ---- display ---------------------------------------------------------------------------
class Stepper:
    """Shows one figure per step and waits on the terminal, or saves without prompting.

    show() returns "next" (go on), "skip" (abandon this event) or "quit"."""

    def __init__(self, save, run, event):
        self.save, self.run, self.event = save, run, event
        self.fig = None
        self.silent = False                 # True for steps before --from-step: compute only
        if not save:
            plt.ion()

    def show(self, step, fig):
        if self.silent:
            plt.close(fig)
            return "next"
        fig.suptitle(f"{self.run}: event {self.event}, step {step}: {STEPS[step - 1]}",
                     x=0.01, ha="left")
        fig.tight_layout()
        if self.save:
            PLOTS.mkdir(exist_ok=True)
            out = PLOTS / f"{self.save}_ev{self.event}_s{step}.png"
            fig.savefig(out, dpi=130)
            print(f"  saved {out}")
            plt.close(fig)
            return "next"
        if self.fig is not None:
            plt.close(self.fig)
        self.fig = fig
        fig.show()
        plt.pause(0.2)
        ans = input("  [Enter] next step, n next event, q quit: ").strip().lower()
        return {"q": "quit", "n": "skip"}.get(ans, "next")

    def close(self):
        if self.fig is not None:
            plt.close(self.fig)
            self.fig = None


def grid64(values, t=None, title=""):
    """8x8 grid in channel order: waveforms (values[c] over t) or one number per channel."""
    if t is None:                                     # map
        fig, ax = plt.subplots(figsize=(6.4, 5.6))
        img = np.full((8, 8), np.nan)
        for c, v in values.items():
            img[c // 8, c % 8] = v
        im = ax.imshow(img, cmap="viridis")
        for c, v in values.items():
            ax.text(c % 8, c // 8, f"{c}\n{v:.1f}", ha="center", va="center", fontsize=7,
                    color="w" if v < np.nanmax(img) / 2 else "k")
        fig.colorbar(im, ax=ax, label=title)
        ax.set(xticks=[], yticks=[])
        return fig
    fig, axs = plt.subplots(8, 8, figsize=(14, 9), sharex=True, sharey=True)
    for c, ax in enumerate(axs.flat):
        if c in values:
            ax.plot(t, values[c], color="C0", lw=0.5)
        ax.axvline(0, color="0.4", ls="--", lw=0.5)
        ax.text(0.03, 0.85, str(c), transform=ax.transAxes, fontsize=7)
        ax.tick_params(labelsize=6)
        ax.grid(lw=0.3)
    return fig


# ---- definitions (the ones JScan will implement) ---------------------------------------
def baseline_window(where, nbase, pre, ns):
    return slice(ns - nbase, ns) if where == "end" else slice(0, min(nbase, pre))


def pulse_search(S, nsm):
    """Smoothed sum, its maximum, height and rising edge (None if above H/2 from the start)."""
    Sh = smooth(S, nsm)
    imax = int(Sh.argmax())
    H = Sh[imax]
    below = np.flatnonzero(Sh[:imax] < H / 2)
    irise = int(below[-1]) if len(below) else None
    return Sh, imax, H, irise


def windows(imax, a, b, ns):
    """Charge window around imax and a control window of the same length in the pulse-free part."""
    lo, hi = max(0, imax - a), min(ns, imax + b + 1)
    n = hi - lo
    before, after = (0, lo), (hi, ns)
    seg = after if (after[1] - after[0]) >= (before[1] - before[0]) else before
    length = seg[1] - seg[0]
    short = length < n
    if seg is after:
        ctrl = slice(max(seg[0], seg[1] - n), seg[1])
    else:
        ctrl = slice(seg[0], min(seg[1], seg[0] + n))
    return slice(lo, hi), ctrl, short


def interp_cm(m, win, edge):
    """m with the window replaced by a line between the means of `edge` samples on each side.

    Returns (mt, L, R, lwin, rwin); lwin or rwin is None when the window touches that end of
    the record, and the other end's mean is then used for both."""
    lo, hi = win.start, win.stop
    lwin = slice(max(0, lo - edge), lo) if lo > 0 else None
    rwin = slice(hi, min(len(m), hi + edge)) if hi < len(m) else None
    L = m[lwin].mean() if lwin is not None else m[rwin].mean()
    R = m[rwin].mean() if rwin is not None else L
    mt = m.copy()
    mt[lo:hi] = np.linspace(L, R, hi - lo)
    return mt, L, R, lwin, rwin


def rms(v):
    return float(np.std(v))


# ---- steps -----------------------------------------------------------------------------
def step1(ev, st):
    print(f"step 1: {ev['run']} event {ev['e']}: trigger_id {ev['tid']}, timestamp {ev['ts'] / 1e3:.3f} us"
          + (f", {ev['dt_prev'] / 1e6:.3f} ms and {ev['ntrig_prev']} board triggers after event {ev['e'] - 1}"
             if ev["dt_prev"] is not None else "")
          + f"; pre-trigger {ev['pre']} samples, {ev['dt']:g} ns/sample, {ev['ns']} samples")
    fig = grid64({c: ev["x"][k] for k, c in enumerate(ev["channels"])}, ev["t"])
    return st.show(1, fig)


def step2(ev, st, out):
    x, t = ev["x"], ev["t"]
    m = np.median(x, axis=0)                     # one value per sample
    r = x - m
    rx, rr = np.mean([rms(v) for v in x]), np.mean([rms(v) for v in r])
    out.update(m=m)
    print(f"step 2: m[i] = median over the {len(x)} channels of x[c,i], one value per sample: "
          f"rms {rms(m):.2f} ADC; channel rms {rx:.2f}, residual (channel minus m) rms {rr:.2f}, "
          f"coherent fraction {1 - rr**2 / rx**2:.2f}; rms of m after 0.1 / 1 / 10 us smoothing: "
          + " / ".join(f"{rms(smooth(m, n)):.2f}" for n in (12, 125, 1250)))
    fig, axs = plt.subplots(3, 1, figsize=(11, 8), sharex=True)
    axs[0].plot(t, x.sum(axis=0), color="0.3", lw=0.6)
    axs[0].set_ylabel("plain sum of channels (ADC)")
    axs[1].plot(t, m, color="C0", lw=0.8)
    axs[1].set_ylabel("m: median across channels (ADC)")
    c0 = ev["channels"][0]
    axs[2].plot(t, x[0], color="0.6", lw=0.5, label=f"channel {c0}")
    axs[2].plot(t, r[0], color="C0", lw=0.5, label=f"channel {c0} minus m")
    axs[2].set(ylabel="ADC", xlabel="ns from trigger")
    for ax in axs:
        ax.axvline(0, color="0.4", ls="--", lw=0.8)
        ax.grid(lw=0.3)
    axs[2].legend(fontsize=8, loc="upper right")
    return st.show(2, fig)


def step3(ev, st, out, args):
    x, t, ns, pre, dt = ev["x"], ev["t"], ev["ns"], ev["pre"], ev["dt"]
    a, b = (int(round(v / dt)) for v in args.window)
    bwin = baseline_window(args.baseline, args.baseline_window, pre, ns)
    S0 = (x - out["m"]).sum(axis=0)
    S0 = S0 - np.median(S0[bwin])
    Sh0, imax0, H0, _ = pulse_search(S0, args.smooth)
    pwin0, cwin0, _ = windows(imax0, a, b, ns)
    n = pwin0.stop - pwin0.start
    out.update(S0=S0, imax0=imax0, pwin0=pwin0, cwin0=cwin0, bwin=bwin)
    print(f"step 3: S0 = sum over channels of (x - m), moving average over {args.smooth} samples: "
          f"maximum {H0:.1f} ADC at {t[imax0]:.0f} ns (sample {imax0}); "
          f"W = [imax - a, imax + b] with a = {a}, b = {b} samples: "
          f"[{t[pwin0.start]:.0f}, {t[pwin0.stop - 1]:.0f}] ns = samples [{pwin0.start}, {pwin0.stop - 1}], "
          f"{n} samples" + (" (clipped to the record)" if n < a + b + 1 else ""))
    fig, ax = plt.subplots(figsize=(11, 3.8))
    ax.plot(t, S0, color="0.7", lw=0.5, label="S0 = sum of (x - m)")
    ax.plot(t, Sh0, color="C0", lw=1.2, label=f"moving average {args.smooth}")
    ax.plot(t[imax0], H0, "o", color="C3", ms=6, label="maximum")
    ax.axvspan(t[pwin0.start], t[pwin0.stop - 1], color="C3", alpha=0.12, lw=0, label="W")
    ax.axvline(0, color="0.4", ls="--", lw=0.8)
    ax.axhline(0, color="0.4", lw=0.5)
    ax.set(xlabel="ns from trigger", ylabel="ADC")
    ax.grid(lw=0.3)
    ax.legend(fontsize=8, loc="upper right")
    return st.show(3, fig)


def step4(ev, st, out, args):
    m, t, pwin0 = out["m"], ev["t"], out["pwin0"]
    mt, L, R, lwin, rwin = interp_cm(m, pwin0, args.edge)
    nch = len(ev["x"])
    out.update(mt=mt)
    span = lambda w: f"samples [{w.start}, {w.stop - 1}]" if w is not None else "none (W touches the record edge)"
    print(f"step 4: L = mean of m over {args.edge} samples before W, {span(lwin)}: {L:.2f} ADC; "
          f"R = mean over {args.edge} samples after W, {span(rwin)}: {R:.2f} ADC; "
          f"mt = m outside W, straight line from L to R inside W; "
          f"jump R - L = {R - L:+.2f} ADC per channel = {nch * (R - L):+.0f} ADC in the sum")
    fig, axs = plt.subplots(2, 1, figsize=(11, 7))
    lo = max(0, (lwin.start if lwin is not None else pwin0.start) - 150)
    hi = min(len(m), (rwin.stop if rwin is not None else pwin0.stop) + 150)
    for ax, sl in zip(axs, (slice(0, len(m)), slice(lo, hi))):
        ax.plot(t[sl], m[sl], color="0.6", lw=0.6, label="m")
        ax.plot(t[sl], mt[sl], color="C0", lw=1.2, label="mt (m with the line across W)")
        ax.axvspan(t[pwin0.start], t[pwin0.stop - 1], color="C3", alpha=0.12, lw=0, label="W")
        for w, lab in ((lwin, "L region"), (rwin, "R region")):
            if w is not None:
                ax.axvspan(t[w.start], t[w.stop - 1], color="C2", alpha=0.3, lw=0, label=lab)
        ax.axvline(0, color="0.4", ls="--", lw=0.8)
        ax.set_ylabel("ADC per channel")
        ax.grid(lw=0.3)
    axs[0].set_title("whole record", fontsize=9, loc="left")
    axs[1].set_title("around W", fontsize=9, loc="left")
    axs[1].set_xlabel("ns from trigger")
    axs[0].legend(fontsize=8, loc="upper right")
    return st.show(4, fig)


def step5(ev, st, out, args):
    x, t, dt = ev["x"], ev["t"], ev["dt"]
    bwin, pwin0, cwin0 = out["bwin"], out["pwin0"], out["cwin0"]
    estimates = {"median": out["m"], "interp": out["mt"], "none": np.zeros_like(out["m"])}
    print(f"step 5: y = x - estimate for every channel, S = sum of y, minus its median over the baseline window; "
          f"per estimate, rms of S in the control window and charge in W:")
    res = {}
    for meth, cm in estimates.items():
        S = (x - cm).sum(axis=0)
        S = S - np.median(S[bwin])
        res[meth] = (S, cm)
        print(f"        {meth:7s} control rms {rms(S[cwin0]):7.1f} ADC   charge {dt * S[pwin0].sum():12.4g} ADC ns")
    S, cm = res[args.cm]
    out.update(S=S, cm=cm)
    print(f"        using {args.cm}")
    fig, ax = plt.subplots(figsize=(11, 3.8))
    ax.plot(t, res["none"][0], color="0.6", lw=0.6, label="no subtraction")
    ax.plot(t, S, color="C0", lw=0.8, label=f"estimate: {args.cm}")
    ax.axvspan(t[pwin0.start], t[pwin0.stop - 1], color="C3", alpha=0.12, lw=0, label="W")
    ax.axvline(0, color="0.4", ls="--", lw=0.8)
    ax.set(xlabel="ns from trigger", ylabel="sum of channels (ADC)")
    ax.grid(lw=0.3)
    ax.legend(fontsize=8, loc="upper right")
    return st.show(5, fig)


def step6(ev, st, out, args):
    S, t, ns, pre = out["S"], ev["t"], ev["ns"], ev["pre"]
    bwin = baseline_window(args.baseline, args.baseline_window, pre, ns)
    ewin = baseline_window("end", args.baseline_window, pre, ns)
    pwin = baseline_window("pre", args.baseline_window, pre, ns)
    B, s1 = float(np.median(S[bwin])), rms(S[bwin])
    lend, lpre = float(np.median(S[ewin])), float(np.median(S[pwin]))
    S = S - B
    out.update(S=S, B=B, s1=s1, bwin=bwin)
    print(f"step 6: baseline ({args.baseline}, {bwin.stop - bwin.start} samples) {B:.1f} ADC, "
          f"per-sample rms {s1:.1f}; end-of-record level {lend:.1f}, pre-trigger level {lpre:.1f}, "
          f"difference {lend - lpre:+.1f}")
    fig, ax = plt.subplots(figsize=(11, 3.8))
    ax.plot(t, S, color="C0", lw=0.7)
    for w, lab, lev in ((ewin, "end window", lend), (pwin, "pre-trigger window", lpre)):
        ax.axvspan(t[w.start], t[w.stop - 1], color="0.5", alpha=0.15, lw=0)
        ax.hlines(lev - B, t[0], t[-1], color="C3" if w is ewin else "C2", lw=0.8, ls=":",
                  label=f"{lab} level {lev - B:+.1f}")
    ax.axhline(0, color="0.4", lw=0.5)
    ax.axvline(0, color="0.4", ls="--", lw=0.8)
    ax.set(xlabel="ns from trigger", ylabel="sum minus baseline (ADC)")
    ax.grid(lw=0.3)
    ax.legend(fontsize=8, loc="upper right")
    return st.show(6, fig)


def step7(ev, st, out, args):
    S, t, ns, dt = out["S"], ev["t"], ev["ns"], ev["dt"]
    a, b = (int(round(v / dt)) for v in args.window)
    Sh, imax, H, irise = pulse_search(S, args.smooth)
    pwin, cwin, short = windows(imax, a, b, ns)
    sN = rms(Sh[cwin])
    out.update(Sh=Sh, imax=imax, H=H, irise=irise, pwin=pwin, cwin=cwin, short=short, sN=sN)
    print(f"step 7: smoothed ({args.smooth} samples) maximum {H:.1f} ADC at {t[imax]:.0f} ns, "
          f"{H / sN:.1f} x control rms {sN:.1f}; "
          + (f"rising edge (half maximum) at {t[irise]:.0f} ns" if irise is not None else "no rising edge"))
    fig, ax = plt.subplots(figsize=(11, 3.8))
    ax.plot(t, S, color="0.7", lw=0.5, label="sum")
    ax.plot(t, Sh, color="C0", lw=1.2, label=f"moving average {args.smooth}")
    ax.plot(t[imax], H, "o", color="C3", ms=6, label="maximum")
    if irise is not None:
        ax.axvline(t[irise], color="C2", lw=1, label="rising edge")
    ax.axvline(0, color="0.4", ls="--", lw=0.8)
    ax.axhline(0, color="0.4", lw=0.5)
    ax.set(xlabel="ns from trigger", ylabel="sum minus baseline (ADC)")
    ax.grid(lw=0.3)
    ax.legend(fontsize=8, loc="upper right")
    return st.show(7, fig)


def step8(ev, st, out, args):
    S, t, dt, s1 = out["S"], ev["t"], ev["dt"], out["s1"]
    pwin, cwin = out["pwin"], out["cwin"]
    n = pwin.stop - pwin.start
    Q, Qc = dt * S[pwin].sum(), dt * S[cwin].sum()
    sQ = s1 * np.sqrt(n) * dt
    out.update(Q=Q, Qc=Qc, sQ=sQ)
    print(f"step 8: charge window [{t[pwin.start]:.0f}, {t[pwin.stop - 1]:.0f}] ns ({n} samples): "
          f"Q {Q:.4g} ADC ns; control window [{t[cwin.start]:.0f}, {t[cwin.stop - 1]:.0f}] ns"
          + (" (shorter than the charge window)" if out["short"] else "")
          + f": {Qc:.4g}; white-noise width {sQ:.4g}; Q/width {Q / sQ:.1f}, control/width {Qc / sQ:.1f}")
    fig, ax = plt.subplots(figsize=(11, 3.8))
    ax.plot(t, S, color="C0", lw=0.7)
    ax.axvspan(t[pwin.start], t[pwin.stop - 1], color="C3", alpha=0.15, lw=0, label="charge window")
    ax.axvspan(t[cwin.start], t[cwin.stop - 1], color="0.5", alpha=0.2, lw=0, label="control window")
    ax.axhline(0, color="0.4", lw=0.5)
    ax.axvline(0, color="0.4", ls="--", lw=0.8)
    ax.set(xlabel="ns from trigger", ylabel="sum minus baseline (ADC)")
    ax.grid(lw=0.3)
    ax.legend(fontsize=8, loc="upper right")
    return st.show(8, fig)


def step9(ev, st, out, args):
    x, dt, chans = ev["x"], ev["dt"], ev["channels"]
    y = x - out["cm"]
    bwin, pwin, cwin = out["bwin"], out["pwin"], out["cwin"]
    n = pwin.stop - pwin.start
    Bc = np.median(y[:, bwin], axis=1)
    q = dt * (y[:, pwin] - Bc[:, None]).sum(axis=1)
    qc = dt * (y[:, cwin] - Bc[:, None]).sum(axis=1)       # same charge in the control window
    sq = rms(qc)                                            # one noise for all channels, drift included
    swhite = np.mean(np.std(y[:, cwin], axis=1)) * np.sqrt(n) * dt
    z = q / sq
    on = z > args.nsigma
    M = int(on.sum())
    out.update(M=M)
    print(f"step 9: per-channel charge noise {sq:.4g} ADC ns from the control window (white-noise "
          f"expectation {swhite:.4g}); multiplicity {M} channels above {args.nsigma:g} sigma, their charge "
          f"{q[on].sum():.4g} vs Q {out['Q']:.4g}; channels: {[int(c) for c in chans[on]]}")
    fig = grid64({int(c): z[k] for k, c in enumerate(chans)}, title="charge significance (sigma)")
    return st.show(9, fig)


def step10(ev, st, out, args):
    Sh, t, pwin, imax, sN = out["Sh"], ev["t"], out["pwin"], out["imax"], out["sN"]
    masked = Sh.copy()
    masked[pwin] = -np.inf
    i2 = int(masked.argmax())
    P2 = masked[i2]
    lo, hi = sorted((imax, i2))
    valley = float(Sh[lo:hi + 1].min())               # lowest point between the two maxima
    pile = P2 > args.nsigma * sN and P2 - valley > args.nsigma * sN
    real = out["Q"] >= args.nsigma * out["sQ"]
    verdict = "empty" if not real else ("real" if out["irise"] is not None else "tail")
    print(f"step 10: {verdict.upper()}" + (" + PILE-UP" if pile else "")
          + f": Q/width {out['Q'] / out['sQ']:.1f} (cut {args.nsigma:g}), "
          + ("rising edge found" if out["irise"] is not None else "no rising edge")
          + f", multiplicity {out['M']}; second maximum {P2:.1f} ADC at {t[i2]:.0f} ns "
          f"({P2 / sN:.1f} x smoothed control rms), valley between the maxima {valley:.1f} ADC, "
          f"prominence {P2 - valley:.1f} ({(P2 - valley) / sN:.1f} x)")
    return "next"


# ---- main ------------------------------------------------------------------------------
def load_event(f, e, meas):
    wf = f["waveforms"]
    A = f.attrs
    x = wf[e].astype(np.float64) - meas[:, None]
    ts, tid = int(f["timestamp_ns"][e]), int(f["trigger_id"][e])
    ev = dict(e=e, x=x, ts=ts, tid=tid, channels=f["channels"][:], ns=wf.shape[2],
              pre=int(A["pretrigger"]), dt=float(A["tsample_ns"]), dt_prev=None, ntrig_prev=None)
    ev["t"] = (np.arange(ev["ns"]) - ev["pre"]) * ev["dt"]
    if e > 0:
        ev["dt_prev"] = ts - int(f["timestamp_ns"][e - 1])
        ev["ntrig_prev"] = (tid - int(f["trigger_id"][e - 1])) % WRAP
    return ev


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("run", help="run name (128K_37V_th19) or path to a .h5 file")
    ap.add_argument("--data", default=str(DATA), help="directory of the run files")
    ap.add_argument("--events", default="0", help="range a-b or list a,b,c (0-based)")
    ap.add_argument("--cm", default="interp", choices=CM_METHODS, help="common-mode estimate")
    ap.add_argument("--baseline", default="end", choices=("end", "pre"), help="baseline window of the sum")
    ap.add_argument("--baseline-window", type=int, default=250, help="length of the baseline window (samples)")
    ap.add_argument("--smooth", type=int, default=64, help="moving-average width of the pulse search (samples)")
    ap.add_argument("--window", default="1000,4000", help="charge window: ns before,after the maximum")
    ap.add_argument("--nsigma", type=float, default=5.0, help="significance threshold")
    ap.add_argument("--edge", type=int, default=32, help="samples averaged at each end of W for the line")
    ap.add_argument("--from-step", type=int, default=1, choices=range(1, 11), help="first step shown")
    ap.add_argument("--save", metavar="NAME", help="no prompts; save figures as src/py/plots/NAME_evE_sN.png")
    args = ap.parse_args()
    args.window = tuple(float(v) for v in args.window.split(","))
    if len(args.window) != 2:
        sys.exit("--window needs two numbers: before,after in ns")

    path = find_run(args.run, Path(args.data))
    events = parse_events(args.events)
    f = h5py.File(path, "r")
    nev = f["waveforms"].shape[0]
    bad = [e for e in events if not 0 <= e < nev]
    if bad:
        sys.exit(f"events out of range 0..{nev - 1}: {bad}")
    meas = np.asarray(f.attrs["baselines"], float)

    for e in events:
        ev = load_event(f, e, meas)
        ev["run"] = path.stem
        st = Stepper(args.save, path.stem, e)
        out = {}
        steps = [lambda: step1(ev, st), lambda: step2(ev, st, out)] + [
            (lambda f: (lambda: f(ev, st, out, args)))(f)
            for f in (step3, step4, step5, step6, step7, step8, step9, step10)]
        ans = "next"
        for k, step in enumerate(steps, start=1):
            st.silent = k < args.from_step
            if st.silent:
                with contextlib.redirect_stdout(io.StringIO()):
                    ans = step()
            else:
                ans = step()
            if ans != "next":
                break
        st.close()
        if ans == "quit":
            break
        print()


if __name__ == "__main__":
    main()
