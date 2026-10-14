#!/usr/bin/env python3
"""spe.py: single-photoelectron pulses in the channels, step by step on one event.

    python3 src/py/spe.py RUN [--events 0-9] [--channel 0] [--region tail|all] [--tail 10]
            [--smooth 3] [--threshold 10] [--window 40,200] [--from-step N] [--save NAME] [--summary]

Interactive: Enter goes to the next step, n to the next event, q quits. See docs/spe_doc.tex.
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
from selection import Stepper, grid64, rms

STEPS = ("read the channel", "smooth and find candidates", "amplitude and charge",
         "all 64 channels", "rates", "accumulated spectra")
QBIN, QMAX = 50.0, 8000.0          # histogram of q: bin width and range (ADC ns)
ABIN, AMAX = 1.0, 80.0             # histogram of A (ADC)


# ---- definitions (the ones the JScan SPE pass will implement) ---------------------------
def search_region(kind, tail_us, ns, dt):
    if kind == "all":
        return slice(0, ns)
    return slice(max(0, ns - int(round(tail_us * 1000 / dt))), ns)


def find_candidates(x, region, nsm, thr, sep):
    """Local maxima of the smoothed x inside region, above thr, at least sep samples apart
    (taken in decreasing order of height). Returns (smoothed, sorted candidate indices)."""
    xs = smooth(x, nsm)
    lo, hi = region.start, region.stop
    seg = xs[lo:hi]
    i = np.arange(1, len(seg) - 1)
    loc = i[(seg[i] > thr) & (seg[i] >= seg[i - 1]) & (seg[i] > seg[i + 1])] + lo
    kept = []
    for j in loc[np.argsort(xs[loc])[::-1]]:
        if all(abs(j - k) >= sep for k in kept):
            kept.append(int(j))
    return xs, sorted(kept)


def measure(x, cands, nsm, a, b, dt):
    """Per candidate: sample of the raw maximum, amplitude A, charge q over [iA - a, iA + b]."""
    out = []
    for i0 in cands:
        lo, hi = max(0, i0 - nsm), min(len(x), i0 + nsm + 1)
        iA = lo + int(x[lo:hi].argmax())
        q = dt * x[max(0, iA - a):min(len(x), iA + b + 1)].sum()
        out.append((iA, float(x[iA]), float(q)))
    return out


def channel_pulses(x, region, args, dt):
    nsm = args.smooth
    a, b = (int(round(v / dt)) for v in args.window)
    xs, cands = find_candidates(x, region, nsm, args.threshold, a + b + 1)
    return xs, cands, measure(x, cands, nsm, a, b, dt)


def mode_of(values, binw, vmax):
    if len(values) == 0:
        return float("nan"), None, None
    h, edges = np.histogram(values, bins=np.arange(0, vmax + binw, binw))
    k = int(h.argmax())
    return float(0.5 * (edges[k] + edges[k + 1])), h, edges


def crosstalk_ratio(values, m):
    v = np.asarray(values)
    low = np.sum(v < 1.5 * m)
    return np.sum((v >= 1.5 * m) & (v < 2.5 * m)) / low if low else float("nan")


# ---- steps -----------------------------------------------------------------------------
def step1(ev, st, out, args):
    x, t, region = out["x"], ev["t"], out["region"]
    print(f"step 1: {ev['run']} event {ev['e']}, channel {args.channel}: baseline {out['b']:.1f} ADC "
          f"(median over the record), rms {rms(x):.2f}; search region {args.region}: "
          f"[{t[region.start]:.0f}, {t[region.stop - 1]:.0f}] ns = samples [{region.start}, {region.stop - 1}], "
          f"{(region.stop - region.start) * ev['dt'] / 1000:.1f} us")
    fig, ax = plt.subplots(figsize=(11, 3.6))
    ax.plot(t, x, color="C0", lw=0.5)
    ax.axvspan(t[region.start], t[region.stop - 1], color="0.5", alpha=0.15, lw=0, label="search region")
    ax.axvline(0, color="0.4", ls="--", lw=0.8)
    ax.set(xlabel="ns from trigger", ylabel=f"channel {args.channel} minus baseline (ADC)")
    ax.grid(lw=0.3)
    ax.legend(fontsize=8, loc="upper right")
    return st.show(1, fig)


def step2(ev, st, out, args):
    x, t, region = out["x"], ev["t"], out["region"]
    xs, cands, meas = channel_pulses(x, region, args, ev["dt"])
    out.update(xs=xs, cands=cands, meas=meas)
    print(f"step 2: moving average over {args.smooth} samples, local maxima above {args.threshold:g} ADC, "
          f"isolation {sum(args.window):g} ns: {len(cands)} candidates"
          + ("; (time ns, smoothed height): " + ", ".join(f"({t[i]:.0f}, {xs[i]:.1f})" for i in cands)
             if cands else ""))
    sl = region
    fig, ax = plt.subplots(figsize=(11, 3.6))
    ax.plot(t[sl], x[sl], color="0.7", lw=0.5, label="raw")
    ax.plot(t[sl], xs[sl], color="C0", lw=0.9, label=f"moving average {args.smooth}")
    ax.axhline(args.threshold, color="C3", lw=0.8, ls=":", label="threshold")
    if cands:
        ax.plot(t[cands], xs[cands], "v", color="C3", ms=7, label="candidates")
    ax.set(xlabel="ns from trigger", ylabel="ADC")
    ax.grid(lw=0.3)
    ax.legend(fontsize=8, loc="upper right")
    return st.show(2, fig)


def step3(ev, st, out, args):
    x, t, dt = out["x"], ev["t"], ev["dt"]
    a, b = (int(round(v / dt)) for v in args.window)
    meas = out["meas"]
    print(f"step 3: charge window [-{args.window[0]:g}, +{args.window[1]:g}] ns around the raw maximum; "
          f"per candidate (time ns, A ADC, q ADC ns, q/A ns):")
    for iA, A, q in meas:
        print(f"        ({t[iA]:8.0f}, {A:5.1f}, {q:8.0f}, {q / A if A else float('nan'):6.1f})")
    n = min(6, len(meas))
    fig, axs = plt.subplots(1, max(n, 1), figsize=(3.2 * max(n, 1) + 1, 3.4), squeeze=False)
    for ax, (iA, A, q) in zip(axs.flat, meas[:n]):
        lo, hi = max(0, iA - a - 12), min(len(x), iA + b + 13)
        ax.plot(t[lo:hi], x[lo:hi], color="C0", lw=0.8)
        ax.axvspan(t[max(0, iA - a)], t[min(len(x) - 1, iA + b)], color="C3", alpha=0.12, lw=0)
        ax.plot(t[iA], A, "o", color="C3", ms=5)
        ax.set_title(f"A {A:.0f}, q {q:.0f}", fontsize=9, loc="left")
        ax.set_xlabel("ns", fontsize=8)
        ax.tick_params(labelsize=7)
        ax.grid(lw=0.3)
    if n == 0:
        axs.flat[0].text(0.5, 0.5, "no candidates", ha="center", transform=axs.flat[0].transAxes)
    axs.flat[0].set_ylabel("ADC", fontsize=8)
    return st.show(3, fig)


def all_channels(ev, out, args):
    """Steps 1 to 3 on every channel; fills out['per'] = {channel: (x, xs, cands, meas)}."""
    if "per" in out:
        return
    per = {}
    for k, c in enumerate(ev["channels"]):
        x = ev["w"][k] - np.median(ev["w"][k])
        xs, cands, meas = channel_pulses(x, out["region"], args, ev["dt"])
        per[int(c)] = (x, xs, cands, meas)
    out["per"] = per


def step4(ev, st, out, args):
    all_channels(ev, out, args)
    per, t, region = out["per"], ev["t"], out["region"]
    counts = {c: len(v[2]) for c, v in per.items()}
    tot = sum(counts.values())
    print(f"step 4: candidates per channel (8 x 8, channel order); total {tot}, mean {tot / len(per):.2f}")
    for r in range(8):
        print("        " + " ".join(f"{counts.get(8 * r + k, 0):3d}" for k in range(8)))
    fig, axs = plt.subplots(8, 8, figsize=(14, 9), sharex=True, sharey=True)
    for c, ax in enumerate(axs.flat):
        if c in per:
            x, xs, cands, _ = per[c]
            ax.plot(t[region], xs[region], color="C0", lw=0.4)
            if cands:
                ax.plot(t[cands], xs[cands], "v", color="C3", ms=4)
        ax.text(0.03, 0.82, str(c), transform=ax.transAxes, fontsize=7)
        ax.tick_params(labelsize=6)
        ax.grid(lw=0.3)
    return st.show(4, fig)


def step5(ev, st, out, args):
    all_channels(ev, out, args)
    per, region, dt = out["per"], out["region"], ev["dt"]
    TR = (region.stop - region.start) * dt * 1e-9
    rate = {c: len(v[2]) / TR for c, v in per.items()}
    vals = np.array(list(rate.values()))
    m, s = vals.mean(), vals.std()
    far = [c for c, r in rate.items() if abs(r - m) > 3 * s] if s else []
    out["TR"] = TR
    print(f"step 5: search time {TR * 1e6:.1f} us per channel; rate mean {m / 1e3:.1f} kHz, "
          f"spread {s / 1e3:.1f} kHz; channels beyond 3 sigma: {far}")
    fig = grid64({c: r / 1e3 for c, r in rate.items()}, title="dark rate (kHz), this event")
    return st.show(5, fig)


def step6(ev, st, out, args, acc):
    all_channels(ev, out, args)
    for c, (_, _, _, meas) in out["per"].items():
        for _, A, q in meas:
            acc["q"][c].append(q)
            acc["A"][c].append(A)
    acc["events"] += 1
    c0 = args.channel
    sets = (("q", f"channel {c0}", acc["q"][c0], QBIN, QMAX, "q (ADC ns)"),
            ("A", f"channel {c0}", acc["A"][c0], ABIN, AMAX, "A (ADC)"),
            ("q", "all channels", sum(acc["q"].values(), []), QBIN, QMAX, "q (ADC ns)"),
            ("A", "all channels", sum(acc["A"].values(), []), ABIN, AMAX, "A (ADC)"))
    print(f"step 6: spectra over {acc['events']} events:")
    fig, axs = plt.subplots(2, 2, figsize=(11, 7))
    for ax, (what, who, vals, binw, vmax, lab) in zip(axs.flat, sets):
        m, h, edges = mode_of(vals, binw, vmax)
        xt = crosstalk_ratio(vals, m) if len(vals) else float("nan")
        print(f"        {what} {who:13s}: {len(vals):6d} pulses, mode {m:7.1f}, "
              f"pulses in [1.5, 2.5] x mode over pulses below: {xt:.3f}")
        if h is not None:
            ax.stairs(h, edges, color="C0", fill=True, alpha=0.6)
            ax.axvline(m, color="C3", lw=1, ls="--", label=f"mode {m:.0f}")
            ax.legend(fontsize=8, loc="upper right")
        ax.set(xlabel=lab, ylabel="pulses", title=f"{what}, {who}, {len(vals)} pulses")
        ax.title.set_fontsize(9)
        ax.grid(lw=0.3)
    return st.show(6, fig)


# ---- main ------------------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("run", help="run name (128K_37V_th19) or path to a .h5 file")
    ap.add_argument("--data", default=str(DATA), help="directory of the run files")
    ap.add_argument("--events", default="0", help="range a-b or list a,b,c (0-based)")
    ap.add_argument("--channel", type=int, default=0, help="channel shown in steps 1-3 and 6")
    ap.add_argument("--region", default="tail", choices=("tail", "all"), help="search region")
    ap.add_argument("--tail", type=float, default=10.0, help="length of the tail region (us)")
    ap.add_argument("--smooth", type=int, default=3, help="moving-average width (samples)")
    ap.add_argument("--threshold", type=float, default=10.0, help="candidate threshold (ADC)")
    ap.add_argument("--window", default="40,200", help="charge window: ns before,after the maximum")
    ap.add_argument("--from-step", type=int, default=1, choices=range(1, 7), help="first step shown")
    ap.add_argument("--save", metavar="NAME", help="no prompts; save figures as src/py/plots/NAME_evE_sN.png")
    ap.add_argument("--summary", action="store_true", help="only the final step 6, after the last event")
    args = ap.parse_args()
    args.window = tuple(float(v) for v in args.window.split(","))
    if len(args.window) != 2:
        sys.exit("--window needs two numbers: before,after in ns")

    path = find_run(args.run, Path(args.data))
    events = parse_events(args.events)
    f = h5py.File(path, "r")
    wf = f["waveforms"]
    nev, nch, ns = wf.shape
    bad = [e for e in events if not 0 <= e < nev]
    if bad:
        sys.exit(f"events out of range 0..{nev - 1}: {bad}")
    channels = f["channels"][:]
    if args.channel not in channels:
        sys.exit(f"channel {args.channel} not in the file")
    dt, pre = float(f.attrs["tsample_ns"]), int(f.attrs["pretrigger"])
    t = (np.arange(ns) - pre) * dt
    region = search_region(args.region, args.tail, ns, dt)
    acc = {"q": {int(c): [] for c in channels}, "A": {int(c): [] for c in channels}, "events": 0}

    ans = "next"
    for e in events:
        w = wf[e].astype(np.float64)
        ev = dict(run=path.stem, e=e, w=w, t=t, dt=dt, channels=channels)
        k = int(np.flatnonzero(channels == args.channel)[0])
        out = dict(region=region, b=float(np.median(w[k])))
        out["x"] = w[k] - out["b"]
        st = Stepper(args.save, path.stem, e, names=STEPS)
        last = e == events[-1]
        steps = [lambda: step1(ev, st, out, args), lambda: step2(ev, st, out, args),
                 lambda: step3(ev, st, out, args), lambda: step4(ev, st, out, args),
                 lambda: step5(ev, st, out, args), lambda: step6(ev, st, out, args, acc)]
        for n, step in enumerate(steps, start=1):
            silent = n < args.from_step or (args.summary and not (n == 6 and last))
            st.silent = silent
            if silent:
                with contextlib.redirect_stdout(io.StringIO()):
                    ans = step()
            else:
                if args.summary and n == 6:                       # final figure of the summary mode
                    st.save = None if not args.save else args.save
                ans = step()
            if ans != "next":
                break
        st.close()
        if ans == "quit":
            break
        if not args.summary:
            print()
    if args.summary and args.save:
        out_path = PLOTS / f"{args.save}_spe.png"
        src = PLOTS / f"{args.save}_ev{events[-1]}_s6.png"
        if src.exists():
            src.rename(out_path)
            print(f"saved {out_path}")


if __name__ == "__main__":
    main()
