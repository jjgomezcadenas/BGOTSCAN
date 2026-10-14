#!/usr/bin/env python3
"""template.py: the shape of the single-photoelectron pulse, step by step.

    python3 src/py/template.py RUN [--events 0-199] [--channel 0] [--region tail|all] [--tail 10]
            [--threshold 12] [--half 8] [--from-step N] [--save NAME] [--summary] [--write FILE]
            [--dead LIST] [--lags 32] [--pad 8] [--write-filter FILE]

Interactive: Enter goes to the next step, n to the next event, q quits. See docs/manual/template_doc.tex.
"""
import argparse
import contextlib
import io
import sys
from pathlib import Path

import h5py
import numpy as np
import matplotlib.pyplot as plt

from scipy.linalg import toeplitz

from bgo.io import DATA, PLOTS, find_run, parse_events
from bgo.signal import rms, clipped_rms, search_region, find_candidates, matched
from bgo.plot import Stepper, grid64

STEPS = ("candidates on the channel", "aligned pulses of the channel", "all channels, accumulated",
         "first template", "refined template", "per-channel check", "noise model and whitened filter")
DEAD = "19,20,21,22,23,25,27,33,35,41,46,47"      # chanstat_doc.tex, 2026-10-14


# ---- definitions -----------------------------------------------------------------------
def aligned(x, cands, H, core=2):
    """Raw samples from H before to H after each candidate, minus a local baseline: the median of the
    window samples more than `core` samples before the peak (a tail cannot reach them). Candidates too
    close to an edge are skipped."""
    out = []
    far = np.arange(-H, H + 1) < -core
    for i in cands:
        if i - H >= 0 and i + H + 1 <= len(x):
            p = x[i - H:i + H + 1]
            out.append(p - np.median(p[far]))
    return np.array(out) if out else np.zeros((0, 2 * H + 1))


def template_stats(T):
    H = (len(T) - 1) // 2
    j = np.arange(-H, H + 1)
    above = np.flatnonzero(T >= 0.5)
    lo, hi = above[0], above[-1]
    # linear interpolation of the half-maximum crossings
    left = j[lo] - (T[lo] - 0.5) / (T[lo] - T[lo - 1]) if lo > 0 else j[lo]
    right = j[hi] + (T[hi] - 0.5) / (T[hi] - T[hi + 1]) if hi < len(T) - 1 else j[hi]
    return dict(fwhm=float(right - left), undershoot=float(T[H + 1:].min()), ssq=float(np.sum(T * T)))


def fmt_template(T):
    return " ".join(f"{v:+.3f}" for v in T)


# ---- steps -----------------------------------------------------------------------------
def step1(ev, st, out, args):
    x, t, region = out["x"], ev["t"], out["region"]
    _, cands = find_candidates(x, region, 1, args.threshold, 2 * args.half)
    out["cands"] = cands
    print(f"step 1: {ev['run']} event {ev['e']}, channel {args.channel}: raw local maxima above "
          f"{args.threshold:g} ADC in the search region [{t[region.start]:.0f}, {t[region.stop - 1]:.0f}] ns, "
          f"at least {2 * args.half} samples apart: {len(cands)} candidates"
          + ("; (time ns, A): " + ", ".join(f"({t[i]:.0f}, {x[i]:.0f})" for i in cands) if cands else ""))
    fig, ax = plt.subplots(figsize=(11, 3.6))
    ax.plot(t[region], x[region], color="C0", lw=0.5)
    ax.axhline(args.threshold, color="C3", lw=0.8, ls=":", label="threshold")
    if cands:
        ax.plot(t[cands], x[cands], "v", color="C3", ms=7, label="candidates")
    ax.set(xlabel="ns from trigger", ylabel=f"channel {args.channel} minus baseline (ADC)")
    ax.grid(lw=0.3)
    ax.legend(fontsize=8, loc="upper right")
    return st.show(1, fig)


def step2(ev, st, out, args):
    H = args.half
    P = aligned(out["x"], out["cands"], H)
    out["P"] = P
    j = np.arange(-H, H + 1)
    print(f"step 2: {len(P)} pulses aligned on their maximum, {2 * H + 1} samples each")
    fig, ax = plt.subplots(figsize=(7, 4))
    for p in P:
        ax.plot(j, p, color="C0", lw=0.6, alpha=0.5)
    if len(P):
        ax.plot(j, P.mean(axis=0), color="C3", lw=2, label=f"mean of {len(P)}")
        ax.legend(fontsize=8)
    ax.set(xlabel="sample offset from the maximum", ylabel="ADC")
    ax.grid(lw=0.3)
    return st.show(2, fig)


def accumulate(ev, out, args, acc):
    """Steps 1 and 2 on all channels of this event; store the pulses and the region waveforms."""
    H, region = args.half, out["region"]
    for k, c in enumerate(ev["channels"]):
        x = ev["w"][k] - np.median(ev["w"][k])
        _, cands = find_candidates(x, region, 1, args.threshold, 2 * H)
        P = aligned(x, cands, H)
        if len(P):
            acc["pulses"][int(c)].append(P)
        acc["waves"][int(c)].append(x[region].astype(np.float32))
    acc["events"] += 1


def all_pulses(acc, key="pulses"):
    per = {c: (np.concatenate(v) if v else np.zeros((0, 2 * acc["H"] + 1))) for c, v in acc[key].items()}
    allp = np.concatenate([p for p in per.values() if len(p)]) if any(len(p) for p in per.values()) else None
    return per, allp


def step3(ev, st, out, args, acc):
    H = args.half
    per, allp = all_pulses(acc)
    counts = {c: len(p) for c, p in per.items()}
    tot = sum(counts.values())
    print(f"step 3: pulses accumulated over {acc['events']} events, per channel (8 x 8); total {tot}")
    for r in range(8):
        print("        " + " ".join(f"{counts.get(8 * r + k, 0):4d}" for k in range(8)))
    j = np.arange(-H, H + 1)
    fig, ax = plt.subplots(figsize=(8, 4.5))
    if allp is not None:
        jj = np.tile(j, len(allp))
        ax.hist2d(jj, allp.ravel(), bins=[np.arange(-H - 0.5, H + 1.5), np.arange(-20, 50, 1)], cmap="Blues")
        ax.plot(j, allp.mean(axis=0), "o-", color="C3", ms=4, lw=1.5, label=f"mean of {len(allp)} pulses")
        ax.legend(fontsize=8, loc="upper right")
    ax.set(xlabel="sample offset from the maximum", ylabel="ADC")
    ax.grid(lw=0.3)
    return st.show(3, fig)


def step4(ev, st, out, args, acc):
    H = args.half
    _, allp = all_pulses(acc)
    if allp is None or len(allp) < 2:
        print("step 4: no pulses accumulated")
        return st.show(4, plt.figure(figsize=(4, 3)))
    mean = allp.mean(axis=0)
    err = allp.std(axis=0) / np.sqrt(len(allp))
    T1 = mean / mean[H]
    s = template_stats(T1)
    acc["T1"] = T1
    j = np.arange(-H, H + 1)
    print(f"step 4: first template from {len(allp)} pulses, peak {mean[H]:.2f} +- {err[H]:.2f} ADC:\n"
          f"        T1 = {fmt_template(T1)}\n"
          f"        FWHM {s['fwhm']:.2f} samples ({8 * s['fwhm']:.0f} ns), undershoot {s['undershoot']:+.3f}, "
          f"sum of squares {s['ssq']:.3f} (matched-filter gain {np.sqrt(s['ssq']):.2f})")
    fig, ax = plt.subplots(figsize=(7, 4))
    ax.errorbar(j, T1, yerr=err / mean[H], fmt="o-", color="C0", ms=4, lw=1.2, capsize=2)
    ax.axhline(0, color="0.4", lw=0.5)
    ax.set(xlabel="sample offset from the maximum", ylabel="template (peak = 1)", title="first template")
    ax.grid(lw=0.3)
    return st.show(4, fig)


def step5(ev, st, out, args, acc):
    H = args.half
    T1 = acc.get("T1")
    if T1 is None:
        print("step 5: no first template")
        return st.show(5, plt.figure(figsize=(4, 3)))
    per2 = {}
    nsig = None
    for c, waves in acc["waves"].items():
        pulses = []
        for x in waves:
            x = x.astype(np.float64)
            y = matched(x, T1)
            sx, sy = clipped_rms(x), clipped_rms(y[H:len(y) - H])
            nsig = args.threshold / sx if sx else 3.0
            _, cands = find_candidates(y, slice(H, len(y) - H), 1, nsig * sy, 2 * H)
            # align on the raw maximum within one sample of the filtered maximum
            iA = [i - 1 + int(x[i - 1:i + 2].argmax()) for i in cands]
            P = aligned(x, iA, H)
            if len(P):
                pulses.append(P)
        per2[c] = np.concatenate(pulses) if pulses else np.zeros((0, 2 * H + 1))
    acc["pulses2"] = per2
    allp = np.concatenate([p for p in per2.values() if len(p)])
    mean = allp.mean(axis=0)
    T2 = mean / mean[H]
    acc["T2"], acc["n2"] = T2, len(allp)
    s1, s2 = template_stats(T1), template_stats(T2)
    j = np.arange(-H, H + 1)
    print(f"step 5: filtered with T1, candidates above {nsig:.1f} sigma of the filtered "
          f"noise, aligned on the raw maximum within one sample: {len(allp)} pulses, peak {mean[H]:.2f} ADC\n"
          f"        T2 = {fmt_template(T2)}\n"
          f"        T2 - T1 = {fmt_template(T2 - T1)}\n"
          f"        FWHM {s2['fwhm']:.2f} samples, undershoot {s2['undershoot']:+.3f}, "
          f"sum of squares {s2['ssq']:.3f} (T1: {s1['ssq']:.3f})")
    fig, ax = plt.subplots(figsize=(7, 4))
    ax.plot(j, T1, "o--", color="0.5", ms=4, lw=1, label="first template")
    ax.plot(j, T2, "o-", color="C0", ms=4, lw=1.5, label="refined template")
    ax.axhline(0, color="0.4", lw=0.5)
    ax.set(xlabel="sample offset from the maximum", ylabel="template (peak = 1)")
    ax.grid(lw=0.3)
    ax.legend(fontsize=8)
    return st.show(5, fig)


def step6(ev, st, out, args, acc):
    H = args.half
    T2, per2 = acc.get("T2"), acc.get("pulses2")
    if T2 is None:
        print("step 6: no refined template")
        return st.show(6, plt.figure(figsize=(4, 3)))
    res = {}
    for c, P in per2.items():
        if len(P) >= 5 and P.mean(axis=0)[H] > 0:
            m = P.mean(axis=0)
            res[c] = float(rms(m / m[H] - T2 + (m / m[H] - T2).mean() * 0))
    med = float(np.median(list(res.values()))) if res else float("nan")
    far = sorted(c for c, r in res.items() if r > 3 * med)
    print(f"step 6: per-channel mean pulse against T2: residual rms over {2 * H + 1} samples, median {med:.3f}; "
          f"channels above 3 x median: {far}; channels with fewer than 5 pulses: "
          f"{sorted(c for c in per2 if c not in res)}")
    fig = grid64(res, title="residual to the template (rms, peak units)")
    return st.show(6, fig)


def noise_model(segments, K):
    """Autocorrelation R[0..K] and averaged periodogram of pure-noise segments (mean removed)."""
    R = np.zeros(K + 1)
    n = 0
    spec = None
    for x in segments:
        x = x - x.mean()
        L = len(x)
        for k in range(K + 1):
            R[k] += np.dot(x[:L - k], x[k:]) / (L - k)
        P = np.abs(np.fft.rfft(x)) ** 2 / L
        spec = P if spec is None else spec + P
        n += 1
    return R / n, spec / n


def whitened_filter(T, R, pad):
    """Best linear unbiased estimate of the height of a pulse T in noise of autocorrelation R:
    w = C^-1 Tp / (Tp' C^-1 Tp), window = T padded with `pad` zeros on each side. Returns (w, Tp, sigma_y)."""
    Tp = np.concatenate([np.zeros(pad), T, np.zeros(pad)])
    L = len(Tp)
    if len(R) < L:
        R = np.concatenate([R, np.zeros(L - len(R))])
    C = toeplitz(R[:L])
    CinvT = np.linalg.solve(C, Tp)
    denom = Tp @ CinvT
    return CinvT / denom, Tp, float(np.sqrt(1.0 / denom))


def step7(ev, st, out, args, acc):
    T = acc.get("T1")
    if T is None:
        print("step 7: no template")
        return st.show(7, plt.figure(figsize=(4, 3)))
    dead = [c for c in args.dead if c in acc["waves"] and acc["waves"][c]]
    segs = [x.astype(np.float64) for c in dead for x in acc["waves"][c]]
    if not segs:
        print("step 7: no dead-channel segments stored")
        return st.show(7, plt.figure(figsize=(4, 3)))
    K, H, pad, dt = args.lags, args.half, args.pad, ev["dt"]
    R, spec = noise_model(segs, K)
    rho = R / R[0]
    sigma = np.sqrt(R[0])
    below = np.flatnonzero(rho < np.exp(-1))
    corr_len = int(below[0]) if len(below) else K
    L = len(segs[0])
    freq = np.fft.rfftfreq(L, d=dt * 1e-9) / 1e6
    cum = np.cumsum(spec) / spec.sum()
    frac = {f: float(cum[np.searchsorted(freq, f)]) for f in (10, 20, 30)}
    # predicted noises
    ssq = np.sum(T * T)
    Cm = toeplitz(R[:len(T)])
    sig_m = float(np.sqrt(T @ Cm @ T) / ssq)
    w, Tp, sig_w = whitened_filter(T, R, pad)
    # measured on the dead channels
    ym = np.concatenate([matched(x, T)[H:len(x) - H] for x in segs[:400]])
    Hp = (len(w) - 1) // 2
    yw = np.concatenate([np.correlate(x, w, mode="valid") for x in segs[:400]])
    mm, mw = clipped_rms(ym), clipped_rms(yw)
    acc.update(w=w, R=R, sig_w=sig_w)
    print(f"step 7: noise from {len(dead)} dead channels, {len(segs)} segments of {L} samples: rms {sigma:.2f} ADC; "
          f"autocorrelation rho(1..8) = " + " ".join(f"{v:+.2f}" for v in rho[1:9])
          + f"; correlation length {corr_len} samples; noise power below 10 / 20 / 30 MHz: "
          + " / ".join(f"{100 * frac[f]:.0f}%" for f in (10, 20, 30)))
    print(f"        predicted noise: single sample {sigma:.2f}, matched filter {sig_m:.2f} (white-noise expectation "
          f"{sigma / np.sqrt(ssq):.2f}), whitened filter ({len(w)} samples) {sig_w:.2f}; "
          f"measured on the dead channels: matched {mm:.2f}, whitened {mw:.2f}; "
          f"gains over a single sample: matched {sigma / mm:.2f}, whitened {sigma / mw:.2f}")
    fig, axs = plt.subplots(1, 3, figsize=(15, 4))
    axs[0].plot(np.arange(K + 1), rho, "o-", color="C0", ms=3)
    axs[0].axhline(0, color="0.4", lw=0.5)
    axs[0].set(xlabel="lag (samples)", ylabel="autocorrelation", title="noise autocorrelation")
    tspec = np.abs(np.fft.rfft(np.concatenate([T, np.zeros(L - len(T))]))) ** 2
    axs[1].plot(freq[1:], spec[1:] / spec[1:].max(), color="C0", lw=0.8, label="noise (dead channels)")
    axs[1].plot(freq[1:], tspec[1:] / tspec[1:].max(), color="C3", lw=1.2, label="template")
    axs[1].set(xlabel="frequency (MHz)", ylabel="power (peak = 1)", title="power spectra", yscale="log",
               ylim=(1e-3, 2))
    axs[1].legend(fontsize=8)
    j = np.arange(-Hp, Hp + 1)
    axs[2].plot(j, np.concatenate([np.zeros(pad), T / ssq, np.zeros(pad)]), "o--", color="0.5", ms=3, label="matched: T / sum T^2")
    axs[2].plot(j, w, "o-", color="C0", ms=3, label="whitened")
    axs[2].axhline(0, color="0.4", lw=0.5)
    axs[2].set(xlabel="sample offset", ylabel="weight", title="filter weights")
    axs[2].legend(fontsize=8)
    for ax in axs:
        ax.grid(lw=0.3)
    return st.show(7, fig)


# ---- main ------------------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("run", help="run name (128K_37V_th19) or path to a .h5 file")
    ap.add_argument("--data", default=str(DATA), help="directory of the run files")
    ap.add_argument("--events", default="0", help="range a-b or list a,b,c (0-based)")
    ap.add_argument("--channel", type=int, default=0, help="channel shown in steps 1 and 2")
    ap.add_argument("--region", default="tail", choices=("tail", "all"), help="search region")
    ap.add_argument("--tail", type=float, default=10.0, help="length of the tail region (us)")
    ap.add_argument("--threshold", type=float, default=12.0, help="candidate threshold on the raw waveform (ADC)")
    ap.add_argument("--half", type=int, default=8, help="half length of the template (samples)")
    ap.add_argument("--from-step", type=int, default=1, choices=range(1, 8), help="first step shown")
    ap.add_argument("--save", metavar="NAME", help="no prompts; save figures as src/py/plots/NAME_evE_sN.png")
    ap.add_argument("--summary", action="store_true", help="steps 1-2 silent; steps 3-6 after the last event")
    ap.add_argument("--write", metavar="FILE", help="write the first template as CSV")
    ap.add_argument("--dead", default=DEAD, help="dead channels (noise sample of step 7), comma list")
    ap.add_argument("--lags", type=int, default=32, help="lags of the noise autocorrelation (samples)")
    ap.add_argument("--pad", type=int, default=8, help="zeros padded on each side of the template (step 7)")
    ap.add_argument("--write-filter", metavar="FILE", help="write the whitened filter weights as CSV")
    args = ap.parse_args()
    args.dead = [int(v) for v in args.dead.split(",") if v.strip()]

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
    acc = {"pulses": {int(c): [] for c in channels}, "waves": {int(c): [] for c in channels},
           "events": 0, "H": args.half}

    ans = "next"
    for e in events:
        w = wf[e].astype(np.float64)
        ev = dict(run=path.stem, e=e, w=w, t=t, dt=dt, channels=channels)
        k = int(np.flatnonzero(channels == args.channel)[0])
        out = dict(region=region)
        out["x"] = w[k] - np.median(w[k])
        st = Stepper(args.save, path.stem, e, names=STEPS)
        last = e == events[-1]
        accumulate(ev, out, args, acc)
        steps = [lambda: step1(ev, st, out, args), lambda: step2(ev, st, out, args),
                 lambda: step3(ev, st, out, args, acc), lambda: step4(ev, st, out, args, acc),
                 lambda: step5(ev, st, out, args, acc), lambda: step6(ev, st, out, args, acc),
                 lambda: step7(ev, st, out, args, acc)]
        for n, step in enumerate(steps, start=1):
            if args.summary and not (last and n >= 3):
                continue
            st.silent = n < args.from_step
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
        if not args.summary:
            print()
    if args.summary and args.save:
        for n in range(3, 8):
            src = PLOTS / f"{args.save}_ev{events[-1]}_s{n}.png"
            if src.exists():
                src.rename(PLOTS / f"{args.save}_s{n}.png")
        print(f"figures saved as {PLOTS}/{args.save}_s3..7.png")
    if args.write and acc.get("T1") is not None:
        H = args.half
        _, allp = all_pulses(acc)
        s = template_stats(acc["T1"])
        with open(args.write, "w") as fo:
            fo.write(f"# SPE template (first template, raw-maximum selection) of {path.stem}: {len(allp)} pulses, "
                     f"threshold {args.threshold:g} ADC, {acc['events']} events, sum of squares {s['ssq']:.4f}, "
                     f"FWHM {s['fwhm']:.2f} samples\n")
            fo.write("offset_samples,value\n")
            for j, v in zip(range(-H, H + 1), acc["T1"]):
                fo.write(f"{j},{v:.5f}\n")
        print(f"template written to {args.write}")


    if args.write_filter and acc.get("w") is not None:
        w = acc["w"]
        Hp = (len(w) - 1) // 2
        with open(args.write_filter, "w") as fo:
            fo.write(f"# whitened matched filter of {path.stem}: noise rms {np.sqrt(acc['R'][0]):.3f} ADC from "
                     f"channels {args.dead}, {args.lags} lags, pad {args.pad}, predicted sigma_y {acc['sig_w']:.3f}\n")
            fo.write("offset_samples,weight\n")
            for j, v in zip(range(-Hp, Hp + 1), w):
                fo.write(f"{j},{v:.6f}\n")
        print(f"filter written to {args.write_filter}")


if __name__ == "__main__":
    main()
