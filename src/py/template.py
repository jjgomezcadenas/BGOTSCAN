#!/usr/bin/env python3
"""template.py: the shape of the single-photoelectron pulse, step by step.

    python3 src/py/template.py RUN [--events 0-199] [--channel 0] [--region tail|all] [--tail 10]
            [--threshold 12] [--half 8] [--from-step N] [--save NAME] [--summary] [--write FILE]

Interactive: Enter goes to the next step, n to the next event, q quits. See docs/template_doc.tex.
"""
import argparse
import contextlib
import io
import sys
from pathlib import Path

import h5py
import numpy as np
import matplotlib.pyplot as plt

from wfview import DATA, PLOTS, find_run, parse_events
from selection import Stepper, grid64, rms
from spe import search_region, find_candidates

STEPS = ("candidates on the channel", "aligned pulses of the channel", "all channels, accumulated",
         "first template", "refined template", "per-channel check")


# ---- definitions -----------------------------------------------------------------------
def clipped_rms(v, nsig=3.0):
    """rms after removing samples beyond nsig times a first rms (pulses excluded from the noise)."""
    s = rms(v)
    return rms(v[np.abs(v - v.mean()) < nsig * s]) if s else 0.0


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


def matched(x, T):
    """Least-squares amplitude of a pulse of shape T at every sample: y_i = sum_j T_j x_{i+j} / sum T^2.
    Returns y aligned with x (edges set to zero)."""
    H = (len(T) - 1) // 2
    y = np.zeros_like(x)
    y[H:len(x) - H] = np.correlate(x, T, mode="valid") / np.sum(T * T)
    return y


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
            P = aligned(x, cands, H)
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
    print(f"step 5: filtered with T1, candidates above {nsig:.1f} sigma of the filtered noise, aligned on the "
          f"filtered maximum: {len(allp)} pulses, peak {mean[H]:.2f} ADC\n"
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
    ap.add_argument("--from-step", type=int, default=1, choices=range(1, 7), help="first step shown")
    ap.add_argument("--save", metavar="NAME", help="no prompts; save figures as src/py/plots/NAME_evE_sN.png")
    ap.add_argument("--summary", action="store_true", help="steps 1-2 silent; steps 3-6 after the last event")
    ap.add_argument("--write", metavar="FILE", help="write the refined template as CSV")
    args = ap.parse_args()

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
                 lambda: step5(ev, st, out, args, acc), lambda: step6(ev, st, out, args, acc)]
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
        for n in range(3, 7):
            src = PLOTS / f"{args.save}_ev{events[-1]}_s{n}.png"
            if src.exists():
                src.rename(PLOTS / f"{args.save}_s{n}.png")
        print(f"figures saved as {PLOTS}/{args.save}_s3..6.png")
    if args.write and acc.get("T2") is not None:
        H = args.half
        s = template_stats(acc["T2"])
        with open(args.write, "w") as fo:
            fo.write(f"# SPE template of {path.stem}: {acc['n2']} pulses, threshold {args.threshold:g} ADC, "
                     f"{acc['events']} events, sum of squares {s['ssq']:.4f}\n")
            fo.write("offset_samples,value\n")
            for j, v in zip(range(-H, H + 1), acc["T2"]):
                fo.write(f"{j},{v:.5f}\n")
        print(f"template written to {args.write}")


if __name__ == "__main__":
    main()
