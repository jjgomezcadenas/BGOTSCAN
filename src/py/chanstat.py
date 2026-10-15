#!/usr/bin/env python3
"""chanstat.py: the 64 channels over many events: noise, dark-pulse candidates, lit fraction, lit charge.

    python3 src/py/chanstat.py RUN [--events 0-499] [--tail 10] [--thresholds 12,16,20] [--smooth 64]
            [--window 1000,4000] [--nsigma 5] [--save NAME] [--csv FILE]

One table and one figure of four 8x8 maps. See docs/manual/chanstat_doc.tex.
"""
import argparse
import csv
import sys
from pathlib import Path

import h5py
import numpy as np
import matplotlib.pyplot as plt

from bgo.io import DATA, find_run, parse_events, read_map, plot_dir
from bgo.signal import smooth, rms, clipped_rms, search_region, find_candidates
from bgo.plot import map8


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("run", help="run name (128K_37V_th19) or path to a .h5 file")
    ap.add_argument("--data", default=str(DATA), help="directory of the run files")
    ap.add_argument("--events", default="0-499", help="range a-b or list a,b,c (0-based)")
    ap.add_argument("--tail", type=float, default=10.0, help="search region: last T us of the record")
    ap.add_argument("--thresholds", default="12,16,20", help="dark-pulse thresholds (ADC)")
    ap.add_argument("--smooth", type=int, default=64, help="moving-average width of the pulse search (samples)")
    ap.add_argument("--window", default="1000,4000", help="charge window: ns before,after the maximum")
    ap.add_argument("--nsigma", type=float, default=5.0, help="lit if charge > nsigma x own charge noise")
    ap.add_argument("--dead", default="19,20,21,22,23,25,27,33,35,41,46,47", help="dead channels, comma list")
    ap.add_argument("--save", metavar="NAME", help="save the figures as src/py/plots/NAME_chanstat.png and NAME_wander.png")
    ap.add_argument("--csv", metavar="FILE", help="write the table to FILE")
    args = ap.parse_args()
    thr = [float(v) for v in args.thresholds.split(",")]
    dead = {int(v) for v in args.dead.split(",") if v.strip()}
    a_ns, b_ns = (float(v) for v in args.window.split(","))

    path = find_run(args.run, Path(args.data))
    events = parse_events(args.events)
    f = h5py.File(path, "r")
    wf = f["waveforms"]
    nev, nch, ns = wf.shape
    bad = [e for e in events if not 0 <= e < nev]
    if bad:
        sys.exit(f"events out of range 0..{nev - 1}: {bad}")
    channels = [int(c) for c in f["channels"][:]]
    dt = float(f.attrs["tsample_ns"])
    region = search_region("tail", args.tail, ns, dt)
    a, b = int(round(a_ns / dt)), int(round(b_ns / dt))
    n = a + b + 1
    cmap = read_map()

    live_rows = np.array([c not in dead for c in channels]); dead_rows = ~live_rows
    wander = []                               # (live median, dead mean) over the whole record, EMPTY events
    rms_acc = np.zeros(nch)
    slow_acc = np.zeros(nch)                 # rms of the 1 us moving average in the region: slow wander
    coh_num = 0.0; coh_den = 0.0             # coherence of the slow wander with the live-channel median
    cand = np.zeros((len(thr), nch), int)
    q_all, qc_all = [], []
    for e in events:
        w = wf[e].astype(np.float64)
        x = w - np.median(w, axis=1, keepdims=True)
        # noise and dark-pulse candidates in the search region
        xs = smooth(x[:, region], 125)
        med = np.median(xs, axis=0)
        for k in range(nch):
            rms_acc[k] += clipped_rms(x[k, region])
            slow_acc[k] += rms(xs[k])
            coh_num += np.sum((xs[k] - xs[k].mean()) * (med - med.mean())); coh_den += np.sum((xs[k] - xs[k].mean()) ** 2)
            for j, A in enumerate(thr):
                cand[j, k] += len(find_candidates(x[k], region, 1, A, 16)[1])
        # scintillation pulse on the sum, charge per channel in the window and in a control window
        S = smooth(x.sum(axis=0), args.smooth)
        imax = int(S.argmax())
        if imax < ns - region.stop + region.start - 250 and len(wander) < 2000:     # pulse early: the tail region is clean
            wander.append((np.median(x[live_rows][:, region], axis=0),
                           x[dead_rows][:, region].mean(axis=0) if dead_rows.any() else np.zeros(region.stop - region.start)))
        lo, hi = max(0, imax - a), min(ns, imax + b + 1)
        ctrl = slice(ns - n, ns) if hi <= ns - n else slice(0, n)
        q_all.append(dt * x[:, lo:hi].sum(axis=1))
        qc_all.append(dt * x[:, ctrl].sum(axis=1))
    N = len(events)
    q, qc = np.array(q_all), np.array(qc_all)
    sigma = rms_acc / N
    slow = slow_acc / N
    sq = qc.std(axis=0)
    lit = q > args.nsigma * sq
    frac = lit.mean(axis=0)
    qlit = np.array([q[lit[:, k], k].mean() if lit[:, k].any() else np.nan for k in range(nch)])

    # ---- table ---------------------------------------------------------------------------
    hdr = (f"{'ch':>3s} {'x':>6s} {'y':>6s} {'rms':>5s} {'slow':>5s} " + " ".join(f"{'n>' + f'{A:g}':>6s}" for A in thr)
           + f" {'sig_q':>7s} {'lit%':>6s} {'q_lit':>9s}")
    print(f"{path.stem}: {N} events, search region last {args.tail:g} us, window {a_ns:g}/{b_ns:g} ns, "
          f"lit if q > {args.nsigma:g} x own charge noise")
    print(hdr)
    rows = []
    for k, c in enumerate(channels):
        x_mm, y_mm = cmap.get(c, (np.nan, np.nan))
        rows.append([c, x_mm, y_mm, sigma[k]] + [int(cand[j, k]) for j in range(len(thr))]
                    + [sq[k], 100 * frac[k], qlit[k]])
        print(f"{c:3d} {x_mm:6.1f} {y_mm:6.1f} {sigma[k]:5.2f} {slow[k]:5.2f} " + " ".join(f"{cand[j, k]:6d}" for j in range(len(thr)))
              + f" {sq[k]:7.0f} {100 * frac[k]:6.1f} {qlit[k]:9.0f}")
    med_rms, med_n, med_f = np.median(sigma), np.median(cand[-1]), np.median(frac)
    flag = lambda cond: [channels[k] for k in range(nch) if cond[k]]
    print(f"\nmedians: rms {med_rms:.2f}, slow (1 us average) {np.median(slow):.2f}, candidates above {thr[-1]:g}: {med_n:.0f}, lit {100 * med_f:.1f}%")
    print(f"slow wander: white-noise expectation {med_rms / np.sqrt(125):.2f}; coherent fraction with the live-channel median {coh_num / coh_den:.2f}")
    print(f"noise more than 20% from the median:            {flag(np.abs(sigma - med_rms) > 0.2 * med_rms)}")
    print(f"candidates above {thr[-1]:g} more than 3x the median: {flag(cand[-1] > 3 * med_n)}")
    print(f"candidates above {thr[-1]:g} below a third of the median: {flag(cand[-1] < med_n / 3)}")
    print(f"lit fraction below a third of the median:       {flag(frac < med_f / 3)}")
    print(f"lit fraction above 3x the median:               {flag(frac > 3 * med_f)}")
    if args.csv:
        with open(args.csv, "w", newline="") as fo:
            wr = csv.writer(fo)
            wr.writerow(["channel", "x_mm", "y_mm", "rms"] + [f"n_above_{A:g}" for A in thr]
                        + ["sigma_q", "lit_fraction", "q_lit"])
            for r in rows:
                wr.writerow([r[0], r[1], r[2], f"{r[3]:.3f}"] + r[4:4 + len(thr)]
                            + [f"{r[-3]:.1f}", f"{r[-2] / 100:.4f}", f"{r[-1]:.1f}"])
        print(f"table written to {args.csv}")

    # ---- the wander in the EMPTY events --------------------------------------------------------
    if wander:
        L = np.array([w[0] for w in wander]); D = np.array([w[1] for w in wander])
        L = L - L.mean(axis=1, keepdims=True); D = D - D.mean(axis=1, keepdims=True)
        nl = int(live_rows.sum()); nd = int(dead_rows.sum())
        white_L = np.sqrt(np.pi / 2) * med_rms / np.sqrt(nl)          # noise of a median of nl white channels
        white_D = med_rms / np.sqrt(max(nd, 1))
        scales = [1, 12, 62, 125, 312, 625]                            # samples: 8 ns to 5 us
        nr = L.shape[1]
        rl = [np.sqrt(np.mean([rms(smooth(v, n)) ** 2 for v in L])) for n in scales]
        rd = [np.sqrt(np.mean([rms(smooth(v, n)) ** 2 for v in D])) for n in scales]
        # correlation between the live median and the dead mean at the 1 us scale
        Ls = np.array([smooth(v, 125) for v in L]); Ds = np.array([smooth(v, 125) for v in D])
        corr = np.sum(Ls * Ds) / np.sqrt(np.sum(Ls ** 2) * np.sum(Ds ** 2))
        # the wander of the sum of the live channels over a 5 us window: rms of the window mean minus
        # the end-of-record mean, times the number of live channels (coherent part only)
        h = nr // 2
        wsum = nl * np.sqrt(np.mean([(v[:h].mean() - v[h:].mean()) ** 2 for v in L]))
        wdead = nd * np.sqrt(np.mean([(v[:h].mean() - v[h:].mean()) ** 2 for v in D])) * nl / max(nd, 1)
        print(f"\nwander in the last {args.tail:g} us of {len(wander)} events whose pulse is early:")
        print(f"  rms of the live-channel median at 8 ns / 0.1 / 0.5 / 1 / 2.5 / 5 us smoothing: "
              + " / ".join(f"{v:.3f}" for v in rl) + f"  (white noise alone: {white_L:.3f} at 8 ns, falling as 1/sqrt(n))")
        print(f"  rms of the dead-channel mean, same scales:                 " + " / ".join(f"{v:.3f}" for v in rd)
              + f"  (white: {white_D:.3f})")
        print(f"  correlation live median vs dead mean at the 1 us scale: {corr:.2f}")
        print(f"  coherent wander of the sum of the {nl} live channels between two adjacent 5 us windows: {wsum:.1f} ADC per sample "
              f"({wsum * 626 * 8:.3g} ADC ns over a window); the same from the dead-channel mean scaled to {nl} channels: {wdead:.1f} "
              f"(white noise alone would give {nl * med_rms * np.sqrt(2 / h) * np.sqrt(np.pi / 2) / np.sqrt(nl):.1f})")
        figw, axs = plt.subplots(1, 2, figsize=(12, 4.2))
        us = [n * 8e-3 for n in scales]
        axs[0].plot(us, rl, "o-", color="C0", label=f"median of the {nl} live channels")
        axs[0].plot(us, rd, "s-", color="C3", label=f"mean of the {nd} dead channels")
        axs[0].plot(us, [white_L / np.sqrt(n) for n in scales], "--", color="C0", lw=0.8, label="white noise, live median")
        axs[0].plot(us, [white_D / np.sqrt(n) for n in scales], "--", color="C3", lw=0.8, label="white noise, dead mean")
        axs[0].set(xscale="log", yscale="log", xlabel="smoothing scale (us)", ylabel="rms (ADC)", title="common mode in the tail region")
        axs[0].grid(lw=0.3, which="both"); axs[0].legend(fontsize=8)
        freq = np.fft.rfftfreq(nr, d=8e-9) / 1e6
        PL = np.mean([np.abs(np.fft.rfft(v)) ** 2 for v in L], axis=0); PD = np.mean([np.abs(np.fft.rfft(v)) ** 2 for v in D], axis=0)
        axs[1].plot(freq[1:], PL[1:], color="C0", lw=0.8, label="live median")
        axs[1].plot(freq[1:], PD[1:], color="C3", lw=0.8, label="dead mean")
        axs[1].set(xscale="log", yscale="log", xlabel="frequency (MHz)", ylabel="power", title="power spectra")
        axs[1].grid(lw=0.3, which="both"); axs[1].legend(fontsize=8)
        figw.suptitle(f"{path.stem}: the wander in the last {args.tail:g} us of {len(wander)} events", x=0.01, ha="left")
        figw.tight_layout()
        if args.save:
            figw.savefig(plot_dir(args.run) / f"{args.save}_wander.png", dpi=130)
    else:
        print("\nno event with an early pulse: no wander monitor")

    # ---- figure --------------------------------------------------------------------------
    fig, axs = plt.subplots(2, 2, figsize=(12, 11))
    map8(axs[0, 0], dict(zip(channels, sigma)), "noise rms in the search region (ADC)")
    map8(axs[0, 1], dict(zip(channels, cand[-1])), f"dark-pulse candidates above {thr[-1]:g} ADC", "{:.0f}")
    map8(axs[1, 0], dict(zip(channels, 100 * frac)), f"lit fraction (%), q > {args.nsigma:g} sigma", "{:.0f}")
    map8(axs[1, 1], dict(zip(channels, qlit / 1e3)), "mean charge when lit (1e3 ADC ns)", "{:.0f}")
    fig.suptitle(f"{path.stem}: {N} events", x=0.01, ha="left")
    fig.tight_layout()
    if args.save:
        out = plot_dir(args.run) / f"{args.save}_chanstat.png"
        fig.savefig(out, dpi=130)
        print(f"saved {out}")
    else:
        plt.show()


if __name__ == "__main__":
    main()
