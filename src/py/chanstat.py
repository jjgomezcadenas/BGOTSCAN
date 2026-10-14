#!/usr/bin/env python3
"""chanstat.py: the 64 channels over many events: noise, dark-pulse candidates, lit fraction, lit charge.

    python3 src/py/chanstat.py RUN [--events 0-499] [--tail 10] [--thresholds 12,16,20] [--smooth 64]
            [--window 1000,4000] [--nsigma 5] [--save NAME] [--csv FILE]

One table and one figure of four 8x8 maps. See docs/chanstat_doc.tex.
"""
import argparse
import csv
import sys
from pathlib import Path

import h5py
import numpy as np
import matplotlib.pyplot as plt

from bgo.io import DATA, PLOTS, find_run, parse_events, read_map
from bgo.signal import smooth, clipped_rms, search_region, find_candidates
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
    ap.add_argument("--save", metavar="NAME", help="save the figure as src/py/plots/NAME_chanstat.png")
    ap.add_argument("--csv", metavar="FILE", help="write the table to FILE")
    args = ap.parse_args()
    thr = [float(v) for v in args.thresholds.split(",")]
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

    rms_acc = np.zeros(nch)
    cand = np.zeros((len(thr), nch), int)
    q_all, qc_all = [], []
    for e in events:
        w = wf[e].astype(np.float64)
        x = w - np.median(w, axis=1, keepdims=True)
        # noise and dark-pulse candidates in the search region
        for k in range(nch):
            rms_acc[k] += clipped_rms(x[k, region])
            for j, A in enumerate(thr):
                cand[j, k] += len(find_candidates(x[k], region, 1, A, 16)[1])
        # scintillation pulse on the sum, charge per channel in the window and in a control window
        S = smooth(x.sum(axis=0), args.smooth)
        imax = int(S.argmax())
        lo, hi = max(0, imax - a), min(ns, imax + b + 1)
        ctrl = slice(ns - n, ns) if hi <= ns - n else slice(0, n)
        q_all.append(dt * x[:, lo:hi].sum(axis=1))
        qc_all.append(dt * x[:, ctrl].sum(axis=1))
    N = len(events)
    q, qc = np.array(q_all), np.array(qc_all)
    sigma = rms_acc / N
    sq = qc.std(axis=0)
    lit = q > args.nsigma * sq
    frac = lit.mean(axis=0)
    qlit = np.array([q[lit[:, k], k].mean() if lit[:, k].any() else np.nan for k in range(nch)])

    # ---- table ---------------------------------------------------------------------------
    hdr = (f"{'ch':>3s} {'x':>6s} {'y':>6s} {'rms':>5s} " + " ".join(f"{'n>' + f'{A:g}':>6s}" for A in thr)
           + f" {'sig_q':>7s} {'lit%':>6s} {'q_lit':>9s}")
    print(f"{path.stem}: {N} events, search region last {args.tail:g} us, window {a_ns:g}/{b_ns:g} ns, "
          f"lit if q > {args.nsigma:g} x own charge noise")
    print(hdr)
    rows = []
    for k, c in enumerate(channels):
        x_mm, y_mm = cmap.get(c, (np.nan, np.nan))
        rows.append([c, x_mm, y_mm, sigma[k]] + [int(cand[j, k]) for j in range(len(thr))]
                    + [sq[k], 100 * frac[k], qlit[k]])
        print(f"{c:3d} {x_mm:6.1f} {y_mm:6.1f} {sigma[k]:5.2f} " + " ".join(f"{cand[j, k]:6d}" for j in range(len(thr)))
              + f" {sq[k]:7.0f} {100 * frac[k]:6.1f} {qlit[k]:9.0f}")
    med_rms, med_n, med_f = np.median(sigma), np.median(cand[-1]), np.median(frac)
    flag = lambda cond: [channels[k] for k in range(nch) if cond[k]]
    print(f"\nmedians: rms {med_rms:.2f}, candidates above {thr[-1]:g}: {med_n:.0f}, lit {100 * med_f:.1f}%")
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

    # ---- figure --------------------------------------------------------------------------
    fig, axs = plt.subplots(2, 2, figsize=(12, 11))
    map8(axs[0, 0], dict(zip(channels, sigma)), "noise rms in the search region (ADC)")
    map8(axs[0, 1], dict(zip(channels, cand[-1])), f"dark-pulse candidates above {thr[-1]:g} ADC", "{:.0f}")
    map8(axs[1, 0], dict(zip(channels, 100 * frac)), f"lit fraction (%), q > {args.nsigma:g} sigma", "{:.0f}")
    map8(axs[1, 1], dict(zip(channels, qlit / 1e3)), "mean charge when lit (1e3 ADC ns)", "{:.0f}")
    fig.suptitle(f"{path.stem}: {N} events", x=0.01, ha="left")
    fig.tight_layout()
    if args.save:
        PLOTS.mkdir(exist_ok=True)
        out = PLOTS / f"{args.save}_chanstat.png"
        fig.savefig(out, dpi=130)
        print(f"saved {out}")
    else:
        plt.show()


if __name__ == "__main__":
    main()
