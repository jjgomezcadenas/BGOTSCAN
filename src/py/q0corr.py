#!/usr/bin/env python3
"""q0corr.py: the seed-channel response correction of offline.py, step by step, on one run.

    python3 src/py/q0corr.py RUN [--channels 24,34,48,3] [--min-events 40] [--save NAME]

Reads pass/RUN_offline.h5 (written by offline.py). Three figures. See docs/manual/q0corr_doc.tex.
"""
import argparse

import h5py
import numpy as np
import matplotlib.pyplot as plt

from bgo.io import PLOTS, PASS
from bgo.fit import fit_photopeak, find_photopeak, peak_bg
from bgo.plot import map8


def seed_positions(vals, seed, sel, channels, nmin):
    """Photopeak position per seed channel: highest maximum of the smoothed 60-bin histogram of vals
    over the selected events whose seed is the channel (channels with at least nmin of them)."""
    top = np.percentile(vals[sel], 99.5)
    out = {}
    for c in channels:
        m = sel & (seed == c)
        if m.sum() >= nmin:
            h, edges = np.histogram(vals[m], bins=60, range=(0, top))
            cc = 0.5 * (edges[1:] + edges[:-1])
            out[c] = find_photopeak(h, cc, 0.25, 3)[0]
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("run", help="run name, e.g. 195K_41V_th20")
    ap.add_argument("--channels", default="24,34,48,3", help="seed channels drawn one by one")
    ap.add_argument("--min-events", type=int, default=40, help="events a channel needs as seed")
    ap.add_argument("--save", metavar="NAME", help="save the figures as src/py/plots/NAME_q0corr_*.png")
    args = ap.parse_args()
    shown = [int(v) for v in args.channels.split(",")]

    path = PASS / f"{args.run}_offline.h5"
    with h5py.File(path, "r") as f:
        q_eq = f["Q_0_eq"][:].astype(float)
        q_pos = f["Q_0_pos"][:].astype(float)
        seed = f["seed"][:].astype(int)
        cls = f["cls"][:]
        pile = f["pileup"][:].astype(bool)
    sel = (cls == 1) & ~pile
    channels = sorted(set(seed[sel]))
    print(f"{path.name}: {int(sel.sum())} REAL events without pile-up")

    # 1. the photopeak position per seed channel, before and after
    mu = seed_positions(q_eq, seed, sel, channels, args.min_events)
    mu2 = seed_positions(q_pos, seed, sel, channels, args.min_events)
    med, med2 = np.median(list(mu.values())), np.median(list(mu2.values()))
    print(f"  {len(mu)} seed channels with {args.min_events}+ events; median position {med:.4g} before, {med2:.4g} after")
    print("  seed channel: events, mu_k / median before, mu'_k / median after")
    for c in sorted(mu, key=mu.get):
        print(f"    ch {c:2d}: {int((sel & (seed == c)).sum()):4d}  {mu[c] / med:5.2f}  {mu2.get(c, np.nan) / med2:5.2f}")
    r1 = np.array(list(mu.values())) / med
    r2 = np.array([mu2[c] for c in mu if c in mu2]) / med2
    print(f"  spread (rms / median): {r1.std():.3f} before, {r2.std():.3f} after; range {r1.min():.2f} to {r1.max():.2f} before, "
          f"{r2.min():.2f} to {r2.max():.2f} after")

    figs = {}
    # 2. the maps
    fig, axs = plt.subplots(1, 2, figsize=(11, 5))
    map8(axs[0], {c: v / med for c, v in mu.items()}, "before: photopeak position by seed channel, Q_0_eq (median = 1)")
    map8(axs[1], {c: v / med2 for c, v in mu2.items()}, "after: the same on Q_0_pos")
    fig.suptitle(f"{args.run}: the factor f_k = mu_k / median, and what is left after dividing by it", x=0.01, ha="left")
    fig.tight_layout(); figs["seedmap"] = fig

    # 3. a few seed-channel populations, before and after
    fig, axs = plt.subplots(1, 2, figsize=(13, 4.5), sharey=True)
    top = 1.7 * med
    for ax, vals, pos, lab in ((axs[0], q_eq, mu, "Q_0_eq"), (axs[1], q_pos, mu2, "Q_0_pos")):
        for i, c in enumerate(shown):
            m = sel & (seed == c)
            if m.sum() == 0:
                continue
            ax.hist(vals[m], bins=50, range=(0, top), histtype="step", color=f"C{i}", lw=1.3,
                    label=f"seed {c}: {int(m.sum())} events" + (f", peak {pos[c] / (med if lab == 'Q_0_eq' else med2):.2f}" if c in pos else ""))
            if c in pos:
                ax.axvline(pos[c], color=f"C{i}", ls="--", lw=1)
        ax.axvline(med if lab == "Q_0_eq" else med2, color="k", lw=1, label="median of the mu_k")
        ax.set(xlabel=f"{lab} (ADC ns)", ylabel="events", title=f"{lab}: events by seed channel"); ax.title.set_fontsize(10)
        ax.grid(lw=0.3); ax.legend(fontsize=7)
    fig.suptitle(f"{args.run}: the same events before (left) and after (right) the division by f_k", x=0.01, ha="left")
    fig.tight_layout(); figs["seeds"] = fig

    # 4. the whole sample
    fig, axs = plt.subplots(1, 2, figsize=(13, 4.5))
    for ax, vals, lab in ((axs[0], q_eq, "Q_0_eq"), (axs[1], q_pos, "Q_0_pos")):
        fi = fit_photopeak(vals[sel])
        ax.stairs(fi["h"], fi["edges"], color="C0", fill=True, alpha=0.6)
        if fi["ok"]:
            x = np.linspace(fi["lo"], fi["hi"], 200)
            ax.plot(x, peak_bg(x, *fi["p"]), color="C3", lw=1.3, label=f"peak {fi['mu']:.4g} +- {fi['err']:.2g}, FWHM {100 * fi['fwhm_rel']:.1f}%")
            ax.legend(fontsize=8)
            print(f"  {lab}: peak {fi['mu']:.4g} +- {fi['err']:.2g}, FWHM {100 * fi['fwhm_rel']:.1f}%")
        ax.set(xlabel=f"{lab} (ADC ns)", ylabel="events", title=lab); ax.title.set_fontsize(10); ax.grid(lw=0.3)
    fig.suptitle(f"{args.run}: the whole sample before and after, {int(sel.sum())} REAL events without pile-up", x=0.01, ha="left")
    fig.tight_layout(); figs["spectra"] = fig

    if args.save:
        PLOTS.mkdir(exist_ok=True)
        for tag, fg in figs.items():
            fg.savefig(PLOTS / f"{args.save}_q0corr_{tag}.png", dpi=130)
        print(f"  figures saved as {PLOTS}/{args.save}_q0corr_{{{','.join(figs)}}}.png")
    else:
        plt.show()


if __name__ == "__main__":
    main()
