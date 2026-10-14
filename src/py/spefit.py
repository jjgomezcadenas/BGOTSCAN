#!/usr/bin/env python3
"""spefit.py: single-photoelectron spectra, fits and equalization factors from the output of the pass.

    python3 src/py/spefit.py RUN|FILE [--channel 48] [--min-pulses 30] [--dead LIST] [--save NAME] [--csv FILE]

Reads pass/RUN_pass.h5 (JScan pass) only. See docs/spe_doc.tex, "spefit.py".
"""
import argparse
import csv

import h5py
import numpy as np
import matplotlib.pyplot as plt

from bgo.io import PLOTS, find_pass
from bgo.fit import ABIN, AMAX, gauss, fit_peak, crosstalk_ratio
from bgo.plot import map8

DEAD = "19,20,21,22,23,25,27,33,35,41,46,47"      # chanstat_doc.tex, 2026-10-14


def draw_fit(ax, f, label=True):
    """Histogram of A with the fitted Gaussian over the fit range."""
    if f["h"] is None:
        return
    ax.stairs(f["h"], f["edges"], color="C0", fill=True, alpha=0.6)
    if np.isfinite(f["A1"]) and "lo" in f:
        a = np.linspace(f["lo"], f["hi"], 100)
        c = 0.5 * (f["edges"][1:] + f["edges"][:-1])
        k = int(np.clip(np.searchsorted(c, f["A1"]), 0, len(f["h"]) - 1))
        ax.plot(a, gauss(a, f["h"][k], f["A1"], f["w1"]), color="C3", lw=1.2,
                label=(f"A1 {f['A1']:.1f}, w {f['w1']:.1f}" + ("" if f["ok"] else " (rejected)")) if label else None)


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("run", help="run name (pass/RUN_pass.h5) or path to a pass file")
    ap.add_argument("--channel", type=int, default=48, help="channel shown alone")
    ap.add_argument("--min-pulses", type=int, default=30, help="channels with fewer candidates are not fitted")
    ap.add_argument("--dead", default=DEAD, help="dead channels, comma list (excluded from the reference)")
    ap.add_argument("--save", metavar="NAME", help="save figures as src/py/plots/NAME_spefit*.png")
    ap.add_argument("--csv", metavar="FILE", help="write the per-channel table")
    args = ap.parse_args()
    dead = {int(v) for v in args.dead.split(",") if v.strip()}

    path = find_pass(args.run)
    with h5py.File(path, "r") as f:
        A = f.attrs
        run, nsig, tail = str(A["run"], "utf8") if isinstance(A["run"], bytes) else str(A["run"]), float(A["nsigma"]), float(A["tail_us"])
        nev = int(A["events_last"] - A["events_first"] + 1)
        template = A["template_file"]
        template = template.decode() if isinstance(template, bytes) else str(template)
        channels = [int(c) for c in A["channels"]]
        spe_ch = f["spe/channel"][:]
        spe_A = f["spe/A"][:]
        sigma = f["spe_sigma"][:]                      # (nev, nch) as seen by h5py
    sy = {c: float(sigma[:, k].mean()) for k, c in enumerate(channels)}
    thr = {c: nsig * sy[c] for c in channels}
    amps = {c: spe_A[spe_ch == c] for c in channels}
    live = [c for c in channels if len(amps[c]) > 0 and c not in dead]
    T_search = nev * tail * 1e-6                       # s per channel

    print(f"{path.name}: run {run}, {nev} events, template {template.split('/')[-1]}, threshold {nsig:g} sigma "
          f"(about {np.mean([thr[c] for c in live]):.1f} ADC), search {tail:g} us per event; "
          f"{len(spe_A)} candidates, {len(live)} live channels, dead excluded: {sorted(dead)}")
    fits = {c: fit_peak(amps[c], thr[c], sy[c]) if len(amps[c]) >= args.min_pulses else None for c in channels}
    allA = np.concatenate([amps[c] for c in live])
    fit_all = fit_peak(allA, float(np.mean([thr[c] for c in live])), float(np.mean([sy[c] for c in live])))
    ok = {c: fits[c] for c in live if fits[c] is not None and fits[c]["ok"]}
    med = float(np.median([fi["A1"] for fi in ok.values()])) if ok else np.nan
    g = {c: fi["A1"] / med for c, fi in ok.items()}

    c0 = args.channel
    for who, fi, vals in ((f"channel {c0}", fits.get(c0), amps.get(c0, [])), ("all live channels", fit_all, allA)):
        if fi is not None and fi["ok"]:
            print(f"  {who:18s}: {fi['n']:6d} pulses, A1 {fi['A1']:.2f} +- {fi['err']:.2f} ADC, width {fi['w1']:.2f}, "
                  f"crosstalk ratio {crosstalk_ratio(vals, fi['A1']):.3f}")
        else:
            print(f"  {who:18s}: {len(vals):6d} pulses, no fit" + (f" ({fi['why']})" if fi is not None else ""))
    print(f"  per channel: n, A1, error, width, g = A1 / median A1 ({med:.2f}), dark rate (kHz)")
    rows = []
    for c in channels:
        n = len(amps[c]); rate = n / T_search / 1e3
        fi = fits[c]
        tag = " dead" if c in dead else ""
        if fi is None:
            print(f"  {c:3d}: {n:6d}  not fitted (fewer than {args.min_pulses} pulses)  rate {rate:6.2f}{tag}")
            rows.append([c, n, "", "", "", "", f"{rate:.3f}"])
        elif fi["ok"] and c in g:
            print(f"  {c:3d}: {n:6d}  {fi['A1']:6.2f} +- {fi['err']:4.2f}  w {fi['w1']:5.2f}  g {g[c]:5.3f}  rate {rate:6.2f}")
            rows.append([c, n, f"{fi['A1']:.3f}", f"{fi['err']:.3f}", f"{fi['w1']:.3f}", f"{g[c]:.4f}", f"{rate:.3f}"])
        else:
            why = fi["why"] if not fi["ok"] else "dead"
            print(f"  {c:3d}: {n:6d}  no fit ({why})" + (f"; mean {fi['A1']:.2f}, width {fi['w1']:.2f}" if np.isfinite(fi["A1"]) else "")
                  + f"  rate {rate:6.2f}{tag}")
            rows.append([c, n, "", "", "", "", f"{rate:.3f}"])
    off = sorted(c for c, v in g.items() if abs(v - 1) > 0.15)
    print(f"  median A1 {med:.2f} ADC over {len(ok)} fitted channels; g more than 15% from 1: {off}")
    if args.csv:
        with open(args.csv, "w", newline="") as fo:
            fo.write(f"# SPE amplitudes of {run} from {path.name}: {nev} events, threshold {nsig:g} sigma; "
                     f"g = A1 / median A1 ({med:.3f}) over the fitted live channels\n")
            wr = csv.writer(fo)
            wr.writerow(["channel", "n", "A1", "err", "width", "g", "rate_khz"])
            wr.writerows(rows)
        print(f"  table written to {args.csv}")

    # ---- figures ---------------------------------------------------------------------
    fig, axs = plt.subplots(2, 2, figsize=(12, 10))
    for ax, (who, fi, vals) in zip(axs[0], ((f"channel {c0}", fits.get(c0), amps.get(c0, [])),
                                            ("all live channels", fit_all, allA))):
        if fi is not None:
            draw_fit(ax, fi)
            if np.isfinite(fi["A1"]):
                ax.legend(fontsize=8)
        ax.set(xlabel="A (ADC)", ylabel="pulses", title=f"{who}, {len(vals)} pulses")
        ax.title.set_fontsize(9)
        ax.grid(lw=0.3)
    if ok:
        map8(axs[1, 0], {c: fi["A1"] for c, fi in ok.items()}, "A1 (ADC)")
        map8(axs[1, 1], g, "g = A1 / median")
    fig.suptitle(f"{run}: SPE from the pass, {nev} events", x=0.01, ha="left")
    fig.tight_layout()
    fig2, axs2 = plt.subplots(8, 8, figsize=(16, 14), sharex=True)
    for c, ax in enumerate(axs2.flat):
        fi = fits.get(c)
        if fi is not None:
            draw_fit(ax, fi, label=False)
        ax.text(0.03, 0.82, f"{c}" + (" dead" if c in dead else "") + (f"\n{fi['A1']:.1f}" if fi is not None and fi["ok"] else ""),
                transform=ax.transAxes, fontsize=7)
        ax.tick_params(labelsize=6)
        ax.set_xlim(0, AMAX)
        ax.grid(lw=0.3)
    fig2.suptitle(f"{run}: per-channel SPE spectra and fits (A in ADC)", x=0.01, ha="left")
    fig2.tight_layout()
    if args.save:
        PLOTS.mkdir(exist_ok=True)
        for tag, fg in (("spefit", fig), ("spefit_channels", fig2)):
            out = PLOTS / f"{args.save}_{tag}.png"
            fg.savefig(out, dpi=130)
            print(f"  saved {out}")
    else:
        plt.show()


if __name__ == "__main__":
    main()
