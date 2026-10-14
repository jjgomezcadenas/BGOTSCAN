#!/usr/bin/env python3
"""wfview.py: look at the summed SiPM waveforms of a BGO temperature-scan run.

    python3 src/py/wfview.py RUN [--events 0-3 | 0,5,7] [--baseline-window 150]
                                 [--smooth N] [--channels] [--average] [--save NAME]

Reads the raw WebDump HDF5 file with h5py. See docs/manual/wfview_doc.tex.
"""
import argparse
import sys
from pathlib import Path

import h5py
import numpy as np
import matplotlib.pyplot as plt

from bgo.io import DATA, PLOTS, find_run, parse_events
from bgo.signal import smooth

CHUNK = 500


def sums(wf, events, meas, nbase):
    """Per event: summed waveform with WebDump baselines and with the pre-trigger baseline."""
    x = wf[events].astype(np.float32)                     # (nev, nch, ns)
    s_meas = (x - meas[None, :, None]).sum(axis=1)
    base = np.median(x[:, :, :nbase], axis=2)             # per event and channel
    s_pre = (x - base[..., None]).sum(axis=1)
    return x, base, s_meas, s_pre


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("run", help="run name (128K_37V_th19) or path to a .h5 file")
    ap.add_argument("--data", default=str(DATA), help="directory of the run files")
    ap.add_argument("--events", default="0-3", help="range a-b or list a,b,c (0-based)")
    ap.add_argument("--baseline-window", type=int, default=150,
                    help="pre-trigger samples used for the baseline (from the start of the record)")
    ap.add_argument("--smooth", type=int, default=0, help="moving-average width in samples (0 = off)")
    ap.add_argument("--channels", action="store_true", help="show the 64 channels of the first event")
    ap.add_argument("--average", action="store_true", help="show the average summed pulse of the run")
    ap.add_argument("--save", metavar="NAME", help="save PNGs under src/py/plots/ with this prefix")
    args = ap.parse_args()

    path = find_run(args.run, Path(args.data))
    events = parse_events(args.events)
    f = h5py.File(path, "r")
    wf = f["waveforms"]
    A = f.attrs
    nev, nch, ns = wf.shape
    dt, pre = float(A["tsample_ns"]), int(A["pretrigger"])
    meas = np.asarray(A["baselines"], float)
    channels = f["channels"][:]
    t = (np.arange(ns) - pre) * dt
    nbase = min(args.baseline_window, pre)
    bad = [e for e in events if not 0 <= e < nev]
    if bad:
        sys.exit(f"events out of range 0..{nev - 1}: {bad}")
    print(f"{path.name}: {nev} events, {nch} channels x {ns} samples, {dt:g} ns/sample, "
          f"pre-trigger {pre} samples; baseline window {nbase} samples"
          + (f"; smoothing {args.smooth} samples" if args.smooth > 1 else ""))

    # ---- summed waveforms of the selected events --------------------------------------
    x, base, s_meas, s_pre = sums(wf, events, meas, nbase)
    s_sm = smooth(s_pre, args.smooth) if args.smooth > 1 else None
    k = len(events)
    ncol = min(k, 2)
    nrow = (k + ncol - 1) // ncol
    fig, axs = plt.subplots(nrow, ncol, figsize=(7 * ncol, 3.2 * nrow), sharex=True, squeeze=False)
    for j, (e, ax) in enumerate(zip(events, axs.flat)):
        imax = int(s_pre[j].argmax())
        integ = s_pre[j].sum() * dt
        print(f"event {e:5d}: baseline sum {base[j].sum():9.1f}  max {s_pre[j][imax]:7.1f} at "
              f"{t[imax]:8.0f} ns  integral {integ:11.4g} ADC ns")
        ax.plot(t, s_meas[j], color="0.7", lw=0.6, label="WebDump baselines")
        ax.plot(t, s_pre[j], color="C0", lw=0.8, label="pre-trigger baseline")
        if s_sm is not None:
            ax.plot(t, s_sm[j], color="C1", lw=1.4, label=f"moving average {args.smooth}")
        ax.axvspan(t[0], t[nbase - 1], color="0.5", alpha=0.12, lw=0)
        ax.axvline(0, color="0.4", ls="--", lw=0.8)
        ax.set_title(f"event {e}: max {s_pre[j][imax]:.0f} at {t[imax]:.0f} ns, "
                     f"integral {integ:.3g} ADC ns", fontsize=9, loc="left")
        ax.grid(lw=0.3)
    for ax in axs.flat[k:]:
        ax.set_visible(False)
    for ax in axs[-1]:
        ax.set_xlabel("ns from trigger")
    for ax in axs[:, 0]:
        ax.set_ylabel("sum of channels (ADC)")
    axs.flat[0].legend(loc="upper right", fontsize=8)
    fig.suptitle(f"{path.stem}: summed waveforms", x=0.01, ha="left")
    fig.tight_layout()
    figs = [("sum", fig)]

    # ---- channels of the first selected event -----------------------------------------
    if args.channels:
        e = events[0]
        y = x[0] - base[0][:, None]
        full = np.full((64, ns), np.nan)
        full[channels] = y
        fig, axs = plt.subplots(8, 8, figsize=(14, 9), sharex=True, sharey=True)
        for ch, ax in enumerate(axs.flat):
            ax.plot(t, full[ch], color="C0", lw=0.5)
            ax.axvline(0, color="0.4", ls="--", lw=0.5)
            ax.text(0.03, 0.88, str(ch), transform=ax.transAxes, fontsize=7)
            ax.tick_params(labelsize=6)
            ax.grid(lw=0.3)
        fig.suptitle(f"{path.stem}: event {e}, channels (ADC, pre-trigger baseline) vs ns from trigger",
                     x=0.01, ha="left")
        fig.tight_layout()
        figs.append((f"ev{e}_channels", fig))

    # ---- average summed pulse of the run ----------------------------------------------
    if args.average:
        acc_meas = np.zeros(ns)
        acc_pre = np.zeros(ns)
        for a in range(0, nev, CHUNK):
            b = min(nev, a + CHUNK)
            _, _, sm, sp = sums(wf, slice(a, b), meas, nbase)
            acc_meas += sm.sum(axis=0, dtype=np.float64)
            acc_pre += sp.sum(axis=0, dtype=np.float64)
        acc_meas /= nev
        acc_pre /= nev
        imax = int(acc_pre.argmax())
        print(f"average pulse: max {acc_pre[imax]:.1f} at {t[imax]:.0f} ns, "
              f"integral {acc_pre.sum() * dt:.4g} ADC ns, pre-trigger level (WebDump baselines) "
              f"{acc_meas[:nbase].mean():+.1f}, end-of-record level {acc_pre[-200:].mean():+.1f}")
        fig, ax = plt.subplots(figsize=(11, 3.6))
        ax.plot(t, acc_meas, color="0.6", lw=0.8, label="WebDump baselines")
        ax.plot(t, acc_pre, color="C0", lw=1.2, label="pre-trigger baseline")
        ax.axvspan(t[0], t[nbase - 1], color="0.5", alpha=0.12, lw=0)
        ax.axvline(0, color="0.4", ls="--", lw=0.8)
        ax.axhline(0, color="0.4", lw=0.5)
        ax.set(title=f"{path.stem}: average summed pulse ({nev} events)",
               xlabel="ns from trigger", ylabel="sum of channels (ADC)")
        ax.grid(lw=0.3)
        ax.legend(loc="upper right", fontsize=8)
        fig.tight_layout()
        figs.append(("average", fig))

    if args.save:
        PLOTS.mkdir(exist_ok=True)
        for tag, fig in figs:
            out = PLOTS / f"{args.save}_{tag}.png"
            fig.savefig(out, dpi=130)
            print(f"saved {out}")
    else:
        plt.show()


if __name__ == "__main__":
    main()
