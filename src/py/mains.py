#!/usr/bin/env python3
"""mains.py: is the slow baseline drift the power line? The baseline of each event against the phase of
the event within the period of the mains (50 Hz: 20 ms), from the pass output.

    python3 src/py/mains.py RUN|FILE [--freq 50] [--save NAME]

If the drift comes from the mains, the baseline traces a sine wave (with harmonics) against the phase;
any other drift gives a cloud. See docs/manual/offline_doc.tex, "mains.py".
"""
import argparse

import h5py
import numpy as np
import matplotlib.pyplot as plt

from bgo.io import PLOTS, find_pass


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("run", help="run name (pass/RUN_pass.h5) or path to a pass file")
    ap.add_argument("--freq", type=float, default=50.0, help="mains frequency (Hz)")
    ap.add_argument("--save", metavar="NAME", help="save the figure as src/py/plots/NAME_mains.png")
    args = ap.parse_args()
    path = find_pass(args.run)
    with h5py.File(path, "r") as f:
        run = f.attrs["run"]
        run = run.decode() if isinstance(run, bytes) else str(run)
        ts = f["events/timestamp_ns"][:] * 1e-9            # s
        B = f["events/baseline"][:].astype(float)            # end-of-record level of the sum (ADC)
    period = 1.0 / args.freq
    phase = (ts % period) / period                            # 0 to 1 within one period
    nb = 40
    edges = np.linspace(0, 1, nb + 1)
    idx = np.clip((phase * nb).astype(int), 0, nb - 1)
    prof = np.array([np.median(B[idx == k]) if (idx == k).any() else np.nan for k in range(nb)])
    perr = np.array([1.25 * B[idx == k].std() / np.sqrt((idx == k).sum()) if (idx == k).sum() > 1 else np.nan for k in range(nb)])
    # amplitude of the first harmonics by a least-squares fit of sines and cosines
    rows = [np.ones_like(phase)]
    for h in (1, 2, 3):
        rows += [np.cos(2 * np.pi * h * phase), np.sin(2 * np.pi * h * phase)]
    Amat = np.array(rows).T
    coef, *_ = np.linalg.lstsq(Amat, B, rcond=None)
    amps = [np.hypot(coef[1 + 2 * k], coef[2 + 2 * k]) for k in range(3)]
    resid = B - Amat @ coef
    print(f"{path.name}: {len(B)} events, baseline rms {B.std():.2f} ADC; folded at {args.freq:g} Hz: amplitude of the "
          f"harmonics 1, 2, 3 = {amps[0]:.2f}, {amps[1]:.2f}, {amps[2]:.2f} ADC; rms after removing them {resid.std():.2f}; "
          f"duration {ts[-1] - ts[0]:.1f} s")
    fig, axs = plt.subplots(1, 2, figsize=(12, 4.2))
    axs[0].plot(phase, B, ".", ms=2, alpha=0.3, color="C0")
    c = 0.5 * (edges[1:] + edges[:-1])
    axs[0].errorbar(c, prof, yerr=perr, fmt="o-", color="C3", ms=4, lw=1.2, label="median per phase bin")
    axs[0].set(xlabel=f"phase within the {args.freq:g} Hz period", ylabel="baseline of the sum (ADC)",
               title="baseline against the mains phase")
    axs[0].legend(fontsize=8); axs[0].grid(lw=0.3)
    axs[1].plot(ts - ts[0], B, ".", ms=2, alpha=0.3, color="C0")
    axs[1].set(xlabel="time in the run (s)", ylabel="baseline of the sum (ADC)", title="baseline against time")
    axs[1].grid(lw=0.3)
    fig.suptitle(f"{run}: the end-of-record baseline of the sum", x=0.01, ha="left")
    fig.tight_layout()
    if args.save:
        PLOTS.mkdir(exist_ok=True)
        out = PLOTS / f"{args.save}_mains.png"
        fig.savefig(out, dpi=130)
        print(f"saved {out}")
    else:
        plt.show()


if __name__ == "__main__":
    main()
