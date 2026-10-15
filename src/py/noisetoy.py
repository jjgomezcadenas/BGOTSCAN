#!/usr/bin/env python3
"""noisetoy.py: what white, independent noise would give in the figures of cmnoise.py and wavefit.py.

    python3 src/py/noisetoy.py [--sigma 3.8] [--samples 1500] [--channels 64] [--events 236] [--common 0] [--save NAME]

Figures: channel (one channel), channels (independent channels), common (with --common > 0). See docs/manual/noisetoy_doc.tex.
"""
import argparse

import numpy as np
import matplotlib.pyplot as plt

from bgo.io import PLOTS
from bgo.signal import smooth

DT = 8.0            # ns


def autocorr(v, nlag):
    v = v - v.mean()
    f = np.fft.rfft(v, 2 * len(v))
    ac = np.fft.irfft(f * np.conj(f))[:nlag]
    return ac / ac[0]


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--sigma", type=float, default=3.8, help="rms per sample (counts)")
    ap.add_argument("--samples", type=int, default=1500, help="samples per record")
    ap.add_argument("--channels", type=int, default=64, help="channels")
    ap.add_argument("--events", type=int, default=236, help="records")
    ap.add_argument("--common", type=float, default=0.0, help="rms of a slow component common to all channels (counts)")
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--save", metavar="NAME", help="save the figures as src/py/plots/NAME_noisetoy_*.png")
    args = ap.parse_args()
    rng = np.random.default_rng(args.seed)
    S, C, E, sig = args.samples, args.channels, args.events, args.sigma
    X = rng.normal(0, sig, (E, C, S))
    freq = np.fft.rfftfreq(S, d=DT * 1e-9) / 1e6
    figs = {}

    # ---- one channel ----------------------------------------------------------------------------------
    x = X[:, 0, :]
    P = np.mean(np.abs(np.fft.rfft(x - x.mean(axis=1, keepdims=True), axis=1)) ** 2, axis=0)
    nlag = int(10000 / DT)
    ac_raw = np.mean([autocorr(v, nlag) for v in x], axis=0)
    ac_sm = np.mean([autocorr(smooth(v, 62), nlag) for v in x], axis=0)
    lag = np.arange(nlag) * DT / 1000
    widths = np.array([1, 2, 4, 8, 16, 31, 62, 125, 250, 500, 625])
    rms_w = [np.sqrt(np.mean([smooth(v, n)[n:S - n].var() for v in x])) for n in widths]
    fig, axs = plt.subplots(1, 3, figsize=(15, 4.2))
    axs[0].plot(freq[1:], P[1:], color="C0", lw=0.6)
    axs[0].axhline(sig ** 2 * S, color="k", lw=0.8, ls="--", label=r"$\sigma^2 N$")
    axs[0].set(xscale="log", yscale="log", xlabel="frequency (MHz)", ylabel="power", title=f"power spectrum, {E} records")
    axs[1].plot(lag, ac_raw, color="C0", lw=1, label="raw samples")
    axs[1].plot(lag, ac_sm, color="C3", lw=1, label="62-sample running mean")
    axs[1].axhline(0, color="k", lw=0.5)
    axs[1].set(xlabel="lag (us)", ylabel="autocorrelation", xlim=(0, 10), title="autocorrelation")
    axs[2].plot(widths * DT / 1000, rms_w, "o", color="C0", label="simulated")
    axs[2].plot(widths * DT / 1000, sig / np.sqrt(widths), color="k", lw=0.8, ls="--", label=r"$\sigma/\sqrt{N}$")
    axs[2].set(xscale="log", yscale="log", xlabel="running-mean width (us)", ylabel="rms (counts)", title="rms of the running mean")
    for ax in axs:
        ax.grid(lw=0.3, which="both"); ax.legend(fontsize=8)
    fig.suptitle(f"white noise, one channel: sigma = {sig:g} counts per sample, {S} samples of {DT:g} ns per record", x=0.01, ha="left")
    fig.tight_layout(); figs["channel"] = fig

    # ---- independent channels -------------------------------------------------------------------------
    def matrix_and_sum(Y, title, name):
        Ys = smooth(Y, 125); Ys = Ys - Ys.mean(axis=2, keepdims=True)
        cov = np.einsum("eci,edi->cd", Ys, Ys) / (E * S)
        d = np.sqrt(np.diag(cov)); Cm = cov / np.outer(d, d)
        off = Cm[~np.eye(C, dtype=bool)]
        ns = np.array([1, 2, 4, 8, 16, 32, 52, 64])
        ns = ns[ns <= C]
        r_raw = [np.sqrt(np.mean(Y[:, :n, :].sum(axis=1).var(axis=1))) for n in ns]
        r_sm = [np.sqrt(np.mean(Ys[:, :n, :].sum(axis=1).var(axis=1))) for n in ns]
        fig, axs = plt.subplots(1, 2, figsize=(13, 5.2))
        im = axs[0].imshow(Cm, vmin=-0.2, vmax=0.6, cmap="RdBu_r")
        axs[0].set(xlabel="channel", ylabel="channel", title=f"correlation of the 125-sample running means; off-diagonal mean {off.mean():+.3f}, rms {off.std():.3f}")
        axs[0].title.set_fontsize(9)
        fig.colorbar(im, ax=axs[0], shrink=0.8)
        axs[1].plot(ns, r_raw, "o", color="C0", label="raw samples")
        axs[1].plot(ns, sig * np.sqrt(ns), color="C0", lw=0.8, ls="--", label=r"$\sigma\sqrt{n}$")
        axs[1].plot(ns, r_sm, "s", color="C3", label="125-sample running means")
        axs[1].plot(ns, r_sm[0] * np.sqrt(ns), color="C3", lw=0.8, ls="--", label=r"$\sqrt{n}$ from one channel")
        axs[1].plot(ns, r_sm[0] * ns, color="0.5", lw=0.8, ls=":", label=r"$n$ (fully common)")
        axs[1].set(xscale="log", yscale="log", xlabel="channels summed, n", ylabel="rms of the sum (counts)", title="rms of the sum of n channels")
        axs[1].grid(lw=0.3, which="both"); axs[1].legend(fontsize=8)
        fig.suptitle(title, x=0.01, ha="left")
        fig.tight_layout(); figs[name] = fig

    matrix_and_sum(X, f"{C} independent channels of white noise, sigma = {sig:g}, {E} records", "channels")
    if args.common > 0:
        slow = smooth(rng.normal(0, 1, (E, 1, S)), 125)
        slow *= args.common / slow.std()
        matrix_and_sum(X + slow, f"{C} channels of white noise plus a common slow component of {args.common:g} counts rms", "common")

    if args.save:
        PLOTS.mkdir(exist_ok=True)
        for name, fg in figs.items():
            fg.savefig(PLOTS / f"{args.save}_noisetoy_{name}.png", dpi=130)
        print(f"figures saved as {PLOTS}/{args.save}_noisetoy_{{{','.join(figs)}}}.png")
    else:
        plt.show()


if __name__ == "__main__":
    main()
