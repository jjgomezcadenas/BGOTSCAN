#!/usr/bin/env python3
"""cmnoise.py: the noise common to the biased channels, in the tail of the events.

    python3 src/py/cmnoise.py RUN [--events 0-299] [--tail 12] [--smooth 125] [--dead LIST]
            [--offline FILE] [--kick 16] [--nsigma 5] [--traces 4] [--save NAME]

Five figures: traces, rms, corr, spectra, kick. See docs/manual/cmnoise_doc.tex.
"""
import argparse
import json
import sys
from pathlib import Path

import h5py
import numpy as np
import matplotlib.pyplot as plt

from bgo.io import DATA, PASS, PLOTS, find_run, parse_events, read_map
from bgo.signal import smooth, rms, clipped_rms

DEAD = "19,20,21,22,23,25,27,33,35,41,46,47"      # chanstat_doc.tex, 2026-10-14


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("run", help="run name (195K_41V_th20) or path to a .h5 file")
    ap.add_argument("--data", default=str(DATA), help="directory of the run files")
    ap.add_argument("--events", default="0-299", help="range a-b or list a,b,c (0-based)")
    ap.add_argument("--tail", type=float, default=12.0, help="the tail: last T us of the record")
    ap.add_argument("--smooth", type=int, default=125, help="running-mean width of the slow component (samples)")
    ap.add_argument("--dead", default=DEAD, help="dead channels, comma list")
    ap.add_argument("--offline", help="per-event file of offline.py (default pass/RUN_offline.h5; none: every event)")
    ap.add_argument("--kick", type=int, default=16, help="half-width of the kick template (samples)")
    ap.add_argument("--nsigma", type=float, default=5.0, help="kick threshold in units of the filtered noise")
    ap.add_argument("--traces", type=int, default=4, help="events drawn in the trace figure")
    ap.add_argument("--save", metavar="NAME", help="save the figures as src/py/plots/NAME_cmnoise_*.png")
    args = ap.parse_args()
    dead = {int(v) for v in args.dead.split(",") if v.strip()}

    path = find_run(args.run, Path(args.data))
    run = path.stem
    events = parse_events(args.events)
    f = h5py.File(path, "r")
    wf = f["waveforms"]
    nev, nch, ns = wf.shape
    bad = [e for e in events if not 0 <= e < nev]
    if bad:
        sys.exit(f"events out of range 0..{nev - 1}: {bad}")
    channels = [int(c) for c in f["channels"][:]]
    dt = float(f.attrs["tsample_ns"])
    pre = int(f.attrs["pretrigger"])
    base = np.asarray(f.attrs["baselines"], float)[:, None]
    ntail = int(round(args.tail * 1000 / dt))
    tail = slice(ns - ntail, ns)
    t_us = np.arange(ns) * dt / 1000
    live = np.array([c not in dead for c in channels]); nl = int(live.sum()); nd = int((~live).sum())
    cmap = read_map()

    # ---- clean events: REAL without pile-up, maximum more than 4 us before the tail ------------------
    off = Path(args.offline) if args.offline else PASS / f"{run}_offline.h5"
    trig_rate = None
    if off.is_file():
        with h5py.File(off, "r") as fo:
            cls, pile, tmax = fo["cls"][:], fo["pileup"][:], fo["t_max"][:]
        clean = (cls == 1) & ~pile & (tmax < (ns - ntail) * dt - 4000)
        kept = [e for e in events if clean[e]]
        js = off.with_suffix(".json")
        if js.is_file():
            trig_rate = json.load(open(js)).get("trigger_rate_hz")
        print(f"{run}: {len(kept)} of {len(events)} events kept (REAL without pile-up, maximum > 4 us before the tail), from {off.name}")
    else:
        kept = list(events)
        print(f"{run}: {len(kept)} events, no offline file ({off}): every event used")
    if len(kept) < 10:
        sys.exit("fewer than 10 clean events")

    # ---- accumulate -----------------------------------------------------------------------------------
    N = args.smooth
    L_raw, D_raw, M_raw = [], [], []            # live mean, dead mean, live median over the tail (per event)
    L_off, D_off, off_ch = [], [], []            # offsets of the record from the run baseline over the tail: live mean, dead mean, per channel
    cov = np.zeros((nch, nch)); nsm_acc = 0
    PL = PD = PC = None
    kick_L = np.zeros(ns); kick_D = np.zeros(ns)
    sigma_ch = []
    for e in kept:
        x = wf[e].astype(np.float64) - base
        xt = x[:, tail]
        L = xt[live].mean(axis=0); D = xt[~live].mean(axis=0); M = np.median(xt[live], axis=0)
        L_off.append(L.mean()); D_off.append(D.mean()); off_ch.append(xt.mean(axis=1))
        L_raw.append(L); D_raw.append(D); M_raw.append(M)
        xs = smooth(xt, N); xs = xs - xs.mean(axis=1, keepdims=True)
        cov += xs @ xs.T; nsm_acc += xs.shape[1]
        sigma_ch.append(np.array([clipped_rms(v) for v in xt]))
        pl = np.abs(np.fft.rfft(L - L.mean())) ** 2; pd_ = np.abs(np.fft.rfft(D - D.mean())) ** 2
        pc = np.abs(np.fft.rfft((L - D) - (L - D).mean())) ** 2
        PL = pl if PL is None else PL + pl; PD = pd_ if PD is None else PD + pd_; PC = pc if PC is None else PC + pc
        kick_L += np.median(x[live], axis=0); kick_D += x[~live].mean(axis=0)
    n = len(kept)
    L_raw, D_raw, M_raw = np.array(L_raw), np.array(D_raw), np.array(M_raw)
    L_off, D_off, off_ch = np.array(L_off), np.array(D_off), np.array(off_ch)
    sigma = np.median(np.array(sigma_ch), axis=0)             # per channel, in the tail
    sig_white = float(np.median(sigma[live]))
    PL /= n; PD /= n; PC /= n
    kick_L /= n; kick_D /= n
    freq = np.fft.rfftfreq(ntail, d=dt * 1e-9) / 1e6

    off_diag = ~np.eye(nch, dtype=bool)

    # ---- how much -------------------------------------------------------------------------------------
    Lc = L_raw - L_raw.mean(axis=1, keepdims=True); Dc = D_raw - D_raw.mean(axis=1, keepdims=True)
    Ls = np.array([smooth(v, N) for v in Lc]); Ds = np.array([smooth(v, N) for v in Dc])
    Ls -= Ls.mean(axis=1, keepdims=True); Ds -= Ds.mean(axis=1, keepdims=True)
    r = dict(L=np.array([rms(v) for v in Lc]), D=np.array([rms(v) for v in Dc]),
             Ls=np.array([rms(v) for v in Ls]), Ds=np.array([rms(v) for v in Ds]))
    white = dict(L=sig_white / np.sqrt(nl), D=sig_white / np.sqrt(nd),
                 Ls=sig_white / np.sqrt(nl * N), Ds=sig_white / np.sqrt(nd * N))
    corr_LD = float(np.sum(Ls * Ds) / np.sqrt(np.sum(Ls ** 2) * np.sum(Ds ** 2)))
    print(f"  tail: last {args.tail:g} us ({ntail} samples); {nl} live and {nd} dead channels; "
          f"white noise per channel {sig_white:.2f} counts (median over the live channels)")
    print(f"  rms over the tail, median over events (white-noise value in brackets):")
    print(f"    live mean  raw {np.median(r['L']):.3f} [{white['L']:.3f}]   smoothed {N} samples {np.median(r['Ls']):.3f} [{white['Ls']:.3f}]")
    print(f"    dead mean  raw {np.median(r['D']):.3f} [{white['D']:.3f}]   smoothed {N} samples {np.median(r['Ds']):.3f} [{white['Ds']:.3f}]")
    print(f"  correlation of the smoothed live and dead means: {corr_LD:.2f}")
    print(f"  offset of the record from the run baseline (mean over the tail), scatter between events: live mean {L_off.std():.2f} counts "
          f"(x {nl} = {nl * L_off.std():.0f} on the sum), dead mean {D_off.std():.2f}; correlation of the two {np.corrcoef(L_off, D_off)[0, 1]:.2f}")
    oc = off_ch - off_ch.mean(axis=0)
    Co = np.corrcoef(oc.T)
    print(f"  per-channel offsets: scatter between events {np.median(oc[:, live].std(axis=0)):.2f} counts (live), {np.median(oc[:, ~live].std(axis=0)):.2f} (dead); "
          f"correlation between channels live-live {Co[np.ix_(live, live)][off_diag[np.ix_(live, live)]].mean():.2f}, "
          f"dead-dead {Co[np.ix_(~live, ~live)][off_diag[np.ix_(~live, ~live)]].mean():.2f}, live-dead {Co[np.ix_(live, ~live)].mean():.2f}")
    # Gaussian shape of the smoothed live mean
    pooled = Ls.ravel()
    k4 = float(np.mean(pooled ** 4) / np.mean(pooled ** 2) ** 2)
    print(f"  pooled smoothed live mean: rms {pooled.std():.3f}, kurtosis {k4:.2f} (Gaussian: 3)")

    # ---- among whom -----------------------------------------------------------------------------------
    cov /= nsm_acc
    d = np.sqrt(np.diag(cov)); C = cov / np.outer(d, d)
    ll = C[np.ix_(live, live)][off_diag[np.ix_(live, live)]].mean()
    dd = C[np.ix_(~live, ~live)][off_diag[np.ix_(~live, ~live)]].mean()
    ld = C[np.ix_(live, ~live)].mean()
    print(f"  correlation of the smoothed channels over the tail: live-live {ll:.2f}, dead-dead {dd:.2f}, live-dead {ld:.2f}")
    # the four 4x4 matrices of the map: row-major 8x8, matrix = 2 * (row // 4) + col // 4
    mat = np.array([2 * ((c // 8) // 4) + (c % 8) // 4 for c in channels])
    within, across = [], []
    for i in range(nch):
        for j in range(i + 1, nch):
            if live[i] and live[j]:
                (within if mat[i] == mat[j] else across).append(C[i, j])
    print(f"  live pairs within a 4x4 matrix {np.mean(within):.2f}, across matrices {np.mean(across):.2f}")
    blocks = np.full((4, 4), np.nan)
    for a in range(4):
        for b in range(4):
            sel_a = live & (mat == a); sel_b = live & (mat == b)
            sub = C[np.ix_(sel_a, sel_b)]
            blocks[a, b] = sub[off_diag[np.ix_(sel_a, sel_b)]].mean() if a == b else sub.mean()
    print("  mean correlation by matrix (rows/columns: top-left, top-right, bottom-left, bottom-right):")
    for a in range(4):
        print("    " + " ".join(f"{blocks[a, b]:5.2f}" for b in range(4)))
    Cl = C[np.ix_(live, live)]
    w, v = np.linalg.eigh(Cl)
    pc1 = v[:, -1] * np.sign(v[:, -1].sum())
    print(f"  first principal component of the live channels: {100 * w[-1] / w.sum():.0f}% of the variance, "
          f"loadings mean {pc1.mean():.3f}, spread {pc1.std():.3f} (uniform: {1 / np.sqrt(nl):.3f}, 0)")

    # ---- kicks ----------------------------------------------------------------------------------------
    K = args.kick
    lo, hi = max(0, pre - 10 * K), min(ns, pre + 10 * K)
    step = float(np.mean(kick_L[pre + 2 * K:pre + 10 * K]) - np.mean(kick_L[lo:pre - 2 * K]))
    # the search uses a boxcar of 2K+1 samples: any excursion of the live median of about that duration
    W = 2 * K + 1
    n_kick = 0; amps = []; filt_example = None; thr_example = None; thr_all = []
    for e in range(n):
        m = M_raw[e] - np.median(M_raw[e])
        fl = smooth(m, W)
        thr = args.nsigma * clipped_rms(fl)
        thr_all.append(thr)
        if filt_example is None:
            filt_example, thr_example = fl, thr
        above = np.flatnonzero(np.abs(fl) > thr)
        i = 0
        while i < len(above):
            j0 = above[i]
            grp = [j0]
            while i + 1 < len(above) and above[i + 1] - grp[-1] <= 2 * K + 1:
                i += 1; grp.append(above[i])
            jbest = max(grp, key=lambda j: abs(fl[j]))
            amps.append(fl[jbest]); n_kick += 1
            i += 1
    amps = np.array(amps)
    rate_us = n_kick / (n * args.tail)
    thr_med = float(np.median(thr_all))
    print(f"  around the trigger, average of the events: the live median steps by {step:+.2f} counts per channel from before to after "
          f"(the pulse in the median channel), the dead mean by {np.mean(kick_D[pre + 2 * K:pre + 10 * K]) - np.mean(kick_D[lo:pre - 2 * K]):+.2f}; no kick shape to use")
    print(f"  kicks in the tail (excursions of the live median above {args.nsigma:g} sigma = {thr_med:.2f} counts per channel, "
          f"{nl * thr_med:.0f} on the sum, over {W} samples = {W * dt:.0f} ns): {n_kick} in {n} events, {n_kick / n:.2f} per event, {rate_us:.3f} per us"
          + (f"; expected from the trigger rate {trig_rate / 1e3:.1f} kHz: {trig_rate * args.tail * 1e-6:.2f} per tail, "
             f"{trig_rate * 1e-6:.3f} per us" if trig_rate else ""))
    if n_kick:
        print(f"    amplitudes: median {np.median(np.abs(amps)):.2f} counts per channel, fraction positive {np.mean(amps > 0):.2f}")

    # ---- figures --------------------------------------------------------------------------------------
    figs = {}
    nt = min(args.traces, n)
    fig, axs = plt.subplots(nt, 1, figsize=(12, 2.6 * nt), sharex=True, squeeze=False)
    ylim = 4 * max(np.median(r["Ls"]), white["Ds"]) + 0.5
    for k in range(nt):
        ax = axs[k, 0]
        ax.plot(t_us[tail], L_raw[k], ".", ms=2, color="C0", alpha=0.35)
        ax.plot(t_us[tail], D_raw[k], ".", ms=2, color="C3", alpha=0.35)
        ax.plot(t_us[tail], smooth(L_raw[k], N), color="C0", lw=1.5, label=f"mean of the {nl} live channels")
        ax.plot(t_us[tail], smooth(D_raw[k], N), color="C3", lw=1.5, label=f"mean of the {nd} dead channels")
        ax.set_ylim(L_raw[k].mean() - ylim, L_raw[k].mean() + ylim)
        ax.set_ylabel("ADC counts"); ax.grid(lw=0.3)
        ax.text(0.01, 0.9, f"event {kept[k]}", transform=ax.transAxes, fontsize=9)
        if k == 0:
            ax.legend(fontsize=8, loc="upper right")
    axs[-1, 0].set_xlabel("time (us)")
    fig.suptitle(f"{run}: the tail of the record, points every 8 ns, lines the {N}-sample running mean", x=0.01, ha="left")
    fig.tight_layout(); figs["traces"] = fig

    fig, axs = plt.subplots(1, 3, figsize=(15, 4.2))
    for key, col, lab in (("L", "C0", "live mean"), ("D", "C3", "dead mean")):
        axs[0].hist(r[key], bins=40, histtype="step", color=col, lw=1.2, label=lab)
        axs[0].axvline(white[key], color=col, ls="--", lw=0.8)
        axs[1].hist(r[key + "s"], bins=40, histtype="step", color=col, lw=1.2, label=lab)
        axs[1].axvline(white[key + "s"], color=col, ls="--", lw=0.8)
    axs[0].set(xlabel="rms over the tail (ADC), 8 ns samples", ylabel="events", title="raw (dashed: white noise)")
    axs[1].set(xlabel=f"rms over the tail (ADC), {N}-sample running mean", ylabel="events", title="smoothed (dashed: white noise)")
    hv, he = np.histogram(pooled, bins=80)
    cc = 0.5 * (he[1:] + he[:-1])
    axs[2].stairs(hv, he, color="C0", fill=True, alpha=0.5, label="smoothed live mean, all samples")
    g = hv.sum() * (he[1] - he[0]) * np.exp(-cc ** 2 / (2 * pooled.var())) / np.sqrt(2 * np.pi * pooled.var())
    axs[2].plot(cc, g, color="k", lw=1, label="Gaussian, same rms")
    axs[2].set(xlabel="ADC", ylabel="samples", title=f"shape (kurtosis {k4:.2f})", yscale="log")
    for ax in axs:
        ax.grid(lw=0.3); ax.legend(fontsize=8)
    fig.suptitle(f"{run}: how much of the tail noise is common, {n} events", x=0.01, ha="left")
    fig.tight_layout(); figs["rms"] = fig

    fig, ax = plt.subplots(figsize=(8, 7))
    im = ax.imshow(C, vmin=-0.2, vmax=0.6, cmap="RdBu_r")
    for k, c in enumerate(channels):
        if not live[k]:
            ax.axhline(k, color="0.5", lw=0.3); ax.axvline(k, color="0.5", lw=0.3)
    for b in (16, 32, 48):
        ax.axhline(b - 0.5, color="k", lw=0.5); ax.axvline(b - 0.5, color="k", lw=0.5)
    ax.set(xlabel="channel", ylabel="channel", title=f"correlation of the {N}-sample running means over the tail (grey: dead channels)")
    ax.title.set_fontsize(9)
    fig.colorbar(im, ax=ax, shrink=0.8)
    fig.suptitle(f"{run}: live-live {ll:.2f}, dead-dead {dd:.2f}, live-dead {ld:.2f}", x=0.01, ha="left")
    fig.tight_layout(); figs["corr"] = fig

    fig, ax = plt.subplots(figsize=(8, 4.5))
    ax.plot(freq[1:], PL[1:], color="C0", lw=0.8, label="live mean")
    ax.plot(freq[1:], PD[1:], color="C3", lw=0.8, label="dead mean")
    ax.plot(freq[1:], PC[1:], color="C2", lw=0.8, label="live minus dead")
    ax.set(xscale="log", yscale="log", xlabel="frequency (MHz)", ylabel="power", title=f"{run}: power spectra over the tail, {n} events")
    ax.title.set_fontsize(10); ax.grid(lw=0.3, which="both"); ax.legend(fontsize=8)
    fig.tight_layout(); figs["spectra"] = fig

    fig, axs = plt.subplots(2, 2, figsize=(13, 8))
    tt = (np.arange(lo, hi) - pre) * dt
    axs[0, 0].plot(tt, kick_L[lo:hi] - np.median(kick_L[lo:hi]), color="C0", lw=1, label="live median")
    axs[0, 0].plot(tt, kick_D[lo:hi] - np.median(kick_D[lo:hi]), color="C3", lw=1, label="dead mean")
    axs[0, 0].axvline(0, color="k", lw=0.5)
    axs[0, 0].set(xlabel="time from the trigger (ns)", ylabel="ADC per channel", title="around the trigger, average of the events")
    axs[0, 1].hist(thr_all, bins=30, color="C0", alpha=0.7)
    axs[0, 1].set(xlabel="ADC per channel", ylabel="events", title=f"{args.nsigma:g} sigma threshold of the {W}-sample running mean, per event")
    axs[1, 0].plot(t_us[tail], filt_example, color="C0", lw=0.7)
    axs[1, 0].axhline(thr_example, color="C3", lw=0.8); axs[1, 0].axhline(-thr_example, color="C3", lw=0.8)
    axs[1, 0].set(xlabel="time (us)", ylabel="ADC per channel", title=f"event {kept[0]}: live median of the tail, {W}-sample running mean, {args.nsigma:g} sigma lines")
    if n_kick:
        axs[1, 1].hist(amps, bins=40, color="C0", alpha=0.7)
    axs[1, 1].set(xlabel="amplitude (ADC per channel)", ylabel="kicks", title=f"{n_kick} kicks in {n} events: {n_kick / n:.2f} per event"
                  + (f", expected {trig_rate * args.tail * 1e-6:.2f}" if trig_rate else ""))
    for ax in axs.flat:
        ax.grid(lw=0.3); ax.title.set_fontsize(9)
    axs[0, 0].legend(fontsize=8)
    fig.suptitle(f"{run}: excursions of the live median in the tail (kicks of the trigger logic?)", x=0.01, ha="left")
    fig.tight_layout(); figs["kick"] = fig

    if args.save:
        PLOTS.mkdir(exist_ok=True)
        for name, fg in figs.items():
            fg.savefig(PLOTS / f"{args.save}_cmnoise_{name}.png", dpi=130)
        print(f"  figures saved as {PLOTS}/{args.save}_cmnoise_{{{','.join(figs)}}}.png")
    else:
        plt.show()


if __name__ == "__main__":
    main()
