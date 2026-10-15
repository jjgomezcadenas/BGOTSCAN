#!/usr/bin/env python3
"""wavefit.py: the slow wave per channel in the tail, and a baseline model that follows it, tested on an empty window.

    python3 src/py/wavefit.py RUN [--events 0-299] [--tail 12] [--window 1,5] [--before 1] [--smooth 62]
            [--band 0.2,1.0] [--dead LIST] [--offline FILE] [--examples 3] [--save NAME]

Six figures: autocorr, spectra, examples, hist, gainmap, sum. See docs/manual/wavefit_doc.tex.
"""
import argparse
import sys
from pathlib import Path

import h5py
import numpy as np
import matplotlib.pyplot as plt

from bgo.io import DATA, PASS, PLOTS, find_run, parse_events
from bgo.signal import smooth, clipped_rms
from bgo.plot import map8

DEAD = "19,20,21,22,23,25,27,33,35,41,46,47"      # chanstat_doc.tex, 2026-10-14


def autocorr(v, nlag):
    """Autocorrelation of v (mean removed) for lags 0..nlag-1, normalised at lag 0."""
    v = v - v.mean()
    n = len(v)
    f = np.fft.rfft(v, 2 * n)
    ac = np.fft.irfft(f * np.conj(f))[:nlag]
    return ac / ac[0]


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("run", help="run name (195K_41V_th20) or path to a .h5 file")
    ap.add_argument("--data", default=str(DATA), help="directory of the run files")
    ap.add_argument("--events", default="0-299", help="range a-b or list a,b,c (0-based)")
    ap.add_argument("--tail", type=float, default=12.0, help="the tail: last T us of the record")
    ap.add_argument("--window", default="1,5", help="empty window: start (us into the tail), length (us)")
    ap.add_argument("--before", type=float, default=1.0, help="us before the window used by the fit")
    ap.add_argument("--after", type=float, default=2.0, help="us right after the window averaged for the 'after' baseline (no slope)")
    ap.add_argument("--smooth", type=int, default=62, help="running-mean width for the autocorrelation (samples)")
    ap.add_argument("--band", default="0.2,1.0", help="MHz band of the wave frequency search")
    ap.add_argument("--dead", default=DEAD, help="dead channels, comma list")
    ap.add_argument("--offline", help="per-event file of offline.py (default pass/RUN_offline.h5; none: every event)")
    ap.add_argument("--examples", type=int, default=3, help="events drawn in the by-eye figure")
    ap.add_argument("--save", metavar="NAME", help="save the figures as src/py/plots/NAME_wavefit_*.png")
    args = ap.parse_args()
    dead = {int(v) for v in args.dead.split(",") if v.strip()}
    w_start, w_len = (float(v) for v in args.window.split(","))
    band = tuple(float(v) for v in args.band.split(","))

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
    base = np.asarray(f.attrs["baselines"], float)[:, None]
    ntail = int(round(args.tail * 1000 / dt))
    t0 = ns - ntail
    tail = slice(t0, ns)
    live = np.array([c not in dead for c in channels]); nl = int(live.sum()); nd = int((~live).sum())
    # the empty window, the fit region and the end-of-record region, as sample indices in the tail
    i_w0 = int(round(w_start * 1000 / dt)); i_w1 = i_w0 + int(round(w_len * 1000 / dt))
    i_f0 = max(0, i_w0 - int(round(args.before * 1000 / dt)))
    if i_w1 + 250 > ntail:
        sys.exit("the window reaches into the end-of-record region: shorten --window or lengthen --tail")
    W = slice(i_w0, i_w1); nW = i_w1 - i_w0
    fit_idx = np.r_[np.arange(i_f0, i_w0), np.arange(i_w1, ntail)]
    E = slice(ntail - 250, ntail)
    tt = np.arange(ntail) * dt * 1e-9                 # seconds, from the start of the tail
    t_us = (t0 + np.arange(ntail)) * dt / 1000

    off = Path(args.offline) if args.offline else PASS / f"{run}_offline.h5"
    if off.is_file():
        with h5py.File(off, "r") as fo:
            cls, pile, tmax = fo["cls"][:], fo["pileup"][:], fo["t_max"][:]
        clean = (cls == 1) & ~pile & (tmax < t0 * dt - 4000)
        kept = [e for e in events if clean[e]]
        print(f"{run}: {len(kept)} of {len(events)} events kept (REAL without pile-up, maximum > 4 us before the tail), from {off.name}")
    else:
        kept = list(events)
        print(f"{run}: {len(kept)} events, no offline file ({off}): every event used")
    n = len(kept)
    if n < 10:
        sys.exit("fewer than 10 clean events")
    print(f"  tail: last {args.tail:g} us; empty window {w_start:g} to {w_start + w_len:g} us into it ({nW} samples); "
          f"fit region: {args.before:g} us before the window and everything after it; end-of-record baseline: last 250 samples")

    # ---- first pass: spectra per channel and the wave frequency ----------------------------------------
    X = np.empty((n, nch, ntail), np.float32)
    for k, e in enumerate(kept):
        X[k] = (wf[e].astype(np.float64) - base)[:, tail]
    freq = np.fft.rfftfreq(ntail, d=dt * 1e-9) / 1e6
    P = np.mean(np.abs(np.fft.rfft(X - X.mean(axis=2, keepdims=True), axis=2)) ** 2, axis=0)   # (nch, nfreq)
    inband = (freq >= band[0]) & (freq <= band[1])
    fpk = np.array([freq[inband][int(P[k, inband].argmax())] for k in range(nch)])
    f0 = float(np.median(fpk[live]))
    print(f"  wave frequency: peak of the spectrum in {band[0]:g} to {band[1]:g} MHz, live channels median {f0:.3f} MHz "
          f"(spread {fpk[live].std():.3f}), dead channels median {np.median(fpk[~live]):.3f} MHz")

    # ---- autocorrelation of the smoothed tail ----------------------------------------------------------
    nlag = int(round(10000 / dt))
    Xs = smooth(X, args.smooth)
    ac = np.zeros((nch, nlag))
    for k in range(n):
        for c in range(nch):
            ac[c] += autocorr(Xs[k, c], nlag)
    ac /= n
    lag_us = np.arange(nlag) * dt / 1000
    ac_cls = {"dead": ac[~live].mean(axis=0), "live": ac[live].mean(axis=0)}
    ac_num = {}
    for name, a in ac_cls.items():
        # skip the white-noise spike (the smoothing width) before looking for the zero crossing
        start = args.smooth
        neg = np.flatnonzero(a[start:] < 0)
        if len(neg):
            iz = start + neg[0]
            imin = iz + int(a[iz:min(nlag, iz + 2 * (iz - start) + 300)].argmin())
            pos = np.flatnonzero(a[imin:] > 0)
            imax2 = imin + int(a[imin:min(nlag, imin + 2 * (imin - iz) + 300)].argmax()) if len(pos) else None
            ac_num[name] = (lag_us[iz], lag_us[imin], a[imin], lag_us[imax2] if imax2 else np.nan, a[imax2] if imax2 else np.nan)
            print(f"  autocorrelation, {name} channels ({args.smooth}-sample smoothing): first zero at {lag_us[iz]:.2f} us, "
                  f"first negative lobe {a[imin]:+.2f} at {lag_us[imin]:.2f} us, second positive lobe "
                  + (f"{a[imax2]:+.2f} at {lag_us[imax2]:.2f} us" if imax2 else "none"))
        else:
            ac_num[name] = None
            print(f"  autocorrelation, {name} channels: no zero crossing within 10 us")

    # ---- the empty window: three baselines -------------------------------------------------------------
    tf = tt[fit_idx]
    A_line = np.c_[np.ones_like(tf), tf]
    A_wave = np.c_[np.ones_like(tf), tf, np.sin(2 * np.pi * f0 * 1e6 * tf), np.cos(2 * np.pi * f0 * 1e6 * tf)]
    tw = tt[W]
    B_line = np.c_[np.ones_like(tw), tw]
    B_wave = np.c_[np.ones_like(tw), tw, np.sin(2 * np.pi * f0 * 1e6 * tw), np.cos(2 * np.pi * f0 * 1e6 * tw)]
    # least squares for all channels and events at once: Y (nfit, n*nch)
    Y = X[:, :, fit_idx].reshape(n * nch, -1).T
    c_line = np.linalg.lstsq(A_line, Y, rcond=None)[0]
    c_wave = np.linalg.lstsq(A_wave, Y, rcond=None)[0]
    XW = X[:, :, W].reshape(n * nch, -1)                             # (n*nch, nW)
    q_E = dt * (XW - X[:, :, E].mean(axis=2).reshape(-1, 1)).sum(axis=1)
    q_line = dt * (XW - (B_line @ c_line).T).sum(axis=1)
    q_wave = dt * (XW - (B_wave @ c_wave).T).sum(axis=1)
    i_a1 = min(ntail, i_w1 + int(round(args.after * 1000 / dt)))
    q_A = dt * (XW - X[:, :, i_w1:i_a1].mean(axis=2).reshape(-1, 1)).sum(axis=1)
    q = {"end": q_E.reshape(n, nch), "after": q_A.reshape(n, nch), "line": q_line.reshape(n, nch), "wave": q_wave.reshape(n, nch)}
    # the pass: the 125-sample running mean of the dead-channel mean subtracted from every channel, then the end-of-record baseline
    cm = smooth(X[:, ~live, :].mean(axis=1), 125)                                   # (n, ntail)
    Xc = X - cm[:, None, :]
    q["end+cm"] = dt * (Xc[:, :, W] - Xc[:, :, E].mean(axis=2, keepdims=True)).sum(axis=2)
    sig_w = np.median([clipped_rms(X[k, c, fit_idx]) for k in range(min(n, 50)) for c in range(nch)])
    white = sig_w * np.sqrt(nW) * dt
    scat = {}
    print(f"  empty-window charge, scatter over channels and events (ADC ns); white-noise floor {white:.0f}:")
    for name, sel_c in (("dead", ~live), ("live", live)):
        scat[name] = {m: clipped_rms(q[m][:, sel_c].ravel()) for m in q}
        s = scat[name]
        print(f"    {name:4s} channels: end-of-record {s['end']:.0f}, mean of the {args.after:g} us after the window {s['after']:.0f} "
              f"(gain {s['end'] / s['after']:.2f}), line {s['line']:.0f} (gain {s['end'] / s['line']:.2f}), "
              f"line + wave at {f0:.2f} MHz {s['wave']:.0f} (gain {s['end'] / s['wave']:.2f}); "
              f"end-of-record after the dead-channel common mode (the pass) {s['end+cm']:.0f}")
    gain_ch = {c: clipped_rms(q["end"][:, k]) / clipped_rms(q["wave"][:, k]) for k, c in enumerate(channels)}
    qs = {m: q[m][:, live].sum(axis=1) for m in q}
    scat["sum"] = {m: clipped_rms(qs[m]) for m in q}
    s = scat["sum"]
    print(f"    sum of the {nl} live channels: end-of-record {s['end']:.0f}, after {s['after']:.0f} (gain {s['end'] / s['after']:.2f}), "
          f"line {s['line']:.0f} (gain {s['end'] / s['line']:.2f}), "
          f"line + wave {s['wave']:.0f} (gain {s['end'] / s['wave']:.2f}); independent channels would give {np.sqrt(nl) * scat['live']['end']:.0f} / "
          f"{np.sqrt(nl) * scat['live']['wave']:.0f}")
    print(f"    sum after the dead-channel common mode (the pass): {s['end+cm']:.0f}; the noise of the {nd} dead channels' mean times {nl}: "
          f"{nl * scat['dead']['end'] / np.sqrt(nd):.0f}, in quadrature with the raw sum {np.hypot(s['end'], nl * scat['dead']['end'] / np.sqrt(nd)):.0f}")

    # ---- figures --------------------------------------------------------------------------------------
    figs = {}
    widths = np.array([1, 2, 4, 8, 16, 31, 62, 125, 250, 500, 625])
    rms_w = {}
    for name, sel_c in (("dead", ~live), ("live", live)):
        sub = X[:min(n, 60)][:, sel_c, :]
        rms_w[name] = [np.sqrt(np.mean([smooth(v, w)[w:ntail - w].var() for ev in sub for v in ev])) for w in widths]
    fig, axs = plt.subplots(1, 2, figsize=(14, 4.5))
    for name, col in (("dead", "C3"), ("live", "C0")):
        axs[0].plot(lag_us, ac_cls[name], color=col, lw=1.3, label=f"{name} channels")
        axs[1].plot(widths * dt / 1000, rms_w[name], "o" if name == "live" else "s", color=col, label=f"{name} channels")
    axs[0].axhline(0, color="k", lw=0.5)
    axs[0].set(xlabel="lag (us)", ylabel="autocorrelation", xlim=(0, 10),
               title=f"autocorrelation of the {args.smooth}-sample running mean, single channels, {n} events")
    axs[1].plot(widths * dt / 1000, sig_w / np.sqrt(widths), color="k", lw=0.8, ls="--", label=r"white noise, $\sigma/\sqrt{N}$")
    axs[1].set(xscale="log", yscale="log", xlabel="running-mean width (us)", ylabel="rms (counts)", title="rms of the running mean of one channel against its width")
    for ax in axs:
        ax.title.set_fontsize(9); ax.grid(lw=0.3, which="both"); ax.legend(fontsize=8)
    fig.suptitle(f"{run}: the slow noise of one channel over the tail", x=0.01, ha="left")
    fig.tight_layout(); figs["autocorr"] = fig

    fig, axs = plt.subplots(1, 2, figsize=(14, 5))
    for k in range(nch):
        axs[0].plot(freq[1:], P[k, 1:], color="C3" if not live[k] else "C0", lw=0.4, alpha=0.6)
    axs[0].axvspan(band[0], band[1], color="0.5", alpha=0.1, lw=0)
    axs[0].set(xscale="log", yscale="log", xlabel="frequency (MHz)", ylabel="power", title="single channels (red: dead, blue: live); grey: search band")
    axs[0].title.set_fontsize(9); axs[0].grid(lw=0.3, which="both")
    map8(axs[1], dict(zip(channels, fpk)), f"peak frequency (MHz), live median {f0:.2f}", "{:.2f}")
    fig.suptitle(f"{run}: the wave in the tail, {n} events", x=0.01, ha="left")
    fig.tight_layout(); figs["spectra"] = fig

    ne = min(args.examples, n)
    c_dead = int(np.flatnonzero(~live)[0]); c_live = int(np.flatnonzero(live)[0])
    fig, axs = plt.subplots(ne, 2, figsize=(14, 2.8 * ne), sharex=True, squeeze=False)
    for r in range(ne):
        for col_i, (c, lab) in enumerate(((c_dead, "dead"), (c_live, "live"))):
            ax = axs[r, col_i]
            x = X[r, c]
            j = r * nch + c
            ax.plot(t_us, smooth(x, 16), color="0.6", lw=0.6, label="waveform, 16-sample running mean")
            ax.axvspan(t_us[i_w0], t_us[i_w1 - 1], color="C2", alpha=0.12, lw=0, label="empty window")
            ax.axhline(x[E].mean(), color="k", lw=0.9, ls="--", label="end-of-record baseline")
            ax.plot(t_us, np.c_[np.ones_like(tt), tt] @ c_line[:, j], color="C1", lw=1, label="line")
            ax.plot(t_us, np.c_[np.ones_like(tt), tt, np.sin(2 * np.pi * f0 * 1e6 * tt), np.cos(2 * np.pi * f0 * 1e6 * tt)] @ c_wave[:, j],
                    color="C3", lw=1.2, label=f"line + wave at {f0:.2f} MHz")
            ax.set_ylim(np.median(x) - 4, np.median(x) + 4)
            ax.set_ylabel("ADC"); ax.grid(lw=0.3)
            ax.text(0.01, 0.88, f"event {kept[r]}, channel {channels[c]} ({lab})", transform=ax.transAxes, fontsize=9)
            if r == 0 and col_i == 0:
                ax.legend(fontsize=7, loc="lower right", ncol=2)
    axs[-1, 0].set_xlabel("time (us)"); axs[-1, 1].set_xlabel("time (us)")
    fig.suptitle(f"{run}: the baselines through the empty window", x=0.01, ha="left")
    fig.tight_layout(); figs["examples"] = fig

    fig, axs = plt.subplots(1, 2, figsize=(14, 4.5))
    for ax, (name, sel_c) in zip(axs, (("dead", ~live), ("live", live))):
        rng = 4 * scat[name]["end"]
        for m, col in (("end+cm", "C0"), ("end", "k"), ("after", "C2"), ("line", "C1"), ("wave", "C3")):
            ax.hist(q[m][:, sel_c].ravel() / 1e3, bins=80, range=(-rng / 1e3, rng / 1e3), histtype="step", color=col, lw=1.2,
                    label=f"{m}: scatter {scat[name][m] / 1e3:.2f}")
        ax.axvline(white / 1e3, color="0.5", lw=0.8, ls=":"); ax.axvline(-white / 1e3, color="0.5", lw=0.8, ls=":")
        ax.set(xlabel="charge of the empty window (1e3 ADC ns)", ylabel="channels x events", title=f"{name} channels; dotted: white-noise floor {white / 1e3:.2f}")
        ax.title.set_fontsize(9); ax.grid(lw=0.3); ax.legend(fontsize=8)
    fig.suptitle(f"{run}: empty-window charge with the three baselines, {n} events", x=0.01, ha="left")
    fig.tight_layout(); figs["hist"] = fig

    fig, ax = plt.subplots(figsize=(7, 6))
    map8(ax, gain_ch, "gain of line + wave over the end-of-record baseline (scatter ratio)", "{:.2f}")
    fig.suptitle(f"{run}", x=0.01, ha="left")
    fig.tight_layout(); figs["gainmap"] = fig

    fig, ax = plt.subplots(figsize=(8, 4.5))
    rng = 4 * scat["sum"]["end"]
    rng = 4 * scat["sum"]["end+cm"]
    for m, col in (("end+cm", "C0"), ("end", "k"), ("after", "C2"), ("line", "C1"), ("wave", "C3")):
        ax.hist(qs[m] / 1e3, bins=60, range=(-rng / 1e3, rng / 1e3), histtype="step", color=col, lw=1.2,
                label=f"{m}: scatter {scat['sum'][m] / 1e3:.1f}")
    ax.set(xlabel="charge of the empty window, sum of the live channels (1e3 ADC ns)", ylabel="events",
           title=f"{run}: the sum of the {nl} live channels, {n} events")
    ax.title.set_fontsize(10); ax.grid(lw=0.3); ax.legend(fontsize=8)
    fig.tight_layout(); figs["sum"] = fig

    if args.save:
        PLOTS.mkdir(exist_ok=True)
        for name, fg in figs.items():
            fg.savefig(PLOTS / f"{args.save}_wavefit_{name}.png", dpi=130)
        print(f"  figures saved as {PLOTS}/{args.save}_wavefit_{{{','.join(figs)}}}.png")
    else:
        plt.show()


if __name__ == "__main__":
    main()
