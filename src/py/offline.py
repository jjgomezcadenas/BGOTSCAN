#!/usr/bin/env python3
"""offline.py: the analysis of one run from the output of the pass: classes, charges, photopeaks, decay.

    python3 src/py/offline.py RUN|FILE [--window A,B] [--smooth 64] [--nsigma 5] [--k 3] [--gain FILE]
            [--dead LIST] [--decay-range A,B] [--decay-model 2]
            [--veto 2] [--save NAME] [--out FILE]

Reads pass/RUN_pass.h5 and the equalization table only. See docs/manual/offline_doc.tex.
"""
import argparse
import csv
import json
import re

import h5py
import numpy as np
import matplotlib.pyplot as plt

from bgo.io import PLOTS, REPO, find_pass, read_map
from bgo.signal import smooth, rms, clipped_rms, pulse_search, windows
from bgo.fit import fit_photopeak, find_photopeak, peak_bg, fit_decay, decay, fit_decay2, decay2
from bgo.plot import map8

DEAD = "19,20,21,22,23,25,27,33,35,41,46,47"      # chanstat_doc.tex, 2026-10-14
GAIN = REPO / "src" / "spe_gain_195K_41V_th20.csv"
CLASSES = ("EMPTY", "REAL", "TAIL", "SHORT")
ESTIMATORS = ("Q", "Q_all", "Q_0", "Q_kept", "Q_all_eq", "Q_0_eq", "Q_kept_eq", "Q_0_pos", "Q_0_posf", "Q_all_posf")


def attr_str(a):
    return a.decode() if isinstance(a, (bytes, np.bytes_)) else str(a)


def read_gain(path):
    """{channel: g} and the median A1 from the header of the table of spefit.py."""
    g, a1 = {}, np.nan
    with open(path) as f:
        for line in f:
            if line.startswith("#"):
                m = re.search(r"median A1 \(([0-9.]+)\)", line)
                if m:
                    a1 = float(m.group(1))
                continue
            break
    with open(path) as f:
        for row in csv.DictReader(l for l in f if not l.startswith("#")):
            if row["g"]:
                g[int(row["channel"])] = float(row["g"])
    return g, a1


def seed_peak(vals, nbins=60, top=None):
    """Photopeak position of one seed-channel population: the Gaussian-on-a-step fit of fit_photopeak when it
    works (100 events or more), else the mode of the smoothed histogram interpolated with a parabola through
    the three bins around its maximum. Returns (mu, 'fit' | 'mode')."""
    if len(vals) >= 100:
        fi = fit_photopeak(vals, nbins=nbins)
        if fi["ok"]:
            return fi["mu"], "fit"
    h, edges = np.histogram(vals, bins=nbins, range=(0, top))
    cc = 0.5 * (edges[1:] + edges[:-1])
    hs = np.convolve(h, np.ones(3) / 3, mode="same")
    k = int(np.clip(hs.argmax(), 1, nbins - 2))
    y0, y1, y2 = hs[k - 1], hs[k], hs[k + 1]
    d = y0 - 2 * y1 + y2
    shift = 0.5 * (y0 - y2) / d if d < 0 else 0.0
    return cc[k] + shift * (cc[1] - cc[0]), "mode"


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("run", help="run name (pass/RUN_pass.h5) or path to a pass file")
    ap.add_argument("--window", help="charge window of the sum: ns before,after the maximum (default: the pass window)")
    ap.add_argument("--smooth", type=int, default=64, help="moving-average width of the pulse search (samples)")
    ap.add_argument("--nsigma", type=float, default=5.0, help="significance of the charge cut, pile-up and multiplicity")
    ap.add_argument("--k", type=float, default=3.0, help="zero-suppression threshold of Q_kept (sigma)")
    ap.add_argument("--gain", default=str(GAIN), help="equalization table (spefit.py)")
    ap.add_argument("--dead", default=DEAD, help="dead channels, comma list")
    ap.add_argument("--decay-range", default="0.5,8", help="us after the maximum for the decay fit")
    ap.add_argument("--decay-model", type=int, default=2, choices=(1, 2), help="one or two exponentials")
    ap.add_argument("--veto", type=float, default=2.0, help="us beyond the charge window where no second maximum is searched")
    ap.add_argument("--save", metavar="NAME", help="save figures as src/py/plots/NAME_*.png")
    ap.add_argument("--out", metavar="FILE", help="per-event output (default pass/RUN_offline.h5)")
    args = ap.parse_args()
    dead = {int(v) for v in args.dead.split(",") if v.strip()}
    t0d, t1d = (float(v) * 1000 for v in args.decay_range.split(","))

    path = find_pass(args.run)
    with h5py.File(path, "r") as f:
        A = {k: f.attrs[k] for k in f.attrs}
        run = attr_str(A["run"])
        dt, pre = float(A["tsample_ns"]), int(A["pretrigger"])
        channels = [int(c) for c in A["channels"]]
        T = np.asarray(A["template"], float)
        win_pass = tuple(float(v) for v in A["window_ns"])
        scale = float(f["sum"].attrs["scale"])
        S_all = f["sum"][:].astype(np.float32) / scale           # (nev, ns)
        qc = f["qc"][:].astype(np.float64)                        # (nev, nch)
        qcc = f["qc_ctrl"][:].astype(np.float64)
        index = f["events/index"][:]
        ev = f["events"]
        # SHORT: the window or its after-window outside the record (flag of the pass); older passes with the adjacent
        # baseline: no room for the full pre-window or the after-window; passes before that: never
        if "short" in ev:
            short = ev["short"][:] > 0
        elif "room_before" in ev:
            short = ~((ev["room_before"][:] > 0) & (ev["room_after"][:] > 0))
        else:
            short = np.zeros(len(index), bool)
        ctrl_ok = ev["ctrl_ok"][:] > 0 if "ctrl_ok" in ev else np.ones(len(index), bool)
        # pre-window samples used in the per-channel baseline, and the full pre-window
        nb_full = int(round(float(A["adj_before_us"]) * 1000 / dt)) if "adj_before_us" in A else 0
        nb_used = ev["nb_used"][:] if "nb_used" in ev else np.full(len(index), nb_full)
    nev, ns = S_all.shape
    t = (np.arange(ns) - pre) * dt
    a_ns, b_ns = (float(v) for v in args.window.split(",")) if args.window else win_pass
    a, b = int(round(a_ns / dt)), int(round(b_ns / dt))
    veto = int(round(args.veto * 1000 / dt))
    g_tab, a1 = read_gain(args.gain)
    q_pe = a1 * T.sum() * dt if np.isfinite(a1) else np.nan

    # ---- per event: pulse, windows, Q ----------------------------------------------------
    imax = np.zeros(nev, int); H = np.zeros(nev); irise = np.full(nev, -1)
    Q = np.zeros(nev); Qc = np.zeros(nev); sN = np.zeros(nev); P2 = np.zeros(nev); V = np.zeros(nev); i2 = np.zeros(nev, int)
    for e in range(nev):
        S = S_all[e].astype(np.float64)
        Sh, im, h, ir = pulse_search(S, args.smooth)
        pwin, cwin, _ = windows(im, a, b, ns)
        imax[e], H[e], irise[e] = im, h, -1 if ir is None else ir
        Q[e] = dt * S[pwin].sum(); Qc[e] = dt * S[cwin].sum()
        sN[e] = rms(Sh[cwin])
        masked = Sh.copy(); masked[max(0, pwin.start - veto):min(ns, pwin.stop + veto)] = -np.inf
        j = int(masked.argmax()); i2[e] = j; P2[e] = masked[j]
        lo, hi = sorted((im, j)); V[e] = Sh[lo:hi + 1].min()
    sigma_Q = clipped_rms(Qc)
    real = H >= args.nsigma * sN                 # the height of the smoothed pulse against the smoothed noise
    cls = np.where(~real, 0, np.where(irise < 0, 2, np.where(short, 3, 1)))     # EMPTY, REAL, TAIL, SHORT
    pile = (P2 > args.nsigma * sN) & (P2 - V > args.nsigma * sN)
    dt_pile = (i2 - imax) * dt

    # ---- per event: charges from the channels ---------------------------------------------
    live = np.array([c not in dead for c in channels])
    g = np.array([g_tab.get(c, 1.0) for c in channels])
    no_g = [c for c, l in zip(channels, live) if l and c not in g_tab]
    okc = ctrl_ok & np.isfinite(qcc).all(axis=1)
    sq = np.array([clipped_rms(qcc[okc, k]) for k in range(len(channels))])
    qm = np.where(live, qc, 0.0)
    pos = qm > 0
    kept = qm > args.k * sq
    sig = qm > args.nsigma * sq
    Q_all = qm.sum(axis=1); Q_0 = (qm * pos).sum(axis=1); Q_kept = (qm * kept).sum(axis=1)
    qe = qm / g
    Q_all_eq = qe.sum(axis=1); Q_0_eq = (qe * pos).sum(axis=1); Q_kept_eq = (qe * kept).sum(axis=1)
    mult = (sig & live).sum(axis=1)
    est = dict(Q=Q, Q_all=Q_all, Q_0=Q_0, Q_kept=Q_kept, Q_all_eq=Q_all_eq, Q_0_eq=Q_0_eq, Q_kept_eq=Q_kept_eq)
    # the same sums over the control charges, with the masks of the event: the noise of each estimator
    qmc = np.where(live, np.nan_to_num(qcc), 0.0); qec = qmc / g
    ctrl = dict(Q=Qc, Q_all=qmc.sum(axis=1), Q_0=(qmc * pos).sum(axis=1), Q_kept=(qmc * kept).sum(axis=1),
                Q_all_eq=qec.sum(axis=1), Q_0_eq=(qec * pos).sum(axis=1), Q_kept_eq=(qec * kept).sum(axis=1))
    # the channel with the largest charge (a proxy for the position of the light) and its distance to the dead block
    cmap = read_map()
    xy = np.array([cmap.get(c, (np.nan, np.nan)) for c in channels])
    seed = np.where(live, qe, -np.inf).argmax(axis=1)
    dead_xy = np.array([cmap[c] for c in dead if c in cmap])
    d_dead = (np.sqrt(((xy[seed][:, None, :] - dead_xy[None, :, :]) ** 2).sum(axis=2)).min(axis=1)
              if len(dead_xy) else np.full(nev, np.inf))

    # ---- per run: classes, photopeaks, average pulse ---------------------------------------
    sel = (cls == 1) & ~pile
    frac = {n: float((cls == i).mean()) for i, n in enumerate(CLASSES)}
    pile_frac = float(pile[cls == 1].mean()) if (cls == 1).any() else np.nan
    # resolution against the distance of the seed channel to the dead block
    print(f"  seed channel distance to the nearest dead channel, Q_0_eq photopeak (REAL without pile-up):")
    for lab, m in (("on the block (< 6 mm)", d_dead < 6), ("next to it (6 to 12 mm)", (d_dead >= 6) & (d_dead < 12)),
                   ("far (>= 12 mm)", d_dead >= 12), ("far (>= 18 mm)", d_dead >= 18)):
        fi = fit_photopeak(est["Q_0_eq"][sel & m])
        print(f"    {lab:26s}: {int((sel & m).sum()):5d} events, " + (f"peak {fi['mu']:.4g}, FWHM {100 * fi['fwhm_rel']:.1f}%"
              if fi["ok"] else f"no fit ({fi['why']})"))
    # photopeak per seed channel: position (mode of the smoothed histogram) and width where the fit works,
    # on the unequalized and on the equalized charge (the seed is the same: the largest equalized charge)
    seed_ch = np.array(channels)[seed]
    n_map = {c: int((sel & (seed_ch == c)).sum()) for c in channels}
    by_seed = {}
    for name in ("Q_0", "Q_0_eq"):
        vals = est[name]
        pos_map, wid_map = {}, {}
        for c in channels:
            m = sel & (seed_ch == c)
            if n_map[c] >= 40:
                h, edges = np.histogram(vals[m], bins=60, range=(0, np.percentile(vals[sel], 99.5)))
                cc = 0.5 * (edges[1:] + edges[:-1])
                mu0, s0 = find_photopeak(h, cc, 0.25, 3)
                pos_map[c] = mu0
                fi = fit_photopeak(vals[m], nbins=60) if n_map[c] >= 100 else None
                wid_map[c] = fi["fwhm_rel"] if fi is not None and fi["ok"] else 2.355 * s0 / mu0
        by_seed[name] = (pos_map, wid_map)
        if pos_map:
            pv = np.array(list(pos_map.values())); wv = np.array(list(wid_map.values()))
            print(f"  {name} photopeak per seed channel ({len(pos_map)} channels with 40+ events): position median {np.median(pv):.4g}, "
                  f"spread (rms/median) {pv.std() / np.median(pv):.2f}, lowest/highest {pv.min() / np.median(pv):.2f}/{pv.max() / np.median(pv):.2f}; "
                  f"width of one seed channel: median {100 * np.median(wv):.0f}%, range {100 * wv.min():.0f} to {100 * wv.max():.0f}%")
    pos_map, wid_map = by_seed["Q_0_eq"]
    if pos_map:
        pv = np.array(list(pos_map.values()))
        pos0, wid0 = by_seed["Q_0"]
        med0 = np.median(list(pos0.values()))
        print("  per seed channel: events, Q_0 peak / median, Q_0_eq peak / median, gain factor, Q_0_eq FWHM")
        for c in sorted(pos_map, key=lambda c: pos_map[c] / np.median(pv)):
            print(f"    ch {c:2d}: {n_map[c]:4d}  {pos0[c] / med0:5.2f}  {pos_map[c] / np.median(pv):5.2f}  g {g_tab.get(c, 1.0):4.2f}  {100 * wid_map[c]:3.0f}%")
    # response correction by seed channel: the charge of the event divided by the photopeak position of its
    # seed channel relative to the median (channels with fewer than 40 events: no correction)
    f_seed = np.array([pos_map.get(c, np.nan) / np.median(pv) if pos_map else np.nan for c in channels])
    f_ev = f_seed[seed]
    no_f = sel & ~np.isfinite(f_ev)
    est["Q_0_pos"] = np.where(np.isfinite(f_ev), est["Q_0_eq"] / np.where(np.isfinite(f_ev), f_ev, 1.0), est["Q_0_eq"])
    ctrl["Q_0_pos"] = np.where(np.isfinite(f_ev), ctrl["Q_0_eq"] / np.where(np.isfinite(f_ev), f_ev, 1.0), ctrl["Q_0_eq"])
    print(f"  Q_0_pos = Q_0_eq / f(seed channel), f from the peaks above; {int(no_f.sum())} of {int(sel.sum())} events "
          f"with a seed channel without factor (left as they are)")
    # the same correction with a fitted (or interpolated) peak per seed channel, on the masked sum and on the
    # unmasked sum, each factor measured on the sum it corrects
    for new, base in (("Q_0_posf", "Q_0_eq"), ("Q_all_posf", "Q_all_eq")):
        vals = est[base]
        top = np.percentile(vals[sel], 99.5)
        pk, how = {}, {"fit": 0, "mode": 0}
        for c in channels:
            m = sel & (seed_ch == c)
            if n_map[c] >= 40:
                pk[c], h_ = seed_peak(vals[m], 60, top)
                how[h_] += 1
        if not pk:
            est[new] = vals.copy(); ctrl[new] = ctrl[base].copy()
            continue
        medp = np.median(list(pk.values()))
        fk = np.array([pk.get(c, np.nan) / medp for c in channels])[seed]
        fin = np.isfinite(fk)
        est[new] = np.where(fin, vals / np.where(fin, fk, 1.0), vals)
        ctrl[new] = np.where(fin, ctrl[base] / np.where(fin, fk, 1.0), ctrl[base])
        pv_ = np.array(list(pk.values())) / medp
        print(f"  {new} = {base} / f'(seed channel), f' from the peak of each seed channel on {base} ({how['fit']} fitted, "
              f"{how['mode']} interpolated modes): spread {pv_.std():.3f}, range {pv_.min():.2f} to {pv_.max():.2f}; "
              f"{int((sel & ~fin).sum())} events without factor")
    fits = {n: fit_photopeak(v[sel]) for n, v in est.items()}
    for name in ("Q_0_eq", "Q_0_pos", "Q_all_posf"):
        vals = est[name]
        print(f"  {name} by multiplicity and time of the maximum:")
        for lab, m in (("multiplicity <= 6", mult <= 6), ("multiplicity 7 to 10", (mult > 6) & (mult <= 10)), ("multiplicity > 10", mult > 10),
                       ("t_max < 0.3 us", t[imax] < 300), ("t_max >= 0.3 us", t[imax] >= 300)):
            fi = fit_photopeak(vals[sel & m])
            print(f"    {lab:22s}: {int((sel & m).sum()):5d} events, " + (f"peak {fi['mu']:.4g}, FWHM {100 * fi['fwhm_rel']:.1f}%"
                  if fi["ok"] else f"no fit ({fi['why']})"))
    # the pre-window: the events with the full pre-window against those with a partial one or none
    pre_groups = (("full pre-window", nb_used >= nb_full), ("partial pre-window", (nb_used > 0) & (nb_used < nb_full)),
                  ("after-window only", nb_used == 0))
    pre_fits = {}
    if nb_full > 0:
        print(f"  pre-window of the per-channel baseline ({nb_full} samples when full), REAL without pile-up:")
        for lab, m in pre_groups:
            pre_fits[lab] = {n: fit_photopeak(est[n][sel & m]) for n in ("Q_0_eq", "Q_0_pos")}
            line = f"    {lab:19s}: {int((sel & m).sum()):5d} events"
            for n, fi in pre_fits[lab].items():
                line += f"; {n} " + (f"peak {fi['mu']:.4g} +- {fi['err']:.2g}, FWHM {100 * fi['fwhm_rel']:.1f}%" if fi["ok"] else f"no fit ({fi['why']})")
            print(line)
    # average pulse of the REAL events without pile-up, aligned on the maximum
    centre = 300
    acc = np.zeros(ns); cnt = np.zeros(ns)
    for e in np.flatnonzero(sel):
        sh = centre - imax[e]
        lo, hi = max(0, sh), min(ns, ns + sh)
        acc[lo:hi] += S_all[e, lo - sh:hi - sh]; cnt[lo:hi] += 1
    pulse = np.where(cnt > 0.5 * sel.sum(), acc / np.maximum(cnt, 1), np.nan)
    tp = (np.arange(ns) - centre) * dt
    ok = np.isfinite(pulse)
    dec = fit_decay2(tp[ok], pulse[ok], t0d, t1d) if args.decay_model == 2 else fit_decay(tp[ok], pulse[ok], t0d, t1d)
    if args.decay_model == 2 and dec["ok"] and (not 0.05 < dec["f1"] < 0.95 or not np.isfinite(dec["err1"] + dec["err2"])
                                               or dec["err2"] == 0 or abs(dec["tau2"] - dec["tau1"]) < 0.1 * dec["tau2"]
                                               or dec["err1"] > 0.5 * dec["tau1"] or dec["err2"] > 0.5 * dec["tau2"]):
        # one component carries all the light, or the two time constants agree within 10%, or an error is larger
        # than half its value: the two exponentials are degenerate
        dec = fit_decay(tp[ok], pulse[ok], t0d, t1d)
        args.decay_model = 1
        print("  decay: the two-exponential fit is degenerate, one exponential reported")
    if args.decay_model == 2 and dec["ok"]:
        dec["tau"], dec["err"] = dec["tau2"], dec["err2"]          # the slow component is "the" decay time
        dec_label = (f"tau1 {dec['tau1'] / 1000:.3f} +- {dec['err1'] / 1000:.3f} us ({100 * dec['f1']:.0f}% of the light), "
                     f"tau2 {dec['tau2'] / 1000:.3f} +- {dec['err2'] / 1000:.3f} us")
    elif dec["ok"]:
        dec_label = f"tau {dec['tau'] / 1000:.3f} +- {dec['err'] / 1000:.3f} us"

    # ---- text ------------------------------------------------------------------------------
    T3 = float(A["T3"]) if "T3" in A else np.nan
    print(f"{path.name}: run {run}, {nev} events, T3 {T3:.1f} K, bias {float(A['V_bias']):g} V; "
          f"window {a_ns:g}/{b_ns:g} ns (pass {win_pass[0]:g}/{win_pass[1]:g}), charge noise sigma_Q {sigma_Q:.4g} ADC ns; "
          f"q_pe {q_pe:.4g} ADC ns from A1 {a1:.2f} (gain table {args.gain.split('/')[-1]})")
    print(f"  classes: " + ", ".join(f"{n} {100 * frac[n]:.1f}%" for n in CLASSES)
          + f"; pile-up in REAL {100 * pile_frac:.1f}%; {sel.sum()} REAL without pile-up used")
    noise = {n: clipped_rms(ctrl[n][sel & okc]) if n in ctrl else np.nan for n in est}      # ADC ns, from the control charges
    # the same over the EMPTY events (no pulse in the record): the noise of the electronics with nothing before it
    empty = (cls == 0) & okc
    noise_e = {n: clipped_rms(ctrl[n][empty]) if n in ctrl and empty.sum() >= 20 else np.nan for n in est}
    print(f"  {'estimator':10s} {'peak ADCns':>11s} {'err':>7s} {'peak PE':>8s} {'PE/keV':>7s} {'sigma':>8s} {'FWHM/peak':>9s} {'n':>6s} "
          f"{'noise':>7s} {'noise/peak':>10s} {'after':>6s} {'nz.EMPTY':>8s} {'after.E':>7s}")
    after_e = {}
    for n, fi in fits.items():
        if fi["ok"]:
            nt = 2.355 * noise[n] / fi["mu"] if np.isfinite(noise[n]) else np.nan
            after = np.sqrt(max(fi["fwhm_rel"] ** 2 - nt ** 2, 0.0)) if np.isfinite(nt) else np.nan
            nte = 2.355 * noise_e[n] / fi["mu"] if np.isfinite(noise_e[n]) else np.nan
            after_e[n] = np.sqrt(max(fi["fwhm_rel"] ** 2 - nte ** 2, 0.0)) if np.isfinite(nte) else np.nan
            print(f"  {n:10s} {fi['mu']:11.0f} {fi['err']:7.0f} {fi['mu'] / q_pe:8.1f} {fi['mu'] / q_pe / 511:7.2f} "
                  f"{fi['s']:8.0f} {100 * fi['fwhm_rel']:8.1f}% {fi['n']:6d} {noise[n]:7.0f} {100 * nt:9.1f}% {100 * after:5.1f}% "
                  f"{noise_e[n]:8.0f} {100 * after_e[n]:6.1f}%")
        else:
            print(f"  {n:10s} no fit ({fi['why']})")
    print(f"  noise: clipped rms of the same sum over the control charges of the selected events (FWHM term 2.355 noise/peak); "
          f"after: the width with it subtracted in quadrature; {int((sel & okc).sum())} events with a valid control window")
    print(f"  nz.EMPTY: the same rms over the {int(empty.sum())} EMPTY events with a valid control window (no pulse in the record; "
          f"error on the rms {100 / np.sqrt(max(2 * empty.sum(), 1)):.0f}%); after.E: the width with that term subtracted")
    # the control sum against the charge of the event: does the pulse feed the tail?
    fq = fits["Q_all_eq"]
    vq, vc = est["Q_all_eq"], ctrl["Q_all_eq"]
    top = np.percentile(vq[sel], 99.5)
    edges = np.linspace(0, top, 13)
    cq, cmean, cerr, crms, cn = [], [], [], [], []
    for lo_, hi_ in zip(edges[:-1], edges[1:]):
        m = sel & okc & (vq >= lo_) & (vq < hi_)
        if m.sum() >= 30:
            v = vc[m]
            cq.append(0.5 * (lo_ + hi_)); cmean.append(v.mean()); cerr.append(v.std() / np.sqrt(m.sum())); crms.append(clipped_rms(v)); cn.append(int(m.sum()))
    cq, cmean, cerr, crms = (np.array(a) for a in (cq, cmean, cerr, crms))
    slope = np.polyfit(cq, cmean, 1, w=1 / cerr)[0] if len(cq) >= 3 else np.nan
    e_mean = vc[empty].mean() if empty.sum() else np.nan
    e_err = vc[empty].std() / np.sqrt(empty.sum()) if empty.sum() else np.nan
    print(f"  control sum (Q_all_eq) against the event charge: mean over EMPTY events {e_mean:.0f} +- {e_err:.0f} ADC ns; "
          f"slope of the mean against the charge {100 * slope:.2f}% ; rms in the bins: " + ", ".join(f"{r:.0f}" for r in crms))
    print(f"  decay ({args.decay_model} exp.): {dec_label} over [{t0d / 1000:g}, {t1d / 1000:g}] us after the maximum"
          if dec["ok"] else f"  decay: no fit ({dec['why']})")
    print(f"  pile-up: second maximum searched outside the window extended by {args.veto:g} us on each side")
    if no_g:
        print(f"  live channels without equalization factor (g = 1): {no_g}")

    # ---- output files --------------------------------------------------------------------
    out = args.out or str(path.parent / f"{run}_offline.h5")
    with h5py.File(out, "w") as f:
        for name, arr in (("index", index), ("cls", cls.astype(np.int8)), ("pileup", pile), ("dt_pile", dt_pile.astype(np.float32)),
                          ("t_max", t[imax].astype(np.float32)), ("H", H.astype(np.float32)),
                          ("t_rise", np.where(irise >= 0, t[np.maximum(irise, 0)], np.nan).astype(np.float32)),
                          ("Q", Q.astype(np.float32)), ("Q_ctrl", Qc.astype(np.float32)), ("mult", mult.astype(np.int16)),
                          ("seed", np.array(channels)[seed].astype(np.int16)), ("d_dead", d_dead.astype(np.float32))):
            f[name] = arr
        for n in ESTIMATORS[1:]:
            f[n] = est[n].astype(np.float32)
        f["pulse"] = pulse.astype(np.float32); f["pulse_t"] = tp.astype(np.float32)
        f.attrs.update(dict(run=run, pass_file=str(path), gain_file=str(args.gain), window_ns=[a_ns, b_ns], smooth=args.smooth,
                            nsigma=args.nsigma, k=args.k, dead=sorted(dead), decay_range_us=[t0d / 1000, t1d / 1000],
                            sigma_Q=sigma_Q, q_pe=q_pe))
    summary = dict(run=run, nev=int(nev), T1=float(A.get("T1", np.nan)), T2=float(A.get("T2", np.nan)), T3=T3,
                   V_bias=float(A["V_bias"]), thr_adc=int(A["thr_adc"]), trigger_rate_hz=float(A["trigger_rate_hz"]),
                   window_ns=[a_ns, b_ns], pass_window_ns=list(win_pass), nsigma=args.nsigma, k=args.k,
                   sigma_Q=float(sigma_Q), q_pe=float(q_pe), A1=float(a1), classes=frac, pileup_frac=pile_frac,
                   n_used=int(sel.sum()), decay_us=dec["tau"] / 1000 if dec["ok"] else None,
                   decay_err_us=dec["err"] / 1000 if dec["ok"] else None, decay_range_us=[t0d / 1000, t1d / 1000],
                   decay_model=args.decay_model, veto_us=args.veto,
                   decay_fast_us=dec.get("tau1", np.nan) / 1000 if dec["ok"] and args.decay_model == 2 else None,
                   decay_fast_err_us=dec.get("err1", np.nan) / 1000 if dec["ok"] and args.decay_model == 2 else None,
                   decay_fast_fraction=dec.get("f1") if dec["ok"] and args.decay_model == 2 else None,
                   no_gain=no_g, nb_full=nb_full,
                   prewindow={lab: {n: dict(mu=fi["mu"], err=fi["err"], fwhm_rel=fi["fwhm_rel"], n=fi["n"], ok=fi["ok"])
                                    for n, fi in d.items()} for lab, d in pre_fits.items()},
                   estimators={n: dict(mu=fi["mu"], err=fi["err"], s=fi["s"], fwhm_rel=fi["fwhm_rel"], n=fi["n"], ok=fi["ok"],
                                       noise=float(noise[n]) if np.isfinite(noise[n]) else None,
                                       noise_empty=float(noise_e[n]) if np.isfinite(noise_e[n]) else None)
                               for n, fi in fits.items()})
    with open(out.replace("_offline.h5", "_offline.json"), "w") as fo:
        json.dump(summary, fo, indent=1, default=float)
    print(f"  written {out} and {out.replace('_offline.h5', '_offline.json')}")

    # ---- figures ---------------------------------------------------------------------------
    figs = {}
    fig, axs = plt.subplots(2, 5, figsize=(19, 7.5))
    for ax, (n, fi) in zip(axs.flat, fits.items()):
        if fi["h"] is not None:
            ax.stairs(fi["h"], fi["edges"], color="C0", fill=True, alpha=0.6)
            if fi["ok"]:
                x = np.linspace(fi["lo"], fi["hi"], 200)
                ax.plot(x, peak_bg(x, *fi["p"]), color="C3", lw=1.3,
                        label=f"peak {fi['mu']:.3g}, FWHM {100 * fi['fwhm_rel']:.1f}%")
                ax.legend(fontsize=7)
        ax.set(title=n, xlabel="ADC ns", ylabel="events")
        ax.title.set_fontsize(9)
        ax.grid(lw=0.3)
    fig.suptitle(f"{run}: charge spectra of the REAL events without pile-up ({sel.sum()})", x=0.01, ha="left")
    fig.tight_layout(); figs["spectra"] = fig

    fig, ax = plt.subplots(figsize=(8, 4.5))
    top = np.percentile(Q, 99.5)
    for i, n in enumerate(CLASSES):
        if (cls == i).any():
            ax.hist(Q[cls == i], bins=200, range=(0, top), histtype="step", lw=1.2, label=f"{n} ({100 * frac[n]:.1f}%)")
    ax.hist(Q[(cls == 1) & pile], bins=200, range=(0, top), histtype="step", lw=1.0, ls="--", label=f"REAL with pile-up ({100 * pile_frac:.1f}% of REAL)")
    ax.set(xlabel="Q (ADC ns)", ylabel="events", title=f"{run}: Q of the sum by class")
    ax.title.set_fontsize(10); ax.grid(lw=0.3); ax.legend(fontsize=8)
    fig.tight_layout(); figs["classes"] = fig

    fig, axs = plt.subplots(1, 2, figsize=(12, 4.2))
    for ax, yscale in zip(axs, ("linear", "log")):
        ax.plot(tp[ok] / 1000, pulse[ok], color="C0", lw=0.8, label=f"mean of {sel.sum()} REAL events")
        if dec["ok"]:
            s2 = dec["sel"]
            model = (decay2(tp[ok][s2] - t0d, dec["A1"], dec["tau1"], dec["A2"], dec["tau2"], dec["C"]) if args.decay_model == 2
                     else decay(tp[ok][s2] - t0d, dec["A"], dec["tau"], dec["C"]))
            ax.plot(tp[ok][s2] / 1000, model, color="C3", lw=1.5, label=dec_label)
        ax.set(xlabel="us from the maximum", ylabel="sum (ADC)", yscale=yscale)
        if yscale == "log":
            ax.set_ylim(max(1.0, np.nanmin(pulse[ok][pulse[ok] > 0]) if (pulse[ok] > 0).any() else 1.0), None)
        ax.grid(lw=0.3); ax.legend(fontsize=8)
    fig.suptitle(f"{run}: average pulse and decay fit", x=0.01, ha="left")
    fig.tight_layout(); figs["pulse"] = fig

    fig, ax = plt.subplots(figsize=(8, 4))
    if pile.any():
        ax.hist(dt_pile[pile] / 1000, bins=60, color="C0", alpha=0.7)
    ax.axvspan(-a_ns / 1000, b_ns / 1000, color="0.5", alpha=0.15, lw=0, label="charge window")
    ax.axvspan(-a_ns / 1000 - args.veto, b_ns / 1000 + args.veto, color="0.5", alpha=0.08, lw=0, label=f"veto {args.veto:g} us")
    ax.set(xlabel="time of the second maximum minus t_max (us)", ylabel="events", title=f"{run}: pile-up, {int(pile.sum())} events")
    ax.title.set_fontsize(10); ax.grid(lw=0.3); ax.legend(fontsize=8)
    fig.tight_layout(); figs["pileup"] = fig

    fig, axs = plt.subplots(1, 2, figsize=(12, 4.2))
    if fq["ok"]:
        for ax in axs:
            ax.axvspan(fq["mu"] - fq["s"], fq["mu"] + fq["s"], color="C3", alpha=0.12, lw=0, label="photopeak +- 1 sigma")
    axs[0].errorbar(cq, cmean, yerr=cerr, fmt="o", color="C0", ms=4, label="REAL without pile-up, by charge")
    axs[0].axhline(e_mean, color="C2", lw=1.2, label=f"EMPTY events: {e_mean:.0f} +- {e_err:.0f}")
    axs[0].axhspan(e_mean - e_err, e_mean + e_err, color="C2", alpha=0.2, lw=0)
    axs[0].set(xlabel="Q_all_eq of the event (ADC ns)", ylabel="mean control sum (ADC ns)", title=f"mean of the control sum; slope {100 * slope:.2f}%")
    axs[1].plot(cq, crms, "o", color="C0", ms=4, label="REAL without pile-up, by charge")
    axs[1].axhline(noise_e["Q_all_eq"], color="C2", lw=1.2, label=f"EMPTY events: {noise_e['Q_all_eq']:.0f}")
    axs[1].axhline(noise["Q_all_eq"], color="C1", lw=1.2, ls="--", label=f"all REAL: {noise['Q_all_eq']:.0f}")
    axs[1].set(xlabel="Q_all_eq of the event (ADC ns)", ylabel="clipped rms of the control sum (ADC ns)", title="width of the control sum", ylim=(0, None))
    for ax in axs:
        ax.title.set_fontsize(10); ax.grid(lw=0.3); ax.legend(fontsize=7)
    fig.suptitle(f"{run}: the control sum of the live channels against the charge of the event", x=0.01, ha="left")
    fig.tight_layout(); figs["ctrlvsq"] = fig

    if pos_map:
        fig, axs = plt.subplots(1, 3, figsize=(16, 5))
        map8(axs[0], {c: v / np.median(pv) for c, v in pos_map.items()}, "photopeak position by seed channel (median = 1)")
        map8(axs[1], {c: 100 * v for c, v in wid_map.items()}, "FWHM by seed channel (%)", "{:.0f}")
        map8(axs[2], {c: v for c, v in n_map.items() if v > 0}, "events per seed channel", "{:.0f}")
        fig.suptitle(f"{run}: Q_0_eq by seed channel, REAL without pile-up", x=0.01, ha="left")
        fig.tight_layout(); figs["seedmap"] = fig
    if pre_fits:
        fig, axs = plt.subplots(1, 2, figsize=(12, 4.2))
        for ax, n in zip(axs, ("Q_0_eq", "Q_0_pos")):
            top = np.percentile(est[n][sel], 99.5)
            for i, (lab, m) in enumerate(pre_groups):
                fi = pre_fits[lab][n]
                if (sel & m).sum() == 0:
                    continue
                label = f"{lab}: {int((sel & m).sum())} events" + (f", peak {fi['mu']:.3g}, FWHM {100 * fi['fwhm_rel']:.1f}%" if fi["ok"] else "")
                if fi["h"] is not None:                              # the histogram the fit was made on
                    ax.stairs(fi["h"], fi["edges"], baseline=None, color=f"C{i}", lw=1.2, label=label)
                else:
                    ax.hist(est[n][sel & m], bins=200, range=(0, top), histtype="step", color=f"C{i}", lw=1.2, label=label)
                if fi["ok"]:
                    x = np.linspace(fi["lo"], fi["hi"], 200)
                    ax.plot(x, peak_bg(x, *fi["p"]), color=f"C{i}", lw=1.0, ls="--")
            ax.set(title=n, xlabel="ADC ns", ylabel="events", xlim=(0, top)); ax.title.set_fontsize(9); ax.grid(lw=0.3); ax.legend(fontsize=7)
        fig.suptitle(f"{run}: the photopeak by pre-window of the baseline, REAL without pile-up", x=0.01, ha="left")
        fig.tight_layout(); figs["prewindow"] = fig
    if args.save:
        PLOTS.mkdir(exist_ok=True)
        for tag, fg in figs.items():
            fg.savefig(PLOTS / f"{args.save}_{tag}.png", dpi=130)
        print(f"  figures saved as {PLOTS}/{args.save}_{{{','.join(figs)}}}.png")
    else:
        plt.show()


if __name__ == "__main__":
    main()
