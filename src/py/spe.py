#!/usr/bin/env python3
"""spe.py: single-photoelectron spectra and equalization factors, with a matched filter.

    python3 src/py/spe.py RUN [--events 0-9] [--channel 0] [--template FILE] [--region tail|all]
            [--tail 10] [--nsigma 4.5] [--from-step N] [--save NAME] [--summary] [--csv FILE]

Interactive: Enter goes to the next step, n to the next event, q quits. See docs/spe_doc.tex.
"""
import argparse
import contextlib
import csv
import io
import json
import sys
from pathlib import Path

import h5py
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

from wfview import DATA, PLOTS, find_run, parse_events, smooth
from selection import Stepper, grid64, rms

STEPS = ("read the channel", "matched filter and candidates", "the pulses",
         "all 64 channels", "rates", "spectra, fit and equalization")
SRC = Path(__file__).resolve().parent.parent
ABIN, AMAX = 1.0, 80.0             # histogram of A (ADC)


# ---- definitions (the ones the JScan SPE pass will implement) ---------------------------
def search_region(kind, tail_us, ns, dt):
    if kind == "all":
        return slice(0, ns)
    return slice(max(0, ns - int(round(tail_us * 1000 / dt))), ns)


def find_candidates(x, region, nsm, thr, sep):
    """Local maxima of the smoothed x inside region, above thr, at least sep samples apart
    (taken in decreasing order of height). Returns (smoothed, sorted candidate indices)."""
    xs = smooth(x, nsm)
    lo, hi = region.start, region.stop
    seg = xs[lo:hi]
    i = np.arange(1, len(seg) - 1)
    loc = i[(seg[i] > thr) & (seg[i] >= seg[i - 1]) & (seg[i] > seg[i + 1])] + lo
    kept = []
    for j in loc[np.argsort(xs[loc])[::-1]]:
        if all(abs(j - k) >= sep for k in kept):
            kept.append(int(j))
    return xs, sorted(kept)


def read_template(path):
    T = []
    with open(path) as f:
        for row in csv.DictReader(l for l in f if not l.startswith("#")):
            T.append(float(row["value"]))
    T = np.array(T)
    if len(T) % 2 == 0 or abs(T[len(T) // 2] - 1) > 1e-6:
        sys.exit(f"template {path}: expected an odd number of values with 1 at the centre")
    return T


def matched(x, T):
    """Least-squares height of a pulse of shape T at every sample: y_i = sum_j T_j x_{i+j} / sum T^2.
    Returns y aligned with x (edges set to zero)."""
    H = (len(T) - 1) // 2
    y = np.zeros_like(x)
    y[H:len(x) - H] = np.correlate(x, T, mode="valid") / np.sum(T * T)
    return y


def clipped_rms(v, nsig=3.0):
    s = rms(v)
    return rms(v[np.abs(v - v.mean()) < nsig * s]) if s else 0.0


def channel_search(x, T, region, nsig):
    """Filtered waveform, its noise, the threshold and the candidates of one channel."""
    H = (len(T) - 1) // 2
    y = matched(x, T)
    inner = slice(max(region.start, H), min(region.stop, len(x) - H))
    sy = clipped_rms(y[inner])
    thr = nsig * sy
    _, cands = find_candidates(y, inner, 1, thr, 2 * H)
    res = [rms(x[i - H:i + H + 1] - y[i] * T) for i in cands]
    return y, sy, thr, cands, res


def gauss(a, n, mu, w):
    return n * np.exp(-0.5 * ((a - mu) / w) ** 2)


def fit_peak(A, thr, sy):
    """Gaussian fit of the single-photoelectron peak of the amplitudes A above the threshold.
    Returns dict(n, A1, err, w1, ok, why, h, edges, lo, hi)."""
    A = np.asarray(A)
    out = dict(n=len(A), A1=np.nan, err=np.nan, w1=np.nan, ok=False, why="", h=None, edges=None)
    if len(A) < 30:
        out["why"] = "fewer than 30 pulses"
        return out
    h, edges = np.histogram(A, bins=np.arange(0, AMAX + ABIN, ABIN))
    c = 0.5 * (edges[1:] + edges[:-1])
    lo = thr + sy
    vis = c >= lo
    if not vis.any() or h[vis].sum() < 30:
        out.update(why="fewer than 30 pulses above threshold", h=h, edges=edges)
        return out
    mode = c[vis][int(h[vis].argmax())]
    hi = 1.6 * mode
    sel = (c >= lo) & (c <= hi)
    out.update(h=h, edges=edges, lo=lo, hi=hi)
    if h[sel].sum() < 30 or sel.sum() < 4:
        out["why"] = "fewer than 30 pulses in the fit range"
        return out
    try:
        p, cov = curve_fit(gauss, c[sel], h[sel], p0=(h[sel].max(), mode, max(sy, 1.0)),
                           bounds=([0, lo, 0.5 * sy], [np.inf, hi, hi]),
                           sigma=np.sqrt(np.maximum(h[sel], 1.0)), absolute_sigma=True, maxfev=5000)
    except RuntimeError:
        out["why"] = "fit did not converge"
        return out
    n, mu, w = p
    w = abs(w)
    err = float(np.sqrt(cov[1, 1])) if np.isfinite(cov[1, 1]) else np.nan
    out.update(A1=float(mu), err=err, w1=float(w))
    if not lo <= mu <= hi:
        out["why"] = "mean outside the fit range"
    elif not 0.5 * sy <= w <= mu / 2:
        out["why"] = "width out of bounds"
    else:
        out["ok"] = True
    return out


def crosstalk_ratio(values, m):
    v = np.asarray(values)
    low = np.sum(v < 1.5 * m)
    return np.sum((v >= 1.5 * m) & (v < 2.5 * m)) / low if low else float("nan")


# ---- steps -----------------------------------------------------------------------------
def step1(ev, st, out, args):
    x, t, region, T = out["x"], ev["t"], out["region"], ev["T"]
    print(f"step 1: {ev['run']} event {ev['e']}, channel {args.channel}: baseline {out['b']:.1f} ADC "
          f"(median over the record), rms {rms(x):.2f}; search region {args.region}: "
          f"[{t[region.start]:.0f}, {t[region.stop - 1]:.0f}] ns = samples [{region.start}, {region.stop - 1}]; "
          f"template {ev['tpath'].name}, {len(T)} samples, sum of squares {np.sum(T * T):.3f}")
    fig, ax = plt.subplots(figsize=(11, 3.6))
    ax.plot(t, x, color="C0", lw=0.5)
    ax.axvspan(t[region.start], t[region.stop - 1], color="0.5", alpha=0.15, lw=0, label="search region")
    ax.axvline(0, color="0.4", ls="--", lw=0.8)
    ax.set(xlabel="ns from trigger", ylabel=f"channel {args.channel} minus baseline (ADC)")
    ax.grid(lw=0.3)
    ax.legend(fontsize=8, loc="upper right")
    return st.show(1, fig)


def step2(ev, st, out, args):
    x, t, region, T = out["x"], ev["t"], out["region"], ev["T"]
    y, sy, thr, cands, res = channel_search(x, T, region, args.nsigma)
    out.update(y=y, sy=sy, thr=thr, cands=cands, res=res)
    print(f"step 2: matched filter: noise of the filtered waveform {sy:.2f} ADC, threshold {args.nsigma:g} sigma = "
          f"{thr:.1f} ADC: {len(cands)} candidates"
          + ("; (time ns, A): " + ", ".join(f"({t[i]:.0f}, {y[i]:.1f})" for i in cands) if cands else ""))
    fig, ax = plt.subplots(figsize=(11, 3.6))
    ax.plot(t[region], x[region], color="0.7", lw=0.5, label="raw")
    ax.plot(t[region], y[region], color="C0", lw=0.9, label="matched filter")
    ax.axhline(thr, color="C3", lw=0.8, ls=":", label="threshold")
    if cands:
        ax.plot(t[cands], y[cands], "v", color="C3", ms=7, label="candidates")
    ax.set(xlabel="ns from trigger", ylabel="ADC")
    ax.grid(lw=0.3)
    ax.legend(fontsize=8, loc="upper right")
    return st.show(2, fig)


def step3(ev, st, out, args):
    x, t, y, T = out["x"], ev["t"], out["y"], ev["T"]
    H = (len(T) - 1) // 2
    cands, res = out["cands"], out["res"]
    print("step 3: per candidate (time ns, A ADC, residual rms of the template fit):")
    for i, r in zip(cands, res):
        print(f"        ({t[i]:8.0f}, {y[i]:5.1f}, {r:4.1f})")
    n = min(6, len(cands))
    fig, axs = plt.subplots(1, max(n, 1), figsize=(3.2 * max(n, 1) + 1, 3.4), squeeze=False)
    j = np.arange(-H, H + 1)
    for ax, i, r in zip(axs.flat, cands[:n], res[:n]):
        lo, hi = i - H - 6, i + H + 7
        ax.plot(t[lo:hi], x[lo:hi], color="C0", lw=0.8, label="raw")
        ax.plot(t[i + j], y[i] * T, "o-", color="C3", ms=3, lw=1, label="A x template")
        ax.set_title(f"A {y[i]:.1f}, residual {r:.1f}", fontsize=9, loc="left")
        ax.set_xlabel("ns", fontsize=8)
        ax.tick_params(labelsize=7)
        ax.grid(lw=0.3)
    if n == 0:
        axs.flat[0].text(0.5, 0.5, "no candidates", ha="center", transform=axs.flat[0].transAxes)
    else:
        axs.flat[0].legend(fontsize=7)
    axs.flat[0].set_ylabel("ADC", fontsize=8)
    return st.show(3, fig)


def all_channels(ev, out, args):
    """Steps 1 and 2 on every channel; fills out['per'] = {channel: (y, sy, thr, cands)}."""
    if "per" in out:
        return
    per = {}
    for k, c in enumerate(ev["channels"]):
        x = ev["w"][k] - np.median(ev["w"][k])
        y, sy, thr, cands, _ = channel_search(x, ev["T"], out["region"], args.nsigma)
        per[int(c)] = (y, sy, thr, cands)
    out["per"] = per


def step4(ev, st, out, args):
    all_channels(ev, out, args)
    per, t, region = out["per"], ev["t"], out["region"]
    counts = {c: len(v[3]) for c, v in per.items()}
    tot = sum(counts.values())
    print(f"step 4: candidates per channel (8 x 8, channel order); total {tot}, mean {tot / len(per):.2f}")
    for r in range(8):
        print("        " + " ".join(f"{counts.get(8 * r + k, 0):3d}" for k in range(8)))
    fig, axs = plt.subplots(8, 8, figsize=(14, 9), sharex=True, sharey=True)
    for c, ax in enumerate(axs.flat):
        if c in per:
            y, sy, thr, cands = per[c]
            ax.plot(t[region], y[region], color="C0", lw=0.4)
            if cands:
                ax.plot(t[cands], y[cands], "v", color="C3", ms=4)
        ax.text(0.03, 0.82, str(c), transform=ax.transAxes, fontsize=7)
        ax.tick_params(labelsize=6)
        ax.grid(lw=0.3)
    return st.show(4, fig)


def step5(ev, st, out, args):
    all_channels(ev, out, args)
    per, region, dt = out["per"], out["region"], ev["dt"]
    TR = (region.stop - region.start) * dt * 1e-9
    rate = {c: len(v[3]) / TR for c, v in per.items()}
    vals = np.array(list(rate.values()))
    m, s = vals.mean(), vals.std()
    far = [c for c, r in rate.items() if abs(r - m) > 3 * s] if s else []
    print(f"step 5: search time {TR * 1e6:.1f} us per channel; rate mean {m / 1e3:.1f} kHz, "
          f"spread {s / 1e3:.1f} kHz; channels beyond 3 sigma: {far}")
    fig = grid64({c: r / 1e3 for c, r in rate.items()}, title="dark rate (kHz), this event")
    return st.show(5, fig)


def step6(ev, st, out, args, acc):
    all_channels(ev, out, args)
    for c, (y, sy, thr, cands) in out["per"].items():
        acc["A"][c].extend(float(y[i]) for i in cands)
        acc["thr"][c].append(thr)
        acc["sy"][c].append(sy)
    acc["events"] += 1
    c0 = args.channel
    live = [c for c in acc["A"] if acc["A"][c]]
    allA = [a for c in live for a in acc["A"][c]]
    thr_all = float(np.mean([np.mean(acc["thr"][c]) for c in live])) if live else 0.0
    sy_all = float(np.mean([np.mean(acc["sy"][c]) for c in live])) if live else 1.0
    fits = {c: fit_peak(acc["A"][c], float(np.mean(acc["thr"][c])), float(np.mean(acc["sy"][c]))) for c in live}
    fit_all = fit_peak(allA, thr_all, sy_all)
    fit_c0 = fits.get(c0) or fit_peak([], 0, 1)
    ok = {c: f for c, f in fits.items() if f["ok"]}
    med = float(np.median([f["A1"] for f in ok.values()])) if ok else np.nan
    g = {c: f["A1"] / med for c, f in ok.items()} if ok else {}

    print(f"step 6: spectra over {acc['events']} events, {len(live)} live channels, threshold about "
          f"{thr_all:.1f} ADC ({args.nsigma:g} sigma, filtered noise {sy_all:.2f})")
    for who, f, vals in ((f"channel {c0}", fit_c0, acc["A"].get(c0, [])), ("all live channels", fit_all, allA)):
        if f["ok"]:
            print(f"        {who:18s}: {f['n']:6d} pulses, A1 {f['A1']:.2f} +- {f['err']:.2f} ADC, width {f['w1']:.2f}, "
                  f"crosstalk ratio {crosstalk_ratio(vals, f['A1']):.3f}")
        else:
            print(f"        {who:18s}: {f['n']:6d} pulses, no fit ({f['why']})")
    print(f"        per channel: n, A1, error, width, g = A1 / median A1 ({med:.2f})")
    rows = []
    for c in sorted(acc["A"]):
        f = fits.get(c)
        if f is None:
            print(f"        {c:3d}: no candidates")
            rows.append([c, 0, "", "", "", ""])
        elif f["ok"]:
            print(f"        {c:3d}: {f['n']:6d}  {f['A1']:6.2f} +- {f['err']:4.2f}  w {f['w1']:5.2f}  g {g[c]:5.3f}")
            rows.append([c, f["n"], f"{f['A1']:.3f}", f"{f['err']:.3f}", f"{f['w1']:.3f}", f"{g[c]:.4f}"])
        else:
            print(f"        {c:3d}: {f['n']:6d}  no fit ({f['why']})"
                  + (f"; mean {f['A1']:.2f}, width {f['w1']:.2f}" if np.isfinite(f["A1"]) else ""))
            rows.append([c, f["n"], "", "", "", ""])
    off = sorted(c for c, v in g.items() if abs(v - 1) > 0.15)
    print(f"        median A1 {med:.2f} ADC over {len(ok)} fitted channels; g more than 15% from 1: {off}")
    acc["rows"] = rows

    fig, axs = plt.subplots(2, 2, figsize=(12, 10))
    for ax, (who, f, vals) in zip(axs[0], ((f"channel {c0}", fit_c0, acc["A"].get(c0, [])),
                                           ("all live channels", fit_all, allA))):
        if f["h"] is not None:
            ax.stairs(f["h"], f["edges"], color="C0", fill=True, alpha=0.6)
            if np.isfinite(f["A1"]) and "lo" in f:
                a = np.linspace(f["lo"], f["hi"], 100)
                ax.plot(a, gauss(a, f["h"].max() if not f["ok"] else
                                 gauss(f["A1"], 1, f["A1"], f["w1"]) * f["h"][int(f["A1"] // ABIN)], f["A1"], f["w1"]),
                        color="C3", lw=1.5, label=f"A1 {f['A1']:.1f}, w {f['w1']:.1f}" + ("" if f["ok"] else " (rejected)"))
                ax.legend(fontsize=8)
        ax.set(xlabel="A (ADC)", ylabel="pulses", title=f"{who}, {f['n']} pulses")
        ax.title.set_fontsize(9)
        ax.grid(lw=0.3)
    if ok:
        grid_into(axs[1, 0], {c: f["A1"] for c, f in ok.items()}, "A1 (ADC)")
        grid_into(axs[1, 1], g, "g = A1 / median")
    return st.show(6, fig)


def grid_into(ax, values, title, fmt="{:.2f}"):
    img = np.full((8, 8), np.nan)
    for c, v in values.items():
        img[c // 8, c % 8] = v
    im = ax.imshow(img, cmap="viridis")
    vmax = np.nanmax(img)
    for c, v in values.items():
        ax.text(c % 8, c // 8, f"{c}\n{fmt.format(v)}", ha="center", va="center", fontsize=6.5,
                color="w" if v < np.nanmin(img) + 0.5 * (vmax - np.nanmin(img)) else "k")
    ax.set(xticks=[], yticks=[], title=title)
    ax.title.set_fontsize(9)
    plt.colorbar(im, ax=ax, fraction=0.046)


# ---- main ------------------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("run", help="run name (128K_37V_th19) or path to a .h5 file")
    ap.add_argument("--data", default=str(DATA), help="directory of the run files")
    ap.add_argument("--events", default="0", help="range a-b or list a,b,c (0-based)")
    ap.add_argument("--channel", type=int, default=0, help="channel shown in steps 1-3 and 6")
    ap.add_argument("--template", metavar="FILE", help="pulse shape (default src/spe_template_RUN.csv)")
    ap.add_argument("--region", default="tail", choices=("tail", "all"), help="search region")
    ap.add_argument("--tail", type=float, default=10.0, help="length of the tail region (us)")
    ap.add_argument("--nsigma", type=float, default=4.5, help="candidate threshold (sigma of the filtered noise)")
    ap.add_argument("--from-step", type=int, default=1, choices=range(1, 7), help="first step shown")
    ap.add_argument("--save", metavar="NAME", help="no prompts; save figures as src/py/plots/NAME_evE_sN.png")
    ap.add_argument("--summary", action="store_true", help="only the final step 6, after the last event")
    ap.add_argument("--csv", metavar="FILE", help="write the per-channel table of step 6")
    ap.add_argument("--dump", metavar="FILE", help="write the candidates of every channel per event to JSON")
    args = ap.parse_args()

    path = find_run(args.run, Path(args.data))
    tpath = Path(args.template) if args.template else SRC / f"spe_template_{path.stem}.csv"
    if not tpath.is_file():
        sys.exit(f"template not found: {tpath} (make it with template.py --write, or give --template)")
    T = read_template(tpath)
    events = parse_events(args.events)
    f = h5py.File(path, "r")
    wf = f["waveforms"]
    nev, nch, ns = wf.shape
    bad = [e for e in events if not 0 <= e < nev]
    if bad:
        sys.exit(f"events out of range 0..{nev - 1}: {bad}")
    channels = f["channels"][:]
    if args.channel not in channels:
        sys.exit(f"channel {args.channel} not in the file")
    dt, pre = float(f.attrs["tsample_ns"]), int(f.attrs["pretrigger"])
    t = (np.arange(ns) - pre) * dt
    region = search_region(args.region, args.tail, ns, dt)
    acc = {"A": {int(c): [] for c in channels}, "thr": {int(c): [] for c in channels},
           "sy": {int(c): [] for c in channels}, "events": 0}

    ans = "next"
    dump = []
    for e in events:
        w = wf[e].astype(np.float64)
        ev = dict(run=path.stem, e=e, w=w, t=t, dt=dt, channels=channels, T=T, tpath=tpath)
        k = int(np.flatnonzero(channels == args.channel)[0])
        out = dict(region=region, b=float(np.median(w[k])))
        out["x"] = w[k] - out["b"]
        st = Stepper(args.save, path.stem, e, names=STEPS)
        last = e == events[-1]
        steps = [lambda: step1(ev, st, out, args), lambda: step2(ev, st, out, args),
                 lambda: step3(ev, st, out, args), lambda: step4(ev, st, out, args),
                 lambda: step5(ev, st, out, args), lambda: step6(ev, st, out, args, acc)]
        for n, step in enumerate(steps, start=1):
            silent = n < args.from_step or (args.summary and not (n == 6 and last))
            st.silent = silent
            if silent:
                with contextlib.redirect_stdout(io.StringIO()):
                    ans = step()
            else:
                ans = step()
            if ans != "next":
                break
        st.close()
        if args.dump and "per" in out:
            dump.append(dict(event=e, channels={str(c): [dict(time=float(t[i]), A=float(y[i])) for i in cands]
                                                for c, (y, sy, thr, cands) in out["per"].items()},
                             sigma={str(c): float(v[1]) for c, v in out["per"].items()}))
        if ans == "quit":
            break
        if not args.summary:
            print()
    if args.dump:
        with open(args.dump, "w") as fo:
            json.dump(dump, fo)
        print(f"dump written to {args.dump} ({len(dump)} events)")
    if args.summary and args.save:
        src = PLOTS / f"{args.save}_ev{events[-1]}_s6.png"
        if src.exists():
            src.rename(PLOTS / f"{args.save}_spe.png")
            print(f"saved {PLOTS / (args.save + '_spe.png')}")
    if args.csv and acc.get("rows"):
        with open(args.csv, "w", newline="") as fo:
            fo.write(f"# SPE amplitudes of {path.stem}: {acc['events']} events, template {tpath.name}, "
                     f"threshold {args.nsigma:g} sigma; g = A1 / median A1 over the fitted channels\n")
            wr = csv.writer(fo)
            wr.writerow(["channel", "n", "A1", "err", "width", "g"])
            wr.writerows(acc["rows"])
        print(f"table written to {args.csv}")


if __name__ == "__main__":
    main()
