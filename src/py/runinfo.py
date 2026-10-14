#!/usr/bin/env python3
"""runinfo.py: run conditions, acquisition rate and true trigger rate of the scan.

    python3 src/py/runinfo.py [--data DIR] [--runs A B ...] [--gap MS] [--save NAME] [--no-figure]

Reads only the attributes, timestamp_ns and trigger_id of each run (never the
waveforms) and temperatures.txt. See docs/manual/runinfo_doc.tex.
"""
import argparse
import ast
import re
import sys
from pathlib import Path

import h5py
import numpy as np
import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parent
REPO = HERE.parent.parent
DATA = REPO / "2026_10_02_bgo_cryo_tempscan" / "data_cryo"
PLOTS = HERE / "plots"
WRAP = 2 ** 24          # trigger_id is a 24-bit counter
NAME_RE = re.compile(r"^(\d+)K_(\d+)(?:d(\d+))?V_th(\d+)$")


def parse_name(stem):
    """'185K_39d5V_th19' -> (185, 39.5, 19); None if the name does not match."""
    m = NAME_RE.match(stem)
    if not m:
        return None
    tnom, vint, vdec, th = m.groups()
    bias = float(vint + "." + (vdec or "0"))
    return int(tnom), bias, int(th)


def read_temperatures(path):
    """Laura's temperatures.txt: lists names, T1, T2, T3 -> {name: (T1, T2, T3)}."""
    if not path.is_file():
        return {}
    vals = {}
    for line in path.read_text().splitlines():
        if "=" in line:
            k, v = line.split("=", 1)
            vals[k.strip()] = ast.literal_eval(v.strip())
    try:
        return {n: (t1, t2, t3) for n, t1, t2, t3
                in zip(vals["names"], vals["T1"], vals["T2"], vals["T3"])}
    except KeyError:
        return {}


def run_info(path, gap_ms):
    """Conditions and bookkeeping of one run, as a dict."""
    with h5py.File(path, "r") as f:
        A = dict(f.attrs)
        t = f["timestamp_ns"][:].astype(np.int64)
        tid = f["trigger_id"][:].astype(np.int64)
    n = len(t)
    dur = (t[-1] - t[0]) / 1e9
    ntrig = int((tid[-1] - tid[0]) % WRAP)
    dt = np.diff(t) / 1e6                                  # ms
    gap = dt > gap_ms
    idx = np.flatnonzero(gap)
    sizes = np.diff(np.concatenate(([0], idx + 1, [n])))   # events per burst
    return dict(
        name=path.stem, start=str(A["start_time"]), events=n,
        majority=int(A["itl_majority"]), duration=dur,
        rate_rec=(n - 1) / dur, ntrig=ntrig, rate_trig=ntrig / dur,
        kept=100 * (n - 1) / ntrig if ntrig else float("nan"),
        nbursts=len(sizes), burst=float(np.median(sizes)),
        gap=float(np.median(dt[gap])) if gap.any() else float("nan"),
        inburst=float(np.median(dt[~gap])) if (~gap).any() else float("nan"),
        record=int(A["record_length"]), pre=int(A["pretrigger"]), ts=float(A["tsample_ns"]),
    )


def print_table(rows):
    hdr = (f"{'run':18s} {'start':8s} {'T1':>6s} {'T2':>6s} {'T3':>6s} {'bias':>5s} {'th':>3s} {'maj':>3s} "
           f"{'dur s':>6s} {'rec Hz':>7s} {'trig Hz':>8s} {'kept %':>7s} {'bursts':>6s} {'ev/b':>5s} "
           f"{'gap ms':>7s} {'in-b ms':>7s}")
    print(hdr)
    print("-" * len(hdr))
    for r in rows:
        tt = lambda v: f"{v:6.1f}" if v is not None else "     -"
        print(f"{r['name']:18s} {r['start'][11:19]:8s} {tt(r['T1'])} {tt(r['T2'])} {tt(r['T3'])} "
              f"{r['bias']:5.1f} {r['th']:3d} {r['majority']:3d} {r['duration']:6.1f} {r['rate_rec']:7.1f} "
              f"{r['rate_trig']:8.0f} {r['kept']:7.2f} {r['nbursts']:6d} {r['burst']:5.0f} "
              f"{r['gap']:7.2f} {r['inburst']:7.2f}")
    r = rows[0]
    print(f"\n{r['events']} events per run, {r['record']} samples x {r['ts']:g} ns, "
          f"pre-trigger {r['pre']} samples; majority {r['majority']}; "
          f"recorded rate = (events-1)/duration, trig Hz from the 24-bit trigger counter.")


def groups_by_setpoint(rows):
    """Runs taken at one temperature setpoint: consecutive runs whose nominal T differ by < 5 K."""
    groups, cur = [], [rows[0]]
    for a, b in zip(rows, rows[1:]):
        if abs(b["tnom"] - a["tnom"]) < 5:
            cur.append(b)
        else:
            groups.append(cur)
            cur = [b]
    groups.append(cur)
    return groups


def figure(rows):
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 4.6))
    groups = groups_by_setpoint(rows)
    rec = np.mean([r["rate_rec"] for r in rows])

    # left: trigger rate vs bias, one colour per setpoint group
    for i, g in enumerate(groups):
        g = sorted(g, key=lambda r: r["bias"])
        x = [r["bias"] for r in g]
        y = [r["rate_trig"] for r in g]
        tlo, thi = min(r["tnom"] for r in g), max(r["tnom"] for r in g)
        lab = f"{tlo} K" if len(g) == 1 else f"{tlo}-{thi} K"
        ax1.plot(x, y, "o-" if len(g) > 1 else "o", color=f"C{i}", ms=6, lw=1.2, label=lab)
        for r in g:
            ax1.annotate(f"{r['T3']:.1f}" if r["T3"] is not None else "", (r["bias"], r["rate_trig"]),
                         textcoords="offset points", xytext=(5, 4), fontsize=7, color="0.3")
    ax1.axhline(rec, color="0.5", lw=1, ls="--")
    ax1.text(ax1.get_xlim()[0], rec * 1.15, f"recorded: {rec:.0f} Hz (readout limit)",
             fontsize=8, color="0.4", va="bottom")
    ax1.set(yscale="log", xlabel="SiPM bias (V)", ylabel="board trigger rate (Hz)",
            title="trigger rate vs bias; labels: T3 (K)")
    ax1.grid(lw=0.3, which="both")
    ax1.legend(fontsize=8, title="setpoint", title_fontsize=8)

    # right: trigger rate vs T3, labels bias; the rising-bias sequence joined
    xs = [r["T3"] for r in rows if r["T3"] is not None]
    ys = [r["rate_trig"] for r in rows if r["T3"] is not None]
    ax2.plot(xs, ys, "o", color="0.3", ms=6)
    for r in rows:
        if r["T3"] is not None:
            ax2.annotate(f"{r['bias']:g} V", (r["T3"], r["rate_trig"]),
                         textcoords="offset points", xytext=(5, 4), fontsize=7, color="0.3")
    seq = [max(g, key=lambda r: r["bias"]) for g in groups]       # highest bias of each setpoint
    ax2.plot([r["T3"] for r in seq], [r["rate_trig"] for r in seq], "-", color="C0", lw=1.2,
             label="highest bias at each setpoint")
    ax2.axhline(rec, color="0.5", lw=1, ls="--")
    ax2.set(yscale="log", xlabel="T3, crystal (K)", ylabel="board trigger rate (Hz)",
            title="trigger rate vs temperature; labels: bias")
    ax2.grid(lw=0.3, which="both")
    ax2.legend(fontsize=8, loc="lower right")
    fig.suptitle("BGO scan 2026-10-02: true trigger rate from the board counter", x=0.01, ha="left")
    fig.tight_layout()
    return fig


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--data", default=str(DATA), help="directory of the run files")
    ap.add_argument("--runs", nargs="*", help="restrict to these runs (default: all *.h5)")
    ap.add_argument("--gap", type=float, default=10.0, help="readout-gap threshold in ms")
    ap.add_argument("--save", metavar="NAME", help="save the figure as src/py/plots/NAME_rates.png")
    ap.add_argument("--no-figure", action="store_true", help="table only")
    args = ap.parse_args()

    data = Path(args.data)
    if args.runs:
        paths = [data / (r if r.endswith(".h5") else r + ".h5") for r in args.runs]
    else:
        paths = sorted(data.glob("*.h5"))
    missing = [p for p in paths if not p.is_file()]
    if missing or not paths:
        sys.exit(f"no run files: {missing or data}")
    temps = read_temperatures(data.parent / "temperatures.txt")

    rows = []
    for p in paths:
        r = run_info(p, args.gap)
        parsed = parse_name(p.stem)
        if parsed is None:
            print(f"warning: {p.name} does not follow TTTK_BBVV_thNN, skipped", file=sys.stderr)
            continue
        r["tnom"], r["bias"], r["th"] = parsed
        r["T1"], r["T2"], r["T3"] = temps.get(p.stem, (None, None, None))
        rows.append(r)
    rows.sort(key=lambda r: r["start"])
    print_table(rows)

    if args.no_figure:
        return
    fig = figure(rows)
    if args.save:
        PLOTS.mkdir(exist_ok=True)
        out = PLOTS / f"{args.save}_rates.png"
        fig.savefig(out, dpi=130)
        print(f"saved {out}")
    else:
        plt.show()


if __name__ == "__main__":
    main()
