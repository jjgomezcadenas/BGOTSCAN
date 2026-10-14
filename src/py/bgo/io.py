"""bgo.io: paths, run lookup, event lists, and the small CSV inputs."""
import csv
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent          # src/py/bgo
REPO = HERE.parent.parent.parent
DATA = REPO / "2026_10_02_bgo_cryo_tempscan" / "data_cryo"
PLOTS = HERE.parent / "plots"                   # src/py/plots, not in git
PASS = REPO / "pass"                            # output of the JScan pass, not in git
MAP = REPO / "src" / "channel_map.csv"


def parse_events(spec):
    if "-" in spec:
        a, b = spec.split("-")
        return list(range(int(a), int(b) + 1))
    return [int(x) for x in spec.split(",") if x]


def find_run(run, data_dir):
    p = Path(run)
    if p.is_file():
        return p
    for cand in (data_dir / run, data_dir / (run + ".h5")):
        if cand.is_file():
            return cand
    sys.exit(f"run not found: {run} (looked in {data_dir})")


def read_map(path=MAP):
    out = {}
    with open(path) as f:
        for row in csv.DictReader(l for l in f if not l.startswith("#")):
            out[int(row["channel"])] = (float(row["x_mm"]), float(row["y_mm"]))
    return out


def read_template(path):
    T = []
    with open(path) as f:
        for row in csv.DictReader(l for l in f if not l.startswith("#")):
            T.append(float(row["value"]))
    T = np.array(T)
    if len(T) % 2 == 0 or abs(T[len(T) // 2] - 1) > 1e-6:
        sys.exit(f"template {path}: expected an odd number of values with 1 at the centre")
    return T


def find_pass(run, pass_dir=PASS):
    """Path of the pass output of a run: a path as given, or pass_dir/RUN_pass.h5."""
    p = Path(run)
    if p.is_file():
        return p
    cand = pass_dir / f"{p.stem.replace('_pass', '')}_pass.h5"
    if cand.is_file():
        return cand
    sys.exit(f"pass output not found: {run} (looked for {cand})")
