#!/usr/bin/env python3
"""make_fixture_real.py: cut a few real events of a run into a small fixture with the same layout.

    python3 src/jscan/test/make_fixture_real.py 197K_39V_th20 --events 0,1,18

Writes src/jscan/test/data/real/<run>.h5 (LZF-compressed, all attributes copied) and the
Python reference dumps next to it, by running selection.py and spe.py on the fixture:
<run>_selection.json and <run>_spe.json. See docs/manual/jscan_doc.tex, "The pass".
"""
import argparse
import os
import subprocess
import sys
from pathlib import Path

import h5py
import numpy as np

HERE = Path(__file__).resolve().parent
REPO = HERE.parent.parent.parent
DATA = REPO / "2026_10_02_bgo_cryo_tempscan" / "data_cryo"
PY = REPO / "src" / "py"


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("run")
    ap.add_argument("--events", default="0,1,18", help="events to copy (0-based, comma list)")
    args = ap.parse_args()
    src = DATA / (args.run + ".h5")
    if not src.is_file():
        sys.exit(f"run not found: {src}")
    events = [int(e) for e in args.events.split(",")]
    (HERE / "data" / "real").mkdir(exist_ok=True)
    out = HERE / "data" / "real" / f"{args.run}.h5"
    with h5py.File(src, "r") as f, h5py.File(out, "w") as g:
        wf = f["waveforms"]
        g.create_dataset("waveforms", data=np.stack([wf[e] for e in events]),
                         chunks=(1, wf.shape[1], wf.shape[2]), compression="lzf")
        g.create_dataset("channels", data=f["channels"][:])
        g.create_dataset("timestamp_ns", data=f["timestamp_ns"][events])
        g.create_dataset("trigger_id", data=f["trigger_id"][events])
        for k, v in f.attrs.items():
            g.attrs[k] = v
        g.attrs["events"] = len(events)
        g.attrs["stop_reason"] = f"fixture: events {events} of {args.run}"
    print(f"wrote {out} ({out.stat().st_size / 1e6:.2f} MB)")
    ev = ",".join(str(i) for i in range(len(events)))
    template = REPO / "src" / "spe_template_195K_41V_th20.csv"
    for script, dump, extra in (("selection.py", f"{args.run}_selection.json", []),
                                ("spe.py", f"{args.run}_spe.json", ["--template", str(template)])):
        # no prompts: selection.py in its --save mode (figures go to src/py/plots), spe.py in --summary
        mode = ["--save", "fixture"] + (["--summary"] if script == "spe.py" else [])
        cmd = [sys.executable, str(PY / script), str(out), "--events", ev, "--dump", str(HERE / "data" / "real" / dump)] + mode + extra
        print(" ".join(cmd))
        subprocess.run(cmd, check=True, stdout=subprocess.DEVNULL, env={"MPLBACKEND": "Agg", **os.environ})
        print(f"wrote {HERE / 'data' / 'real' / dump}")


if __name__ == "__main__":
    main()
