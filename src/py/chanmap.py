#!/usr/bin/env python3
"""chanmap.py: the 8x8 array of SiPMs on the face of the crystal, with the dead channels and the ringers marked.

    python3 src/py/chanmap.py [--dead LIST] [--ringers LIST] [--save NAME]

One figure. See docs/manual/chanmap_doc.tex.
"""
import argparse

import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

from bgo.io import PLOTS, read_map

DEAD = "19,20,21,22,23,25,27,33,35,41,46,47"      # chanstat_doc.tex, 2026-10-14
RINGERS = "11,12,48,49"                            # spe_doc.tex, 2026-10-14


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--dead", default=DEAD, help="dead channels, comma list")
    ap.add_argument("--ringers", default=RINGERS, help="high-gain channels, comma list")
    ap.add_argument("--save", metavar="NAME", help="save the figure as src/py/plots/NAME_chanmap.png")
    args = ap.parse_args()
    dead = {int(v) for v in args.dead.split(",") if v.strip()}
    ringers = {int(v) for v in args.ringers.split(",") if v.strip()}
    cmap = read_map()

    fig, ax = plt.subplots(figsize=(7, 7))
    for c, (x, y) in cmap.items():
        col = "0.55" if c in dead else ("#f4b183" if c in ringers else "#cfe2f3")
        ax.add_patch(Rectangle((x - 3, y - 3), 6, 6, facecolor=col, edgecolor="k", lw=0.6))
        ax.text(x, y, str(c), ha="center", va="center", fontsize=10, color="w" if c in dead else "k")
    ax.set(xlim=(-25.5, 25.5), ylim=(-25.5, 25.5), xlabel="x (mm)", ylabel="y (mm)", aspect="equal")
    ax.set_title(f"the 64 SiPMs seen from the crystal\ngrey: {len(dead)} dead channels, {100 * len(dead) / 64:.0f}% of the face; "
                 f"orange: {len(ringers)} high-gain channels", fontsize=9)
    ax.grid(lw=0.3)
    fig.tight_layout()
    if args.save:
        PLOTS.mkdir(exist_ok=True)
        out = PLOTS / f"{args.save}_chanmap.png"
        fig.savefig(out, dpi=130)
        print(f"saved {out}")
    else:
        plt.show()


if __name__ == "__main__":
    main()
