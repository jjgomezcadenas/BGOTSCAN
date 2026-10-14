"""bgo.plot: the step-by-step display (Stepper) and the 8x8 channel maps."""
import numpy as np
import matplotlib.pyplot as plt

from .io import PLOTS


class Stepper:
    """Shows one figure per step and waits on the terminal, or saves without prompting.

    show() returns "next" (go on), "skip" (abandon this event), "quit", or ("k", value) when
    ask_k is set and a number was typed."""

    def __init__(self, save, run, event, names):
        self.save, self.run, self.event, self.names = save, run, event, names
        self.fig = None
        self.silent = False                 # True for steps before --from-step: compute only
        if not save:
            plt.ion()

    def show(self, step, fig, ask_k=False):
        if self.silent:
            plt.close(fig)
            return "next"
        fig.suptitle(f"{self.run}: event {self.event}, step {step}: {self.names[step - 1]}",
                     x=0.01, ha="left")
        fig.tight_layout()
        if self.save:
            PLOTS.mkdir(exist_ok=True)
            out = PLOTS / f"{self.save}_ev{self.event}_s{step}.png"
            fig.savefig(out, dpi=130)
            print(f"  saved {out}")
            plt.close(fig)
            return "next"
        if self.fig is not None:
            plt.close(self.fig)
        self.fig = fig
        fig.show()
        plt.pause(0.2)
        if ask_k:
            ans = input("  [Enter] next step, <number> redo with that k, n next event, q quit: ").strip().lower()
            try:
                return ("k", float(ans))
            except ValueError:
                pass
        else:
            ans = input("  [Enter] next step, n next event, q quit: ").strip().lower()
        return {"q": "quit", "n": "skip"}.get(ans, "next")

    def close(self):
        if self.fig is not None:
            plt.close(self.fig)
            self.fig = None


def grid64(values, t=None, title=""):
    """8x8 grid in channel order: waveforms (values[c] over t) or one number per channel."""
    if t is None:                                     # map
        fig, ax = plt.subplots(figsize=(6.4, 5.6))
        img = np.full((8, 8), np.nan)
        for c, v in values.items():
            img[c // 8, c % 8] = v
        im = ax.imshow(img, cmap="viridis")
        for c, v in values.items():
            ax.text(c % 8, c // 8, f"{c}\n{v:.1f}", ha="center", va="center", fontsize=7,
                    color="w" if v < np.nanmax(img) / 2 else "k")
        fig.colorbar(im, ax=ax, label=title)
        ax.set(xticks=[], yticks=[])
        return fig
    fig, axs = plt.subplots(8, 8, figsize=(14, 9), sharex=True, sharey=True)
    for c, ax in enumerate(axs.flat):
        if c in values:
            ax.plot(t, values[c], color="C0", lw=0.5)
        ax.axvline(0, color="0.4", ls="--", lw=0.5)
        ax.text(0.03, 0.85, str(c), transform=ax.transAxes, fontsize=7)
        ax.tick_params(labelsize=6)
        ax.grid(lw=0.3)
    return fig


def map8(ax, values, title, fmt="{:.2f}"):
    """8x8 map of one number per channel, in channel order, value printed in each cell."""
    img = np.full((8, 8), np.nan)
    for c, v in values.items():
        img[c // 8, c % 8] = v
    im = ax.imshow(img, cmap="viridis")
    finite = np.isfinite(img)
    vmin, vmax = (np.nanmin(img), np.nanmax(img)) if finite.any() else (0, 1)
    for c, v in values.items():
        ax.text(c % 8, c // 8, f"{c}\n{fmt.format(v)}", ha="center", va="center", fontsize=6.5,
                color="w" if v < vmin + 0.5 * (vmax - vmin) else "k")
    ax.set(xticks=[], yticks=[], title=title)
    ax.title.set_fontsize(9)
    plt.colorbar(im, ax=ax, fraction=0.046)
