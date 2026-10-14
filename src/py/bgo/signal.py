"""bgo.signal: the waveform definitions shared by the scripts and mirrored by the JScan pass.

The physics of each function is in docs/selection_doc.tex (steps 2 to 9) and docs/spe_doc.tex
(steps 1 to 3); docs/bgo_doc.tex lists them with their signatures.
"""
import numpy as np


def smooth(x, n):
    """Moving average of width n along the last axis (edges use what is available)."""
    if n <= 1:
        return x
    k = np.ones(n)
    num = np.apply_along_axis(lambda v: np.convolve(v, k, mode="same"), -1, x)
    den = np.convolve(np.ones(x.shape[-1]), k, mode="same")
    return num / den


def rms(v):
    return float(np.std(v))


def clipped_rms(v, nsig=3.0):
    s = rms(v)
    return rms(v[np.abs(v - v.mean()) < nsig * s]) if s else 0.0


def baseline_window(where, nbase, pre, ns):
    return slice(ns - nbase, ns) if where == "end" else slice(0, min(nbase, pre))


def pulse_search(S, nsm):
    """Smoothed sum, its maximum, height and rising edge (None if above H/2 from the start)."""
    Sh = smooth(S, nsm)
    imax = int(Sh.argmax())
    H = Sh[imax]
    below = np.flatnonzero(Sh[:imax] < H / 2)
    irise = int(below[-1]) if len(below) else None
    return Sh, imax, H, irise


def windows(imax, a, b, ns):
    """Charge window around imax and a control window of the same length in the pulse-free part."""
    lo, hi = max(0, imax - a), min(ns, imax + b + 1)
    n = hi - lo
    before, after = (0, lo), (hi, ns)
    seg = after if (after[1] - after[0]) >= (before[1] - before[0]) else before
    length = seg[1] - seg[0]
    short = length < n
    if seg is after:
        ctrl = slice(max(seg[0], seg[1] - n), seg[1])
    else:
        ctrl = slice(seg[0], min(seg[1], seg[0] + n))
    return slice(lo, hi), ctrl, short


def interp_cm(m, win, edge):
    """m with the window replaced by a line between the means of `edge` samples on each side.

    Returns (mt, L, R, lwin, rwin); lwin or rwin is None when the window touches that end of
    the record, and the other end's mean is then used for both."""
    lo, hi = win.start, win.stop
    lwin = slice(max(0, lo - edge), lo) if lo > 0 else None
    rwin = slice(hi, min(len(m), hi + edge)) if hi < len(m) else None
    L = m[lwin].mean() if lwin is not None else m[rwin].mean()
    R = m[rwin].mean() if rwin is not None else L
    mt = m.copy()
    mt[lo:hi] = np.linspace(L, R, hi - lo)
    return mt, L, R, lwin, rwin


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


def matched(x, T):
    """Least-squares height of a pulse of shape T at every sample: y_i = sum_j T_j x_{i+j} / sum T^2.
    Returns y aligned with x (edges set to zero)."""
    H = (len(T) - 1) // 2
    y = np.zeros_like(x)
    y[H:len(x) - H] = np.correlate(x, T, mode="valid") / np.sum(T * T)
    return y


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
