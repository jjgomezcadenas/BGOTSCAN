"""bgo.fit: the single-photoelectron peak fit (docs/spe_doc.tex, step 6)."""
import numpy as np
from scipy.optimize import curve_fit

ABIN, AMAX = 1.0, 80.0             # histogram of the filtered amplitude A (ADC)


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
