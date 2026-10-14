"""bgo.fit: the single-photoelectron peak fit (docs/manual/spe_doc.tex, step 6), the photopeak fit
and the decay fit (docs/manual/offline_doc.tex)."""
import numpy as np
from scipy.optimize import curve_fit
from scipy.special import erfc

ABIN, AMAX = 1.0, 80.0             # histogram of the filtered amplitude A (ADC)


def gauss(a, n, mu, w):
    return n * np.exp(-0.5 * ((a - mu) / w) ** 2)


def fit_peak(A, thr, sy):
    """Gaussian fit of the single-photoelectron peak of the amplitudes A above the threshold,
    from one sigma above the threshold to 1.6 times the mode (docs/manual/spe_doc.tex, step 6).
    Returns dict(n, A1, err, w1, ok, why, h, edges, lo, hi)."""
    A = np.asarray(A)
    out = dict(n=len(A), A1=np.nan, err=np.nan, w1=np.nan, ok=False, why="", h=None, edges=None)
    if len(A) < 30:
        out["why"] = "fewer than 30 pulses"
        return out
    h, edges = np.histogram(A, bins=np.arange(0, AMAX + ABIN, ABIN))
    c = 0.5 * (edges[1:] + edges[:-1])
    # one sigma above the threshold: starting at the cut itself (tried) lets the Gaussian stretch
    # over the truncation and the noise spike of the first bins, and fits fewer channels, worse
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


# ---- photopeak and decay (docs/manual/offline_doc.tex) ---------------------------------
def peak_bg(x, n, mu, s, c0, c1):
    """Gaussian on a Compton step (erfc) and a constant."""
    return n * np.exp(-0.5 * ((x - mu) / s) ** 2) + 0.5 * c0 * erfc((x - mu) / (s * np.sqrt(2))) + c1


def find_photopeak(h, c, frac=0.25, nsm=5):
    """Rightmost local maximum of the smoothed histogram above `frac` of its maximum: (centre, sigma),
    sigma from the width of the smoothed peak at half its height (FWHM / 2.355)."""
    hs = np.convolve(h, np.ones(nsm) / nsm, mode="same")
    i = np.arange(1, len(hs) - 1)
    loc = i[(hs[i] >= hs[i - 1]) & (hs[i] > hs[i + 1]) & (hs[i] > frac * hs.max())]
    k = int(loc[-1]) if len(loc) else int(hs.argmax())
    half = 0.5 * hs[k]
    lo = k
    while lo > 0 and hs[lo] > half:
        lo -= 1
    hi = k
    while hi < len(hs) - 1 and hs[hi] > half:
        hi += 1
    width = (c[hi] - c[lo]) / 2.355
    return c[k], max(width, 2 * (c[1] - c[0]))


def fit_photopeak(values, nbins=200, frac=0.25, nsig=2.5):
    """Histogram of `values` (0 to the 99.5th percentile), photopeak located and fitted twice with
    peak_bg over mu +- nsig s. Returns dict(mu, err, s, fwhm_rel, n, ok, why, h, edges, lo, hi, p)."""
    v = np.asarray(values, float)
    v = v[np.isfinite(v)]
    out = dict(mu=np.nan, err=np.nan, s=np.nan, fwhm_rel=np.nan, n=len(v), ok=False, why="", h=None,
               edges=None, lo=np.nan, hi=np.nan, p=None)
    if len(v) < 100:
        out["why"] = "fewer than 100 events"
        return out
    top = np.percentile(v, 99.5)
    h, edges = np.histogram(v, bins=nbins, range=(0, top))
    c = 0.5 * (edges[1:] + edges[:-1])
    out.update(h=h, edges=edges)
    mu, s = find_photopeak(h, c, frac)
    p = None
    for _ in range(2):
        sel = (c > mu - nsig * s) & (c < mu + nsig * s)
        if sel.sum() < 6:
            out["why"] = "fit range too narrow"
            return out
        k = int(np.clip(np.searchsorted(c, mu), 0, len(h) - 1))
        p0 = (max(h[k], 1.0), mu, s, max(h[sel][0] - h[sel][-1], 0.0), max(h[sel][-1], 0.0))
        try:
            p, cov = curve_fit(peak_bg, c[sel], h[sel], p0=p0,
                               bounds=([0, mu - s, 0.5 * s, 0, 0], [np.inf, mu + s, 2 * s, np.inf, np.inf]),
                               sigma=np.sqrt(np.maximum(h[sel], 1.0)), absolute_sigma=True, maxfev=20000)
        except RuntimeError:
            out["why"] = "fit did not converge"
            return out
        mu, s = p[1], abs(p[2])
    err = float(np.sqrt(cov[1, 1]))
    out.update(mu=float(mu), err=err, s=float(s), fwhm_rel=float(2.355 * s / mu),
               lo=float(mu - nsig * s), hi=float(mu + nsig * s), p=p, ok=np.isfinite(err) and err < 0.2 * mu)
    if not out["ok"]:
        out["why"] = "fit on a bound"
    return out


def decay(t, A, tau, C):
    return A * np.exp(-t / tau) + C


def decay2(t, A1, tau1, A2, tau2, C):
    return A1 * np.exp(-t / tau1) + A2 * np.exp(-t / tau2) + C


def fit_decay2(t, y, t0, t1):
    """Two exponentials A1 exp(-t/tau1) + A2 exp(-t/tau2) + C fitted to y(t) for t0 <= t <= t1,
    tau1 < tau2. Returns dict(tau1, err1, tau2, err2, A1, A2, C, f1, ok, why, sel), with f1 the
    fraction of the light (area) in the fast component."""
    sel = (t >= t0) & (t <= t1)
    out = dict(tau1=np.nan, err1=np.nan, tau2=np.nan, err2=np.nan, A1=np.nan, A2=np.nan, C=np.nan,
               f1=np.nan, ok=False, why="", sel=sel)
    if sel.sum() < 20:
        out["why"] = "fewer than 20 samples in the range"
        return out
    tt, yy = t[sel] - t0, y[sel]
    span = t1 - t0
    try:
        p, cov = curve_fit(decay2, tt, yy, p0=(0.5 * max(yy[0], 1.0), 0.1 * span, 0.5 * max(yy[0], 1.0), 0.5 * span, 0.0),
                           bounds=([0, 1e-3 * span, 0, 1e-3 * span, -np.inf], [np.inf, 10 * span, np.inf, 100 * span, np.inf]),
                           maxfev=50000)
    except RuntimeError:
        out["why"] = "fit did not converge"
        return out
    A1, tau1, A2, tau2, C = p
    e1, e2 = np.sqrt(cov[1, 1]), np.sqrt(cov[3, 3])
    if tau1 > tau2:
        A1, tau1, A2, tau2, e1, e2 = A2, tau2, A1, tau1, e2, e1
    area1, area2 = A1 * tau1, A2 * tau2
    out.update(tau1=float(tau1), err1=float(e1), tau2=float(tau2), err2=float(e2), A1=float(A1), A2=float(A2),
               C=float(C), f1=float(area1 / (area1 + area2)) if area1 + area2 > 0 else np.nan, ok=True)
    return out


def fit_decay(t, y, t0, t1):
    """Single exponential A exp(-t/tau) + C fitted to y(t) for t0 <= t <= t1 (same units as t).
    Returns dict(tau, err, A, C, ok, why, sel)."""
    sel = (t >= t0) & (t <= t1)
    out = dict(tau=np.nan, err=np.nan, A=np.nan, C=np.nan, ok=False, why="", sel=sel)
    if sel.sum() < 10:
        out["why"] = "fewer than 10 samples in the range"
        return out
    tt, yy = t[sel] - t0, y[sel]
    try:
        p, cov = curve_fit(decay, tt, yy, p0=(max(yy[0], 1.0), 0.3 * (t1 - t0), 0.0),
                           bounds=([0, 1e-3 * (t1 - t0), -np.inf], [np.inf, 100 * (t1 - t0), np.inf]), maxfev=20000)
    except RuntimeError:
        out["why"] = "fit did not converge"
        return out
    out.update(tau=float(p[1]), err=float(np.sqrt(cov[1, 1])), A=float(p[0]), C=float(p[2]), ok=True)
    return out
