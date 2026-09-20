"""Fixed shell measures with logits scaled toward one half.

The measures stay fixed across all regions. An upper envelope accelerates
the existing adaptive band relaxation; it does not choose a new probability
after observing a state. Numerical values are binary64 diagnostics.
"""
import math

import numpy as np

import activation_occupation as general


def density_roots(counts, block, scale, band_count=None):
    if (not counts or not math.isfinite(scale) or scale < 0
            or any(not isinstance(w, int) or not 0 < w <= block
                   or not isinstance(n, int) or n <= 0 for w, n in counts.items())):
        raise ValueError('positive integer spectrum and finite nonnegative scale required')
    bands = general.bands_for(counts, block, band_count or len(counts))
    roots, active, inactive = [], [], []
    for band in bands:
        if band == [block]:
            roots.append(math.log(counts[block])/block)
            active.append(0.); inactive.append(-math.inf)
            continue
        p = (min(band)+max(band))/(2*block)
        eta = scale*(math.log(p)-math.log1p(-p))
        lp = -float(np.logaddexp(0., -eta))
        ln = -float(np.logaddexp(0., eta))
        roots.append(max(math.log(counts[w])-math.log(math.comb(block, w))
                         - w*lp-(block-w)*ln for w in band)/block)
        active.append(lp); inactive.append(ln)
    return bands, np.array(roots), np.array(active), np.array(inactive)


class Envelope:
    """Upper hull of a_g + b_g x for x >= 0, evaluated in log space."""
    def __init__(self, roots, active, inactive):
        roots, active, inactive = map(np.asarray, (roots, active, inactive))
        if not (roots.ndim == 1 and len(roots) and roots.shape == active.shape == inactive.shape):
            raise ValueError('nonempty matching measure arrays required')
        a = np.exp(roots+inactive); b = np.exp(roots+active)
        if not (np.isfinite(a).all() and np.isfinite(b).all()):
            raise ValueError('nonfinite envelope coefficients')
        hull = []; starts = []
        for j in sorted(range(len(a)), key=lambda j: (b[j], a[j])):
            if hull and b[j] == b[hull[-1]]:
                hull.pop(); starts.pop()
            start = -math.inf
            while hull:
                previous = hull[-1]
                start = (a[previous]-a[j])/(b[j]-b[previous])
                if start > starts[-1]: break
                hull.pop(); starts.pop()
            if not hull: start = -math.inf
            hull.append(j); starts.append(start)
        # Keep the line active at zero, followed by positive breakpoints.
        first = max(j for j, start in enumerate(starts) if start <= 0.)
        hull = hull[first:]; starts = starts[first:]
        self.breaks = np.array([math.log(x) for x in starts[1:]])
        self.log_a = (roots+inactive)[hull]
        self.log_b = (roots+active)[hull]

    def apply(self, left, right):
        left, right = np.broadcast_arrays(left, right)
        # Both-zero entries have log ratio zero and remain exactly zero.
        ratio = np.zeros_like(left)
        np.subtract(right, left, out=ratio, where=~(np.isneginf(left)&np.isneginf(right)))
        index = np.searchsorted(self.breaks, ratio)
        return np.logaddexp(self.log_a[index]+left, self.log_b[index]+right)


def occupation_bounds(regions, counts, block, length, cutoff, lam, scale, band_count=None):
    bands, roots, active, inactive = density_roots(counts, block, scale, band_count)
    envelope = Envelope(roots, active, inactive)
    current = regions.copy()
    result = []
    for q in range(1, len(regions)):
        current = envelope.apply(current[:-1], current[1:])
        value = (math.log(math.comb(length, q))+q*math.log(len(bands))
                 +cutoff*lam+general.terminal_log(current[0], block))
        trivial = math.log(math.comb(length, q))+q*math.log(sum(counts.values()))
        result.append(min(value, trivial))
    return np.array(result)
