"""Split type boxes at zero faces and by relative count uncertainty."""
import math

import numpy as np


def split_box(lower, upper, corners):
    lower = corners.min(axis=0).astype(np.int64)
    upper = corners.max(axis=0).astype(np.int64)
    widths = upper-lower
    if not np.any(widths):
        return []
    # Multinomial curvature is strongest near an absent category. Give
    # those coordinates priority over a wider but nearly fixed large count.
    score = widths.astype(float)**2/(lower+.25*widths+1)
    coordinate = int(np.argmax(score))
    lo, hi = int(lower[coordinate]), int(upper[coordinate])
    midpoint = 0 if lo == 0 else min(hi-1, math.isqrt(lo*hi))
    left = upper.copy(); left[coordinate] = midpoint
    right = lower.copy(); right[coordinate] = midpoint+1
    return [(lower.copy(), left), (right, upper.copy())]
