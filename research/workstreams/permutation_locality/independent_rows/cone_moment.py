"""Propagate valid eleven-coordinate uppers with mature-subset clamping.

Actual tilted measures obey L48<=L56<=M and pointwise density C<=M.
Clamping upper bounds by these inequalities is safe after each region.
The transfer must be a nonnegative valid operator for these coordinates.
"""
from math import log
from pathlib import Path
import sys
import numpy as np
from flint import arb, arb_mat

sys.path.insert(0,str(Path(__file__).resolve().parent.parent))
from occupancy_memory import Z, M, C
from mature_tail import L48, L56, TAIL_TERMINAL
from group_rank_one_verify import up


def clamp(values):
    result = list(values)
    result[L56] = min(result[L56],result[M])
    result[L48] = min(result[L48],result[L56])
    result[C] = min(result[C],result[M])
    return result


def moment(region, steps=256, terminal=TAIL_TERMINAL):
    if region.nrows()!=11 or region.ncols()!=11 or not isinstance(steps,int) or steps<0:
        raise ValueError('an eleven-coordinate operator and nonnegative step count are required')
    region = arb_mat([[up(region[i,j]) for j in range(11)] for i in range(11)])
    if any(region[i,j]<0 or not region[i,j].is_finite() for i in range(11) for j in range(11)):
        raise ValueError('finite nonnegative operator required')
    values = [arb(int(i==Z)) for i in range(11)]
    for _ in range(steps):
        updated = arb_mat([values])*region
        values = clamp([up(updated[0,j]) for j in range(11)])
    return up(sum((values[j] for j,flag in enumerate(terminal) if flag),arb(0)))


def log_moment(region, steps=256, terminal=TAIL_TERMINAL):
    """Scaled binary64 witness objective only, never a certificate."""
    region = np.asarray(region,dtype=float)
    if region.shape!=(11,11) or not np.isfinite(region).all() or (region<0).any():
        raise ValueError('finite nonnegative proposal operator required')
    if not isinstance(steps,int) or steps<0:
        raise ValueError('nonnegative integer step count required')
    values = np.zeros(11)
    values[Z] = 1
    scale = 0.
    for _ in range(steps):
        values = np.asarray(clamp(values@region))
        peak = max(values)
        if peak==0:
            return -np.inf
        values /= peak
        scale += log(peak)
    total = sum(values[j] for j,flag in enumerate(terminal) if flag)
    return scale+log(total) if total else -np.inf
