"""Search-only direct count caps with floating candidate selection.

Every MGF and optional atom bound is rebuilt with directed arithmetic.
Floating scores select one such bound for each count; selection errors may
weaken the result but cannot invalidate an individually checked candidate.
The binomial baseline is unchanged. Final replay uses regional_count's
ordinary implementation, not this optional search helper.
"""
from fractions import Fraction as Q
from numbers import Integral
from pathlib import Path
import sys

import numpy as np
from flint import arb

PARENT = Path(__file__).resolve().parents[1]
_import_path = sys.path[:]
try:
    sys.path.insert(0, str(PARENT))
    sys.path.insert(0, str(PARENT/'gf16_packets'))
    import regional_count as rc
finally:
    # Legacy dependencies prepend further research directories themselves.
    # Do not let importing this optional helper redirect later imports such
    # as packed_mixing's dense_cover to a same-named legacy module.
    sys.path[:] = _import_path
    del _import_path


def _positive(value, what):
    if not value.is_finite() or not value > 0:
        raise ArithmeticError(f'finite positive {what} required')
    return value


def _select(moments, tilts, groups):
    """Return proposal indices only; no floating value enters a bound."""
    counts = np.arange(groups+1, dtype=float)
    best = np.full(groups+1, np.inf)
    choices = np.zeros(groups+1, dtype=np.int64)
    for index, (moment, tilt) in enumerate(zip(moments, tilts)):
        try:
            # The moment itself may overflow float while its log is small.
            intercept, slope = float(moment.log()), float(tilt)
        except (OverflowError, ValueError):
            continue
        if not np.isfinite(intercept) or not np.isfinite(slope):
            continue
        with np.errstate(over='ignore', invalid='ignore'):
            scores = intercept-slope*counts
        improved = np.isfinite(scores) & (scores < best)
        choices[improved] = index
        best[improved] = scores[improved]
    # If every score at a count was nonfinite, candidate zero is still a
    # checked bound. Float failures can affect tightness, never validity.
    return choices


def count_mass_caps(features, active, cell, interval, minimum_groups, groups, cap, witnesses):
    """Mirror regional_count.count_mass_caps, evaluating selected atoms only.

    As in the original interface, the caller supplies a valid baseline
    density-ratio cap. Every supplied MGF witness is validated, even when
    the floating selector never chooses it.
    """
    if not cap > 0:
        raise ValueError('positive baseline ratio cap required')
    cap = _positive(rc.aq(cap) if type(cap) in (int, Q) else arb(cap), 'baseline cap')
    rc.mgf_upper(features, active, cell, interval, minimum_groups, groups,
        {'tilt': '0', 'dual': ['0', '0', '0']})
    bounds = [min(arb(1), rc.up(cap*mass))
        for mass in rc.binomial_interval_upper(groups, cell)]
    moments, tilts, ratios = [], [], []
    for witness in witnesses:
        moment = rc.mgf_upper(features, active, cell, interval,
            minimum_groups, groups, witness).exp()
        enabled = witness.get('tilted_atom', False)
        if type(enabled) is not bool:
            raise ValueError('boolean tilted-atom option required')
        tilt = Q(witness['tilt'])
        if enabled:
            lower = (rc.tilted_variance_lower(features, active, cell, interval,
                minimum_groups, groups, witness)
                if 'tilted_variance_dual' in witness else arb(0))
            moment *= rc.tilted_atom_upper(features, interval[0], groups, tilt, lower)
        moments.append(_positive(moment, 'MGF/atom bound'))
        tilts.append(tilt)
        ratios.append(_positive((-rc.aq(tilt)).exp(), 'count recurrence ratio'))
    if moments:
        choices = _select(moments, tilts, groups)
        if (len(choices) != groups+1
                or any(not isinstance(v, Integral) or isinstance(v, bool)
                       or not 0 <= v < len(moments) for v in choices)):
            raise ValueError('one valid integer candidate index per count required')
        j = 0
        while j <= groups:
            selected = int(choices[j])
            # Reset at every run, using only fresh Arb arithmetic. In
            # particular, neither a floating intercept nor a previous
            # run's accumulated value supplies this starting bound.
            moment = _positive(moments[selected]*(-rc.aq(tilts[selected])*j).exp(),
                'selected starting atom bound')
            ratio = ratios[selected]
            while True:
                bounds[j] = min(bounds[j], rc.up(moment))
                j += 1
                if j > groups or choices[j] != selected:
                    break
                moment *= ratio
    for bound in bounds:
        _positive(bound, 'returned count bound')
    return bounds
