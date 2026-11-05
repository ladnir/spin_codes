"""FLOAT cyclic-run count logs for fixed-bank proposal screening.

For a nonconstant binary cycle of length n, weight r, and c runs, the
exact count is n/c * binom(r-1,c-1) * binom(n-r-1,c-1). Evaluate this
identity with log-factorials instead of constructing large integers.
These values are not directed enclosures. Directed replayers continue
to use the exact integer count from hd_shifted_pairs_count.
"""
from functools import lru_cache
from math import lgamma, log

import numpy as np

from hd_shifted_swap_potential import full_counts


@lru_cache(4)
def _log_factorials(length):
    values = np.asarray([lgamma(k+1)/log(2) for k in range(length+1)])
    values.setflags(write=False)
    return values


def log_run_counts_float(pairs, weight):
    """Return log2 count upper-screen weights indexed by both full counts.

    Impossible run counts receive -inf. The two constant cycles each
    have count one. No table is cached across weights; only the small
    log-factorial vector is retained.
    """
    full = np.asarray(tuple(full_counts(pairs, weight)), dtype=np.int64)
    length = 2*pairs
    by_runs = np.full(weight+1, -np.inf)
    if weight in (0, length):
        by_runs[0] = 0.
    else:
        runs = np.arange(1, min(weight, length-weight)+1)
        fact = _log_factorials(length)
        by_runs[runs] = (np.log2(length/runs)
                        +fact[weight-1]-fact[runs-1]-fact[weight-runs]
                        +fact[length-weight-1]-fact[runs-1]-fact[length-weight-runs])
    result = by_runs[weight-np.add.outer(full, full)]
    result.setflags(write=False)
    return result


def run_log_screen_float(pairs, weight, log_regions, log_ratio, *, repetitions=256):
    """The existing common-potential run sum with FLOAT log-factorial counts."""
    if not isinstance(repetitions, int) or repetitions < 2 or repetitions % 2:
        raise ValueError('Positive even regional repetition count required')
    entries = np.asarray([np.asarray(log_regions[f]).ravel()
                          for f in full_counts(pairs, weight)])
    a, b, c, d = entries.T
    logs = repetitions/2*np.maximum(np.logaddexp2(a, b+log_ratio),
                                    np.logaddexp2(c-log_ratio, d))
    terms = log_run_counts_float(pairs, weight)+np.add.outer(logs, logs)
    return float(np.logaddexp2.reduce(terms.ravel())+max(0., -log_ratio))
