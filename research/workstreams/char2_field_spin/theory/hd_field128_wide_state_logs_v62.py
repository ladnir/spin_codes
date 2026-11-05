"""Fast FLOAT evaluation of the SAME V60 packet caps, not new local bounds.

Public local_logs(z, state) accepts states 32, 36, 40 and 0 < float(z) < 1.
The first call per state caches exact rational Bernstein coefficients and
V57's exact clipped factorial moments / tails. Warm calls use positive
log-sums only; they never evaluate an exact full packet at the new mark.
An exact zero coefficient is -inf, never a numerical cutoff.

This is a proposal provider, NOT directed arithmetic or a certificate.
V60's construction and once-global principal-event risk are unchanged.
Preparation checks the pinned closure before and after exact compilation.
For a search, retain source_pins() at entry and check_source_pins() on exit;
warm evaluations deliberately do not hash source files in the hot path.
"""
from dataclasses import dataclass
from fractions import Fraction as F
from functools import lru_cache
from hashlib import sha256
from math import comb, isfinite, log, log1p, log2
from pathlib import Path

import numpy as np
import hd_field128_wide_state_v60 as exact

HERE = Path(__file__).resolve().parent
NAME = Path(__file__).name
PROVIDER_SHA256 = '939ba095547900a6bc234f01426a48073f6097d0cfbfe16e398ddaf7a23813a0'
STATES = (32, 36, 40)
NEG_INF = float('-inf')


def source_pins():
    path = Path(exact.__file__).resolve()
    if path.parent != HERE or sha256(path.read_bytes()).hexdigest() != PROVIDER_SHA256:
        raise RuntimeError('Frozen V60 provider changed')
    pins = exact.source_pins()
    if pins.get(path.name) != PROVIDER_SHA256:
        raise RuntimeError('V60 closure identity mismatch')
    digest = sha256(Path(__file__).read_bytes()).hexdigest()
    if NAME in pins and pins[NAME] != digest:
        raise RuntimeError('Conflicting V62 identity')
    pins[NAME] = digest
    return dict(sorted(pins.items()))


def check_source_pins(expected):
    if source_pins() != expected:
        raise RuntimeError('FLOAT local-cap source closure changed')


def _log_exact(value):
    """Integer logarithms avoid converting a tiny positive rational to zero."""
    if isinstance(value, (float, bool)):
        raise ValueError('Exact nonnegative coefficient required')
    value = F(value)
    if value < 0:
        raise ValueError('Negative coefficient')
    return log2(value.numerator)-log2(value.denominator) if value else NEG_INF


def _mark(z):
    if isinstance(z, (bool, np.bool_)):
        raise ValueError('Boolean is not a mark')
    try:
        z = float(z)
    except (TypeError, ValueError, OverflowError) as error:
        raise ValueError('A finite FLOAT mark is required') from error
    if not isfinite(z) or not 0 < z < 1:
        raise ValueError('Expected finite 0 < z < 1')
    return z


@dataclass(frozen=True)
class _Table:
    n: int
    state: int
    logq: float
    logd: float
    logrho: float
    j: np.ndarray
    k: np.ndarray
    quiet: np.ndarray
    birth: np.ndarray
    live: np.ndarray
    original: np.ndarray
    clipped: np.ndarray
    tails: np.ndarray


def _compile(q, n, state):
    """Private generic geometry supports independent tiny exact controls.

    The public provider fixes q=2^128, n=64 and the three declared states.
    Every coefficient below is rational and independent of the mark z.
    """
    exact._geometry(0, q, n, state)
    shape = (n+1, n+1)
    quiet, birth, live, original, clipped, tails = [
        np.full(shape, NEG_INF) for unused in range(6)]
    for j in range(n+1):
        for c in range(j+1):
            quiet[j, c] = _log_exact(F(comb(j, c), q**(c+min(state, j-c))))
        for c in range(j):
            birth[j, c] = _log_exact(F(comb(j, c), q**c)
                                     *(1-F(1, q**min(state, j-c))))
        inactive = n-j
        for e in range(min(state, inactive+1)):
            live[j, e] = _log_exact(F(comb(inactive, e), q**e)
                                    *(1-F(1, q**(state-e)))/(1-F(1, q**state)))
        bounds = exact._factorial_bounds(j, q, n, state)
        maximum, clipped_exact, tail_exact = exact.tail.factorial_envelope(bounds)
        if len(bounds) != n+1 or bounds[0] != 1-F(1, q**min(state, j)):
            raise ArithmeticError('Original q**j mass changed')
        if len(tail_exact) != maximum+1:
            raise ArithmeticError('Unexpected exact tail geometry')
        original[j] = [_log_exact(x) for x in bounds]
        clipped[j] = [_log_exact(x) for x in clipped_exact]
        tails[j, :maximum+1] = [_log_exact(x) for x in tail_exact]
    j = np.arange(n+1, dtype=float)[:, None]
    k = np.arange(n+1, dtype=float)[None, :]
    for array in (j, k, quiet, birth, live, original, clipped, tails):
        array.setflags(write=False)
    d = q**state-1
    return _Table(n, state, log2(q), _log_exact(d),
                  _log_exact(F(d, (q-1)**state)), j, k,
                  quiet, birth, live, original, clipped, tails)


@lru_cache(maxsize=3)
def _prepared(state):
    # Validation precedes the cache in local_logs; bool must not alias state.
    pins = source_pins()
    try:
        return _compile(exact.FIELD_ORDER, exact.WIDTH, state)
    finally:
        check_source_pins(pins)


def _sum(log_terms):
    return np.logaddexp2.reduce(log_terms, axis=1)


def _evaluate(table, z):
    """Positive log-domain evaluation, with the original V60 minima intact."""
    z = _mark(z)
    lz, lw = log2(z), log1p(-z)/log(2)
    lt = float(np.logaddexp2(lz, lw-table.logq))
    j, k, n = table.j, table.k, table.n
    active_basis = (j-k)*lz+k*lw
    quiet = _sum(table.quiet+active_basis)
    birth = _sum(table.birth+active_basis)
    live = j[:, 0]*lt+_sum(table.live+(n-j-k)*lz+k*lw)

    bernstein_basis = (n-k)*lz+k*lw
    original = _sum(table.original+bernstein_basis)
    clipped = _sum(table.clipped+bernstein_basis)
    tail = np.logaddexp2(n*lz+table.original[:, 0],
                        lw+_sum(table.tails[:, 1:]+(n-k[:, 1:])*lz))
    emission = (table.original[:, 0]
                +np.maximum(0, n-table.state+1-j[:, 0])*lz)
    ret = np.minimum.reduce((birth, emission, table.logd+live,
                             original, clipped, tail))
    rows = np.stack((quiet, table.logrho+birth,
                     ret-table.logd, table.logrho+live), axis=1)
    rows[0, 3] = min(rows[0, 3], (n-table.state+1)*lz)
    if np.any(np.isnan(rows)) or np.any(np.isposinf(rows)):
        raise ArithmeticError('Nonfinite positive log evaluation')
    rows = rows.reshape(n+1, 2, 2)
    rows.setflags(write=False)
    return rows


def local_logs(z, state):
    """Read-only (65,2,2) FLOAT log2 packet caps, before alpha.

    Return is divided by q**state-1 exactly once in logarithmic arithmetic;
    refresh density applies only to birth/live. No setup risk is included.
    """
    if type(state) is not int or state not in STATES:
        raise ValueError('Public state must be 32, 36, or 40')
    z = _mark(z)
    return _evaluate(_prepared(state), z)


def metadata():
    return dict(schema='wide-state-positive-log-proposal-v62', states=list(STATES),
        mathematical_provider_sha256=PROVIDER_SHA256,
        cached_exact_coefficients=True, exact_zeros_preserved=True,
        return_refinement='same V60 basics and V57 original/clipped/tail minimum',
        source_checks='preparation entry/exit; caller pins search entry/exit',
        construction_changed=False, setup_risk_included=False,
        arithmetic='binary64 positive log-sums; not outward rounded',
        directed=False, whole_code_certificate=False, source_sha256=source_pins())
