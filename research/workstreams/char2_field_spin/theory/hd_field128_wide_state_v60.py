"""Exact local caps for three NEW wider-state, random-global-gamma ensembles.

Width64; A interpolates128..127+s to0..63, C=A^T, s=32/36/40.
No CA=0, involution, translated-feedback rank, or old numerical cap is used.
All local bounds condition on every principal I+gamma C_D A_D being
invertible. Its failure probability is charged once outside a code union.
This provider never generates coefficients, runs a screen, or writes files.
"""
from fractions import Fraction as F
from functools import lru_cache
from hashlib import sha256
from math import comb, log2
from pathlib import Path

import hd_field128_full_rs_state_screen as maps
import hd_field128_return_tail_moments_v57 as tail

HERE = Path(__file__).resolve().parent
NAME = Path(__file__).name
STATES = (32, 36, 40)
FIELD_ORDER, WIDTH = 1 << 128, 64
FROZEN = {
    'hd_field128_full_rs_state_screen.py':
        '742f7c8bfd8fed7eb277bddb396736577701b8ae1f02177f97a183b304e1025d',
    'hd_field128_return_tail_moments_v57.py':
        '8c1ed22d4ccd3a7033bbd384f23e46291785c0c1bafc6f332a277b08aec61bef',
    '../HD_FIELD128_WIDE_STATE.md':
        'f93c13b8808584a2302ed524318b1688cb2967a3aa71fc42374ae84770b91348',
}


def verify_frozen():
    for name, digest in FROZEN.items():
        if sha256((HERE/name).read_bytes()).hexdigest() != digest:
            raise RuntimeError('Frozen pure helper changed: '+name)
    tail.verify_frozen()


def source_pins():
    verify_frozen()
    pins = {}
    for closure in (maps.source_pins(), tail.source_pins(), FROZEN,
                    {NAME: sha256(Path(__file__).read_bytes()).hexdigest()}):
        for name, digest in closure.items():
            if name in pins and pins[name] != digest:
                raise RuntimeError('Conflicting source identity: '+name)
            pins[name] = digest
    verify_frozen()
    return dict(sorted(pins.items()))


def check_source_pins(expected):
    if source_pins() != expected:
        raise RuntimeError('Wider-state source closure changed')


def _state(state):
    if type(state) is not int or state not in STATES:
        raise ValueError('Public state must be32,36,40')
    return state


def _geometry(j, q, n, state):
    # Generic tiny-oracle parameters; the caller supplies MDS maps.
    if (any(type(x) is not int for x in (j, q, n, state))
            or not 1 <= state < n <= 64 or not 0 <= j <= n
            or q < n or q & (q-1)):
        raise ValueError('Invalid exact MDS moment parameters')


def _mark(z):
    if isinstance(z, (float, bool)):
        raise ValueError('An exact rational mark is required')
    z = F(z)
    if not 0 < z < 1:
        raise ValueError('Expected0<z<1')
    return z


def joint_probability(j, c, e, q, n, state):
    """Dimension-only joint-zero upper mass, normalized by original q**j."""
    _geometry(j, q, n, state)
    if (type(c) is not int or type(e) is not int
            or not 0 <= c <= j or not 0 <= e <= n-j):
        raise ValueError('Invalid disjoint zero-set sizes')
    t = j-c
    u = min(state, t)
    b = max(0, min(e, state)+u-state)
    return F((q**(u-b)-1)*q**(t-u), q**j)


@lru_cache(maxsize=256)
def _factorial_bounds(j, q, n, state):
    _geometry(j, q, n, state)
    bounds = [F(0)]*(n+1)
    for c in range(j):
        for e in range(min(state, n-j+1)):
            bounds[c+e] += (comb(j, c)*comb(n-j, e)
                            *joint_probability(j, c, e, q, n, state))
    expected = 1-F(1, q**min(state, j))
    if bounds[0] != expected:
        raise ArithmeticError('Original nonzero-syndrome mass changed')
    return tuple(bounds)


@lru_cache(maxsize=512)
def _moments(j, z, q, n, state, tail_refine):
    """Generic exact formulas; the caller supplies principal-G and MDS ranks."""
    _geometry(j, q, n, state)
    z = _mark(z)
    if type(tail_refine) is not bool:
        raise ValueError('Expected explicit Boolean refinement flag')
    w, t, u = 1-z, z+(1-z)/q, min(state, j)
    quiet = sum((F(comb(j, c), q**(c+min(state, j-c)))
                 *z**(j-c)*w**c for c in range(j+1)), F(0))
    birth = sum((F(comb(j, c), q**c)*z**(j-c)*w**c
                 *(1-F(1, q**min(state, j-c))) for c in range(j)), F(0))
    inactive = n-j
    excess = sum((F(comb(inactive, e), q**e)*z**(inactive-e)*w**e
                  *(1-F(1, q**(state-e)))
                  for e in range(min(state, inactive+1))), F(0))
    live = t**j*excess/(1-F(1, q**state))
    mass = 1-F(1, q**u)
    # Principal invertibility makes active Y_S uniform and preserves its
    # syndrome-zero kernel, proving ret<=birth even when CA!=0.
    ret = min(birth, mass*z**max(0, n-state+1-j), (q**state-1)*live)
    bounds = _factorial_bounds(j, q, n, state)
    if tail_refine:
        joint = tail.evaluate_envelope(bounds, z)['refined']
    else:
        joint = sum((h*z**(n-k)*w**k for k, h in enumerate(bounds)), F(0))
    ret = min(ret, joint)
    if not z**n*mass <= ret or any(v < 0 for v in (quiet, birth, ret, live)):
        raise ArithmeticError('Exact local moment invariant failed')
    return quiet, birth, ret, live


def exact_moments(j, z, state, *, tail_refine=True):
    """Exact (quiet,birth,return,live); return has NOT yet been divided by d."""
    state, z = _state(state), _mark(z)
    _geometry(j, FIELD_ORDER, WIDTH, state)
    verify_frozen()
    try:
        return _moments(j, z, FIELD_ORDER, WIDTH, state, tail_refine)
    finally:
        verify_frozen()


def exact_packet(j, z, state, *, tail_refine=True):
    """Row-major00,01,10,11 before alpha; same atom-cone normalization."""
    state, z = _state(state), _mark(z)
    quiet, birth, ret, live = exact_moments(j, z, state, tail_refine=tail_refine)
    d = FIELD_ORDER**state-1
    density = F(d, (FIELD_ORDER-1)**state)
    return (quiet, density*birth, ret/d,
            min(density*live, z**(WIDTH-state+1)) if j == 0 else density*live)


def local_logs(z, state, *, tail_refine=True):
    """FLOAT logs sourced solely from this provider's exact rational entries."""
    import numpy as np
    rows = []
    for j in range(WIDTH+1):
        row = exact_packet(j, z, state, tail_refine=tail_refine)
        rows.append([log2(v.numerator)-log2(v.denominator) if v else -np.inf
                     for v in row])
    result = np.asarray(rows, dtype=float).reshape((WIDTH+1, 2, 2))
    result.setflags(write=False)
    return result


@lru_cache(maxsize=3)
def fixed_maps(state):
    """GF256-byte maps; no numerical cap is imported from the map helper."""
    state = _state(state)
    verify_frozen()
    try:
        a = maps.lagrange(range(64), range(128, 128+state))
        h = maps.lagrange(range(state, 2*state), range(state))
        return a, tuple(zip(*a)), h
    finally:
        verify_frozen()


def construction(state):
    state = _state(state)
    return f'toeplitz-fullfield-rs64-{state}-A128-C128T-global-gamma16-k20-v60'


def principal_bad_gamma_upper(state):
    state = _state(state)
    return F(state*(1 << WIDTH), FIELD_ORDER-(1 << 64))


def metadata(state):
    state = _state(state)
    return dict(schema='wide-state-exact-local-v60', construction=construction(state),
        field_order=FIELD_ORDER, inner_width=WIDTH, state_dimension=state,
        emission_nodes=[128, 127+state], output_nodes=[0, 63], feedback='A^T',
        refresh=f'D2 H{state} D1; interpolate0..{state-1} to{state}..{2*state-1}',
        refresh_diagonals='fresh independent nonzero full-field entries',
        gamma_distribution='one uniform scalar in GF(2^128) minus GF(2^64)',
        gamma_independent_of_other_setup=True,
        conditional_on=f'I{state}+gamma C_D A_D invertible for every coordinate set D',
        conditional_pointwise_in_every_good_gamma=True,
        principal_bad_gamma_upper=str(principal_bad_gamma_upper(state)),
        setup_risk_accounting='add once outside the full conditional code union',
        setup_risk_included_in_packet=False, CA_zero_assumed=False,
        involution_assumed=False, translated_rank_bounds_used=False,
        original_input_normalization='q**j', return_state_normalization=f'q**{state}-1 once',
        refresh_density_applies_to=['birth', 'live'], initial_state='zero',
        final_state='unrestricted', state_persists_across_regions=True, pair_swaps=True,
        outer_measure='invertible-toeplitz', local_arithmetic='exact positive rational',
        old_numerical_bounds_reused=False, source_sha256=source_pins(),
        directed_certificate=False, all_profiles_replayed=False, whole_code_certificate=False)
