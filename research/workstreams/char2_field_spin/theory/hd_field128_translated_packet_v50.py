"""Exact local caps for A128 emission and A64-transpose feedback, state32.

These caps are pointwise in every degree16 gamma satisfying the global
principal-invertibility event. They do not include its failure probability.
Only the return cap changes relative to the generic state32 local moments.
No FLOAT logarithm, fitted bound, coefficient table or selected-weight claim
is imported. This module performs no production run or artifact write.
"""
from fractions import Fraction as F
from functools import lru_cache
from hashlib import sha256
from math import comb
from pathlib import Path

import hd_field128_full_rs_state_screen as generic
import hd_field128_reversed_return_v44 as reversal

HERE = Path(__file__).resolve().parent
NAME = Path(__file__).name
FIELD_ORDER, WIDTH, STATE, GAMMA_DEGREE = 1 << 128, 64, 32, 16
CONSTRUCTION = 'toeplitz-fullfield-rs64-32-A128-C64T-global-gamma16-paired-k20-v50'
FROZEN_SHA256 = {
    'hd_field128_full_rs_state_screen.py':
        '742f7c8bfd8fed7eb277bddb396736577701b8ae1f02177f97a183b304e1025d',
    'hd_field128_reversed_return_v44.py':
        '485685eec3829eb01d4821b0c87cfc98cc3fc27b9ecd8328f4a827c05444341f',
    '../HD_FIELD128_TRANSLATED_FEEDBACK_RETURN.md':
        '4ad4619a30f7b1d8f9e27276fbb76b303e6e8ae3c414e23690e269ace46526e5',
}


def verify_frozen():
    """Byte-check the fixed mathematical providers and translated-map proof."""
    for name, expected in FROZEN_SHA256.items():
        if sha256((HERE/name).read_bytes()).hexdigest() != expected:
            raise RuntimeError('Frozen local-cap dependency changed: '+name)


def source_pins():
    """Complete frozen import closures plus this new provider and its proof."""
    verify_frozen()
    pins = {}
    for closure in (generic.source_pins(), reversal.source_pins(), FROZEN_SHA256,
                    {NAME: sha256(Path(__file__).read_bytes()).hexdigest()}):
        for name, digest in closure.items():
            if name in pins and pins[name] != digest:
                raise RuntimeError('Conflicting source pin: '+name)
            pins[name] = digest
    verify_frozen()
    return dict(sorted(pins.items()))


def check_source_pins(expected):
    """The future directed bank must call this at entry and every phase exit."""
    if source_pins() != expected:
        raise RuntimeError('Translated packet source closure changed')


def _geometry(j, c=0, e=0):
    if (any(type(x) is not int for x in (j, c, e))
            or not 0 <= c <= j <= WIDTH or not 0 <= e <= WIDTH-j):
        raise ValueError('Expected original occupancy and disjoint zero-set sizes')


def _mark(z):
    if isinstance(z, (float, bool)):
        raise ValueError('An exact rational mark is required, not binary64')
    z = F(z)
    if not 0 < z < 1:
        raise ValueError('Expected an exact mark strictly between zero and one')
    return z


def probability_from_ranks(j, c, u, b, q):
    """Count nonzero syndromes with a proved residual rank lower bound b.

    There are j-c free input coordinates, syndrome rank u, and at most
    q**(u-b)-1 admissible nonzero syndromes. Each has q**(j-c-u)
    preimages. The c active coordinates are forced, not newly sampled.
    Divide by the ORIGINAL q**j. This counting interface alone makes no
    assertion that a supplied rank bound is valid for a construction.
    """
    if (any(type(x) is not int for x in (j, c, u, b, q))
            or not 0 <= c <= j or not 0 <= b <= u <= j-c
            or q < 2 or q & (q-1)):
        raise ValueError('Invalid field cardinality or exact rank/count data')
    return F((q**(u-b)-1)*q**(j-c-u), q**j)


def rank_floor(j, c, e):
    """Frozen V44 occupancy-only bound, only at width64/state32/degree16."""
    _geometry(j, c, e)
    return reversal.rank_floor(j, c, e, state=STATE, degree=GAMMA_DEGREE)


def return_probability(j, c, e):
    """Conditional upper probability for active/inactive zeros and CX != 0."""
    _geometry(j, c, e)
    u = min(STATE, j-c)
    return probability_from_ranks(j, c, u, rank_floor(j, c, e), FIELD_ORDER)


@lru_cache(maxsize=65)
def _return_coefficients(j):
    # Group only positive exact terms of equal z degree. These coefficients
    # are independent of the mark, so repeated z witnesses share the work.
    coefficients = [F(0)]*(WIDTH+1)
    for c in range(j):
        active_count = comb(j, c)
        for e in range(min(STATE, WIDTH-j+1)):
            probability = return_probability(j, c, e)
            if probability:
                coefficients[c+e] += active_count*comb(WIDTH-j, e)*probability
    return tuple(coefficients)


def _return_cap(j, z):
    # Horner is inappropriate in the (z,1-z) Bernstein basis: keep this
    # positive sum. No cancellation or approximate subtraction is used.
    zp, wp, w = [F(1)], [F(1)], 1-z
    for unused in range(WIDTH):
        zp.append(zp[-1]*z)
        wp.append(wp[-1]*w)
    return sum((coefficient*zp[WIDTH-k]*wp[k]
                for k, coefficient in enumerate(_return_coefficients(j))
                if coefficient), F(0))


def return_cap(j, z):
    """Exact positive V44 zero-expansion cap before division by q**32-1."""
    _geometry(j)
    z = _mark(z)
    verify_frozen()
    try:
        return _return_cap(j, z)
    finally:
        verify_frozen()


def exact_moments(j, z):
    """Return exact rational (quiet, birth, return_upper, live).

    The generic moments use only the MDS emission/code/kernel properties
    shared by the translated construction. The V44 cap uses its additional
    principal condition and reversal lemma. No setup risk enters a moment.
    """
    _geometry(j)
    z = _mark(z)
    verify_frozen()
    try:
        quiet, birth, basic_return, live = generic.exact_moments(
            j, z, FIELD_ORDER, WIDTH, STATE, False)
        ret = min(basic_return, _return_cap(j, z))
        result = quiet, birth, ret, live
        if any(not isinstance(value, F) or value < 0 for value in result):
            raise ArithmeticError('Expected nonnegative exact local moments')
        return result
    finally:
        verify_frozen()


def exact_packet(j, z):
    """Row-major (00,01,10,11) caps, before any alpha or packet pairing.

    The refresh density multiplies birth/live, never return. Return is
    divided by the uniform-nonzero-state cardinality exactly once.
    """
    z = _mark(z)
    quiet, birth, ret, live = exact_moments(j, z)
    d = FIELD_ORDER**STATE-1
    density = F(d, (FIELD_ORDER-1)**STATE)
    return (quiet, density*birth, ret/d,
            min(density*live, z**(WIDTH-STATE+1)) if j == 0 else density*live)


def principal_bad_gamma_upper():
    """One unconditional setup risk for gamma sampled once from Gamma16."""
    return F(8*comb(WIDTH, STATE), FIELD_ORDER-(1 << 64))


def metadata():
    return dict(schema='translated-state32-exact-local-packet-v50',
        construction=CONSTRUCTION, field_order=FIELD_ORDER, inner_width=WIDTH,
        state_dimension=STATE, local_arithmetic='exact positive rational',
        emission_nodes=[128, 159], feedback_nodes=[64, 95], output_nodes=[0, 63],
        feedback='gamma times transpose of interpolation from64..95 to0..63',
        refresh='D2 H32 D1; H32 interpolates0..31 to32..63',
        refresh_diagonals='fresh independent nonzero full-field entries',
        CA_zero=True, initial_state='zero', final_state='unrestricted',
        pair_swaps=True, state_persists_across_regions=True,
        outer_measure='invertible-toeplitz',
        gamma_base_field_order=256, gamma_extension_degree=GAMMA_DEGREE,
        gamma_distribution='one uniform scalar in GF(2^128) minus GF(2^64)',
        gamma_independent_of_other_setup=True,
        conditional_on='I32+gamma C_D A_D invertible for every coordinate set D',
        conditional_pointwise_in_every_good_gamma=True,
        principal_bad_gamma_upper=str(principal_bad_gamma_upper()),
        setup_risk_accounting='add once outside the full conditional code union',
        setup_risk_included_in_packet=False, additional_gamma_exclusions=False,
        original_input_normalization='q**j; reversal changes ranks, not the input law',
        return_state_normalization='divide by q**32-1 once',
        refresh_density_applies_to=['birth', 'live'],
        source_sha256=source_pins(), old_numerical_bounds_reused=False,
        distance_threshold_included=False, directed_certificate=False,
        all_profiles_replayed=False, whole_code_certificate=False)
