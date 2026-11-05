"""Raw five-cone packet caps for the proposed coupled second-order state.

Y=X+AU; (U',V')=(V,M(U+V+CX)), with CA=0 and fresh MDS refresh M.
State order Z,L,R,B,P. Z/L/R/B use normalized point/flat product references;
P is arbitrary mass on(nonzero U,0), NOT flat density. No fractional powers,
regional averaging, state potential, or production gate is implemented.

Generic rational caps support exact tiny tests. Production FLOAT tables
reuse the audited old R->Z/L row, but NEVER its empty inverse-refresh
sharpening for B->B. New supremum bounds deliberately expose their slack.
"""
from fractions import Fraction as F
from functools import lru_cache
from itertools import product
from math import comb, isfinite

import numpy as np


STATE_ORDER = ('Z', 'L', 'R', 'B', 'P')
Z, L, R, B, P = range(5)


@lru_cache(None)
def kernel_probability(occupancy, rank, q):
    """Exact full-support kernel probability of a rank-min(rank,l) MDS map."""
    if occupancy == 0:
        return F(1)
    if occupancy <= rank:
        return F(0)
    count = sum((-1)**(occupancy-j)*comb(occupancy, j)*q**max(j-rank, 0)
                for j in range(occupancy+1))
    return F(count, (q-1)**occupancy)


def _validate(occupancy, z, q, length, rank):
    occupancy = tuple(occupancy)
    if (not occupancy or any(type(x) is not int or not 0 <= x <= length for x in occupancy)
            or any(type(x) is not int for x in (q, length, rank)) or q < 2
            or not 1 <= rank < length <= q or not 0 < z <= 1):
        raise ValueError('Valid MDS block occupancies, field size, dimensions, and mark required')
    return occupancy


@lru_cache(None)
def support_floors(occupancy, length, rank):
    """Deterministic emission floors for arbitrary live U and U=CX!=0.

    In a returning block CY=CX!=0 also forces at least one output. Blocks
    with1..rank inputs have mandatory nonzero feedback; larger blocks may
    cancel. These are necessary constraints only, not an exact enumerator.
    """
    distance = length-rank+1
    arbitrary, returning = [], []
    for active in product((0, 1), repeat=len(occupancy)):
        if not any(active):
            continue
        arbitrary.append(sum(max(0, distance-l) if bit else l
                             for l, bit in zip(occupancy, active)))
        if all((bit == 0 if l == 0 else bit == 1 if l <= rank else True)
               for l, bit in zip(occupancy, active)):
            returning.append(sum(max(1, distance-l) if bit else l
                                 for l, bit in zip(occupancy, active)))
    return min(arbitrary), min(returning) if returning else None


@lru_cache(None)
def _block_data(z, q, length, rank, exact):
    number = F if exact else float
    one, zero = number(1), number(0)
    base = one/q**rank
    distance = length-rank+1
    h = z+(one-z)/(q-1)
    result = []
    for l in range(length+1):
        outside = length-l
        # E_U z^wt(AU outside S)=q^-rank+excess. The expression is
        # positive and exact for uniform U; no cancellation subtraction.
        excess = sum((number(comb(outside, c))*z**(outside-c)*(one-z)**c
                      *(one/q**c-base) for c in range(min(rank, outside+1))), zero)
        result.append((z**l, h**l*z**max(0, distance-l), excess))
    return base, tuple(result)


def _components(occupancy, z, q, length, rank, *, exact):
    number = F if exact else float
    one, zero = number(1), number(0)
    blocks, weight = len(occupancy), sum(occupancy)
    state = blocks*rank
    dimension = q**state-1
    beta = number(dimension)/(q-1)**state
    quiet = one
    base, block = _block_data(z, q, length, rank, exact)
    uniform_total, uniform_nonzero = one, zero
    zero_state, max_live = one, zero
    for l in occupancy:
        quiet *= number(kernel_probability(l, rank, q))
        empty_state_mark, live_state_mark, excess = block[l]
        # Positive recurrence for product(base+excess)-base^blocks.
        uniform_nonzero = uniform_nonzero*base+uniform_total*excess
        uniform_total *= base+excess
        # Sup over any nonempty set of live state blocks. Actual codeword
        # constraints can only decrease this deliberately generic bound.
        max_live = max(max_live*max(empty_state_mark, live_state_mark),
                       zero_state*live_state_mark)
        zero_state *= empty_state_mark
    j = min(one, (z+one/(q-1))**weight*uniform_nonzero/(one-base**blocks))
    general_floor, return_floor = support_floors(occupancy, length, rank)
    # I+AC is invertible since(CA)=0. An information set of its restricted
    # columns gives the additional generic affine-image marked-mass bound.
    k = min(one, z**general_floor, (z+one/(q-1))**weight)
    nz = max(zero, one-quiet)
    if return_floor is None:
        return_mark, h = zero, zero
    else:
        return_mark = z**return_floor
        rho = sum(min(rank, l) for l in occupancy)
        h = min(nz, one/(q-1)**rho)*return_mark
    return dict(weight=weight, dimension=dimension, beta=beta, quiet=quiet,
                nonzero_feedback=nz, uniform_live_emission=j,
                sup_shifted_return=k, sup_point_return=h,
                sup_point_live=min(one, max_live), return_mark=return_mark,
                generic_weight_floor=general_floor, return_weight_floor=return_floor,
                rank=sum(min(rank, l) for l in occupancy))


def _assemble(occupancy, z, q, length, rank, *, exact):
    number = F if exact else float
    value = _components(occupancy, z, q, length, rank, exact=exact)
    result = [number(0)]*25
    beta, dimension = value['beta'], value['dimension']
    emitted, quiet, nz = z**value['weight'], value['quiet'], value['nonzero_feedback']
    result[5*Z+Z] = quiet*emitted
    result[5*Z+L] = beta*nz*emitted
    result[5*L+B] = beta*emitted
    result[5*L+P] = nz*emitted/dimension
    result[5*R+Z] = nz*value['return_mark']/dimension
    result[5*R+L] = beta*value['uniform_live_emission']
    result[5*B+R] = value['sup_shifted_return']/dimension
    result[5*B+B] = beta*value['uniform_live_emission']
    result[5*P+Z] = value['sup_point_return']
    result[5*P+L] = beta*value['sup_point_live']
    return tuple(result), value


def packet_caps(occupancy, z, *, q=256, length=16, rank=4):
    """Exact rational flat25 cap; assumes block-MDS A/C with CA=0.

    Refresh must have every fixed-nonzero-input output atom <=(q-1)^-state.
    This assumption is verified by the established production refresh lemma
    and by uniform nonzero scalar refresh in the q4 scalar-state fixture.
    """
    z = F(z)
    occupancy = _validate(occupancy, z, q, length, rank)
    return _assemble(occupancy, z, q, length, rank, exact=True)[0]


def packet_details(occupancy, z, *, q=256, length=16, rank=4):
    z = F(z)
    occupancy = _validate(occupancy, z, q, length, rank)
    return _assemble(occupancy, z, q, length, rank, exact=True)[1]


def raw_packet_table(z):
    """Production FLOAT5x5 caps keyed by sorted four-block occupancies.

    Uses the existing alpha1 old-state packet row only at R->Z/L, where
    that map is identical. No old bound is relabelled as a new certificate.
    Values are unscaled density/mass coefficients in STATE_ORDER.
    """
    if isinstance(z, bool) or not isfinite(z) or not 0 < z <= 1:
        raise ValueError('Finite mark in(0,1] required')
    from hd_xor_multigroup import nonzero_matrices, profile_arrays
    occupancies, _ = profile_arrays()
    _, _, old_caps = nonzero_matrices(float(z))
    table = {}
    for occupancy, old in zip(occupancies, old_caps):
        key = tuple(map(int, occupancy))
        values, _ = _assemble(key, float(z), 256, 16, 4, exact=False)
        matrix = np.array(values).reshape(5, 5)
        matrix[R, Z] = min(matrix[R, Z], old[1, 0])
        matrix[R, L] = min(matrix[R, L], old[1, 1])
        if np.any(~np.isfinite(matrix)) or np.any(matrix < 0):
            raise ArithmeticError('Invalid raw five-cone cap')
        table[key] = matrix
    return table


def metadata():
    return dict(state_order=list(STATE_ORDER), recurrence='Y=X+AU; (Unew,Vnew)=(V,M(U+V+CX))',
        field_order=256, packet_width=64, state_symbols=32, refresh_symbols=16,
        regional_alpha=1, point_cone_is_arbitrary_mass=True,
        product_cones_assert_actual_independence=False, state_scaling_applied=False,
        new_entries_use_occupancy_only_supremum_bounds=True,
        old_packet_row_reused_only='R->Z,L', inverse_refresh_empty_cap_used_for_B=False,
        parameter_fits=0, directed=False, whole_code_certificate=False)
