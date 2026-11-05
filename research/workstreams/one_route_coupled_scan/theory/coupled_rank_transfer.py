"""Exact rational evaluation of a relaxed local rank-transfer bound.

This is not a transition probability or a whole-code certificate. For fixed
input, current-state, output, and next-state support masks I, H, J, Hprime,
the module retains the charged-coefficient factor tau**(-p) from PROOF.md,
but drops the preliminary deterministic rank gain d. Every retained output
mask contributes y**popcount(J). The rank-bound interpretation requires
1 <= tau <= 256; positive tau outside that interval is only formal algebra.

Matrices have current-state masks as rows and next-state masks as columns.
Cycle tuples list coordinates in forward rotation order: state coordinate
cycle[k] moves to cycle[(k+1) % width]. Physical pairs always remain
(0,1), (2,3), ... regardless of that state wiring.

Only 32 symbolic (I,H) rows are cached. Row evaluators avoid constructing a
full matrix; the largest supported matrix has 256 by 256 rational entries.
This module runs no calculations at import and provides no command-line job.
"""

from __future__ import annotations

from fractions import Fraction
from functools import lru_cache
from math import comb
from typing import Iterator, TypeAlias


Rational: TypeAlias = int | str | Fraction
Term: TypeAlias = tuple[int, int, int]
SymbolicRow: TypeAlias = tuple[tuple[int, tuple[Term, ...]], ...]
RationalRow: TypeAlias = tuple[Fraction, ...]
RationalMatrix: TypeAlias = tuple[RationalRow, ...]

SUPPORTED_WIDTHS = (2, 4, 8)
SYMBOLIC_CACHE_ROWS = 32


def _width(width: int) -> int:
    if type(width) is not int or width not in SUPPORTED_WIDTHS:
        raise ValueError("width must be one of 2, 4, or 8")
    return width


def _cycle(width: int, cycle: tuple[int, ...] | None) -> tuple[int, ...]:
    _width(width)
    if cycle is None:
        return tuple(range(width))
    result = tuple(cycle)
    if (
        len(result) != width
        or any(type(coordinate) is not int for coordinate in result)
        or set(result) != set(range(width))
    ):
        raise ValueError("cycle must list each coordinate exactly once")
    return result


def _mask(width: int, mask: int, name: str) -> int:
    if type(mask) is not int or not 0 <= mask < (1 << width):
        raise ValueError(f"{name} must be an integer in [0, 2**width)")
    return mask


def _fraction(value: Rational, name: str, *, positive: bool) -> Fraction:
    if isinstance(value, bool) or not isinstance(value, (int, str, Fraction)):
        raise TypeError(f"{name} must be an int, Fraction, or rational string")
    result = Fraction(value)
    if (positive and result <= 0) or (not positive and result < 0):
        qualifier = "positive" if positive else "nonnegative"
        raise ValueError(f"{name} must be {qualifier}")
    return result


def _subsets(mask: int) -> Iterator[int]:
    subset = mask
    while True:
        yield subset
        if subset == 0:
            return
        subset = (subset - 1) & mask


def _rotate(mask: int, cycle: tuple[int, ...]) -> int:
    result = 0
    for index, source in enumerate(cycle):
        if mask & (1 << source):
            result |= 1 << cycle[(index + 1) % len(cycle)]
    return result


def rotate_state_mask(
    width: int, state_mask: int, cycle: tuple[int, ...] | None = None
) -> int:
    """Apply the specified state-coordinate cycle to a support mask."""
    wiring = _cycle(width, cycle)
    return _rotate(_mask(width, state_mask, "state_mask"), wiring)


def pair_antipodal_cycle(width: int) -> tuple[int, ...]:
    """Put the members of each physical pair half a cycle apart."""
    _width(width)
    return tuple(range(0, width, 2)) + tuple(range(1, width, 2))


@lru_cache(maxsize=SYMBOLIC_CACHE_ROWS)
def _symbolic_transitions(
    width: int, input_mask: int, state_mask: int, cycle: tuple[int, ...]
) -> SymbolicRow:
    rotated = _rotate(state_mask, cycle)
    mandatory_output = input_mask ^ state_mask
    optional_output = input_mask & state_mask
    counts: dict[int, dict[tuple[int, int], int]] = {}

    # I xor H <= J <= I union H.
    for optional in _subsets(optional_output):
        output_mask = mandatory_output | optional
        output_weight = output_mask.bit_count()
        mandatory_next = output_mask ^ rotated
        optional_next = output_mask & rotated

        # J xor R(H) <= Hprime <= J union R(H).
        for extra in _subsets(optional_next):
            next_mask = mandatory_next | extra
            if input_mask == 0 and state_mask != 0 and next_mask == 0:
                # R+D is invertible for the specified widths and alphabet.
                continue
            charged = (state_mask & ~(output_mask & next_mask)).bit_count()
            bucket = counts.setdefault(next_mask, {})
            key = (output_weight, charged)
            bucket[key] = bucket.get(key, 0) + 1

    return tuple(
        (
            next_mask,
            tuple(
                (output_weight, charged, multiplicity)
                for (output_weight, charged), multiplicity in sorted(bucket.items())
            ),
        )
        for next_mask, bucket in sorted(counts.items())
    )


def symbolic_transitions(
    width: int,
    input_mask: int,
    state_mask: int,
    cycle: tuple[int, ...] | None = None,
) -> dict[int, tuple[Term, ...]]:
    """Return Hprime -> (output weight, charged count, multiplicity) terms.

    Multiplicity counts output-support masks with the same exponents. It
    does not count field values or coefficient assignments. A fresh mapping
    protects the immutable bounded cache from caller mutation.
    """
    wiring = _cycle(width, cycle)
    source = _mask(width, input_mask, "input_mask")
    state = _mask(width, state_mask, "state_mask")
    return dict(_symbolic_transitions(width, source, state, wiring))


def clear_symbolic_cache() -> None:
    """Release all cached symbolic rows."""
    _symbolic_transitions.cache_clear()


def lookahead_charge(
    width: int,
    state_mask: int,
    output_mask: int,
    next_mask: int,
    after_next_mask: int,
    cycle: tuple[int, ...] | None = None,
) -> int:
    """Count additional fresh charges from RANK_CONTINUATION.md section 6.

    This count applies only before the final step, with the full route and
    input forms fixed independently of inner coefficients. It is disjoint
    from the ordinary charge and does not add encoder operations.
    """
    wiring = _cycle(width, cycle)
    state = _mask(width, state_mask, "state_mask")
    output = _mask(width, output_mask, "output_mask")
    following = _mask(width, next_mask, "next_mask")
    after = _mask(width, after_next_mask, "after_next_mask")
    return (_rotate(state & output & following, wiring) & ~(following | after)).bit_count()


def lookahead_transition_terms(
    width: int,
    input_mask: int,
    state_mask: int,
    next_mask: int,
    after_next_mask: int,
    cycle: tuple[int, ...] | None = None,
) -> tuple[Term, ...]:
    """Terms for one augmented edge (H,H1)->(H1,H2), not a dense matrix.

    Returned exponents include ordinary plus lookahead charges. Final-step
    weights must instead use the ordinary local matrix with a free terminal
    state. No future masks are used to condition coefficient probabilities.
    """
    wiring = _cycle(width, cycle)
    source = _mask(width, input_mask, "input_mask")
    state = _mask(width, state_mask, "state_mask")
    following = _mask(width, next_mask, "next_mask")
    after = _mask(width, after_next_mask, "after_next_mask")
    if source == 0 and state != 0 and following == 0:
        return ()
    rotated = _rotate(state, wiring)
    counts: dict[tuple[int, int], int] = {}
    for extra in _subsets(source & state):
        output = (source ^ state) | extra
        if ((output ^ rotated) & ~following) or (following & ~(output | rotated)):
            continue
        ordinary = (state & ~(output & following)).bit_count()
        lagged = (_rotate(state & output & following, wiring) & ~(following | after)).bit_count()
        key = (output.bit_count(), ordinary + lagged)
        counts[key] = counts.get(key, 0) + 1
    return tuple((weight, charged, count) for (weight, charged), count in sorted(counts.items()))


def symbolic_cache_info():
    """Return the standard bounded-lru cache statistics."""
    return _symbolic_transitions.cache_info()


def _monomials(width: int, y: Rational, tau: Rational) -> tuple[RationalRow, ...]:
    output_marker = _fraction(y, "y", positive=True)
    rank_marker = _fraction(tau, "tau", positive=True)
    output_powers = [Fraction(1)]
    rank_powers = [Fraction(1)]
    inverse_rank = 1 / rank_marker
    for _ in range(width):
        output_powers.append(output_powers[-1] * output_marker)
        rank_powers.append(rank_powers[-1] * inverse_rank)
    return tuple(
        tuple(output_power * rank_power for rank_power in rank_powers)
        for output_power in output_powers
    )


def _evaluate_row(
    width: int,
    input_mask: int,
    state_mask: int,
    cycle: tuple[int, ...],
    monomials: tuple[RationalRow, ...],
) -> RationalRow:
    result = [Fraction(0)] * (1 << width)
    for next_mask, terms in _symbolic_transitions(width, input_mask, state_mask, cycle):
        result[next_mask] = sum(
            (multiplicity * monomials[out_weight][charged]
             for out_weight, charged, multiplicity in terms),
            Fraction(0),
        )
    return tuple(result)


def local_row(
    width: int,
    input_mask: int,
    state_mask: int,
    y: Rational,
    tau: Rational,
    cycle: tuple[int, ...] | None = None,
) -> RationalRow:
    """Evaluate one fixed-input relaxed rank-transfer row exactly."""
    wiring = _cycle(width, cycle)
    source = _mask(width, input_mask, "input_mask")
    state = _mask(width, state_mask, "state_mask")
    return _evaluate_row(width, source, state, wiring, _monomials(width, y, tau))


def local_matrix(
    width: int,
    input_mask: int,
    y: Rational,
    tau: Rational,
    cycle: tuple[int, ...] | None = None,
) -> RationalMatrix:
    """Evaluate the matrix for one fixed input-support mask.

    The matrix sums support branches and rank weights; its rows need not
    sum to one. The deterministic preliminary rank gain is omitted.
    """
    wiring = _cycle(width, cycle)
    source = _mask(width, input_mask, "input_mask")
    monomials = _monomials(width, y, tau)
    return tuple(
        _evaluate_row(width, source, state, wiring, monomials)
        for state in range(1 << width)
    )


def physical_pair_weights(
    width: int, x0: Rational, x1: Rational, x2: Rational
) -> RationalRow:
    """Evaluate input-mask weights for independent physical pair factors.

    A pair contributes x0 for 00, x1/2 for each of 01 and 10, and x2 for
    11. The singleton halves model the independent fair pair swap. These
    are generating weights, not an assertion that actual pair occupancies
    are independent in a region.
    """
    _width(width)
    pair_weights = (
        _fraction(x0, "x0", positive=False),
        _fraction(x1, "x1", positive=False) / 2,
        _fraction(x2, "x2", positive=False),
    )
    values = []
    for mask in range(1 << width):
        value = Fraction(1)
        for offset in range(0, width, 2):
            pair = (mask >> offset) & 3
            value *= pair_weights[pair.bit_count()]
        values.append(value)
    return tuple(values)


def _paired_row(
    width: int,
    state_mask: int,
    cycle: tuple[int, ...],
    weighted_monomials: tuple[tuple[int, tuple[RationalRow, ...]], ...],
) -> RationalRow:
    result = [Fraction(0)] * (1 << width)
    for input_mask, monomials in weighted_monomials:
        for next_mask, terms in _symbolic_transitions(width, input_mask, state_mask, cycle):
            result[next_mask] += sum(
                (multiplicity * monomials[out_weight][charged]
                 for out_weight, charged, multiplicity in terms),
                Fraction(0),
            )
    return tuple(result)


def _weighted_monomials(
    width: int,
    y: Rational,
    tau: Rational,
    x0: Rational,
    x1: Rational,
    x2: Rational,
) -> tuple[tuple[int, tuple[RationalRow, ...]], ...]:
    monomials = _monomials(width, y, tau)
    input_weights = physical_pair_weights(width, x0, x1, x2)
    # Equal pair-composition weights share one small table. There are at
    # most choose(width/2 + 2, 2) distinct weights before coincidences.
    shared: dict[Fraction, tuple[RationalRow, ...]] = {}
    result = []
    for mask, weight in enumerate(input_weights):
        if weight == 0:
            continue
        if weight not in shared:
            shared[weight] = tuple(
                tuple(weight * value for value in row) for row in monomials
            )
        result.append((mask, shared[weight]))
    return tuple(result)


def paired_route_row(
    width: int,
    state_mask: int,
    y: Rational,
    tau: Rational,
    x0: Rational,
    x1: Rational,
    x2: Rational,
    cycle: tuple[int, ...] | None = None,
) -> RationalRow:
    """Evaluate one row of the physical-pair generating matrix."""
    wiring = _cycle(width, cycle)
    state = _mask(width, state_mask, "state_mask")
    weighted = _weighted_monomials(width, y, tau, x0, x1, x2)
    return _paired_row(width, state, wiring, weighted)


def paired_route_matrix(
    width: int,
    y: Rational,
    tau: Rational,
    x0: Rational,
    x1: Rational,
    x2: Rational,
    cycle: tuple[int, ...] | None = None,
) -> RationalMatrix:
    """Evaluate the local physical-pair generating matrix exactly.

    A complete regional argument must extract actual pair-profile counts
    and normalize their placements. This function does neither. It also
    omits outer shortened dimensions and the global message-rank budget.
    Width eight is a finite but potentially substantial calculation; use
    paired_route_row when a complete matrix is unnecessary.
    """
    wiring = _cycle(width, cycle)
    weighted = _weighted_monomials(width, y, tau, x0, x1, x2)
    return tuple(
        _paired_row(width, state, wiring, weighted)
        for state in range(1 << width)
    )


def route_profile_normalizer(M: int, n0: int, n1: int, n2: int) -> int:
    """Return M!/(n0! n1! n2!) for a specified pair-category profile.

    Here n0+n1+n2=M. This counts category placements, not the ordering of
    distinct pair labels within categories. Singleton orientation weights
    are supplied separately by physical_pair_weights.
    """
    if any(type(value) is not int or value < 0 for value in (M, n0, n1, n2)):
        raise ValueError("M, n0, n1, and n2 must be nonnegative integers")
    if n0 + n1 + n2 != M:
        raise ValueError("pair-profile counts must sum to M")
    return comb(M, n0) * comb(M - n0, n1)
