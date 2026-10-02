"""Exact expected support counts for independent random-systematic outers.

Each of ``rows`` independent constituents maps u in F_2^n to (u, R u),
where R is sampled uniformly from all binary n-by-n matrices. The setup is
sampled once; counts below average over setup, not over fresh matrices per
message. The systematic half makes every realization injective.

A group comprises one word from each constituent. Its packet support is the
number of nonzero columns of the rows-by-2n matrix of those words. For a
fixed set of h nonzero messages, let a be their systematic union support.
There are

    binom(n,a) * sum_j (-1)^j binom(h,j) (2^(h-j)-1)^a

such message tuples. Inclusion-exclusion enforces that all h rows are
nonzero. For every fixed tuple, the parity halves are independent uniform
n-bit vectors. Their union support has numerator

    binom(n,b) * (2^h-1)^b

and denominator 2^(h*n). Convolution and the binom(rows,h) active-row
choices give the exact expected group counts. No spectrum concentration
or independence between different messages is used.

If each output column subsequently receives an independent uniform GL(rows,2)
map, its label is uniform nonzero conditional on its being active. Labels
are conditionally independent across columns. For rows=4, independent uniform
nonzero GF(16) multipliers give the same fixed-message law. This GF(16)
statement does not apply unchanged to other row counts. A later independent
uniform column shuffle makes support uniform conditional on its size.

These are outer-count identities, not a distance certificate for SPIN.
"""
from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
from math import comb


def _positive_integer(value: int, name: str) -> None:
    if isinstance(value, bool) or not isinstance(value, int) or value < 1:
        raise ValueError(f"{name} must be a positive integer")


def _active_rows(value: int) -> None:
    if isinstance(value, bool) or not isinstance(value, int) or value < 0:
        raise ValueError("active_rows must be a nonnegative integer")


def systematic_union_counts(n: int, active_rows: int) -> tuple[int, ...]:
    """Count ordered nonzero row messages by the size of their union support.

    Row identities are fixed. The sum is (2^n-1)^active_rows. With zero
    active rows, the sole empty tuple has support zero.
    """
    _positive_integer(n, "n")
    _active_rows(active_rows)
    h = active_rows
    result = tuple(
        comb(n, a) * sum(
            (-1) ** j * comb(h, j) * ((1 << (h - j)) - 1) ** a
            for j in range(h + 1)
        )
        for a in range(n + 1)
    )
    # Python's 0**0=1 counts the single empty assignment in inclusion-exclusion.
    assert all(value >= 0 for value in result)
    assert sum(result) == ((1 << n) - 1) ** h
    return result


def active_row_support_numerators(n: int, active_rows: int) -> tuple[int, ...]:
    """Expected support counts for a fixed active-row set, times 2^(h*n).

    The returned integer at index v counts pairs (message tuple, parity tuple)
    whose combined support is v. Every message row is nonzero; parity rows
    may be zero. Division by 2^(active_rows*n) averages its uniform parity
    tuple. It does not divide by the number of messages.
    """
    systematic = systematic_union_counts(n, active_rows)
    parity = tuple(comb(n, b) * ((1 << active_rows) - 1) ** b
                   for b in range(n + 1))
    result = [0] * (2 * n + 1)
    for a, count in enumerate(systematic):
        if not count:
            continue
        for b, choices in enumerate(parity):
            result[a + b] += count * choices
    assert sum(result) == ((1 << n) - 1) ** active_rows * (1 << (active_rows * n))
    return tuple(result)


@dataclass(frozen=True)
class SupportLayer:
    """One active-row class, including all choices of its row identities.

    ``numerators[v] / denominator`` is an expected number of nonzero group
    messages at support v. The denominator is deliberately unreduced.
    """
    active_rows: int
    numerators: tuple[int, ...]
    denominator: int


def group_support_layers(n: int = 128, rows: int = 4) -> tuple[SupportLayer, ...]:
    """Return exact nonzero group counts, separated by active-row count.

    Each h-layer has denominator 2^(h*n). Its sum after division is
    binom(rows,h)*(2^n-1)^h. The all-zero group is excluded.
    """
    _positive_integer(n, "n")
    _positive_integer(rows, "rows")
    return tuple(
        SupportLayer(h,
                     tuple(comb(rows, h) * value
                           for value in active_row_support_numerators(n, h)),
                     1 << (h * n))
        for h in range(1, rows + 1)
    )


def expected_group_support_counts(n: int = 128, rows: int = 4) -> tuple[Fraction, ...]:
    """Exact averaged group enumerator, indexed by packet support 0..2n.

    Fractions are counts, not probabilities or upper bounds on individual
    sampled spectra. Their sum is 2^(rows*n)-1. Independent setup between
    groups permits products of these averaged measures in a first moment.
    """
    layers = group_support_layers(n, rows)
    result = tuple(sum((Fraction(layer.numerators[v], layer.denominator)
                        for layer in layers), Fraction(0))
                   for v in range(2 * n + 1))
    assert result[0] == 0
    assert sum(result) == (1 << (rows * n)) - 1
    return result


def expected_row_weight_counts(n: int = 128) -> tuple[Fraction, ...]:
    """Independent one-row identity: E[A_w]=(C(2n,w)-C(n,w))/2^n.

    C(n,w) is zero when w>n. Subtracting it excludes the zero systematic
    message and every parity-only vector, which this encoder cannot produce.
    """
    _positive_integer(n, "n")
    return tuple(Fraction(comb(2 * n, w) - (comb(n, w) if w <= n else 0), 1 << n)
                 for w in range(2 * n + 1))
