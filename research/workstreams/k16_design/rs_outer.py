"""Exact expected packet-support counts for the proposed RS[8,4] outer.

Four parallel RS[8,4] codes over GF(256) form one binary [256,128] group.
At each of the eight symbol positions, concatenate the four row symbols
into a 32-bit block. Their block-support enumerator is that of an MDS[8,4]
code over an alphabet of size Q=256^4=2^32. Indeed, the base-field generator
matrix remains MDS over a degree-four extension; grouping the four rows in
an extension-field basis gives the asserted code. Systematic versus
evaluation encoding does not change this complete codeword enumerator.

Independently at each symbol position, sample a uniform GL(32,2) map once
at setup. For each fixed nonzero block, its image is uniform among Q-1
nonzero bit strings. Splitting that image into eight four-bit packets gives
the packet-support probability generating function

    P(z)/(Q-1), where P(z)=(1+15*z)^8-1.

Consequently, if A_h counts messages with h active MDS symbols, the exact
expected nonzero group enumerator is sum_h A_h*(P(z)/(Q-1))^h. This is an
expectation over fixed setup maps, not the spectrum of every sampled map.
It uses no independence between different messages. Independent setup
between groups permits products of these expected counts in first moments.

Conditional on the full packet support, active packet labels are independent
uniform nonzero four-bit strings. An additional independent uniform shuffle
of all 64 packet positions makes the support uniform conditional on its
size. These shuffles and the later regional routing are not included in the
counts below. No inner-code or whole-code distance certificate is claimed.
"""
from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
from math import comb


FIELD_SIZE = 1 << 32
RS_LENGTH = 8
RS_DIMENSION = 4
PACKET_BITS = 4
PACKETS_PER_SYMBOL = 8
GROUP_DIMENSION = 128
GROUP_OUTPUT_BITS = 256
REGION_COUNT = 64
MESSAGE_BITS = 1 << 16
GROUP_COUNT = 512
PHYSICAL_T = 64
PHYSICAL_STEPS_PER_REGION = 32
MACRO_T = 128
MACROS_PER_REGION = 16


def _positive_integer(value: int, name: str) -> None:
    if isinstance(value, bool) or not isinstance(value, int) or value < 1:
        raise ValueError(f"{name} must be a positive integer")


def _is_prime_power(q: int) -> bool:
    """Exact trial division; intended for small fields and powers of two."""
    divisor = 2
    while divisor * divisor <= q:
        if q % divisor == 0:
            while q % divisor == 0:
                q //= divisor
            return q == 1
        divisor = 3 if divisor == 2 else divisor + 2
    return True


def mds_symbol_weight_counts(q: int = FIELD_SIZE, n: int = RS_LENGTH,
                             k: int = RS_DIMENSION) -> tuple[int, ...]:
    """MDS codeword counts indexed by symbol weight 0..n, including zero.

    Require a prime-power field size and 1 <= k <= n <= q. The length
    restriction is sufficient for an ordinary evaluation RS construction;
    extended MDS lengths are deliberately outside this interface. For
    d=n-k+1 and w>=d the formula used is

      A_w = C(n,w) sum_{j=0}^{w-d} (-1)^j C(w,j) (q^(w-d+1-j)-1).

    It follows by inclusion-exclusion on nonzero coordinates in a fixed
    support: shortening to any w coordinates has dimension w-d+1 when
    w>=d, and is the zero code otherwise. The zero codeword is counted once
    at weight zero; the sum at positive weights is q^k-1.
    """
    for value, name in ((q, "q"), (n, "n"), (k, "k")):
        _positive_integer(value, name)
    if q < 2 or not _is_prime_power(q):
        raise ValueError("q must be a prime power at least two")
    if not 1 <= k <= n <= q:
        raise ValueError("require 1 <= k <= n <= q for an evaluation RS code")
    d = n - k + 1
    counts = [0] * (n + 1)
    counts[0] = 1
    for w in range(d, n + 1):
        counts[w] = comb(n, w) * sum(
            (-1) ** j * comb(w, j) * (q ** (w - d + 1 - j) - 1)
            for j in range(w - d + 1)
        )
    assert all(value >= 0 for value in counts)
    assert sum(counts) == q ** k
    return tuple(counts)


def packet_support_numerators(packet_bits: int = PACKET_BITS,
                              packets_per_symbol: int = PACKETS_PER_SYMBOL) -> tuple[int, ...]:
    """Count nonzero bit strings by the number of active fixed-width packets.

    Division by 2^(packet_bits*packets_per_symbol)-1 gives the support law
    of one uniform nonzero block. The numerator at zero is zero.
    """
    _positive_integer(packet_bits, "packet_bits")
    _positive_integer(packets_per_symbol, "packets_per_symbol")
    labels = (1 << packet_bits) - 1
    result = (0,) + tuple(comb(packets_per_symbol, w) * labels ** w
                          for w in range(1, packets_per_symbol + 1))
    assert sum(result) == (1 << (packet_bits * packets_per_symbol)) - 1
    return result


def _convolve(left: tuple[int, ...], right: tuple[int, ...]) -> tuple[int, ...]:
    result = [0] * (len(left) + len(right) - 1)
    for a, left_count in enumerate(left):
        if left_count:
            for b, right_count in enumerate(right):
                if right_count:
                    result[a + b] += left_count * right_count
    return tuple(result)


@dataclass(frozen=True)
class SupportLayer:
    """One symbol-weight class, already including its A_h messages.

    ``numerators[v] / denominator`` is the expected number of nonzero
    group messages at packet support v in this class. Numerators include
    all support positions. The denominator (Q-1)^h is not reduced.
    """
    symbol_weight: int
    numerators: tuple[int, ...]
    denominator: int


def group_support_layers(n: int = RS_LENGTH, k: int = RS_DIMENSION,
                         packet_bits: int = PACKET_BITS,
                         packets_per_symbol: int = PACKETS_PER_SYMBOL) -> tuple[SupportLayer, ...]:
    """Exact nonzero expected counts separated by active MDS-symbol count.

    The alphabet has size Q=2^(packet_bits*packets_per_symbol). Each layer
    has n*packets_per_symbol+1 entries and mass A_h. The zero message is
    excluded. Independent uniform GL(packet_bits*packets_per_symbol,2)
    setup maps at all symbol positions are the required transport premise.
    """
    packet_counts = packet_support_numerators(packet_bits, packets_per_symbol)
    q = 1 << (packet_bits * packets_per_symbol)
    symbol_counts = mds_symbol_weight_counts(q, n, k)
    length = n * packets_per_symbol + 1
    power = (1,)
    layers = []
    for h in range(1, n + 1):
        power = _convolve(power, packet_counts)
        if not symbol_counts[h]:
            continue
        numerators = tuple(symbol_counts[h] * value for value in power)
        numerators += (0,) * (length - len(numerators))
        denominator = (q - 1) ** h
        assert sum(numerators) == symbol_counts[h] * denominator
        layers.append(SupportLayer(h, numerators, denominator))
    return tuple(layers)


def expected_group_support_counts(n: int = RS_LENGTH, k: int = RS_DIMENSION,
                                  packet_bits: int = PACKET_BITS,
                                  packets_per_symbol: int = PACKETS_PER_SYMBOL) -> tuple[Fraction, ...]:
    """Exact expected nonzero-message counts indexed by packet-support size.

    Defaults return 65 Fractions indexed 0..64. Entry zero is zero and the
    exact sum is 2^128-1. In general the length is n*packets_per_symbol+1
    and the sum is 2^(k*packet_bits*packets_per_symbol)-1. These are counts,
    not probabilities, deterministic spectrum bounds, or distance bounds.
    """
    layers = group_support_layers(n, k, packet_bits, packets_per_symbol)
    result = tuple(sum((Fraction(layer.numerators[v], layer.denominator)
                        for layer in layers), Fraction(0))
                   for v in range(n * packets_per_symbol + 1))
    assert result[0] == 0
    assert sum(result) == (1 << (k * packet_bits * packets_per_symbol)) - 1
    return result
