"""Exact binary weight polynomial of the unrandomized RS8 base row.

Polynomial coefficients enumerate the same RS[8,4] code as the public
systematic evaluation interface, since evaluation at 0..3 is bijective.
This is proof-side arithmetic, not an encoder benchmark.
"""

from collections import Counter
from fractions import Fraction
from math import comb


def multiply16(left: int, right: int) -> int:
    result = 0
    for _ in range(4):
        if right & 1:
            result ^= left
        right >>= 1
        left <<= 1
        if left & 16:
            left ^= 0x13
    return result


def row_spectrum() -> list[int]:
    counts: Counter[int] = Counter()
    for message in range(1 << 16):
        coefficients = [(message >> (4 * j)) & 15 for j in range(4)]
        weight = 0
        for x in range(8):
            value = coefficients[3]
            for c in coefficients[2::-1]:
                value = multiply16(value, x) ^ c
            weight += value.bit_count()
        counts[weight] += 1
    result = [counts[w] for w in range(33)]
    assert sum(result) == 1 << 16
    assert result[0] == 1 and min(w for w in range(1, 33) if result[w]) == 6
    assert all(result[w] == result[32 - w] for w in range(33))
    return result


def randomized_rs8_binary_tail(cutoff: int) -> Fraction:
    """Expected number of nonzero RS8/GF32 local words of weight <= cutoff."""
    q = 1 << 32
    ah = {
        h: comb(8, h) * sum(
            (-1) ** j * comb(h, j) * (q ** (h - 4 - j) - 1)
            for j in range(h - 4)
        )
        for h in range(5, 9)
    }
    assert sum(ah.values()) == (1 << 128) - 1
    packet = [Fraction(0)] + [Fraction(comb(32, w), q - 1) for w in range(1, 33)]
    power = [Fraction(1)]
    tail = Fraction(0)
    for h in range(1, 9):
        product = [Fraction(0)] * min(cutoff + 1, len(power) + 32)
        for i, a in enumerate(power):
            if a:
                for j in range(1, min(32, cutoff - i) + 1):
                    product[i + j] += a * packet[j]
        power = product
        if h in ah:
            tail += ah[h] * sum(power)
    return tail


if __name__ == "__main__":
    counts = row_spectrum()
    print("RS[8,4] binary row spectrum:")
    print({w: count for w, count in enumerate(counts) if count})
    print("Eight unrandomized rows: [256,128], binary minimum distance 6")
    tail = randomized_rs8_binary_tail(24)
    print("Randomized RS8 expected nonzero words of weight <= 24:", float(tail))
    assert tail < Fraction(1, 90000)
