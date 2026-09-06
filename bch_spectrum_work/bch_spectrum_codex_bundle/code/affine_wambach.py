"""Reproduce the affine-stabilizer checks for Wambach's explicit minimum words."""

MODULUS = 0x14D  # x^8 + x^6 + x^3 + x^2 + 1
ALPHA = 0x02

P_WORD_OCTAL = """
00000 01000 02000 00400 00000 00000 00000 00040 00014
24404 00214 01002 26602 40020 05000 16443 04446
"""

Q_WORD_OCTAL = """
00400 00000 00200 00000 00000 00001 00200 00000 04200
55102 40005 52220 14001 22104 50304 40010 50403
"""


def gf_mul(a: int, b: int) -> int:
    r = 0
    while b:
        if b & 1:
            r ^= a
        b >>= 1
        a <<= 1
        if a & 0x100:
            a ^= MODULUS
    return r & 0xFF


def gf_pow(a: int, e: int) -> int:
    r = 1
    while e:
        if e & 1:
            r = gf_mul(r, a)
        a = gf_mul(a, a)
        e >>= 1
    return r


def octal_bits(
    text: str,
    *,
    reverse_digits: bool = False,
    reverse_bits_in_digit: bool = False,
) -> list[int]:
    digits = "".join(text.split())
    if len(digits) != 85:
        raise ValueError(f"expected 85 octal digits = 255 bits, got {len(digits)}")
    if reverse_digits:
        digits = digits[::-1]
    bits: list[int] = []
    for c in digits:
        v = int(c, 8)
        b = [(v >> 2) & 1, (v >> 1) & 1, v & 1]
        if reverse_bits_in_digit:
            b.reverse()
        bits.extend(b)
    assert len(bits) == 255
    return bits


def punctured_weight(text: str) -> int:
    return sum(octal_bits(text))


def extended_support(
    text: str,
    beta_exponent: int,
    *,
    reverse_digits: bool = False,
    reverse_bits_in_digit: bool = False,
) -> frozenset[int]:
    """Map cyclic positions i to beta^i and add extension coordinate 0.

    Wambach uses beta=alpha^13 for the [255,131,37] representation and
    beta=alpha^7 for the [255,123,39] representation.

    The original report's printed-octal bit orientation can be represented by
    digit/within-digit reversals.  For the stabilizer claim used in the handoff,
    all four natural choices give the same result: trivial stabilizer.
    """
    bits = octal_bits(
        text,
        reverse_digits=reverse_digits,
        reverse_bits_in_digit=reverse_bits_in_digit,
    )
    beta = gf_pow(ALPHA, beta_exponent)
    x = 1
    support = []
    for bit in bits:
        if bit:
            support.append(x)
        x = gf_mul(x, beta)
    if sum(bits) & 1:
        support.append(0)  # parity extension
    return frozenset(support)


def affine_image(S: frozenset[int], a: int, b: int) -> frozenset[int]:
    return frozenset(gf_mul(a, x) ^ b for x in S)


def affine_stabilizer(S: frozenset[int]) -> list[tuple[int, int]]:
    stab = []
    for a in range(1, 256):
        for b in range(256):
            if affine_image(S, a, b) == S:
                stab.append((a, b))
    return stab


def verify_seed(text: str, beta_exponent: int, expected_punctured_weight: int) -> None:
    assert punctured_weight(text) == expected_punctured_weight
    expected_extended_weight = expected_punctured_weight + (expected_punctured_weight & 1)
    for rd in (False, True):
        for rb in (False, True):
            S = extended_support(
                text,
                beta_exponent,
                reverse_digits=rd,
                reverse_bits_in_digit=rb,
            )
            assert len(S) == expected_extended_weight
            stab = affine_stabilizer(S)
            assert stab == [(1, 0)], (rd, rb, stab)


def main() -> None:
    verify_seed(P_WORD_OCTAL, 13, 37)
    verify_seed(Q_WORD_OCTAL, 7, 39)
    print("Wambach P seed: punctured wt 37, extended wt 38, affine stabilizer 1")
    print("Wambach Q seed: punctured wt 39, extended wt 40, affine stabilizer 1")
    print("Each extended seed therefore has AGL(1,256) orbit size 256*255 = 65280.")


if __name__ == "__main__":
    main()
