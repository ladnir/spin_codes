"""Exact verification of the BCH(255) quotient used by the spectrum model.

The narrow-sense designed-distance-37 and designed-distance-39 generators differ
by one primitive irreducible factor of degree eight.  Consequently their quotient
is an eight-dimensional cyclic module on which a cyclic shift has one orbit on
all 255 nonzero elements.
"""

from affine_wambach import ALPHA, gf_mul, gf_pow


LENGTH = 255


def cyclotomic_coset(exponent: int, modulus: int = LENGTH) -> tuple[int, ...]:
    """Return the binary cyclotomic coset of ``exponent`` modulo ``modulus``."""
    first = exponent % modulus
    value = first
    result = []
    while value not in result:
        result.append(value)
        value = (2 * value) % modulus
    assert value == first
    return tuple(result)


def defining_cosets(designed_distance: int) -> set[frozenset[int]]:
    """Return distinct root cosets for a primitive narrow-sense BCH code."""
    return {
        frozenset(cyclotomic_coset(exponent))
        for exponent in range(1, designed_distance)
    }


def binary_poly_degree(poly: int) -> int:
    """Degree of a GF(2) polynomial stored with coefficient i in bit i."""
    if poly <= 0:
        raise ValueError("polynomial must be nonzero")
    return poly.bit_length() - 1


def binary_poly_mul(left: int, right: int) -> int:
    """Multiply GF(2) polynomials stored as bit sets."""
    product = 0
    while right:
        if right & 1:
            product ^= left
        left <<= 1
        right >>= 1
    return product


def binary_poly_divmod(dividend: int, divisor: int) -> tuple[int, int]:
    """Divide GF(2) polynomials stored as bit sets."""
    divisor_degree = binary_poly_degree(divisor)
    quotient = 0
    remainder = dividend
    while remainder and binary_poly_degree(remainder) >= divisor_degree:
        shift = binary_poly_degree(remainder) - divisor_degree
        quotient ^= 1 << shift
        remainder ^= divisor << shift
    return quotient, remainder


def minimal_polynomial(exponent: int) -> int:
    """Compute the minimal polynomial of alpha**exponent over GF(2)."""
    # Coefficients are temporarily elements of GF(256), in ascending order.
    coefficients = [1]
    for conjugate_exponent in cyclotomic_coset(exponent):
        root = gf_pow(ALPHA, conjugate_exponent)
        updated = [0] * (len(coefficients) + 1)
        for degree, coefficient in enumerate(coefficients):
            updated[degree] ^= gf_mul(coefficient, root)
            updated[degree + 1] ^= coefficient
        coefficients = updated

    if any(coefficient not in (0, 1) for coefficient in coefficients):
        raise ArithmeticError("Frobenius product did not descend to GF(2)")
    return sum(coefficient << degree for degree, coefficient in enumerate(coefficients))


def generator_polynomial(designed_distance: int) -> int:
    """Construct the primitive narrow-sense BCH generator over GF(2)."""
    generator = 1
    for coset in sorted(defining_cosets(designed_distance), key=min):
        generator = binary_poly_mul(generator, minimal_polynomial(min(coset)))
    return generator


def multiply_by_x_mod(value: int, modulus: int) -> int:
    """Multiply a residue by x modulo a monic binary polynomial."""
    degree = binary_poly_degree(modulus)
    high_bit = (value >> (degree - 1)) & 1
    value = (value << 1) & ((1 << degree) - 1)
    if high_bit:
        value ^= modulus & ((1 << degree) - 1)
    return value


def verify_quotient_algebra() -> dict[str, int]:
    """Verify the P/Q dimensions, generator factor, and shift orbit exactly."""
    p_cosets = defining_cosets(37)
    q_cosets = defining_cosets(39)
    new_coset = frozenset({37, 74, 148, 41, 82, 164, 73, 146})

    assert p_cosets < q_cosets
    assert q_cosets - p_cosets == {new_coset}

    p_generator = generator_polynomial(37)
    q_generator = generator_polynomial(39)
    quotient_factor = minimal_polynomial(37)

    p_generator_degree = binary_poly_degree(p_generator)
    q_generator_degree = binary_poly_degree(q_generator)
    assert LENGTH - p_generator_degree == 131
    assert LENGTH - q_generator_degree == 123
    assert q_generator == binary_poly_mul(p_generator, quotient_factor)
    assert binary_poly_degree(quotient_factor) == 8

    xn_minus_one = (1 << LENGTH) | 1  # minus and plus coincide over GF(2)
    assert binary_poly_divmod(xn_minus_one, p_generator)[1] == 0
    assert binary_poly_divmod(xn_minus_one, q_generator)[1] == 0

    # P/Q is the cyclic module GF(2)[x]/(quotient_factor).  Starting from 1,
    # multiplication by x must visit every nonzero eight-bit residue exactly once.
    seen = set()
    value = 1
    while value not in seen:
        seen.add(value)
        value = multiply_by_x_mod(value, quotient_factor)
    assert value == 1
    assert seen == set(range(1, 1 << 8))

    return {
        "p_dimension": LENGTH - p_generator_degree,
        "q_dimension": LENGTH - q_generator_degree,
        "quotient_dimension": binary_poly_degree(quotient_factor),
        "nonzero_shift_orbit": len(seen),
        "quotient_factor": quotient_factor,
    }


def polynomial_text(poly: int) -> str:
    """Format a binary polynomial for a compact reproducibility report."""
    terms = []
    for degree in range(binary_poly_degree(poly), -1, -1):
        if not ((poly >> degree) & 1):
            continue
        if degree == 0:
            terms.append("1")
        elif degree == 1:
            terms.append("x")
        else:
            terms.append(f"x^{degree}")
    return " + ".join(terms)


def main() -> None:
    result = verify_quotient_algebra()
    print("PASS: BCH(255,131,37) generator degree 124")
    print("PASS: BCH(255,123,39) generator degree 132")
    print("PASS: Q <= P and dim(P/Q) = 8")
    print(f"New quotient factor: {polynomial_text(result['quotient_factor'])}")
    print("PASS: cyclic shift has one orbit of length 255 on P/Q minus zero")


if __name__ == "__main__":
    main()
