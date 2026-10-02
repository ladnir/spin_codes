"""Exact RS[16,8] maps over polynomial-basis GF16, modulo x^4+x+1.

Messages give the polynomial's values at 0..7; codewords evaluate the same
degree-at-most-seven polynomial at 0..15. The parity block is XOR convolution
on the additive group {0,..,7}. It is symmetric and its square is the identity:
evaluating on the other coset twice is translation by 8 twice.

Write the group-algebra generators as g_i=1+y_i. In characteristic two,
y_i^2=0. A superset-zeta transform changes from the g basis to the y basis,
where convolution becomes multiplication using only disjoint subsets. The
fixed parity circuit needs 19 nontrivial GF16 products rather than 64; its
two basis changes need 24 XORs. These are scalar algebraic operation counts,
not measurements or predictions of a SIMD implementation's instruction count.

Binary transpose uses the adjoints of the four-by-four multiplication maps.
Ordinary field multiplication is not generally self-adjoint in this basis.
"""


def multiply(a, b):
    """Multiply two GF16 values represented by four polynomial coefficients."""
    result = 0
    for _ in range(4):
        if b & 1:
            result ^= a
        b >>= 1
        a <<= 1
        if a & 16:
            a ^= 0x13
    return result


def inverse(a):
    if not 0 < a < 16:
        raise ValueError('nonzero GF16 element required')
    result = 1
    for _ in range(14):
        result = multiply(result, a)
    return result


def systematic_generator():
    """Lagrange interpolation on 0..7, evaluated at the ordered points 0..15."""
    matrix = []
    for x in range(16):
        row = []
        for j in range(8):
            numerator, denominator = 1, 1
            for i in range(8):
                if i != j:
                    numerator = multiply(numerator, x ^ i)
                    denominator = multiply(denominator, j ^ i)
            row.append(multiply(numerator, inverse(denominator)))
        matrix.append(tuple(row))
    return tuple(matrix)


GENERATOR = systematic_generator()
PARITY_CONVOLUTION = GENERATOR[8]


def superset_zeta(values):
    result = list(values)
    if len(result) != 8:
        raise ValueError('eight GF16 elements required')
    for bit in (1, 2, 4):
        for mask in range(8):
            if not mask & bit:
                result[mask] ^= result[mask | bit]
    return result


SQUARE_ZERO_COEFFICIENTS = tuple(superset_zeta(PARITY_CONVOLUTION))


def adjoint_multiply(coefficient, value):
    """Adjoint for the ordinary dot product of the four polynomial-basis bits."""
    return sum(((multiply(coefficient, 1 << bit) & value).bit_count() & 1) << bit
               for bit in range(4))


def encode_row(values):
    if len(values) != 8:
        raise ValueError('eight GF16 message symbols required')
    result = []
    for row in GENERATOR:
        value = 0
        for a, b in zip(row, values):
            value ^= multiply(a, b)
        result.append(value)
    return tuple(result)


def parity_factored(values, transpose=False):
    """Apply the fixed parity map or its binary adjoint using the same circuit.

    The field matrix is symmetric. Replacing each multiplication map in this
    linear circuit by its binary adjoint therefore gives its binary transpose;
    this statement would not hold for an arbitrary nonsymmetric field matrix.
    """
    action = adjoint_multiply if transpose else multiply
    transformed = superset_zeta(values)
    result = [0] * 8
    for output in range(8):
        subset = output
        while True:
            coefficient = SQUARE_ZERO_COEFFICIENTS[subset]
            if coefficient == 1:
                result[output] ^= transformed[output ^ subset]
            elif coefficient:
                result[output] ^= action(coefficient, transformed[output ^ subset])
            if subset == 0:
                break
            subset = (subset - 1) & output
    return tuple(superset_zeta(result))


def transpose_row(values):
    if len(values) != 16:
        raise ValueError('sixteen GF16 codeword symbols required')
    parity = parity_factored(values[8:], transpose=True)
    return tuple(a ^ b for a, b in zip(values[:8], parity))


def check_parity():
    """Check every binary basis input and return exact scalar circuit counts."""
    assert all(GENERATOR[i][j] == int(i == j)
               for i in range(8) for j in range(8))
    assert all(GENERATOR[8 + i][j] == PARITY_CONVOLUTION[i ^ j]
               for i in range(8) for j in range(8))
    assert SQUARE_ZERO_COEFFICIENTS[0] == 1
    for bit in range(32):
        values = [0] * 8
        values[bit // 4] = 1 << (bit % 4)
        reference = encode_row(values)[8:]
        assert parity_factored(values) == reference
        assert parity_factored(reference) == tuple(values)
    return dict(field_modulus='0x13', evaluation_points=list(range(16)),
                parity_convolution=PARITY_CONVOLUTION,
                square_zero_coefficients=SQUARE_ZERO_COEFFICIENTS,
                direct_nontrivial_products=sum(c not in (0, 1)
                    for row in GENERATOR[8:] for c in row),
                factored_nontrivial_products=sum(1 << (3 - mask.bit_count())
                    for mask, c in enumerate(SQUARE_ZERO_COEFFICIENTS)
                    if c not in (0, 1)),
                change_of_basis_xors=24, binary_basis_checks=32)
