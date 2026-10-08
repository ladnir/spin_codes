"""Exact small-field model for a prospective characteristic-two SPIN study.

This is a mathematical baseline, not a proposed production construction or a
distance certificate. In particular, the one-state inner has a state-gap
distance obstruction. The GV helpers give asymptotic reference points only.
"""

from dataclasses import dataclass
from fractions import Fraction
import math


@dataclass(frozen=True)
class BinaryField:
    """Polynomial-basis arithmetic; the caller supplies an irreducible modulus."""

    degree: int
    modulus: int

    def __post_init__(self):
        if self.degree < 1 or self.modulus.bit_length() != self.degree + 1:
            raise ValueError("modulus must be monic of the stated degree")
        if not self.modulus & 1:
            raise ValueError("modulus must have nonzero constant term")

    @property
    def order(self):
        return 1 << self.degree

    def _check(self, value):
        if not isinstance(value, int) or not 0 <= value < self.order:
            raise ValueError("field element is outside the polynomial basis")

    def add(self, left, right):
        self._check(left)
        self._check(right)
        return left ^ right

    def mul(self, left, right):
        self._check(left)
        self._check(right)
        result = 0
        while right:
            if right & 1:
                result ^= left
            right >>= 1
            left <<= 1
            if left & self.order:
                left ^= self.modulus
        return result

    def power(self, value, exponent):
        self._check(value)
        if not isinstance(exponent, int) or exponent < 0:
            raise ValueError("exponent must be a nonnegative integer")
        result = 1
        while exponent:
            if exponent & 1:
                result = self.mul(result, value)
            exponent >>= 1
            value = self.mul(value, value)
        return result

    def inverse(self, value):
        self._check(value)
        if value == 0:
            raise ZeroDivisionError("zero has no multiplicative inverse")
        return self.power(value, self.order - 2)

    def div(self, numerator, denominator):
        return self.mul(numerator, self.inverse(denominator))


GF4 = BinaryField(2, 0b111)
GF256 = BinaryField(8, 0x11B)


def marked_support_transition(order, active_input):
    """Return exact (constant, z) coefficients of the two-state transition.

    Rows describe the incoming state (zero, uniform nonzero); columns describe
    the next state in the same order. The marker z records nonzero output.
    The input is fixed before independent nonzero coefficients a,b are sampled.
    With p=1/(order-1), f=p+(1-p)z, the active-input matrix is
    [[0,z],[p*f,(1-p)*f]]. The inactive-input matrix is [[1,0],[0,z]].
    """
    if not isinstance(order, int) or order < 2 or order & (order - 1):
        raise ValueError("order must be a power of two")
    zero = (Fraction(0), Fraction(0))
    one = (Fraction(1), Fraction(0))
    marker = (Fraction(0), Fraction(1))
    if not active_input:
        return ((one, zero), (zero, marker))
    p = Fraction(1, order - 1)
    return ((zero, marker), ((p * p, p * (1 - p)),
                            (p * (1 - p), (1 - p) * (1 - p))))


def _check_vectors(field, values, *coefficients):
    for value in values:
        field._check(value)
    for vector in coefficients:
        if len(vector) != len(values):
            raise ValueError("coefficient vectors must match the input length")
        for value in vector:
            field._check(value)
            if value == 0:
                raise ValueError("all coefficients must be nonzero")


def recursive_encode(field, values, a, b):
    """Compute y_i=x_i+q_i; q_(i+1)=a_i*q_i+b_i*x_i, with q_0=0."""
    _check_vectors(field, values, a, b)
    state = 0
    result = []
    for value, ai, bi in zip(values, a, b):
        result.append(value ^ state)
        state = field.mul(ai, state) ^ field.mul(bi, value)
    return result


def recursive_decode(field, outputs, a, b):
    """Invert the unit lower-triangular recursive map, without a matrix inverse."""
    _check_vectors(field, outputs, a, b)
    state = 0
    result = []
    for output, ai, bi in zip(outputs, a, b):
        value = output ^ state
        result.append(value)
        state = field.mul(ai, state) ^ field.mul(bi, value)
    return result


def recursive_transpose(field, values, a, b):
    """Apply the field-linear transpose, not the binary-basis adjoint."""
    _check_vectors(field, values, a, b)
    state = 0
    result = [0] * len(values)
    for i in range(len(values) - 1, -1, -1):
        result[i] = values[i] ^ field.mul(b[i], state)
        state = values[i] ^ field.mul(a[i], state)
    return result


def independent_scales(field, a, b):
    """Convert recurrence coefficients to delta_i and gamma_i.

    delta_0=1, delta_(i+1)=a_i*delta_i, gamma_i=b_i/delta_(i+1).
    The last gamma is unused by encoding. delta_0 is fixed, not random.
    Uniform independent nonzero a,b induce uniform independent nonzero
    delta_1,...,delta_n and gamma_0,...,gamma_(n-1); delta_n is omitted here.
    """
    _check_vectors(field, [0] * len(a), a, b)
    delta, gamma = [], []
    scale = 1
    for ai, bi in zip(a, b):
        delta.append(scale)
        scale = field.mul(scale, ai)
        gamma.append(field.div(bi, scale))
    return delta, gamma


def scaled_prefix_encode(field, values, delta, gamma):
    """Compute y_i=x_i+delta_i*sum_(j<i) gamma_j*x_j."""
    _check_vectors(field, values, delta, gamma)
    result = []
    prefix = 0
    for value, di, gi in zip(values, delta, gamma):
        result.append(value ^ field.mul(di, prefix))
        prefix ^= field.mul(gi, value)
    return result


def scaled_suffix_transpose(field, values, delta, gamma):
    """Compute z_i+gamma_i*sum_(j>i) delta_j*z_j."""
    _check_vectors(field, values, delta, gamma)
    result = [0] * len(values)
    suffix = 0
    for i in range(len(values) - 1, -1, -1):
        result[i] = values[i] ^ field.mul(gamma[i], suffix)
        suffix ^= field.mul(delta[i], values[i])
    return result


def field_dot(field, left, right):
    if len(left) != len(right):
        raise ValueError("dot-product lengths differ")
    result = 0
    for a, b in zip(left, right):
        result ^= field.mul(a, b)
    return result


def qary_entropy(delta, m):
    """Asymptotic q-ary entropy H_q(delta), where q=2**m; not a certificate."""
    if not isinstance(m, int) or m < 1:
        raise ValueError("m must be a positive integer")
    if not 0 <= delta <= 1:
        raise ValueError("delta must be in [0,1]")
    ln2 = math.log(2.0)
    binary = 0.0
    if delta:
        binary -= delta * math.log(delta) / ln2
    if delta < 1:
        binary -= (1 - delta) * math.log1p(-delta) / ln2
    log_nonzero = m + math.log1p(-math.ldexp(1.0, -m)) / ln2
    return (binary + delta * log_nonzero) / m


def gv_distance(rate, m):
    """Asymptotic GV reference root H_(2**m)(delta)=1-rate."""
    if not 0 <= rate <= 1:
        raise ValueError("rate must be in [0,1]")
    qary_entropy(0.0, m)  # Validate m before using it in the endpoint.
    low, high = 0.0, 1.0 - math.ldexp(1.0, -m)
    if rate == 0:
        return high
    if rate == 1:
        return 0.0
    for _ in range(100):
        middle = (low + high) / 2
        if qary_entropy(middle, m) < 1 - rate:
            low = middle
        else:
            high = middle
    return (low + high) / 2


def gv_rate(delta, m):
    """Asymptotic reference 1-H_q(delta), not a finite-length rate guarantee."""
    qary_entropy(0.0, m)  # Validate m before evaluating the domain endpoint.
    if not 0 <= delta <= 1.0 - math.ldexp(1.0, -m):
        raise ValueError("delta is outside the increasing entropy branch")
    return 1.0 - qary_entropy(delta, m)
