"""Exact small-field audits of cheaper coordinate-scalar distributions.

The entropy figures are floating-point design screens, not certificates.
The optional quiet receipt uses the existing directed quiet-ball count.
"""

from collections import Counter
from fractions import Fraction
from itertools import product
from math import log2

from hd_blockwise_shared_transfer import multiply, power


def evaluation_code(field, length, dimension):
    """Tiny evaluation RS code; exhaustive audit helper, not an encoder."""
    if not 1 <= dimension <= length <= field:
        raise ValueError('Invalid evaluation code parameters')
    words = []
    for coefficients in product(range(field), repeat=dimension):
        word = []
        for node in range(length):
            value = 0
            for coefficient in reversed(coefficients):
                value = multiply(value, node, field) ^ coefficient
            word.append(value)
        words.append(tuple(word))
    return tuple(words)


def randomized_measure(code, law, field):
    """Exact nonzero expected code measure under independent scalar laws."""
    if (not code or any(not 1 <= c < field or p <= 0 for c, p in law.items())
            or sum(law.values()) != 1):
        raise ValueError('Need a normalized nonzero scalar law')
    length = len(code[0])
    answer = Counter()
    for scalars in product(law, repeat=length):
        probability = Fraction(1)
        for scalar in scalars:
            probability *= law[scalar]
        for word in code:
            if any(word):
                image = tuple(multiply(c, x, field) for c, x in zip(scalars, word))
                answer[image] += probability
    return answer


def code_image_law(code, subgroup, field, *, fix_first=False):
    """Exact distribution on entire code images, including the zero word."""
    length = len(code[0])
    answer = Counter()
    count = len(subgroup)**(length-int(fix_first))
    for tail in product(subgroup, repeat=length-int(fix_first)):
        scalars = (1,)+tail if fix_first else tail
        image = tuple(sorted(tuple(multiply(c, x, field)
                                   for c, x in zip(scalars, word)) for word in code))
        answer[image] += Fraction(1, count)
    return answer


def entropy_screen(field=65536, delta=Fraction(2, 5)):
    """Leading ordinary-ball exponent per symbol at outer rate one half."""
    d = float(delta)
    if not 0 < d < 1 or field < 2:
        raise ValueError('Need an interior distance and a field size')
    binary_entropy = -d*log2(d)-(1-d)*log2(1-d)
    ball = binary_entropy+d*log2(field-1)
    return dict(field=field, distance=str(delta), rigorous=False,
                required_log2_scalar_set_size=2*ball,
                required_scalar_set_size=2**(2*ball),
                exponents_per_symbol={str(m):ball-log2(m)/2
                                      for m in (255, 257, 21845, 32768, 65535)},
                gf256_loss_per_active_row=128*log2(65535/255))


def base_quiet_dense_receipt(first=7519, precision=192):
    """GF256* coefficients: a directed quiet tail, not all-code distance."""
    from flint import arb, ctx
    from hd_blockwise_shared_quiet_compose import quiet_ball_log2

    if not 1 <= first <= 8192:
        raise ValueError('Invalid first active-row count')
    previous = ctx.prec
    ctx.prec = precision
    try:
        # quiet_ball_log2 already divides the count by the 255-element orbit.
        ball = quiet_ball_log2(precision)
        charge = 128*arb(255).log()/arb(2).log()
        endpoint = ball-first*charge
        if not endpoint < -40:
            raise ArithmeticError('The requested quiet tail did not close')
        return dict(field=256, extension_field=65536, first=first, last=8192,
                    projective_quiet_ball_log2=str(ball),
                    failure_log2_upper=str(endpoint), directed=True,
                    nonquiet_covered=False, whole_code_certificate=False)
    finally:
        ctx.prec = previous


if __name__ == '__main__':
    import json
    print(json.dumps(entropy_screen(), indent=2))
    print(json.dumps(base_quiet_dense_receipt(), indent=2))
