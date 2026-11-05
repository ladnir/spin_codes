"""One-product MDS pair candidate and its restricted positive-atom bound.

Sample alpha uniformly from GF(q) excluding0,1 and set
t=alpha*(a+b), output(a+t,b+t).  The map is an involution and both
singletons expand.  Its double-input law depends on whether a=b; it is
NOT a support-only Markov kernel.  For product markers in[0,1], the
unequal-input law is a universal pointwise majorant.  This helper makes
no multi-layer or native-performance claim.
"""
from fractions import Fraction

from hd_outer_pair_layers import field_multiply


def _parameters(mask, q):
    if type(q) is not int or q < 4 or q & (q-1) or type(mask) is not int or not 0 <= mask < 4:
        raise ValueError('A binary-field order q>=4 and a two-coordinate mask required')


def apply_pair(a, b, alpha, *, modulus=0x11b):
    q = 1 << (modulus.bit_length()-1)
    if not 0 <= a < q or not 0 <= b < q or not 1 < alpha < q:
        raise ValueError('Field inputs and alpha outside{0,1} required')
    shared = field_multiply(alpha, a ^ b, modulus=modulus)
    return a ^ shared, b ^ shared


def support_law(mask, *, equal_nonzero=None, q=256):
    """Exact output law; double inputs REQUIRE the additional equality flag."""
    _parameters(mask, q)
    if mask == 0:
        return (Fraction(1), Fraction(0), Fraction(0), Fraction(0))
    if mask != 3:
        return (Fraction(0), Fraction(0), Fraction(0), Fraction(1))
    if type(equal_nonzero) is not bool:
        raise ValueError('The exact double-input law requires its equality flag')
    if equal_nonzero:
        return (Fraction(0), Fraction(0), Fraction(0), Fraction(1))
    return (Fraction(0), Fraction(1, q-2), Fraction(1, q-2), Fraction(q-4, q-2))


def atom_majorant(mask, first, second, *, q=256):
    """Valid for a single layer and product atoms both in[0,1].

    The difference unequal-minus-equal is
    (first+second-2*first*second)/(q-2), which is nonnegative here.
    This is not coefficientwise domination of probability distributions.
    """
    _parameters(mask, q)
    first, second = Fraction(first), Fraction(second)
    if not 0 <= first <= 1 or not 0 <= second <= 1:
        raise ValueError('The normalized atom majorant requires both markers in[0,1]')
    law = support_law(mask, equal_nonzero=False, q=q)
    return law[0]+law[1]*first+law[2]*second+law[3]*first*second


def positive_mask_envelope(mask, *, q=256):
    """Unnormalized entrywise envelope, safe for ANY nonnegative output function.

    Double-input mass is1+2/(q-2); retain that surcharge if this looser
    interface is used instead of the product-atom majorant.
    """
    _parameters(mask, q)
    if mask != 3:
        return support_law(mask, q=q)
    return (Fraction(0), Fraction(1, q-2), Fraction(1, q-2), Fraction(1))


def primitive_counts():
    return dict(variable_products=1, fixed_products=0, binary_xors=3,
                alpha_choices=254, field_order=256, final_diagonal_not_included=True,
                input_gauge_transfer_proved=False, native_performance_claim=False)
