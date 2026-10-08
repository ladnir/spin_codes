"""Directed class bounds for the GF65536 randomized-outer SPIN ensemble.

These are whole-class bounds, not a complete distance certificate.  Scalar
randomizers are independent uniform nonzero field elements at every outer
coordinate.  Routing is regional, not a uniform global permutation.
"""

from dataclasses import dataclass
from fractions import Fraction
from functools import lru_cache
from math import comb


@dataclass(frozen=True)
class Parameters:
    field: int = 65536
    rows: int = 8192
    outer_n: int = 256
    outer_k: int = 128
    width: int = 64
    state: int = 4
    limit: int = 838860

    @property
    def length(self):
        return self.rows*self.outer_n

    @property
    def packets(self):
        return self.length//self.width

    @property
    def beta(self):
        return Fraction(self.field**self.outer_n,
                        (self.field-1)**(self.outer_n-self.outer_k))


def prescribed_mds_count(weight, length, dimension, field):
    """Exact number of MDS words with one prescribed exact support."""
    if not 0 <= weight <= length or not 0 <= dimension <= length:
        raise ValueError('Invalid MDS parameters')
    if weight == 0:
        return 1
    distance = length-dimension+1
    if weight < distance:
        return 0
    return (field-1)*sum((-1)**j*comb(weight-1, j)*field**(weight-distance-j)
                         for j in range(weight-distance+1))


@lru_cache(None)
def quiet_weight_polynomial(width=64, state=4, field=65536):
    return tuple(comb(width, w)*prescribed_mds_count(w, width, width-state, field)
                 for w in range(width+1))


def independent_nonzero_syndrome_zero(weight, state=4, field=65536):
    """Pr[C X=0] for fixed support and independent nonzero amplitudes.

    Any <=state columns of C are independent.  A restriction to more columns
    checks an MDS code; the returned law is BEFORE output conditioning.
    """
    if weight < 0:
        raise ValueError('Negative support size')
    if weight == 0:
        return Fraction(1)
    if weight <= state:
        return Fraction(0)
    return Fraction(prescribed_mds_count(weight, weight, weight-state, field),
                    (field-1)**weight)


def outer_pointwise_counts(length=256, dimension=128, field=65536):
    """Expected multiplicity at each fixed nonzero word, by its weight."""
    return (Fraction(0),)+tuple(
        Fraction(prescribed_mds_count(w, length, dimension, field), (field-1)**w)
        for w in range(1, length+1))


def _arb_fraction(value):
    from flint import arb
    return arb(value.numerator)/arb(value.denominator)


def _log2(value):
    from flint import arb
    return value.log()/arb(2).log()


def ball_log2(parameters=Parameters(), *, quiet=False,
              mark=None):
    """Directed log2 of a tilted ball COUNT, including the zero word.

    The quiet ball constrains all packets except the final packet.  It is a
    superset of the quiet encoder outputs at every fixed route and mixer.
    """
    p = parameters
    if mark is None:
        mark = Fraction(133 if quiet else 2,
                        (128 if quiet else 3)*(p.field-1))
    if not 0 < mark < 1:
        raise ValueError('Need a rational mark strictly between zero and one')
    z = _arb_fraction(mark)
    last = (1+(p.field-1)*z)**p.width
    if quiet:
        polynomial = quiet_weight_polynomial(p.width, p.state, p.field)
        # The short polynomial is exact rational before the logarithm.
        value = sum((a*mark**w for w, a in enumerate(polynomial)), Fraction(0))
        return ((p.packets-1)*_log2(_arb_fraction(value))
                +_log2(last)-p.limit*_log2(z))
    return p.length*_log2(1+(p.field-1)*z)-p.limit*_log2(z)


def dense_class_log2(active_rows, parameters=Parameters(), *, quiet=False,
                     mark=None, ball=None):
    """Bound for all words with AT LEAST r active rows, over all locations.

    A fixed output has one inverse input at each fixed route and mixer.  Its
    active-row set is determined, so there is no binomial subset union.  Divide
    nonzero witness count by field-1 using projective scaling.
    """
    from flint import arb
    p = parameters
    if not 1 <= active_rows <= p.rows:
        raise ValueError('Invalid active-row count')
    if ball is None:
        ball = ball_log2(p, quiet=quiet, mark=mark)
    return ball-((p.outer_n-p.outer_k)*active_rows+1)*_log2(arb(p.field-1))


def dense_tail_receipt(parameters=Parameters(), *, quiet=False, failure_bits=40):
    """Find and certify an entire terminal range, with directed comparisons."""
    p = parameters
    ball = ball_log2(p, quiet=quiet)
    left, right = 1, p.rows
    if not dense_class_log2(right, p, ball=ball) < -failure_bits:
        raise ArithmeticError('No certified dense tail at this target')
    while left < right:
        middle = (left+right)//2
        endpoint = dense_class_log2(middle, p, ball=ball)
        if endpoint < -failure_bits:
            right = middle
        elif endpoint > -failure_bits:
            left = middle+1
        else:
            raise ArithmeticError('Increase Arb precision for the threshold test')
    endpoint = dense_class_log2(left, p, ball=ball)
    previous = dense_class_log2(left-1, p, ball=ball)
    if not endpoint < -failure_bits or not previous > -failure_bits:
        raise ArithmeticError('Threshold enclosure failed')
    return dict(class_kind='quiet' if quiet else 'all_words',
                active_rows_first=left, active_rows_last=p.rows,
                ball_log2_upper=str(ball), union_log2_upper=str(endpoint),
                previous_tail_envelope_log2=str(previous),
                all_rows_log2_upper=str(dense_class_log2(p.rows, p, ball=ball)),
                mark=str(Fraction(133 if quiet else 2,
                                  (128 if quiet else 3)*(p.field-1))),
                requested_failure_bits=failure_bits,
                directed_acceptance=True, whole_code_certificate=False)


def quiet_rank_marker_log2(active_rows, marker, alpha=Fraction(1),
                          parameters=Parameters()):
    """Route-conditioned quiet bound at weight mark z=1.

    Uses q^-min(j,s) for j potential input slots and a Bernoulli coefficient
    bound for each region.  The actual route remains uniform without
    replacement; its binomial normalization is charged outer_n times.
    """
    from flint import arb
    p = parameters
    if not 1 <= active_rows < p.rows or not 0 < marker < 1 or not 0 < alpha <= 1:
        raise ValueError('Need interior occupancy/marker and 0<alpha<=1')
    marker, alpha = Fraction(marker), Fraction(alpha)
    ap, aa = _arb_fraction(marker), _arb_fraction(alpha)
    total = arb(0)
    for j in range(p.width+1):
        rank_factor = arb(p.field)**(-aa*min(j, p.state))
        total += comb(p.width, j)*ap**j*(1-ap)**(p.width-j)*rank_factor
    log_choose = _log2(arb(comb(p.rows, active_rows)))
    log_binomial_mass = (log_choose+active_rows*_log2(ap)
                         +(p.rows-active_rows)*_log2(1-ap))
    return (log_choose+aa*(active_rows*_log2(_arb_fraction(p.beta))
                          -_log2(arb(p.field-1)))
            +(p.packets-1)*_log2(total)-p.outer_n*log_binomial_mass)


def quiet_rank_witnesses(parameters=Parameters(), stop=1023):
    """Float screening proposes rational witnesses; it certifies nothing."""
    from math import lgamma, log, log2
    p = parameters
    coefficients = tuple(comb(p.width, j) for j in range(p.width+1))
    markers = sorted(set([Fraction(j, 128) for j in range(1, 128)]
                         +[Fraction(1, 1 << j) for j in range(8, 18)]))
    candidates = []
    log_beta = p.outer_n*log2(p.field)-(p.outer_n-p.outer_k)*log2(p.field-1)
    for numerator in range(1, 33):
        alpha = Fraction(numerator, 32)
        a = float(alpha)
        ranks = tuple(p.field**(-a*min(j, p.state)) for j in range(p.width+1))
        for marker in markers:
            m = float(marker)
            value = sum(coefficients[j]*m**j*(1-m)**(p.width-j)*ranks[j]
                        for j in range(p.width+1))
            slope = a*log_beta-p.outer_n*log2(m)+p.outer_n*log2(1-m)
            intercept = ((p.packets-1)*log2(value)
                         -p.outer_n*p.rows*log2(1-m)-a*log2(p.field-1))
            candidates.append((slope, intercept, marker, alpha))
    witnesses = []
    for rows in range(1, stop+1):
        log_choose = (lgamma(p.rows+1)-lgamma(rows+1)-lgamma(p.rows-rows+1))/log(2)
        best = min(candidates, key=lambda item: item[0]*rows+item[1])
        witnesses.append((rows, best[2], best[3],
                          best[0]*rows+best[1]-(p.outer_n-1)*log_choose))
    return tuple(witnesses)


def quiet_sparse_receipt(parameters=Parameters(), failure_bits=40, stop=1023):
    """Certify a contiguous sparse quiet class using rational witnesses.

    Float screening only chooses the mark and fractional exponent.  Every
    accepted occupancy and the final union use directed Arb arithmetic.
    """
    from flint import arb
    total, weakest, last = arb(0), None, 0
    witnesses = quiet_rank_witnesses(parameters, stop)
    for rows, marker, alpha, _ in witnesses:
        bound = quiet_rank_marker_log2(rows, marker, alpha, parameters)
        candidate = total+arb(2)**bound
        if not candidate < arb(2)**(-failure_bits):
            break
        total, last = candidate, rows
        if weakest is None or bound > weakest[0]:
            weakest = (bound, rows, marker, alpha)
    if not last:
        raise ArithmeticError('No sparse quiet occupancy certified')
    return dict(class_kind='quiet', active_rows_first=1, active_rows_last=last,
                union_log2_upper=str(_log2(total)),
                weakest_active_rows=weakest[1], weakest_log2_upper=str(weakest[0]),
                weakest_marker=str(weakest[2]), weakest_alpha=str(weakest[3]),
                weight_ignored=True, final_packet_unconstrained=True,
                directed_acceptance=True, whole_code_certificate=False)


def main():
    import json
    from flint import ctx
    ctx.prec = 192
    print(json.dumps(dense_tail_receipt(), indent=2))
    print(json.dumps(dense_tail_receipt(quiet=True), indent=2))
    print(json.dumps(quiet_sparse_receipt(), indent=2))
    print(json.dumps(dict(nonzero_syndrome_zero_at_weight_5=str(
        independent_nonzero_syndrome_zero(5))), indent=2))


if __name__ == '__main__':
    main()
