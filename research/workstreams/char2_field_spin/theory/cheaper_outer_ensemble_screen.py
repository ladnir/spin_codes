"""Exact necessary gates and local laws for cheaper outer ensembles.

No native code, timing, or whole-code certificate. All acceptance comparisons
below use integers/Fractions. The GL postmix gate concerns the OUTER distance
only; its matrices are shared by the L rows unless explicitly changed.
"""

from fractions import Fraction as F
from math import comb
import json


def window_gate(rows, support_regions, subcode_dimension, state=8,
                total_regions=256, threshold=None, alignment=16,
                discard_tail_rows=0):
    """Any supported outer subcode yields a deterministic confinement gate.

    Choose E aligned slots in one window per supported region. Shortening
    removes at most h*L-E dimensions. State-zero exits remove at most s*h.
    """
    if not (rows > 0 and 0 < subcode_dimension <= support_regions <= total_regions
            and alignment > 0 and rows % alignment == 0 and state >= 0
            and 0 <= discard_tail_rows < rows):
        raise ValueError('invalid aligned geometry')
    if threshold is None:
        threshold = (rows*total_regions+9)//10-1
    discarded_dimension = subcode_dimension*discard_tail_rows
    need = ((support_regions-subcode_dimension)*rows+discarded_dimension
            +state*support_regions+1)
    budget = alignment*((need+alignment-1)//alignment)
    possible = budget <= min(threshold, support_regions*rows)
    return dict(rows=rows,support_regions=support_regions,dimension=subcode_dimension,
                budget=budget,threshold=threshold,
                discarded_tail_rows=discard_tail_rows,
                shortened_dimension_lower=(budget-(support_regions-subcode_dimension)*rows
                                            -discarded_dimension),
                exit_rank_upper=state*support_regions,deterministically_bad=possible)


def diagonal_binary_postmix_lower(q=256, length=32, support=13):
    """Bonferroni lower for the fixed binary (I+J8) times random diagonal mix.

    Use just one component pair and ratios of distinct monic irreducible
    quadratics. All arithmetic is exact. For the main parameters, each event
    gives a weight26 word in eight postmixed GRS[32,22] components.

    For support<=4, distinct ratio events on the SAME support need not be
    disjoint; the last intersection term deliberately retains them. This
    permits tiny-field exhaustive checks of the pair upper bound.
    """
    if not (q >= 2 and 1 <= support <= length <= q):
        raise ValueError('invalid evaluation geometry')
    irreducibles = q*(q-1)//2
    ratios = irreducibles*(irreducibles-1)
    supports = comb(length,support)
    single = F(supports*ratios,(q-1)**(support-1))
    pair_upper = F(0)
    for shared in range(max(0,2*support-length),support+1):
        ordered_support_pairs = (supports*comb(support,shared)
                                 *comb(length-support,support-shared))
        joint = F(1,(q-1)**(2*support-2 if shared==0 else 2*support-shared-1))
        # A repeated event is not a Bonferroni pair.
        same_ratio_pairs = ratios if shared < support else 0
        # Degree-four cross multiplication rules out five shared points.
        different_ratio_pairs = ratios*(ratios-1) if shared <= 4 else 0
        pair_upper += F(ordered_support_pairs*(same_ratio_pairs+different_ratio_pairs),2)*joint
    return dict(first_moment=single,pair_upper=pair_upper,
                failure_lower=max(F(0),single-pair_upper),irreducibles=irreducibles,
                projective_ratios=ratios)


def support_count(length, dimension, q, weight):
    """Nonzero MDS words with one specified support (not all supports)."""
    parity = length-dimension
    if not parity < weight <= length:
        return 0
    return (q-1)*sum((-1)**j*comb(weight-1,j)*q**(weight-parity-1-j)
                     for j in range(weight-parity))


def multiply_truncated(left, right, degree):
    out = [0]*(min(degree,len(left)+len(right)-2)+1)
    for i,a in enumerate(left):
        for j,b in enumerate(right[:max(0,degree-i+1)]):
            out[i+j] += a*b
    return out


def gl_postmix_lowweight_upper(cutoff=26, q=256, components=8,
                             length=32, dimension=22):
    """Projective first moment for outer distance after independent GL(c,q).

    The base is c copies of the SAME MDS code. A nonzero column becomes
    uniform among q^c-1 vectors, independently at each of the length columns.
    """
    alphabet = q**components
    numerator = [0]+[comb(components,w)*(q-1)**w for w in range(1,components+1)]
    polynomial = [1]
    answer = F(0)
    for weight in range(1,min(length,cutoff)+1):
        polynomial = multiply_truncated(polynomial,numerator,cutoff)
        count = support_count(length,dimension,alphabet,weight)
        answer += F(comb(length,weight)*count*sum(polynomial),
                    (q-1)*(alphabet-1)**weight)
    return answer


def gl_postmix_kernel_mean(allowed, q=256, components=8, dimension=22):
    """Exact E[q^dim(shortening)] for ONE postmixed outer row.

    allowed[j] is the number of allowed scalar coordinates in column j.
    This expectation must NOT be raised to L for a shared outer setup.
    """
    length = len(allowed)
    if any(not 0 <= t <= components for t in allowed):
        raise ValueError('invalid allowed-coordinate count')
    alphabet = q**components
    elementary = [F(1)]
    for t in allowed:
        probability = F(q**t-1,alphabet-1)
        elementary = [elementary[0]]+[elementary[i]+probability*elementary[i-1]
                       for i in range(1,len(elementary))]+[probability*elementary[-1]]
    return 1+sum(support_count(length,dimension,alphabet,u)*elementary[u]
                 for u in range(1,length+1))


def sparse_zero_probability(weight, q, density):
    """A sum of weight independent sparse-uniform GF(q) entries."""
    if weight < 0 or q < 2 or not 0 <= density <= 1:
        raise ValueError('invalid sparse law')
    return F(1,q)+F(q-1,q)*(1-F(q, q-1)*density)**weight


def sparse_kernel_mean(allowed_source, disallowed_parity, q, density):
    """Exact E[q^dim ker(P_restricted)] for independent sparse entries."""
    return sum(comb(allowed_source,a)*(q-1)**a
               *sparse_zero_probability(a,q,density)**disallowed_parity
               for a in range(allowed_source+1))


def sparse_basis_failure_lower(density, source=176, parity=80, outer_cutoff=26):
    """Exact probability some systematic generator column has weight<=cutoff.

    All parity entries have independent Bernoulli nonzero indicators. The
    deterministic window gate converts this outer event to whole-code failure.
    """
    if not 0 <= density <= 1 or not 1 <= outer_cutoff <= parity+1:
        raise ValueError('invalid independent-entry sparse law')
    per_column = sum(comb(parity,j)*density**j*(1-density)**(parity-j)
                     for j in range(outer_cutoff))
    return 1-(1-per_column)**source


def shear_pair_law(mask, q):
    """Exact support law of u=x+b*y; v=y+a*u, a,b uniform nonzero.

    bit0 is the first wire. The law is independent of the nonzero input values.
    """
    r = F(1,q-1)
    if mask == 0:
        return {0:F(1)}
    if mask == 1:
        return {3:F(1)}
    if mask == 2:
        return {1:r,3:1-r}
    if mask == 3:
        return {2:r,1:(1-r)*r,3:(1-r)**2}
    raise ValueError('two-wire mask required')


def butterfly_support_law(mask, q=256, stages=3):
    """Exact eight-wire support kernel for independent two-shear switches."""
    if not 0 <= mask < 256 or not 0 <= stages <= 3:
        raise ValueError('invalid butterfly')
    law = {mask:F(1)}
    for bit in range(stages):
        for left in range(8):
            right = left ^ (1 << bit)
            if left >= right:
                continue
            following = {}
            clear = ~((1 << left)|(1 << right))
            for previous,probability in law.items():
                pair = ((previous >> left)&1)|(((previous >> right)&1) << 1)
                for output,weight in shear_pair_law(pair,q).items():
                    result = (previous&clear)|((output&1) << left)|(((output >> 1)&1) << right)
                    following[result] = following.get(result,F(0))+probability*weight
            law = following
    return law


def integer_bits(upper):
    if not 0 < upper < 1:
        raise ValueError('not a nontrivial upper bound')
    bits = upper.denominator.bit_length()-upper.numerator.bit_length()
    while upper > F(2)**(-bits):
        bits -= 1
    while upper <= F(2)**(-bits-1):
        bits += 1
    return bits


if __name__ == '__main__':
    gates = [window_gate(6016,h,k) for h,k in ((16,11),(32,22),(64,44),(128,88),(26,1),(27,1))]
    full = gl_postmix_lowweight_upper()
    diagonal = diagonal_binary_postmix_lower()
    print(json.dumps(dict(window_gates=gates,
        useful_k_window_gates=[window_gate(6016,h,k,discard_tail_rows=59)
                              for h,k in ((16,11),(32,22),(64,44),(26,1))],
        full_gl8_outer_weight_at_most26_upper_bits=integer_bits(full),
        binary_diagonal_failure_lower_display=float(diagonal['failure_lower']),
        binary_diagonal_failure_lower_gt_2neg40=diagonal['failure_lower'] > F(1,2**40),
        independent_sparse_basis_failure_gt_2neg40_at_density073=(
            sparse_basis_failure_lower(F(73,100)) > F(1,2**40)),
        full_gl8_outer_only=True,whole_code=False,
        postmix_products_per_four_rows=dict(full_gl8=2048,two_shear_three_stages=768,
                                           obstructed_one_product_three_stages=384))))
