"""Class-specific mature-tail arrivals, preserving the overlapping-tail cone.

Only one coupled column form is replaced at a time. The replacement is
no larger for actual masses M >= L56 >= L48 >= 0. Propagated coordinate
uppers need not be nested, so numerical matrix products need not improve.
No construction parameter or existing source module is changed.
"""
from fractions import Fraction as Q
from itertools import combinations, product
from types import SimpleNamespace
from pathlib import Path
import sys

from flint import arb, arb_mat, ctx

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from occupancy_memory import M, C
from mature_tail import L48, L56
from group_rank_one_verify import up


def exact_fraction(value):
    """The existing outward coefficient must be a point, not an interval."""
    if not value.is_exact():
        raise ValueError('exact outward endpoint coefficient required')
    fraction = value.fmpq()
    return Q(int(fraction.numerator),int(fraction.denominator))


def exact_dyadic(value):
    """Preserve coupled differences exactly, even at restored lower precision."""
    value = Q(value)
    denominator = value.denominator
    assert value >= 0 and denominator & (denominator-1) == 0
    precision = ctx.prec
    try:
        ctx.prec = max(precision, value.numerator.bit_length()+2)
        result = arb(value.numerator)*arb(2)**(-(denominator.bit_length()-1))
        assert result.is_exact() and exact_fraction(result) == value
        return result
    finally:
        ctx.prec = precision


def split_class_bounds(bounds):
    """Return exact dyadic (a,b,c), jointly bounding all source classes."""
    assert 48 in bounds and 56 in bounds and any(v >= 64 for v in bounds)
    assert all(v in (48,56) or v >= 64 for v in bounds)
    a = max(bound for v,bound in bounds.items() if v >= 64)
    ab = max(a,bounds[56])
    abc = max(ab,bounds[48])
    return a,ab-a,abc-ab


def refine(t, weights, data, powers, rounds=2):
    """Refine tail-arrival columns after the stated transvection transform.

    weights are nonzero packet weights in distinct four-bit windows.
    data uses shape_inner.prepare_shared's inputs/tails/pairs fields.
    powers[k] must upper-bound exp(-lambda*k). The function does not
    infer uniformity after weighting paths by their output.
    """
    weights = tuple(sorted(weights))
    j, total = len(weights),sum(weights)
    if not isinstance(rounds,int) or rounds < 1:
        raise ValueError('a positive integer transvection count is required')
    if not j:
        return t*1
    if not 1 <= j <= 32 or any(not 1 <= b <= 4 for b in weights):
        raise ValueError('one to 32 nonzero packet weights are required')
    if t.nrows() != 11 or t.ncols() != 11:
        raise ValueError('the existing eleven-coordinate envelope is required')
    levels = sorted(data.inputs[3])
    assert 48 in levels and 56 in levels and all(v in (48,56) or v >= 64 for v in levels)
    lazy = arb(2)**(-rounds)
    result = t*1
    for target,cutoff in ((L48,48),(L56,56)):
        # Current shape_inner gives this arrival solely from mature mass
        # and its two subsets, not from its separate pointwise density cap.
        if t[C,target] != 0:
            raise ValueError('density-dependent mature-tail arrival is not supported')
        try:
            old = tuple(exact_fraction(t[source,target]) for source in (M,L56,L48))
        except ValueError:
            # Keeping an existing valid column is safer than mixing its
            # interval coefficients with exact coupled differences.
            continue
        assert min(old) >= 0
        if j == 1:
            probabilities = data.tails[weights[0]][cutoff][3]
        elif j == 2:
            probabilities = data.pairs[weights][cutoff][2]
        else:
            # Conditioning the other packets translates the source state;
            # its original expansion class cannot constrain that translated
            # state. Retain only the arbitrary-state pair bound here.
            pairs = set(combinations(weights,2))
            probability = min(Q(1), min(data.pairs[p][cutoff][0] for p in pairs)
                              *Q(32*31,(34-j)*(33-j)))
            probabilities = {v:probability for v in levels}
        bounds = {}
        for v in levels:
            probability = Q(probabilities[v])
            assert 0 <= probability <= 1
            value = up(lazy*powers[abs(v-total)]*arb(probability.numerator)/probability.denominator)
            old_class = old[0]+(old[1] if v <= 56 else 0)+(old[2] if v <= 48 else 0)
            bounds[v] = min(old_class,exact_fraction(value))
        new = split_class_bounds(bounds)
        # These cumulative inequalities, rather than entrywise coefficient
        # minima, express domination on the nested-subset cone.
        assert all(sum(new[:n]) <= sum(old[:n]) for n in (1,2,3))
        for source,value in zip((M,L56,L48),new):
            result[source,target] = exact_dyadic(value)
        assert tuple(exact_fraction(result[source,target]) for source in (M,L56,L48)) == new
    return result


def self_test():
    """Exact cone tests and independent one-/two-window toy enumeration."""
    checks = 0
    # Arbitrary dyadic coefficient forms and classwise valid upper bounds.
    for old in product((Q(0),Q(1,8),Q(1,2)),repeat=3):
        old_bounds = {64:old[0],72:old[0],56:sum(old[:2]),48:sum(old)}
        for fractions in product((Q(0),Q(1,3),Q(1)),repeat=4):
            bounds = {v:bound*f for (v,bound),f in zip(old_bounds.items(),fractions)}
            new = split_class_bounds(bounds)
            assert all(sum(new[:n]) <= sum(old[:n]) for n in (1,2,3))
            assert new[0] >= bounds[64] and new[0] >= bounds[72]
            assert sum(new[:2]) >= bounds[56] and sum(new) >= bounds[48]
            for m48,m56,mhigh in ((1,0,0),(0,1,0),(0,0,1),(1,2,3)):
                incoming = (m48+m56+mhigh,m48+m56,m48)
                assert sum(a*b for a,b in zip(new,incoming)) <= sum(a*b for a,b in zip(old,incoming))
                checks += 1
    # A separate two-bit state space embedded in128 output coordinates.
    # Its three nonzero images have weights48,56,64, matching the source
    # classes needed by the production tail interface.
    first = (1 << 48)-1
    second = ((1 << 20)-1)|(((1 << 36)-1) << 48)
    images = (0,first,second,first ^ second)
    assert [x.bit_count() for x in images] == [0,48,56,64]
    columns = (1,2,3,1,2,3,1,2,3,1,2,3,1,3,2,1)
    single = [(columns[4*w+b],1 << (4*w+b),w) for w in range(4) for b in range(4)]
    cases = {(1,):[(s,x) for s,x,_ in single],
             (1,1):[(sa ^ sb,xa ^ xb) for sa,xa,wa in single for sb,xb,wb in single if wa != wb]}
    precision = ctx.prec
    try:
        ctx.prec = max(precision,768)
        z = Q(7,8)
        powers = [exact_dyadic(z**w) for w in range(145)]
        for shape,entries in cases.items():
            probabilities = {}
            actual = {}
            for cutoff in (48,56):
                probabilities[cutoff] = {}
                actual[cutoff] = {}
                for state in (1,2,3):
                    v = images[state].bit_count()
                    selected = [(s,x) for s,x in entries if 0 < images[state ^ s].bit_count() <= cutoff]
                    probabilities[cutoff][v] = Q(len(selected),len(entries))
                    actual[cutoff][v] = sum((z**(images[state] ^ x).bit_count() for _,x in selected),Q(0))/(4*len(entries))
            records = {cut:(max(p.values()),None,p,p) for cut,p in probabilities.items()}
            data = SimpleNamespace(inputs=(None,None,None,{48:1,56:1,64:1}),
                                   tails={1:records},pairs={(1,1):records})
            before = arb_mat(11,11)
            for target,cutoff in ((L48,48),(L56,56)):
                upper = z**(48-sum(shape))*max(probabilities[cutoff].values())/4
                before[M,target] = up(arb(upper.numerator)/upper.denominator)
                before[L56,target] = exact_dyadic(Q(1,64))
                before[L48,target] = exact_dyadic(Q(1,128))
            after = refine(before,shape,data,powers)
            for target,cutoff in ((L48,48),(L56,56)):
                old = [exact_fraction(before[source,target]) for source in (M,L56,L48)]
                new = [exact_fraction(after[source,target]) for source in (M,L56,L48)]
                for n48,n56,n64 in product(range(4),repeat=3):
                    incoming = (n48+n56+n64,n48+n56,n48)
                    observed = sum(n*actual[cutoff][v] for n,v in zip((n48,n56,n64),(48,56,64)))
                    new_value = sum(a*b for a,b in zip(new,incoming))
                    old_value = sum(a*b for a,b in zip(old,incoming))
                    assert observed <= new_value <= old_value
                    checks += 1
            assert all(after[i,k] == before[i,k] for i in range(11) for k in range(11)
                       if (i,k) not in {(source,target) for source in (M,L56,L48) for target in (L48,L56)})
            unchanged = refine(before,(),data,powers)
            assert unchanged == before
            for invalid_rounds in (0,-1,Q(3,2),2.0):
                try:
                    refine(before,shape,data,powers,rounds=invalid_rounds)
                except ValueError:
                    checks += 1
                else:
                    raise AssertionError('nonpositive or noninteger rounds accepted')
    finally:
        ctx.prec = precision
    print('Class-tail exact coupled-cone and direct toy-window checks:',checks,'passed',flush=True)


if __name__ == '__main__':
    self_test()
