"""Toy-only exact checks for typed packet averaging and safe fallback."""
from fractions import Fraction as Q
from itertools import product
from math import comb, factorial, prod
from types import SimpleNamespace
import unittest

from flint import arb, ctx, fmpq, fmpq_mat

from typed_inner import build_typed_operators
from typed_placement import typed_placement


def as_fraction(value):
    value=value.fmpq() if isinstance(value,arb) else value
    return Q(int(value.p),int(value.q))


def as_fmpq(value):
    return fmpq(value.numerator,value.denominator)


def counts_at(active):
    return [counts for counts in product(range(active+1),repeat=4) if sum(counts)==active]


def toy_model(cut=3, slots=4):
    all_shapes={}
    for active in range(slots+1):
        all_shapes[active]={counts:fmpq_mat([[1+counts[0]+2*counts[2],1+counts[1]],
                                             [counts[3],2+3*counts[1]+counts[3]]])
                            for counts in counts_at(active)}
    universal=[]
    for active in range(slots+1):
        values=list(all_shapes[active].values())
        universal.append(fmpq_mat([[max(value[i,j] for value in values) for j in range(2)]
                                   for i in range(2)]))
    model=SimpleNamespace(cut=cut,precision=192,
                          shapes={k:v for k,v in all_shapes.items() if k<=cut},fallback=universal)
    return model,all_shapes,universal


def exact_average(shapes, law1, law2, j1, j2):
    result=fmpq_mat(2,2)
    total=Q(0)
    for weights in product(range(5),repeat=j1+j2):
        mass=prod(law1[b] for b in weights[:j1])*prod(law2[b] for b in weights[j1:])
        counts=tuple(weights.count(b) for b in range(1,5))
        total+=mass
        if mass:
            result+=shapes[sum(counts)][counts]*as_fmpq(mass)
    assert total==1
    return result


def homogeneous_average(shapes, universal, law, selected, *, coarse=False):
    zero=law[0]
    result=fmpq_mat(2,2)
    if zero==1:
        return shapes[0][(0,0,0,0)]
    theta=[p/(1-zero) for p in law[1:]]
    for active in range(selected+1):
        activation=Q(comb(selected,active))*zero**(selected-active)*(1-zero)**active
        if coarse:
            mixed=universal[active]
        else:
            mixed=fmpq_mat(2,2)
            for counts,value in shapes[active].items():
                mass=Q(factorial(active),prod(factorial(n) for n in counts))*prod(p**n for p,n in zip(theta,counts))
                mixed+=value*as_fmpq(mass)
        result+=mixed*as_fmpq(activation)
    return result


class TypedInnerTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        ctx.prec=192

    def assert_dominates(self, actual, expected, *, equal=False):
        for i in range(expected.nrows()):
            for j in range(expected.ncols()):
                self.assertTrue(actual[i,j].is_exact())
                if equal:
                    self.assertEqual(as_fraction(actual[i,j]),as_fraction(expected[i,j]))
                else:
                    self.assertGreaterEqual(as_fraction(actual[i,j]),as_fraction(expected[i,j]))

    def test_exact_two_law_enumeration(self):
        model,shapes,_=toy_model(cut=4)
        law1=(Q(1,3),Q(1,6),Q(1,6),Q(1,6),Q(1,6))
        law2=(Q(1,7),Q(2,7),Q(0),Q(1,7),Q(3,7))
        result=build_typed_operators(model,law1,law2,max_counts=(2,2),slots=4)
        self.assertEqual(set(result),set(product(range(3),repeat=2)))
        self.assertNotEqual(shapes[1][(1,0,0,0)]*shapes[1][(0,1,0,0)],
                            shapes[1][(0,1,0,0)]*shapes[1][(1,0,0,0)])
        for (j1,j2),value in result.items():
            self.assert_dominates(value,exact_average(shapes,law1,law2,j1,j2))

    def test_deterministic_limits(self):
        model,shapes,_=toy_model(cut=4)
        for first,second in product(range(5),repeat=2):
            law1=tuple(Q(int(b==first)) for b in range(5))
            law2=tuple(Q(int(b==second)) for b in range(5))
            result=build_typed_operators(model,law1,law2,max_counts=(2,2),slots=4)
            for (j1,j2),value in result.items():
                counts=tuple(j1*int(first==b)+j2*int(second==b) for b in range(1,5))
                self.assert_dominates(value,shapes[sum(counts)][counts],equal=True)

    def test_equal_laws_match_homogeneous_thinning(self):
        model,shapes,universal=toy_model(cut=2)
        law=(Q(1,2),Q(1,4),Q(0),Q(1,4),Q(0))
        result=build_typed_operators(model,law,law,max_counts=(3,3),slots=4)
        for (j1,j2),value in result.items():
            selected=j1+j2
            expected=homogeneous_average(shapes,universal,law,selected,coarse=selected>model.cut)
            self.assert_dominates(value,expected,equal=True)

    def test_universal_fallback_dominates_exact_shape_average(self):
        model,shapes,_=toy_model(cut=1)
        law1=(Q(1,2),Q(1,4),Q(1,4),Q(0),Q(0))
        law2=(Q(1,4),Q(0),Q(0),Q(1,4),Q(1,2))
        result=build_typed_operators(model,law1,law2,max_counts=(2,2),slots=4)
        strict=False
        for (j1,j2),value in result.items():
            expected=exact_average(shapes,law1,law2,j1,j2)
            self.assert_dominates(value,expected)
            if j1+j2>model.cut:
                strict|=any(as_fraction(value[i,j])>as_fraction(expected[i,j])
                            for i in range(2) for j in range(2))
        self.assertTrue(strict)
        all_zero=(Q(1),Q(0),Q(0),Q(0),Q(0))
        zero_result=build_typed_operators(model,all_zero,all_zero,max_counts=(4,4),slots=4)
        for value in zero_result.values():
            self.assert_dominates(value,shapes[0][(0,0,0,0)],equal=True)

    def test_composition_with_typed_placement(self):
        model,shapes,_=toy_model(cut=2,slots=2)
        law1=(Q(1,2),Q(1,2),Q(0),Q(0),Q(0))
        law2=(Q(1,4),Q(0),Q(0),Q(1,4),Q(1,2))
        local=build_typed_operators(model,law1,law2,max_counts=(1,1),slots=2)
        actual=typed_placement(local,(1,1),epochs=2,slots=2)
        total=fmpq_mat(2,2)
        placements=0
        for labels in product((0,1,2),repeat=4):
            if labels.count(1)!=1 or labels.count(2)!=1:
                continue
            placements+=1
            for first,second in product(range(5),repeat=2):
                mass=law1[first]*law2[second]
                weights=[first if label==1 else second if label==2 else 0 for label in labels]
                result=fmpq_mat([[1,0],[0,1]])
                for epoch in range(2):
                    local_weights=weights[epoch*2:(epoch+1)*2]
                    counts=tuple(local_weights.count(b) for b in range(1,5))
                    result=result*shapes[sum(counts)][counts]
                total+=result*as_fmpq(mass)
        self.assertEqual(placements,12)
        self.assert_dominates(actual,total/placements)

    def test_validation_and_precision_restoration(self):
        model,_,_=toy_model()
        law=(Q(1),Q(0),Q(0),Q(0),Q(0))
        previous=ctx.prec
        try:
            ctx.prec=128
            result=build_typed_operators(model,law,law,max_counts=(0,0),slots=4)
            self.assertEqual(set(result),{(0,0)})
            self.assertEqual(ctx.prec,128)
        finally:
            ctx.prec=previous
        for invalid in ((1,0,0,0),(2,0,0,0,0),(2,-1,0,0,0)):
            with self.assertRaises(ValueError):
                build_typed_operators(model,invalid,law)
        with self.assertRaises(ValueError):
            build_typed_operators(model,law,law,max_counts=(-1,1))
        broken,_,_=toy_model()
        del broken.shapes[1][(1,0,0,0)]
        with self.assertRaises(ValueError):
            build_typed_operators(broken,law,law,max_counts=(1,1))
        broken,_,_=toy_model(cut=0,slots=0)
        with self.assertRaises(ValueError):
            build_typed_operators(broken,law,law,max_counts=(1,1))


if __name__=='__main__':
    unittest.main()
