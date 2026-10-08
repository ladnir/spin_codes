"""Exact tests of q-ary blockwise rank-one refresh and SPIN chronology."""

from fractions import Fraction as F
from itertools import product
import unittest

import hd_blockwise_inner as proof


class BlockwiseInnerTests(unittest.TestCase):
    def test_all_small_rank_one_maps_are_invertible(self):
        for q,state in ((2,3),(4,2),(8,2)):
            count = 0
            for u,v in proof.allowed_pairs(q,state):
                pairing = proof.dot(v,u,q)
                denominator = 1^pairing
                self.assertNotEqual(denominator,0)
                reverse_u = proof.scale(proof.inverse(denominator,q),u,q)
                for index in range(state):
                    basis = tuple(int(j==index) for j in range(state))
                    transformed = proof.apply_rank_one(basis,u,v,q)
                    self.assertEqual(proof.apply_rank_one(transformed,reverse_u,v,q),basis)
                count += 1
            self.assertEqual(count,(q**state-1)*q**(state-1)*(q-1))

    def test_exact_image_law(self):
        for q,state in ((2,3),(4,2),(4,3),(8,2)):
            vectors = tuple(product(range(q),repeat=state))
            sources = vectors if q<=4 and state<=2 else (
                (0,)*state,(1,)+(0,)*(state-1),tuple((i+1)%q for i in range(state)))
            for source in sources:
                actual = proof.enumerated_image(source,q)
                expected = {target:p for target in vectors
                            if (p:=proof.theoretical_image(source,target,q))}
                self.assertEqual(actual,expected)

    def test_strict_transvection_needs_projective_correction(self):
        q,source = 4,(1,0)
        strict = proof.enumerated_image(source,q,strict=True)
        self.assertEqual(strict.get((2,0),0),0)
        uniform = F(1,q**2-1)
        self.assertGreater((1-F(1,q))*uniform,0)
        for strict in (False,True):
            actual = proof.enumerated_image(source,q,strict=strict,symmetrize=True)
            expected = {target:p for target in product(range(q),repeat=2)
                        if (p:=proof.theoretical_image(source,target,q,scalar_symmetrized=True))}
            self.assertEqual(actual,expected)

    def test_composed_law_and_feedback_chronology(self):
        q,state = 4,2
        vectors = tuple(product(range(q),repeat=state))
        source = (1,2)
        once = proof.enumerated_image(source,q)
        twice = {target:sum(p*proof.theoretical_image(middle,target,q)
                           for middle,p in once.items()) for target in vectors}
        self.assertEqual(twice,{target:proof.theoretical_image(source,target,q,rounds=2)
                               for target in vectors})
        for rounds in (0,1,2,3):
            weights = proof.law_weights(q,state,rounds)
            self.assertEqual(sum(weights),1)
            self.assertTrue(all(w>=0 for w in weights))
            for syndrome in vectors:
                law = proof.step_law(source,syndrome,q,rounds)
                self.assertEqual(sum(law.values()),1)
                self.assertEqual(law.get((0,0),0),proof.return_probability(source,syndrome,q,rounds))
                if not any(syndrome):
                    self.assertNotIn((0,0),law)

    def test_binary_recovery_and_nonuniform_qary_lazy_branch(self):
        for state in (1,2,3):
            for rounds in (1,2,3):
                point,line,refresh = proof.law_weights(2,state,rounds)
                # Over GF2 a nonzero line is one point, so the first two
                # operators coincide and recover the binary IMT law.
                self.assertEqual(point+line,F(1,2)**rounds)
                self.assertEqual(refresh,1-F(1,2)**rounds)
        q,state=4,2
        same = proof.theoretical_image((1,0),(1,0),q)
        collinear = proof.theoretical_image((1,0),(2,0),q)
        outside = proof.theoretical_image((1,0),(0,1),q)
        self.assertEqual((same,collinear,outside),(F(4,15),F(1,15),F(1,20)))
        self.assertNotEqual(collinear,outside)
        self.assertNotEqual(same,F(1,q)+(1-F(1,q))*F(1,q**state-1))

    def test_uniform_fiber_sampler(self):
        for q,state in ((2,3),(4,2),(8,2)):
            for u in product(range(q),repeat=state):
                if not any(u):
                    continue
                samples = tuple(proof.sample_fiber(u,q))
                target = {v for v in product(range(q),repeat=state)
                          if proof.dot(v,u,q)!=1}
                self.assertEqual(len(samples),len(set(samples)))
                self.assertEqual(set(samples),target)

    def test_literal_blockwise_inverse_and_adjoint(self):
        q=4
        # A,C deliberately do NOT satisfy CA=0 or C=A^T.
        emission = ((1,2),(3,1),(2,2))
        feedback = ((1,0,3),(2,1,1))
        pairs = tuple(proof.allowed_pairs(q,2))
        # The middle identity step models a prescribed periodic refresh.
        schedule = ((pairs[17],pairs[53]),(),(pairs[119],))
        count = 0
        for index in range(9):
            for coefficient in range(1,q):
                flat = tuple(coefficient if j==index else 0 for j in range(9))
                inputs = tuple(flat[3*i:3*i+3] for i in range(3))
                output = proof.encode(inputs,emission,feedback,schedule,q)
                self.assertEqual(proof.decode(output,emission,feedback,schedule,q),inputs)
                for j in range(9):
                    dual_flat = tuple(int(k==j) for k in range(9))
                    dual = tuple(dual_flat[3*i:3*i+3] for i in range(3))
                    transposed = proof.adjoint(dual,emission,feedback,schedule,q)
                    self.assertEqual(proof.dot(tuple(x for row in output for x in row),dual_flat,q),
                                     proof.dot(flat,tuple(x for row in transposed for x in row),q))
                    count += 1
        self.assertEqual(count,243)

    def test_single_vector_law_is_not_a_quotient_rank_law(self):
        # With state forms the coordinate basis and target equal to old
        # state, M-I=uv^T has rank0 only for v=0 and rank1 otherwise.
        q,state,tau=4,3,F(3)
        pairs=tuple(proof.allowed_pairs(q,state))
        actual=sum((F(1) if not any(v) else 1/tau) for _,v in pairs)/len(pairs)
        zero_probability=F(1,q**(state-1)*(q-1))
        self.assertEqual(actual,zero_probability+(1-zero_probability)/tau)
        self.assertGreater(actual,tau**-state)


if __name__=='__main__':
    unittest.main()
