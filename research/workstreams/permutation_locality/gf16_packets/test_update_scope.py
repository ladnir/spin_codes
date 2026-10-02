import unittest
from unittest.mock import patch
from types import SimpleNamespace
from fractions import Fraction as Q
from flint import arb,arb_mat
import occupancy_rank
import occupancy_birth_classes
import single_group


class InnerUpdateScopeTests(unittest.TestCase):
    def test_exact_sampled_transvection_composition(self):
        # Enumerate the actual (u, v) sampler, independently of the
        # lazy/refresh formula used by the transfer-operator tests.
        for bits in (2,3,4):
            size=1<<bits
            samples=[(u,v) for u in range(1,size) for v in range(size)
                     if not (u&v).bit_count()%2]
            transition=[[Q(0) for _ in range(size)] for _ in range(size)]
            for a in range(size):
                for u,v in samples:
                    b=a^(u if (a&v).bit_count()%2 else 0)
                    transition[a][b]+=Q(1,len(samples))
            product=[[Q(int(a==b)) for b in range(size)] for a in range(size)]
            for updates in range(1,9):
                product=[[sum(product[a][c]*transition[c][b] for c in range(size))
                          for b in range(size)] for a in range(size)]
                alpha=Q(1,2**updates)
                for a in range(size):
                    for b in range(size):
                        expected=(Q(int(b==0)) if not a else
                                  alpha*int(a==b)+(1-alpha)*Q(int(b!=0),size-1))
                        self.assertEqual(product[a][b],expected)

    def test_operator_factories_build_the_requested_inner(self):
        matrix=arb_mat([[int(i==j) for j in range(3)] for i in range(3)])
        for updates in (2,3,4):
            args=SimpleNamespace(groups=1,tilts=['.001'],precision=192,exact_feedback=True,updates=updates)
            for module,actual_name in ((occupancy_rank,'rank_return.actual'),(occupancy_birth_classes,'actual')):
                with patch.object(module,'outward',return_value=[matrix]), \
                        patch('occupancy_model.placement',return_value=[matrix]), \
                        patch(module.__name__+'.'+actual_name,return_value={}) as actual:
                    result=module.build_operators(args)
                    actual.assert_called_once_with(updates)
                    self.assertEqual(result['.001','1'][0],[matrix])

    def test_single_group_forwards_update_count(self):
        matrix=arb_mat([[1,0,0],[0,1,0],[0,0,1]])
        sparse=single_group.sparse_cover.sparse
        for updates in (2,3,4):
            with patch.object(sparse,'authenticated_caps',return_value={}), \
                    patch.object(sparse,'weighted_cdf_upper',return_value=[0]*256+[1]), \
                    patch.object(sparse,'integer_cdf',side_effect=lambda v:v), \
                    patch.object(occupancy_rank,'build_operators',return_value={('.001','1'):([matrix,matrix],[])}) as build, \
                    patch.object(single_group,'support_moments',return_value=[arb(2)**-120]*257):
                result=single_group.run(1,['.001'],192,updates=updates)
                self.assertGreater(result,0)
                self.assertEqual(build.call_args.args[0].updates,updates)

    def test_unsupported_update_counts_rejected(self):
        for updates in (0,1,5,True,'3'):
            with self.assertRaises(ValueError):single_group.run(1,['.001'],192,updates=updates)
            args=SimpleNamespace(updates=updates,exact_feedback=True)
            for module in (occupancy_rank,occupancy_birth_classes):
                with self.assertRaises(ValueError):module.build_operators(args)


if __name__=='__main__':unittest.main()
