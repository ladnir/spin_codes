from fractions import Fraction as Q
from types import SimpleNamespace
import unittest
from unittest.mock import patch

from flint import arb, arb_mat, ctx
import shared_relaxed_region_cache as cache_module
import birth_classes


def toy_outward(model, cell, witness):
    polynomial = placement(local_operators(model, witness))
    answer = (polynomial[0][0,0] + cache_module.regional.aq(cell[0])
              + cache_module.regional.aq(Q(witness['parameters'][1])))
    polynomial[0][0,0] = arb(999)
    return answer, witness


class RegionCacheTests(unittest.TestCase):
    def setUp(self):
        self.old_precision = ctx.prec
        ctx.prec = 256
        self.model = SimpleNamespace(data=dict(windows=32, bits=19, updates=2))
        self.witness = dict(parameters=['1/10','0','0'], regional_exact_zero=True)

    def tearDown(self):
        ctx.prec = self.old_precision

    def polynomial(self):
        return [arb_mat([[j+1]]) for j in range(5)]

    def test_hits_retain_data_and_return_unaliased_matrices(self):
        cache = cache_module.RegionalCache()
        with patch.object(cache_module.regional.sc,'G',4), \
                patch.object(cache_module.regional,'local_operators',return_value='local') as local, \
                patch.object(cache_module.regional,'placement',side_effect=lambda _:self.polynomial()):
            first = cache.polynomial(self.model,self.witness)
            first[0][0,0] = arb(99)
            second = cache.polynomial(self.model,self.witness)
            self.assertEqual(second[0][0,0],1)
            self.assertEqual(local.call_count,1)
            self.assertEqual((cache.hits,cache.misses),(1,1))
            self.assertIs(next(iter(cache.entries.values()))[0],self.model.data)

    def test_every_local_selector_precision_tilt_and_data_affect_key(self):
        base = cache_module.local_key(self.model,self.witness)
        variants = [dict(self.witness,regional_exact_zero=False),
            dict(self.witness,regional_joint_return_through=3),
            dict(self.witness,regional_lazy_density_through=6),
            dict(self.witness,regional_feedback_classes_from=3,regional_feedback_classes_through=32),
            dict(self.witness,regional_feedback_classes_from=4,regional_feedback_classes_through=31),
            dict(self.witness,regional_feedback_classes_from=3,regional_feedback_classes_through=32,
                 regional_feedback_uniform_classes=True),
            dict(self.witness,regional_feedback_classes_from=3,regional_feedback_classes_through=32,
                 regional_feedback_uniform_classes=True,regional_feedback_uniform_replace=True),
            dict(self.witness,parameters=['1/11','0','0'])]
        keys = [cache_module.local_key(self.model,v) for v in variants]
        self.assertEqual(len(set([base,*keys])),len(keys)+1)
        ctx.prec = 384
        self.assertNotEqual(cache_module.local_key(self.model,self.witness),base)
        ctx.prec = 256
        self.assertNotEqual(cache_module.local_key(SimpleNamespace(data=dict(self.model.data)),self.witness),base)
        nonlocal_change = dict(self.witness,parameters=['1/10','23','17'],
            regional_direct_counts=True,regional_count_parts=[{'unrelated':'proposal'}])
        self.assertEqual(cache_module.local_key(self.model,nonlocal_change),base)

    def test_reject_unknown_selector_and_invalid_local_values(self):
        for changes in (dict(regional_future_kernel=True),dict(regional_exact_zero=1),
                dict(regional_joint_return_through=5),dict(regional_lazy_density_through=True),
                dict(regional_feedback_uniform_replace=True),dict(parameters=['0','0','0'])):
            with self.subTest(changes=changes),self.assertRaises(ValueError):
                cache_module.local_key(self.model,dict(self.witness,**changes))

    def test_original_body_rechecks_cell_and_duals_and_cannot_mutate_cache(self):
        cache = cache_module.RegionalCache()
        with patch.object(cache_module.regional.sc,'G',4), \
                patch.object(cache_module.regional,'local_operators',return_value='local') as local, \
                patch.object(cache_module.regional,'placement',side_effect=lambda _:self.polynomial()):
            first,_ = cache.evaluate(toy_outward,self.model,(Q(1,8),Q(1,4)),self.witness)
            other = dict(self.witness,parameters=['1/10','2','0'])
            second,_ = cache.evaluate(toy_outward,self.model,(Q(1,4),Q(1,3)),other)
            self.assertEqual(first,arb(9)/8)
            self.assertEqual(second,arb(13)/4)
            self.assertEqual(local.call_count,1)
            self.assertEqual((cache.hits,cache.misses),(1,1))

    def test_context_restores_original_even_on_exception(self):
        original = cache_module.regional.outward
        with self.assertRaises(RuntimeError):
            with cache_module.install() as cache:
                self.assertIsNot(cache_module.regional.outward,original)
                with self.assertRaises(ValueError):
                    with cache_module.install(): pass
                raise RuntimeError('stop')
        self.assertIs(cache_module.regional.outward,original)

    def test_original_regional_body_agrees_with_uncached_on_two_cells(self):
        regional = cache_module.regional
        sc = regional.sc
        local = [arb_mat([[1,0],[0,1]]), arb_mat([[1,1],[0,1]]), arb_mat([[1,0],[1,1]])]
        ordered = regional.placement(local,2)
        components = [(str(i),Q(i+1),sc.probabilities(p),int(i!=0))
                      for i,p in enumerate(map(Q,('0','1/4','3/4','1')))]
        with patch.object(sc,'G',4),patch.object(sc,'REGIONS',2),patch.object(sc,'PACKETS',8):
            model = sc.Model(components,{'windows':32},0,q_min=2,tilt=Q(2,3),
                inner=birth_classes,variance_shuffle=True,variance_bins=2,regional_count=True)
            cell = (Q(1,4),Q(3,4)); child = (Q(1,4),Q(1,2))
            witness = dict(tilt=str(model.tilt),parameters=['1','0','0'],
                variance_dual=list(map(str,model.shuffle_dual(cell))),regional_direct_counts=True,
                variance_partition=[dict(interval=['0','1/8'],dual=['0','0','0']),
                                    dict(interval=['1/8','1/4'],dual=['0','0','0'])])
            with patch.object(regional,'local_operators',return_value=local) as build, \
                    patch.object(regional,'placement',return_value=ordered):
                first,checked = regional.outward(model,cell,witness)
                second,_ = regional.outward(model,child,checked)
                build.reset_mock()
                with cache_module.install() as cache:
                    self.assertEqual(regional.outward(model,cell,checked)[0],first)
                    self.assertEqual(regional.outward(model,child,checked)[0],second)
                    self.assertEqual((cache.hits,cache.misses),(1,1))
                self.assertEqual(build.call_count,1)

    def test_precision_mutation_and_wrong_geometry_are_rejected(self):
        def changed(_):
            ctx.prec = 192
            return self.polynomial()
        for placement in (changed,lambda _:[arb_mat([[1]])]):
            ctx.prec = 256
            with patch.object(cache_module.regional.sc,'G',4), \
                    patch.object(cache_module.regional,'local_operators',return_value='local'), \
                    patch.object(cache_module.regional,'placement',side_effect=placement), \
                    self.assertRaises(ArithmeticError):
                cache_module.RegionalCache().polynomial(self.model,self.witness)

    def test_lru_capacity_and_invalid_capacity(self):
        cache = cache_module.RegionalCache(1)
        with patch.object(cache_module.regional.sc,'G',4), \
                patch.object(cache_module.regional,'local_operators',return_value='local') as local, \
                patch.object(cache_module.regional,'placement',side_effect=lambda _:self.polynomial()):
            for tilt in ('1/10','1/11','1/10'):
                cache.polynomial(self.model,dict(self.witness,parameters=[tilt,'0','0']))
            self.assertEqual(len(cache.entries),1)
            self.assertEqual(local.call_count,3)
        for capacity in (0,33,True,1.5):
            with self.assertRaises(ValueError):cache_module.RegionalCache(capacity)


if __name__ == '__main__':
    unittest.main()
