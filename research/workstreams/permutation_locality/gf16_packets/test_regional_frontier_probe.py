"""Diagnostics must reconstruct the selected bound before ablating it."""
from fractions import Fraction as Q
from types import SimpleNamespace
import unittest
from unittest.mock import patch
import numpy as np
from flint import arb,arb_mat
import regional_frontier_probe as probe


class RegionalFrontierTests(unittest.TestCase):
    def test_direct_mass_and_reference_ratio_have_distinct_normalizers(self):
        model=SimpleNamespace(tilt=Q(1,2),features=[Q(0),Q(1,2)],active=[0,1],q_min=1,
            weights=lambda *args:([Q(3,4),Q(1,2)],None),
            family=lambda *args:(None,None,np.log([1.,2.])))
        cell=(Q(1,4),Q(1,3));parts=[((Q(0),Q(1,4)),(Q(0),Q(0),Q(0)))]
        checked={'regional_count_parts':[{'mgf_witnesses':[]}]}
        with patch.object(probe.sc,'G',4),patch.object(probe.variance,'factor',return_value=Q(2)), \
                patch.object(probe.regional,'count_mass_caps',return_value=[arb(1)/2]*5) as direct, \
                patch.object(probe.regional,'count_ratios',return_value=[arb(2)]*5) as ratios:
            witness={'variance_dual':[],'regional_direct_counts':True}
            rows,without,scale=probe.weighted_parts(model,cell,witness,parts,checked)
            self.assertEqual(scale,1);self.assertEqual(without,[])
            direct.assert_called_once();ratios.assert_not_called()
            np.testing.assert_allclose(rows[0][0],np.arange(5)*np.log(2)-np.log(2))
            witness['regional_direct_counts']=False
            rows,without,scale=probe.weighted_parts(model,cell,witness,parts,checked)
            self.assertEqual(scale,Q(5,4));self.assertEqual(len(without),1)
            ratios.assert_called_once()
            np.testing.assert_allclose(rows[0][0]-without[0][0],np.log(2))

    def test_probe_reconstructs_saved_local_refinements_and_checks_score(self):
        witness={'tilt':'1/2','parameters':['1/8','0','0'],'regional_count_parts':[{}],
            'regional_direct_counts':True,'regional_lazy_density_through':6,'regional_joint_return_through':3}
        row={'cell':['1/4','1/3'],'witness':witness,'proposal':7.,'activity':'2/5'}
        data={'bits':3};source={'threshold':10,'minimum_groups':1}
        model=SimpleNamespace(threshold=10)
        with patch.object(probe.sc,'actual_components',return_value=[]), \
                patch.object(probe.sc,'Model',return_value=model), \
                patch.object(probe.regional,'prepare_witness',return_value=([],witness)), \
                patch.object(probe.regional,'local_operators',return_value=[arb_mat(np.eye(3,dtype=int).tolist())]) as local, \
                patch.object(probe.occupancy,'outward',return_value=[arb_mat(np.eye(3,dtype=int).tolist())]), \
                patch.object(probe,'weighted_parts',return_value=([([],0)],[],Q(1))), \
                patch.object(probe,'evaluate',return_value=7.):
            result=probe.probe(source,row,data)
            local.assert_called_once_with(model,witness)
            self.assertEqual(result['scores']['baseline'],7.)
            self.assertNotIn('no_lazy_returns',result['scores'])
            self.assertNotIn('count_ratios_replaced_by_one',result['scores'])
            self.assertIn('no_returns',result['scores'])
            with self.assertRaisesRegex(ArithmeticError,'does not match'):
                probe.probe(source,dict(row,proposal=8.),data)

    def test_alternative_allocations_restart_from_unrefined_operators(self):
        witness={'tilt':'1/2','parameters':['1/8','0','0'],'regional_count_parts':[{}],
            'regional_direct_counts':True,'regional_lazy_density_through':6}
        row={'cell':['1/4','1/3'],'witness':witness,'proposal':7.,'activity':'2/5'}
        census=object();caps=object()
        data={'bits':3,'joint_only':True,'joint_return_census':census,'lazy_density_through':[8]}
        source={'threshold':10,'minimum_groups':1};model=SimpleNamespace(threshold=10)
        original=[arb_mat(np.eye(3,dtype=int).tolist())]
        refined=[original[0]*2];joint=[original[0]*3]
        with patch.object(probe.sc,'actual_components',return_value=[]), \
                patch.object(probe.sc,'Model',return_value=model), \
                patch.object(probe.regional,'prepare_witness',return_value=([],witness)), \
                patch.object(probe.regional,'local_operators',return_value=refined), \
                patch.object(probe.occupancy,'outward',return_value=original), \
                patch.object(probe,'weighted_parts',return_value=([([],0)],[],Q(1))), \
                patch.object(probe,'evaluate',return_value=7.), \
                patch('fixed_joint_return_probe.refine',return_value=joint) as returns, \
                patch('lazy_density_probe.density_caps',return_value=caps), \
                patch('lazy_density_probe.candidate',return_value=joint) as density:
            result=probe.probe(source,row,data)
            self.assertIs(returns.call_args.args[1],original)
            self.assertIs(returns.call_args.args[2],census)
            self.assertIs(density.call_args_list[0].args[1],original)
            self.assertIs(density.call_args_list[1].args[1],joint)
            self.assertIn('joint_and_lazy_through_8',result['scores'])


if __name__=='__main__':unittest.main()
