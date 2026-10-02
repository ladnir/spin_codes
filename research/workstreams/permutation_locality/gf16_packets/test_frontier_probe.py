from fractions import Fraction as Q
from types import SimpleNamespace
from unittest.mock import patch
import unittest
import numpy as np
from scipy.special import logsumexp
import frontier_probe as diagnostic
import scalar_cover as sc


class FrontierDiagnosticTests(unittest.TestCase):
    def test_all_nonzero_classes_are_modified_without_mutation(self):
        matrix=np.arange(1,65,dtype=float).reshape(8,8)/64
        original=matrix.copy()
        cases=diagnostic.counterfactual_matrices(matrix,{'bits':3})
        np.testing.assert_array_equal(matrix,original)
        np.testing.assert_array_equal(cases['no_returns'][1:,0],np.zeros(7))
        np.testing.assert_array_equal(cases['no_returns'][0],matrix[0])
        np.testing.assert_array_equal(cases['no_lazy_returns'][1:,0],
            np.minimum(matrix[1:,0],matrix[1:,2]/7))
        np.testing.assert_array_equal(cases['no_lazy_mass_to_arbitrary'][1:,1],np.zeros(7))

    def test_variance_partition_decomposition(self):
        model=SimpleNamespace(features=[Q(1,4),Q(3,4)],active=[False,True],q_min=3,
            family=lambda tilt:(None,None,np.log([2.,3.])))
        cell=(Q(1,2),Q(1,2))
        witness={'tilt':'1/8','parameters':['1/10','0','0'],'variance_dual':['0','0']}
        with patch.object(model,'comparison_loss',create=True,return_value=Q(2)):
            count,weighted=diagnostic.outer_logs(model,cell,witness)
            marked,_=diagnostic.outer_logs(model,cell,witness,activity_penalty=Q(1,2))
        self.assertAlmostEqual(count,sc.G*np.log(5),places=10)
        self.assertAlmostEqual(marked,sc.G*np.log(2+3*np.exp(-.5)),places=10)
        self.assertAlmostEqual(weighted-count,sc.REGIONS*np.log(2),places=10)
        witness['variance_partition']=[{'interval':['0','1/8'],'dual':['0','0','0']},
            {'interval':['1/8','1/4'],'dual':['0','0','0']}]
        with patch('variance_partition.factor',side_effect=[Q(2),Q(3)]):
            count,weighted=diagnostic.outer_logs(model,cell,witness)
        self.assertAlmostEqual(count,sc.G*np.log(5)+np.log(2),places=10)
        expected=sc.G*np.log(5)+logsumexp(sc.REGIONS*np.log([2,3]))
        self.assertAlmostEqual(weighted,expected,places=10)


if __name__=='__main__':unittest.main()
