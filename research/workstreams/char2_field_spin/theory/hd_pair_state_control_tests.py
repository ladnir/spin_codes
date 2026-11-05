import unittest
import numpy as np
from hd_pair_state_control import branch_data,joint_branches,combine_branches


class StateControlTests(unittest.TestCase):
    def test_all_labelled_source_counts(self):
        for mixer in ('H0','Hbin'):
            for arity,count in ((1,16),(2,120)):
                data=branch_data(arity,mixer=mixer)
                self.assertEqual(len(data['sources']),count)
                self.assertEqual(data['log_counts'].shape,(count,(2 if arity==1 else 5)**2))
                for source in data['sources']:
                    terms=joint_branches(source,mixer=mixer)
                    self.assertEqual(sum(n for _,_,n in terms),256**arity)
                    self.assertEqual(sum(n for w,_,n in terms if w==1),255*arity)
    def test_marked_branch_minkowski_at_r1(self):
        matrices=np.arange(1.,33.).reshape(2,4,2,2)
        counts=np.arange(1.,9.).reshape(2,4)
        activity=np.array(((0,1,1,2),(0,1,1,2)))
        weight=np.array(((0,1,2,3),(0,1,2,3)))
        got=combine_branches(np.log2(matrices),np.log2(counts),activity,weight,1,2,3)
        expected=np.sum(matrices*(counts*3.**activity*2.**weight)[...,None,None],axis=1)
        np.testing.assert_allclose(2.**got,expected,rtol=1e-14)
    def test_exact_zero_entries_stay_zero(self):
        matrices=np.full((1,2,2,2),-np.inf)
        matrices[:,:,0,0]=0.
        got=combine_branches(matrices,np.zeros((1,2)),np.zeros((1,2)),np.zeros((1,2)),18,1,1)
        self.assertEqual(got[0,0,0],18.)
        self.assertTrue(np.isneginf(got[0,1,1]))


if __name__=='__main__':unittest.main()
