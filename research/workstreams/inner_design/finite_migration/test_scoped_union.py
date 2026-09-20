import copy
import math
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

import parameter_scoped_union as subject


class ScopedUnionTests(unittest.TestCase):
    def fixtures(self):
        geometry=[128,64,20,20]
        inner={'test_map':1}
        prefix=dict(status='BINARY64_IMT_RANGE_BRIDGED_FULL_PARAMETER_BOUND',geometry=geometry,
                    inner=inner,through=320,covered_occupancies=[1,16384],
                    inputs=dict(q1=['q1','q1v'],compositions=['small','smallv'],
                                ranges=[['range1','range1v'],['range2','range2v'],['range3','range3v']]))
        dense=dict(status=subject.scoped.STATUS,geometry=geometry,inner=inner,
                   dense=dict(occupation_min=321,occupation_max=16384,margin_bits=60.,
                              log_union_upper=-60*math.log(2)))
        q1=dict(maps={'t64_s20':inner},rows=[dict(b=128,t=64,s=20,exponent=20,q1_margin_bits=32.)])
        small=dict(geometry=geometry,inner=inner,results=[dict(q=q,log_bound=-1000.) for q in (2,3,4)])
        pairs={'prefix':(prefix,dict(geometry=geometry,covered_occupancies=[1,16384]),{}),
               'dense':(dense,dict(geometry=geometry,covered_occupancies=[321,16384],margin_bits=60.),{}),
               'q1':(q1,dict(geometries_checked=130),{}),
               'small':(small,dict(geometry=geometry,covered_occupancies=[2,4]),{})}
        for i,(first,last) in enumerate(((5,64),(65,256),(257,320)),1):
            pairs['range'+str(i)]=(dict(geometry=geometry,inner=inner,first=first,last=last,
                                      log_bounds=[-1000.]*(last-first+1)),
                                  dict(geometry=geometry,covered_occupancies=[first,last]),{})
        return pairs

    def invoke(self,pairs):
        bank=subject.bridge.bank
        with tempfile.TemporaryDirectory() as folder, \
             patch.object(subject.bridge.union,'pair',side_effect=lambda a,b,status:copy.deepcopy(pairs[Path(a).name])), \
             patch.object(bank.base.grid.ladder.candidate,'sources',return_value={}), \
             patch.object(bank.model.base,'write_new') as write:
            subject.combine([Path('prefix'),Path('prefixv')],[Path('dense'),Path('densev')],Path(folder)/'result.json')
            return write.call_args.args[1]

    def test_exact_components_are_combined(self):
        result=self.invoke(self.fixtures())
        expected=-math.log2(2.**-32+2.**-60)
        self.assertAlmostEqual(result['full_margin_bits'],expected,places=12)
        self.assertEqual(result['covered_occupancies'],[1,16384])
        self.assertTrue(result['full_bound_useful'])
        self.assertFalse(result['full_distance_proved'])

    def test_gap_rejected(self):
        pairs=self.fixtures()
        pairs['range2'][0]['first']=66
        pairs['range2'][1]['covered_occupancies'][0]=66
        with self.assertRaises(AssertionError):
            self.invoke(pairs)

    def test_suffix_overlap_rejected(self):
        pairs=self.fixtures()
        pairs['dense'][0]['dense']['occupation_min']=320
        pairs['dense'][1]['covered_occupancies'][0]=320
        with self.assertRaises(AssertionError):
            self.invoke(pairs)


if __name__=='__main__':
    unittest.main()
