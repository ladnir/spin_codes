import copy
from pathlib import Path
from types import SimpleNamespace
import unittest
from unittest.mock import patch
from contextlib import nullcontext

import parameter_scoped_cover as subject
import type_box_coverage


class ScopedCoverTests(unittest.TestCase):
    def record(self):
        box = dict(lower=[0,0,0],upper=[16384-321,16384,16384],witness={},own_log_bound=0.)
        count = type_box_coverage.lattice_count(box['lower'],box['upper'],16384)
        return dict(status=subject.STATUS,geometry=[128,64,20,20],inner={'test_map':1},
                    dense=dict(occupation_min=321,occupation_max=16384,selected_boxes=[box],
                               integer_types_checked=str(count),log_union_upper=0.,margin_bits=0.))

    def test_geometry_and_partition_rejections(self):
        record = self.record()
        moments = SimpleNamespace(original=SimpleNamespace(record=record['inner']))
        self.assertGreater(subject.check(record,moments),0)
        for field,value in (('geometry',[128,64,19,20]),('inner',{'test_map':2})):
            changed = copy.deepcopy(record)
            changed[field]=value
            with self.assertRaises(AssertionError):
                subject.check(changed,moments)
        changed=copy.deepcopy(record)
        changed['dense']['selected_boxes'] *= 2
        with self.assertRaises(ValueError):
            subject.check(changed,moments)
        changed=copy.deepcopy(record)
        changed['dense']['occupation_min']=322
        with self.assertRaises(ValueError):
            subject.check(changed,moments)

    def test_replay_rejects_changed_numeric_claim(self):
        record=self.record()
        moments=SimpleNamespace(original=SimpleNamespace(record=record['inner']),bound=lambda box:1.)
        with patch.object(subject.bank.model.base,'read',return_value=record), \
             patch.object(subject.bank.model,'authenticate'), \
             patch.object(subject.bank.previous.maps,'use',side_effect=nullcontext), \
             patch.object(subject.bridge,'Moments',return_value=moments), \
             patch.object(subject.bank.model.base,'write_new') as write:
            with self.assertRaises(AssertionError):
                subject.verify(Path('test-receipt.json'),Path('unused-output.json'))
            write.assert_not_called()


if __name__=='__main__':
    unittest.main()
