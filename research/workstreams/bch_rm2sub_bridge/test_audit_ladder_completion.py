from fractions import Fraction as F
import unittest

import audit_ladder_completion as audit


class CompletionGateTests(unittest.TestCase):
    def example(self):
        spec = audit.identity.instance(20)
        a,b,d = F(2)**-41,F(2)**-42,F(2)**-42
        saved = dict(instance=spec,status='LADDER_FULL_DISTANCE_SETUP_CERTIFIED_40_BITS',complete_coverage=True,
                     covered_intervals=[[1,spec['rows']]],sparse_interval=[1,2],dense_interval=[3,spec['rows']],
                     q1_upper=audit.base.encode(a),sparse_rest_upper=audit.base.encode(b),dense_upper=audit.base.encode(d),
                     union_upper=audit.base.encode(a+b+d))
        return saved,{1:a,2:b},d

    def test_exact_target_is_allowed(self):
        saved,values,d = self.example()
        self.assertEqual(audit.full_fields(saved,20,values,d),F(2)**-40)

    def test_missing_occupancy_rejected(self):
        saved,values,d = self.example()
        del values[2]
        with self.assertRaises(ValueError):
            audit.full_fields(saved,20,values,d)

    def test_cross_size_rejected(self):
        saved,values,d = self.example()
        with self.assertRaises(ValueError):
            audit.full_fields(saved,22,values,d)

    def test_changed_cutoff_rejected(self):
        saved,values,d = self.example()
        saved['instance'] = dict(saved['instance'],cutoff=saved['instance']['cutoff']-1)
        with self.assertRaises(ValueError):
            audit.full_fields(saved,20,values,d)

    def test_weak_or_changed_union_rejected(self):
        saved,values,d = self.example()
        with self.assertRaises(ValueError):
            audit.full_fields(saved,20,values,2*d)


if __name__ == '__main__':
    unittest.main()
