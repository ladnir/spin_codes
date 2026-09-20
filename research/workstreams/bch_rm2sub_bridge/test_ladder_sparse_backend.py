import unittest

import certificate_search_backend as old
import ladder_sparse_backend as ladder


class SparseLadderTests(unittest.TestCase):
    def test_old_range_result_unchanged(self):
        spec = ladder.identity.instance(16)
        task = dict(instance=spec,kind='range',interval=[2,3],tilt=-40)
        self.assertEqual(ladder.compute(task),old.compute(task))

    def test_new_sizes_producer_and_replay(self):
        for m in (22,24):
            spec = ladder.identity.instance(m)
            task = dict(instance=spec,kind='range',interval=[2,2],tilt=-70-7*(m-20))
            result = ladder.compute(task)
            replay = ladder.compute(task,512,result)
            self.assertEqual([r['occupation'] for r in result['bounds']],[2])
            self.assertLessEqual(ladder.base.decode(replay['bounds'][0]['upper']),
                                 ladder.base.decode(result['bounds'][0]['upper']))

    def test_wrong_instance_rejected(self):
        spec = dict(ladder.identity.instance(22),rows=8192)
        with self.assertRaises(ValueError):
            ladder.compute(dict(instance=spec,kind='range',interval=[2,2],tilt=-80))


if __name__ == '__main__':
    unittest.main()
