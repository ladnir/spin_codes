from fractions import Fraction as F
import unittest

import search_sparse_only as search


class SparseOnlyTests(unittest.TestCase):
    def test_q1_first(self):
        spec = dict(rows=2048, message_exponent=18)
        self.assertEqual(search.propose(spec, {}, set(), 479)['kind'], 'q1')

    def test_sparse_only_not_dense_endpoint(self):
        spec = dict(rows=2048, message_exponent=18)
        best = {1: F(2)**-52}
        tried = set()
        for _ in range(10):
            task = search.propose(spec, best, tried, 479)
            self.assertEqual(task['kind'], 'range')
            self.assertLessEqual(task['interval'][1], 479)
            tried.add(search.old.task_key(task))

    def test_complete_sparse_stops_without_claiming_full(self):
        spec = dict(rows=2048, message_exponent=18)
        best = {q: F(2)**-80 for q in range(1, 480)}
        self.assertIsNone(search.propose(spec, best, set(), 479))
        self.assertEqual(search.core.summarize(best, 2048, 40)['status'], 'UNRESOLVED')

    def test_timeout_gets_one_fresh_retry(self):
        task = dict(kind='range', interval=[457, 479], tilt=-10)
        attempt = dict(job='old_job', task=task, outcome=dict(status='TIMEOUT'))
        retry = search.interrupted_retry([attempt], set())
        self.assertEqual(retry['retry_of'], 'old_job')
        self.assertEqual(retry['interval'], task['interval'])
        self.assertIsNone(search.interrupted_retry([attempt], {search.old.task_key(retry)}))


if __name__ == '__main__':
    unittest.main()
