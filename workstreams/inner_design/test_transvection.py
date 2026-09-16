from collections import Counter
from fractions import Fraction as F
import unittest

import transvection as tv


class TransvectionTests(unittest.TestCase):
    def test_exact_single_state_kernel_and_involution(self):
        for s in range(2, 6):
            m = (1 << s)-1
            total = m*(1 << (s-1))
            for q in range(1, m+1):
                counts = tv.transition_counts(s, q)
                self.assertNotIn(0, counts)
                for y in range(1, m+1):
                    self.assertEqual(F(counts[y], total), F(y == q,2)+F(1,2*m))
            for u in range(1,m+1):
                for v in range(m+1):
                    if tv.dot(u,v): continue
                    for q in range(m+1):
                        self.assertEqual(tv.update(tv.update(q,u,v),u,v),q)

    def test_setup_sampler_is_uniform(self):
        for s in range(2,7):
            for u in range(1,1<<s):
                counts = Counter(tv.projected_v(u,w) for w in range(1<<s))
                self.assertEqual(set(counts.values()),{2})
                self.assertEqual(set(counts),{v for v in range(1<<s) if not tv.dot(u,v)})

    def test_round_products(self):
        for s in range(2,5):
            m = (1<<s)-1
            kernel = [[F(c, m*(1<<(s-1))) for _,c in sorted(tv.transition_counts(s,q).items())]
                      for q in range(1,m+1)]
            law = [F(y==1) for y in range(1,m+1)]
            for r in range(1,6):
                law = [sum(law[q]*kernel[q][y] for q in range(m)) for y in range(m)]
                self.assertEqual(law,[F(y==1,1<<r)+F((1<<r)-1,(1<<r)*m) for y in range(1,m+1)])


if __name__ == '__main__':
    unittest.main()
