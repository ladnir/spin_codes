import itertools
import random
import unittest
from fractions import Fraction as F

import cheaper_outer_ensemble_screen as screen


def mul(a,b):
    result = 0
    while b:
        if b&1:
            result ^= a
        a <<= 1
        if a&4:
            a ^= 7
        b >>= 1
    return result


class CheaperOuterTests(unittest.TestCase):
    def test_window_thresholds(self):
        for h,k in ((16,11),(32,22),(64,44),(26,1)):
            result = screen.window_gate(6016,h,k)
            self.assertTrue(result['deterministically_bad'])
            self.assertGreater(result['shortened_dimension_lower'],result['exit_rank_upper'])
            self.assertLessEqual(result['budget'],154009)
            useful = screen.window_gate(6016,h,k,discard_tail_rows=59)
            self.assertTrue(useful['deterministically_bad'])
            self.assertGreater(useful['shortened_dimension_lower'],useful['exit_rank_upper'])
        self.assertLessEqual((6016-59)*176,2**20)
        self.assertFalse(screen.window_gate(6016,128,88)['deterministically_bad'])
        self.assertFalse(screen.window_gate(6016,27,1)['deterministically_bad'])

    def test_realized_binary_confinement(self):
        # Four independent [3,2] parity rows; the first3 of12 regions carry them.
        # Eight retained slots leave >3 dimensions before the three exit checks.
        rng = random.Random(7723)
        lengths = (3,3,2)
        for _ in range(24):
            forms = [0]*48
            for region in range(3):
                permutation = list(range(4))
                rng.shuffle(permutation)
                for slot,row in enumerate(permutation):
                    forms[4*region+slot] = ((1 << (2*row)) if region==0 else
                        (1 << (2*row+1)) if region==1 else (3 << (2*row)))
            allowed = {4*r+j for r,w in enumerate(lengths) for j in range(w)}
            equations = [form for i,form in enumerate(forms) if i not in allowed]
            # Inner state follows x; output is x+previous state.
            equations += [forms[4*r+w-1] for r,w in enumerate(lengths)]
            witness = next(x for x in range(1,256)
                           if all((x&equation).bit_count()%2==0 for equation in equations))
            state = 0
            output = []
            for form in forms:
                x = (form&witness).bit_count()%2
                output.append(x^state)
                state = x
            self.assertTrue(any(output))
            self.assertLessEqual(sum(output),8)
            self.assertTrue(all(not output[i] for i in range(48) if i not in allowed))

    def test_one_product_common_flag(self):
        rng = random.Random(881)
        for symbol in range(1,4):
            for _ in range(80):
                values = [symbol]*8
                for bit in range(3):
                    for left in range(8):
                        right = left^(1 << bit)
                        if left < right:
                            u = values[left]^values[right]
                            values[right] ^= mul(rng.randrange(4),u)
                            values[left] = u
                self.assertEqual(values,[0]*7+[symbol])

    def test_two_shear_support_law_exhaustive_gf4(self):
        for x,y in itertools.product(range(4),repeat=2):
            actual = {}
            for a,b in itertools.product(range(1,4),repeat=2):
                u = x^mul(b,y)
                v = y^mul(a,u)
                mask = int(bool(u))+2*int(bool(v))
                actual[mask] = actual.get(mask,F(0))+F(1,9)
            self.assertEqual(actual,screen.shear_pair_law(int(bool(x))+2*int(bool(y)),4))
        for mask in (0,1,2,3,15,85,127,255):
            law = screen.butterfly_support_law(mask,4)
            self.assertEqual(sum(law.values()),1)
            self.assertEqual(0 in law,mask==0)

    def test_sparse_rank_mean_exhaustive_gf4(self):
        for density in (F(1,2),F(3,4)):
            expected = F(0)
            for entries in itertools.product(range(4),repeat=4):
                probability = F(1)
                for value in entries:
                    probability *= density/3 if value else 1-density
                count = 0
                for x,y in itertools.product(range(4),repeat=2):
                    count += (mul(entries[0],x)^mul(entries[1],y)==0
                              and mul(entries[2],x)^mul(entries[3],y)==0)
                expected += probability*count
            self.assertEqual(expected,screen.sparse_kernel_mean(2,2,4,density))
        self.assertEqual(screen.sparse_basis_failure_lower(F(0)),1)
        self.assertEqual(screen.sparse_basis_failure_lower(F(1)),0)
        self.assertGreater(screen.sparse_basis_failure_lower(F(73,100)),F(1,2**40))
        self.assertLess(screen.sparse_basis_failure_lower(F(74,100)),F(1,2**40))

    def test_gl_kernel_and_outer_gate(self):
        self.assertEqual(screen.gl_postmix_kernel_mean([0]*32),1)
        self.assertEqual(screen.gl_postmix_kernel_mean([8]*32),256**176)
        self.assertEqual(screen.gl_postmix_kernel_mean([8]*11+[0]*21),256**8)
        self.assertGreaterEqual(screen.integer_bits(screen.gl_postmix_lowweight_upper()),40)

    def test_diagonal_binary_postmix_pair_upper_exhaustive_gf4(self):
        # Monic irreducible quadratics have no roots in the evaluation field.
        polynomials = [(a,b) for a,b in itertools.product(range(4),repeat=2)
                       if all(mul(x,x)^mul(a,x)^b for x in range(4))]
        self.assertEqual(len(polynomials),6)
        inverse = {a:next(b for b in range(1,4) if mul(a,b)==1) for a in range(1,4)}
        ratios = []
        for f,g in itertools.permutations(polynomials,2):
            values = []
            for x in range(4):
                numerator = mul(x,x)^mul(f[0],x)^f[1]
                denominator = mul(x,x)^mul(g[0],x)^g[1]
                values.append(mul(numerator,inverse[denominator]))
            ratios.append(values)
        supports = list(itertools.combinations(range(4),3))
        first,pairs,union = F(0),F(0),F(0)
        for diagonal in itertools.product(range(1,4),repeat=4):
            event_count = 0
            for support in supports:
                for ratio in ratios:
                    scales = {mul(diagonal[x],inverse[ratio[x]]) for x in support}
                    event_count += len(scales)==1
            first += F(event_count,81)
            pairs += F(event_count*(event_count-1),162)
            union += F(event_count>0,81)
        bound = screen.diagonal_binary_postmix_lower(q=4,length=4,support=3)
        self.assertEqual(first,bound['first_moment'])
        self.assertLessEqual(pairs,bound['pair_upper'])
        self.assertLessEqual(bound['failure_lower'],union)

    def test_diagonal_binary_postmix_actual_threshold(self):
        bound = screen.diagonal_binary_postmix_lower()
        self.assertEqual(bound['irreducibles'],32640)
        self.assertGreater(bound['failure_lower'],F(1,2**40))
        self.assertLess(bound['failure_lower'],F(1,2**38))
        self.assertLess(bound['pair_upper'],bound['first_moment'])


if __name__ == '__main__':
    unittest.main()
