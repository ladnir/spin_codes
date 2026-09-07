"""Exact optimization and reproducibility checks for the BCH ROI probes."""
from fractions import Fraction as F
import itertools
import unittest
import bridge as base
import low_shell_roi as roi
import low_shell_constraint_probe as probe


class LowShellRoiTest(unittest.TestCase):
    def test_knapsack_vertices(self):
        cost = {0:F(2),1:F(3),2:F(5)}
        caps = {0:F(4),1:F(3),2:F(2)}
        for values in itertools.product((0,1,4), repeat=3):
            value = {i:F(v) for i,v in enumerate(values)}
            for budget in map(F,(0,1,7,12,20,40)):
                score, x = roi.knapsack(cost,value,caps,budget)
                vertices = [F(0)]
                for free in cost:
                    others = [i for i in cost if i != free]
                    for sides in itertools.product((0,1), repeat=2):
                        y = {i:caps[i]*s for i,s in zip(others,sides)}
                        left = budget-sum(cost[i]*y[i] for i in others)
                        if left >= 0:
                            y[free] = min(caps[free],left/cost[free])
                            vertices.append(sum(value[i]*y[i] for i in cost))
                self.assertEqual(score,max(vertices))

    def test_exact_roi_replay(self):
        saved = base.read(base.HERE/'generated/low_shell_roi_v1.json')
        self.assertEqual(roi.build(),saved)
        self.assertLess(saved['fixed_coefficient_reoptimization_gain_ceiling_bits'],.003)
        self.assertGreater(saved['primal_shells'][0]['low_shell_score_fraction'],.999)
        for weights in ([38],[38,40],[38,40,42]):
            rows = [r for r in saved['sensitivity'] if r['tightened_weights']==weights]
            bounds = [base.decode(r['conditional_full_upper']) for r in rows]
            self.assertEqual(bounds,sorted(bounds,reverse=True))
        zero_tests = {(r['code'],r['weight']):r['hypothetical_zero_excludes_primal'] for r in saved['moment_screen']}
        self.assertTrue(zero_tests['Pdual',30])
        self.assertFalse(zero_tests['Pdual',32])

    def test_candidate_receipt_replay(self):
        folders = sorted((base.HERE/'generated/low_shell_constraints_v1').glob('*/audit.json'))
        self.assertTrue(folders)
        for path in folders:
            with self.subTest(candidate=path.parent.name):
                model,scales,objective,norm,lp,_,_ = probe.prepare(path.parent.name)
                self.assertEqual((path.parent/'h_38.lp').read_text(),lp)
                result = probe.audit(path.parent,model,scales,objective,norm)
                self.assertEqual(result,base.read(path))
                self.assertEqual(result['status'],'EXACT_LP_IMPLICATION_ONLY_ADDITIONAL_CODE_FACT_NOT_PROVED')


if __name__ == '__main__':
    unittest.main()
