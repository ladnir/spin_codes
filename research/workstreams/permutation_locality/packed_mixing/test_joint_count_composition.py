"""Small synthetic composition tests, never production certificate evidence."""
import copy
from contextlib import ExitStack
from fractions import Fraction as Q
from pathlib import Path
import sys
from types import SimpleNamespace
import unittest
from unittest.mock import patch

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
sys.path.insert(0, str(Path(__file__).resolve().parents[1]/'gf16_packets'))

from flint import arb, arb_mat, ctx
import joint_count_dual as dual
import joint_count_region as integration
import regional_count as rc
import scalar_cover as sc
import variance_partition as variance


class JointCompositionTests(unittest.TestCase):
    def setUp(self):
        old_precision = ctx.prec
        self.addCleanup(setattr, ctx, 'prec', old_precision)
        ctx.prec = 256
        self.cell = (Q(1, 2), Q(1, 2))
        self.model = SimpleNamespace(tilt=Q(3, 16), threshold=0, q_min=2,
            features=(Q(1, 4), Q(3, 4)), active=(1, 1),
            family=lambda tilt: ((Q(1, 2), Q(1, 2)), None, None))
        self.region = [arb_mat([[1, 0], [0, 1]]),
            arb_mat([[rc.aq(Q(1, 2)), rc.aq(Q(1, 4))],
                     [rc.aq(Q(1, 3)), rc.aq(Q(1, 4))]]),
            arb_mat([[rc.aq(Q(1, 8)), rc.aq(Q(1, 2))],
                     [rc.aq(Q(1, 4)), rc.aq(Q(1, 8))]])]
        # This is the exact count law for two independent Bernoulli variables
        # with probabilities 1/4 and 3/4; it is normalized BEFORE tau^-j.
        self.mass = (Q(3, 16), Q(5, 8), Q(3, 16))
        self.witness = dict(regional_direct_counts=True,
            parameters=['1/2', '0', '0'], variance_dual=['0', '0'],
            variance_partition=[
                dict(interval=['0', '1/8'], dual=['0', '0', '0']),
                dict(interval=['1/8', '1/4'], dual=['0', '0', '0'])])
        self.witness['regional_count_parts'] = [dict(mgf_witnesses=[
            dict(tilt='-1/2', tilted_atom=True),
            dict(tilt='1/2', tilted_atom=True)]) for _ in range(2)]
        self.mgf_calls = []
        stack = self.enterContext(ExitStack())
        stack.enter_context(patch.object(sc, 'G', 2))
        stack.enter_context(patch.object(sc, 'REGIONS', 3))
        stack.enter_context(patch.object(rc, 'prepare_witness', self.prepare))
        stack.enter_context(patch.object(rc, 'local_operators', return_value=['synthetic']))
        stack.enter_context(patch.object(rc, 'placement', return_value=self.region))
        stack.enter_context(patch.object(variance, 'factor', return_value=Q(1)))
        stack.enter_context(patch.object(rc, 'count_mass_caps', return_value=[arb(1)]*3))
        stack.enter_context(patch.object(rc, 'mgf_upper', self.raw_mgf))

    def prepare(self, model, cell, witness):
        # Keep actual partition validation while mocking only the expensive
        # production geometry and local operators.
        return variance.validate(cell, witness['variance_partition']), copy.deepcopy(witness)

    def raw_mgf(self, features, active, cell, interval, q_min, groups, witness):
        self.mgf_calls.append((tuple(cell), tuple(interval), q_min, groups, witness['tilt']))
        tilt = Q(witness['tilt'])
        return rc.up(sum((rc.aq(p)*rc.aq(tilt*j).exp()
                          for j, p in enumerate(self.mass)), arb(0)).log())

    def actual(self):
        matrix = sum((rc.aq(p)*rc.aq(self.model.tilt)**(-j)*region
                      for j, (p, region) in enumerate(zip(self.mass, self.region))), arb_mat(2, 2))
        power = matrix**3
        return sum((power[0, j] for j in range(2)), arb(0))

    def test_joint_composition_and_fresh_saved_dual_replay(self):
        upper, checked = integration.outward(self.model, self.cell, self.witness,
                                             selection=[1], moments=True)
        self.assertGreaterEqual(upper, self.actual())
        self.assertEqual(len(checked['parts']), 2)
        self.assertEqual([p['joint'] for p in checked['parts']], [False, True])
        self.assertIsNone(checked['joint_duals'][0])
        self.assertEqual({call[1] for call in self.mgf_calls}, {(Q(1, 8), Q(1, 4))})
        with patch.object(dual, 'linprog', side_effect=AssertionError('replay must not optimize')):
            replayed, replay_checked = integration.outward(self.model, self.cell,
                checked['base_witness'], selection=[1], moments=True,
                replay=checked['joint_duals'])
        self.assertEqual(replayed, upper)
        self.assertEqual(replay_checked['joint_duals'], checked['joint_duals'])
        baseline, _ = integration.outward(self.model, self.cell, self.witness, selection=[])
        self.assertLess(upper, baseline)

    def test_baseline_still_includes_every_variance_part(self):
        upper, checked = integration.outward(self.model, self.cell, self.witness, selection=[])
        matrix = sum((rc.aq(self.model.tilt)**(-j)*region
                      for j, region in enumerate(self.region)), arb_mat(2, 2))
        power = matrix**3
        exact = 2*sum((power[0, j] for j in range(2)), arb(0))
        self.assertGreaterEqual(upper, exact.lower())
        self.assertLess(float(upper-exact), 1e-50)
        self.assertEqual(checked['joint_duals'], [None, None])
        self.assertEqual(self.mgf_calls, [])

    def test_bad_dual_and_replay_geometry_fail_closed(self):
        _, checked = integration.outward(self.model, self.cell, self.witness,
                                         selection=[1], moments=True)
        damaged = copy.deepcopy(checked['joint_duals'])
        damaged[1][0][0]['beta'][0] = '-1'
        with self.assertRaises(ValueError):
            integration.outward(self.model, self.cell, self.witness,
                selection=[1], moments=True, replay=damaged)
        for wrong in (checked['joint_duals'][:1], [None, None],
                      [None, checked['joint_duals'][1][:1]]):
            with self.assertRaises(ValueError):
                integration.outward(self.model, self.cell, self.witness,
                    selection=[1], moments=True, replay=wrong)

    def test_missing_variance_interval_is_rejected(self):
        damaged = copy.deepcopy(self.witness)
        damaged['variance_partition'] = damaged['variance_partition'][1:]
        with self.assertRaises(ValueError):
            integration.outward(self.model, self.cell, damaged, selection=[])


if __name__ == '__main__':
    unittest.main()
