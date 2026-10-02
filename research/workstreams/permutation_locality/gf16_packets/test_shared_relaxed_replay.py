from fractions import Fraction as Q
from types import SimpleNamespace
import unittest
from unittest.mock import patch

from flint import arb, ctx
import shared_relaxed_replay as replay


class SharedFreshReplayTests(unittest.TestCase):
    def test_worker_count_rejected_before_model_build(self):
        with patch.object(replay.proof, 'build_model') as build:
            for count in (0, 9, True, 1.5):
                with self.assertRaises(ValueError):
                    replay.replay_dense({}, workers=count, logdir='unused')
            build.assert_not_called()

    def model(self):
        return SimpleNamespace(root=(Q(0), Q(1)), threshold=104857,
            components=[(Q(1), Q(0), 0), (Q(100), Q(1, 2), 1)])

    def leaves(self):
        return {p: dict(witness=dict(parameters=['1/10', '0', '0'])) for p in ('0', '1')}

    def mapping(self, model, corrupt=None):
        def mapping(fn, jobs):
            for path, cell, witness in jobs:
                row = dict(path=path, cell=list(map(str, cell)), precision=384,
                    component_sha256=replay.search.component_digest(model.components),
                    upper=[1, -50], output_tilt=witness['parameters'][0])
                if corrupt:
                    row.update(corrupt)
                yield row
        return mapping

    def test_fresh_sum_and_full_partition(self):
        model = self.model()
        result = replay.replay(model, self.leaves(), self.mapping(model), 384)
        self.assertEqual(replay.proof.dyadic(result['upper']), Q(2)**-49)
        self.assertEqual(result['threshold'], 104857)
        self.assertEqual(result['leaves'], 2)
        self.assertEqual([r['path'] for r in result['checked']], ['0', '1'])
        with self.assertRaises(ValueError):
            replay.replay(model, {'0': self.leaves()['0']}, self.mapping(model), 384)

    def test_scope_precision_and_comparison_mismatch_rejected(self):
        model = self.model()
        for corrupt in (dict(path='bad'), dict(cell=['0', '1']), dict(precision=256),
                        dict(component_sha256='bad'), dict(output_tilt='1/20'), dict(upper=[0, 0])):
            with self.assertRaises(ArithmeticError):
                replay.replay(model, self.leaves(), self.mapping(model, corrupt), 384)

    def test_missing_or_extra_rows_rejected(self):
        model = self.model()
        with self.assertRaises(ArithmeticError):
            replay.replay(model, self.leaves(), lambda f, j: [], 384)
        def extra(fn, jobs):
            yield from self.mapping(model)(fn, jobs)
            yield dict(extra=True)
        with self.assertRaises(ArithmeticError):
            replay.replay(model, self.leaves(), extra, 384)

    def test_worker_calls_outward_at_requested_precision(self):
        class Model:
            def outward(self, cell, witness):
                self.seen = (cell, witness, ctx.prec)
                return arb(2)**-70
        model = Model()
        prior = ctx.prec
        try:
            witness = dict(parameters=['1/10', '0', '0'])
            row = replay.evaluate(model, ('', (Q(0), Q(1)), witness), 384, 'digest')
            self.assertEqual(model.seen, ((Q(0), Q(1)), witness, 384))
            self.assertEqual(row['upper'], [1, -70])
        finally:
            ctx.prec = prior

    def test_worker_rejects_precision_mutation_or_nonfinite_result(self):
        class Model:
            def outward(self, *args):
                ctx.prec = 128
                return arb(2)**-70
        prior = ctx.prec
        try:
            with self.assertRaises(ArithmeticError):
                replay.evaluate(Model(), ('', (Q(0), Q(1)), {}), 384, 'digest')
        finally:
            ctx.prec = prior
        for value in (arb(0), arb(-1), arb('inf'), arb('nan')):
            with self.assertRaises(ArithmeticError):
                replay.evaluate(SimpleNamespace(outward=lambda *a: value),
                    ('', (Q(0), Q(1)), {}), 384, 'digest')


if __name__ == '__main__':
    unittest.main()
