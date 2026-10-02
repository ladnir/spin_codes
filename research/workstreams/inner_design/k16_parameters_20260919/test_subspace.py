from fractions import Fraction as F
import unittest
import numpy as np
from flint import ctx
import model
import subspace
import refine


class Subspaces(unittest.TestCase):
    def test_nullspace(self):
        rng=__import__('random').Random(1)
        for width in range(2,10):
            for _ in range(20):
                constraints=[rng.randrange(1,1 << width) for _ in range(width//2)]
                basis=subspace.nullspace(constraints,width)
                rank=model.base.independent.search.rank(constraints)
                self.assertEqual(len(basis),width-rank)
                self.assertEqual(model.base.independent.search.rank(basis),len(basis))
                self.assertTrue(all(not ((a&b).bit_count()&1) for a in basis for b in constraints))

    def test_search_and_exact_enumeration(self):
        for s in (12,13):
            record=subspace.search(s)
            self.assertIsNotNone(record)
            self.assertEqual(record['feedback_columns'],model.study.core.inner(64,s)['feedback_columns'])
            e=subspace.Engine(record)
            self.assertGreaterEqual(min(e.spectrum),24)
            self.assertLessEqual(max(e.spectrum),40)
            weights=model.base.independent.search.wm.all_weights(e.a_columns,s)
            self.assertEqual(sum(weights==0),1)
            self.assertEqual(int(weights[1:].min()),24)
            self.assertEqual(int(weights[1:].max()),40)
            # Exact map embedding, not merely an equal weight histogram.
            parent=model.study.core.inner(64,15)['expansion_columns']
            for q in range(1 << s):
                old=0
                for j,v in enumerate(record['parent_basis']):
                    if q >> j & 1: old^=v
                self.assertEqual([(c&q).bit_count()&1 for c in e.a_columns],
                                 [(c&old).bit_count()&1 for c in parent])

    def test_control_adapter(self):
        ctx.prec=256
        e=subspace.Engine(model.study.core.inner(64,12))
        old=refine.Engine(64,12)
        self.assertEqual(e.epoch(F(-2)),old.epoch(F(-2)))
        self.assertEqual(e.region(F(-2),2),old.region(F(-2),2))


if __name__=='__main__': unittest.main()
