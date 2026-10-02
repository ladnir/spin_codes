from fractions import Fraction as Q
import unittest
from flint import arb,ctx
import kernel
import selected
from boundary_probe import composition,CASES,parameters,breakdown
from mixture import group_components,G


class BoundaryProbeTests(unittest.TestCase):
    def test_named_cases_have_exact_active_counts(self):
        rows=[(Q(1),p) for p in (Q(0),Q(1),Q(1,2),Q(2,5),Q(3,5))]
        components=group_components(rows)
        for case in CASES.values():
            counts=composition(components,case)
            self.assertEqual(sum(counts),G)
            self.assertEqual(sum(c*r[3] for c,r in zip(counts,components)),sum(case.values()))
        for case in ({'0000':1},{'0003':G+1},{'9999':1},{'0003':-1}):
            with self.assertRaises(ValueError):composition(components,case)

    def test_quantized_witness_is_positive(self):
        witness=parameters(dict(tilt=[.12,-1.,-2.,-3.,-4.]))
        self.assertEqual(witness[0],Q(3,25))
        self.assertTrue(all(v>0 for v in witness))

    def test_selected_density_replay_matches_float_decomposition(self):
        ctx.prec=192
        rows=[(Q(1),p) for p in (Q(0),Q(1),Q(1,2),Q(2,5),Q(3,5))]
        components=group_components(rows)
        counts=composition(components,{'0003':64,'1111':1})
        data=kernel.prepare(list(range(8)),[1,2,4,3]*32,3)
        witness=[Q(1,8),Q(1,8),Q(1,4),Q(1,2),Q(1)]
        for anchors in (None,(1980,50,10,6,2),(2000,48,0,0,0)):
            upper=selected.outward(data,components,counts,20971,witness,density_anchors=anchors)
            self.assertTrue(upper>0)
            parts=breakdown(data,components,counts,20971,witness,anchors)
            self.assertAlmostEqual(sum(parts.values()),float(upper.log()/arb(2).log()),places=6)
        exact=selected.outward(data,components,counts,20971,witness,exact_shuffle=True)
        self.assertAlmostEqual(sum(breakdown(data,components,counts,20971,witness,exact_shuffle=True).values()),
                               float(exact.log()/arb(2).log()),places=6)
        counts=composition(components,{'0003':32,'0034':32,'1111':1})
        posterior=selected.outward(data,components,counts,20971,witness,posterior_shuffle=True)
        self.assertAlmostEqual(sum(breakdown(data,components,counts,20971,witness,posterior_shuffle=True).values()),
                               float(posterior.log()/arb(2).log()),places=6)
        with self.assertRaises(ValueError):selected.outward(data,components,counts,20971,witness,posterior_shuffle=True,exact_shuffle=True)
        for anchors in ((1,2),(-1,0,0,0,0)):
            with self.assertRaises(ValueError):selected.outward(data,components,counts,20971,witness,density_anchors=anchors)
        with self.assertRaises(ValueError):selected.outward(data,components,counts,-1,witness)


if __name__=='__main__':unittest.main()
