"""Exact toy coset transforms and replay of the bounded payoff probe."""
import math
import unittest
import bridge as base
import kernel_coset_payoff as probe


class KernelCosetPayoffTest(unittest.TestCase):
    def test_toy_coset_transform(self):
        # P=repetition[4,1], D=even[4,3], all three nonzero P-cosets
        # have two words of weight two. Here E=D dual=P and H=u.
        p=[1,0,0,0,1];u=[0,0,2,0,0]
        def transform(values):
            return [sum(values[w]*sum((-1)**i*math.comb(w,i)*math.comb(4-w,j-i)
                        for i in range(max(0,j-(4-w)),min(w,j)+1)) for w in range(5)) for j in range(5)]
        kp,ku=transform(p),transform(u)
        self.assertEqual([a+3*b for a,b in zip(kp,ku)],[8*x for x in p])
        self.assertEqual([a-b for a,b in zip(kp,ku)],[8*x for x in u])

    def test_exact_structure(self):
        import pdual_low_weight as low
        p_rows,_,values=low.data()
        self.assertTrue(all(low.syndrome(w,values)==0 for w in p_rows))
        facts=probe.structure()
        self.assertEqual(facts['D_dimension'],139)
        self.assertEqual(facts['E_dimension'],117)
        self.assertEqual(facts['nonzero_P_coset_orbit_size'],255)

    def test_probe_replay(self):
        folder=base.HERE/'generated/kernel_coset_payoff_v1'
        model,scales,objective,norm,lp,_,_,facts=probe.prepare()
        self.assertEqual(lp,(folder/'h_38.lp').read_text())
        # If the bounded solver is inconclusive, export reconstruction is the
        # only result. Never infer a certificate from a solution file's presence.
        if (folder/'audit.json').exists():
            saved=base.read(folder/'audit.json')
            for name,digest in saved['source_sha256'].items():self.assertEqual(base.sha(base.ROOT/name),digest)
            result=probe.probe.audit(folder,model,scales,objective,norm)
            self.assertEqual(result['paired_upper'],saved['paired_upper'])
            self.assertEqual(result['fixed_coefficients_gain_ceiling_bits'],saved['fixed_coefficients_gain_ceiling_bits'])
            self.assertEqual(facts,saved['structure'])


if __name__=='__main__':unittest.main()
