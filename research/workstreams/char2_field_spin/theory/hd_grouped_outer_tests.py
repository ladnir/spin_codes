import unittest
from fractions import Fraction
from hd_grouped_outer import atom_cap,envelope,toy_mixer_counts,toy_group_counts,field_multiplier_cap,gf4_mul


class GroupedOuterTests(unittest.TestCase):
    def test_same_full_group_envelope(self):
        old=Fraction(65536**256,65535**128)
        for g in (2,4,8,16,32,64):
            self.assertEqual(envelope(65536,g,256//g,128//g),old)

    def test_every_toy_input_output_mixer_atom(self):
        counts,den=toy_mixer_counts()
        self.assertEqual(counts[0][0],den)
        for x in range(1,16):
            self.assertEqual(counts[x][0],0)
            self.assertEqual(sum(counts[x]),den)
            self.assertLessEqual(max(counts[x])*9,den)

    def test_all_group_outputs_and_all_messages(self):
        counts,den=toy_group_counts()
        self.assertEqual(counts[0],0)
        self.assertEqual(sum(counts),255*den)
        self.assertEqual(Fraction(max(counts),den),atom_cap(4,2,4,2))

    def test_one_diagonal_is_not_a_drop_in(self):
        counts,den=toy_group_counts(False)
        self.assertGreater(Fraction(max(counts),den),atom_cap(4,2,4,2))

    def test_gauge_normalization_preserves_every_atom(self):
        original,d0=toy_mixer_counts()
        normalized,d1=toy_mixer_counts(normalized=True)
        for a,b in zip(original,normalized):
            for x,y in zip(a,b):self.assertEqual(x*d1,y*d0)

    def test_extension_multiplier_uniform_nonzero_images(self):
        # GF16 as GF4[v]/(v^2+v+2); Tr_GF4/GF2(2)=1.
        self.assertEqual(2^gf4_mul(2,2),1)
        def mul(x,y):
            a,b=x&3,x>>2;c,d=y&3,y>>2
            p,q,r=gf4_mul(a,c),gf4_mul(b,d),gf4_mul(a^b,c^d)
            return (p^gf4_mul(2,q))+4*(r^p)
        for x in range(1,16):self.assertEqual({mul(x,y) for y in range(1,16)},set(range(1,16)))
        self.assertLess(field_multiplier_cap(65536,2,128,64),atom_cap(65536,2,128,64))

    def test_invalid_dimensions(self):
        for args in ((1,4,64,32),(4,0,4,2),(4,2,4,4)):
            with self.assertRaises(ValueError):atom_cap(*args)

    def test_subfield_screen_density_and_old_field_control(self):
        from math import log2
        import numpy as np
        from hd_grouped_subfield_screen import coefficients,bound,density_bits,LENGTH
        from hd_sumrank_profile_pairs import screen
        self.assertAlmostEqual(density_bits(4,2),log2(15/9),places=14)
        raw=coefficients(.9,.1,q=65536,local_density=False)
        old=screen(128,.9,.1,raw,width=64,forbidden_weight=LENGTH*7//20)
        value,_=bound(128,.9,.1,raw,q=65536,local_density=False)
        self.assertAlmostEqual(value,old['log2_upper_screen']+(LENGTH//64)*density_bits(65536,16),places=7)
        weighted=coefficients(.9,.1,q=65536,local_density=True)
        np.testing.assert_allclose(weighted[:,:,0],raw[:,:,0],atol=1e-12)
        np.testing.assert_allclose(weighted[:,:,1],raw[:,:,1]+.1*density_bits(65536,16),atol=1e-12)
        self.assertLessEqual(bound(128,.9,.1,weighted,q=65536)[0],value)

    def test_subfield_screen_rejects_unported_operator_variants(self):
        from hd_common_refresh_screen import matrices
        for q in (4,1024):
            with self.assertRaises(ValueError):matrices(.9,field_order=q)
        with self.assertRaises(ValueError):matrices(.9,field_order=256,exact_small=1)
        with self.assertRaises(ValueError):matrices(.9,field_order=256,scattered=True)

    def test_inverse_refresh_weight_law_every_toy_input(self):
        from hd_grouped_quiet import weight_transition
        counts,den=toy_mixer_counts();transition=weight_transition(4,2)
        weight=lambda x:bool(x&3)+bool(x>>2)
        for x in range(1,16):
            for b in (1,2):
                self.assertEqual(Fraction(sum(counts[x][y] for y in range(1,16) if weight(y)==b),den),transition[weight(x)-1][b-1])
        for q in (4,256,65536):
            s=2 if q==4 else 16
            matrix=weight_transition(q,s)
            from math import comb
            pi=[Fraction(comb(s,a)*(q-1)**a,q**s-1) for a in range(1,s+1)]
            for b in range(s):self.assertEqual(sum(pi[a]*matrix[a][b] for a in range(s)),pi[b])

    def test_inverse_quiet_cap_matches_all_toy_output_atoms(self):
        from hd_grouped_quiet import quiet_cap,quiet_screen
        counts,den=toy_mixer_counts();z=Fraction(2,3)
        # Identity emission: one coordinate per h=1 block, so the weight
        # floor is exact. Enumerate the full inverse-weighted operator.
        actual=[sum(Fraction(counts[x][y],den)*z**(bool(x&3)+bool(x>>2)) for x in range(1,16)) for y in range(1,16)]
        cap=quiet_cap(z,q=4,state=2,width=2,h=1)
        self.assertEqual(max(actual),cap)
        self.assertLessEqual(cap,1)
        self.assertAlmostEqual(float(cap),quiet_screen(float(z),q=4,state=2,width=2,h=1),places=14)

    def test_inverse_quiet_screen_only_changes_zero_degree_live_entry(self):
        import numpy as np
        from hd_grouped_subfield_screen import coefficients
        old=coefficients(.998,.5,q=256,h=4)
        new=coefficients(.998,.5,q=256,h=4,quiet_inverse=True)
        self.assertLess(new[0,1,1],old[0,1,1])
        old[0,1,1]=new[0,1,1]
        np.testing.assert_array_equal(new,old)

    def test_directed_custom_source_preserves_old_receipt(self):
        from flint import arb,ctx
        from hd_sumrank_paired_certificate import receipt
        from hd_common_refresh_packets import grouped_packets
        z,a,p=Fraction(999,1000),Fraction(1,2),Fraction(1,1000)
        options=dict(paired=False,feedthrough='common2',threshold=1000000)
        with ctx.workprec(256):
            old=receipt(1,2,z,a,p,**options)
            custom=receipt(1,2,z,a,p,packet_source=lambda z,maximum:grouped_packets(z,maximum=maximum,h=2),**options)
            self.assertEqual(old['certified_intervals'],custom['certified_intervals'])
            self.assertTrue(arb(old['union_log2_upper']).overlaps(arb(custom['union_log2_upper'])))
        with self.assertRaises(ValueError):receipt(1,2,z,a,p,field_order=256,**options)

    def test_pinned_subfield_sparse_witness_is_directed(self):
        from flint import arb,ctx
        from hd_grouped_subfield_certificate import evaluate
        from hd_grouped_subfield_witnesses import WITNESSES
        with ctx.workprec(256):
            item=evaluate(WITNESSES[0],256)
            self.assertEqual(item['field_order'],256)
            self.assertEqual(item['certified_intervals'],[[1,1]])
            self.assertEqual(item['mantissa_bits'],16)
            self.assertTrue(arb(item['union_log2_upper'])<-113)
            self.assertFalse(item['whole_code_certificate'])


if __name__=='__main__':unittest.main()
