import unittest
from collections import Counter
from fractions import Fraction as F
from itertools import combinations

from hd_blockwise_maps import matmul,matvec
from hd_blockwise_shared_maps import (matrices,coset_matrices,evaluation,inverse,minor_count,receipt,
    identity,branch_lower,bernoulli_branch_fraction,fixed_subset_branch_fraction,
    occupancy_profile_probability,prescribed_mds_count)


class SharedMapTests(unittest.TestCase):
    def test_exact_matrix_receipt(self):
        r=receipt()
        self.assertEqual(r['F_first_row'],[52,172,158,17,209,32,250,29])
        self.assertEqual(r['H'],((245,29,9,1),(29,0,1,0),(9,1,0,0),(1,0,0,0)))
        self.assertEqual(r['rank_F_minus_I'],32)
        self.assertEqual(r['rank_stacked_C_F_minus_I'],34)
        self.assertFalse(r['whole_certificate'])

    def test_every_tile_minor(self):
        full,*_=matrices()
        local=tuple(tuple(row[:8]) for row in full[:8])
        self.assertEqual(minor_count(local),12869)

    def test_outside_packet_translations_preserve_exact_identities(self):
        for beta in (16,32,64,128,248):
            full,a,c,t,h=matrices(beta)
            self.assertEqual(matmul(full,full),identity(64))
            self.assertEqual(matmul(full,a),matmul(a,t))
            self.assertEqual(matmul(c,full),matmul(t,c))
            self.assertEqual(matmul(h,tuple(zip(*a))),c)
            self.assertFalse(any(x for row in matmul(c,a) for x in row))
            for offset in range(0,64,8):
                self.assertEqual(tuple(tuple(row[offset:offset+8]) for row in full[offset:offset+8]),
                                 tuple(tuple(row[:8]) for row in full[:8]))
        full,*_=matrices(64)
        self.assertEqual(minor_count(tuple(tuple(row[:8]) for row in full[:8])),12869)
        for beta in (0,1,7,256):
            with self.assertRaises(ValueError):
                matrices(beta)

    def test_arbitrary_disjoint_cosets_and_inverse_family(self):
        inputs=tuple(range(0,64,8))
        for outputs in (tuple(range(64,128,8)), (128,72,240,96,184,152,80,216)):
            full,a,c=coset_matrices(inputs,outputs)
            backward,back_a,back_c=coset_matrices(outputs,inputs)
            self.assertEqual(backward,full)
            self.assertEqual(matmul(full,full),identity(64))
            self.assertEqual(matmul(full,a),back_a)
            self.assertEqual(matmul(c,full),back_c)
            for start in range(0,64,8):
                local_c=tuple(row[start:start+8] for row in c)
                self.assertFalse(any(x for row in matmul(local_c,a[start:start+8]) for x in row))
            self.assertFalse(any(x for row in matmul(back_c,back_a) for x in row))
        for invalid in (inputs, (64,)*8, tuple(range(65,129,8)), (64,72)):
            with self.assertRaises(ValueError):
                coset_matrices(inputs,invalid)

    def test_shared_circuit_every_basis_vector(self):
        full,a,c,_,_=matrices()
        for source in range(68):
            x=tuple(int(j==source) for j in range(64))
            state=tuple(int(64+j==source) for j in range(4))
            emitted=[]
            syndrome=[0]*4
            for offset in range(0,64,8):
                z=list(matvec(inverse(evaluation(offset)),x[offset:offset+8]))
                for j in range(4):
                    syndrome[j]^=z[4+j]
                    z[j]^=state[j]
                emitted.extend(matvec(evaluation(offset^8),z))
            self.assertEqual(tuple(emitted),tuple(u^v for u,v in zip(matvec(full,x),matvec(a,state))))
            self.assertEqual(tuple(syndrome),matvec(c,x))

    def test_exact_region_histograms(self):
        for group,rows in [(2,8),(4,8),(8,16)]:
            for active in range(min(rows,5)+1):
                histograms=Counter()
                weight_sum=0
                for subset in combinations(range(rows),active):
                    occupied=set(subset)
                    counts=[0]*(group+1)
                    for begin in range(0,rows,group):
                        counts[sum(i in occupied for i in range(begin,begin+group))]+=1
                    histograms[tuple(counts)]+=1
                    weight_sum+=sum(n*branch_lower(w,group) for w,n in enumerate(counts))
                total=sum(histograms.values())
                for counts,frequency in histograms.items():
                    self.assertEqual(occupancy_profile_probability(rows,active,counts),F(frequency,total))
                self.assertEqual(fixed_subset_branch_fraction(rows,active,group),F(weight_sum,total*rows))

    def test_branch_diagnostics_and_mds_count(self):
        self.assertEqual(bernoulli_branch_fraction(F(1,4),2),F(13,32))
        self.assertLess(bernoulli_branch_fraction(F(1,8),4),F(2,5))
        self.assertGreater(bernoulli_branch_fraction(F(1,8),8),F(3,5))
        self.assertEqual(prescribed_mds_count(8),0)
        self.assertEqual(prescribed_mds_count(9),65535)
        # Sum all support classes of the systematic [16,8] MDS code.
        from math import comb
        self.assertEqual(sum(comb(16,w)*prescribed_mds_count(w) for w in range(17)),65536**8)


if __name__=='__main__':
    unittest.main()
