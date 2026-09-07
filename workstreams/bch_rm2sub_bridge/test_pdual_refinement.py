"""Partition assumptions and adversarial checks for conditional rank witnesses."""
import copy
from pathlib import Path
import unittest
import bridge as base
import pdual_rank_refinement as single
import pdual_multi_pivot as multi
import pdual_case_audit as partition


class PdualRefinementTest(unittest.TestCase):
    def test_full_partition_and_kernel(self):
        saved=base.read(base.HERE/'generated/pdual_partial_cover_v1.json')
        for name,digest in saved['source_sha256'].items():self.assertEqual(base.sha(base.ROOT/name),digest)
        result=partition.audit(saved['forest'])
        self.assertEqual(result,saved['audit'])
        self.assertFalse(result['full_Pdual_distance32_proved'])
        self.assertEqual(result['case_minimum_ranks']['15'],31)
        self.assertEqual(result['extended_F7_kernel']['dimension'],117)
        self.assertEqual(result['extended_F7_kernel']['minimum_distance_lower'],32)
        forged=copy.deepcopy(saved['forest']);del forged['15']['nonzero']
        with self.assertRaises(AssertionError):partition.audit(forged)
        forged=copy.deepcopy(saved['forest']);forged['15']['split']=29
        with self.assertRaises(AssertionError):partition.audit(forged)

    def test_partition_assumptions(self):
        for parent,split in ((7,15),(15,23)):
            zeros,known=single.assumptions(parent,[],[])
            z0,k0=single.assumptions(parent,[split],[])
            z1,k1=single.assumptions(parent,[],[split])
            self.assertEqual(z0,zeros|single.orbit(split))
            self.assertEqual(k0,known)
            self.assertEqual(z1,zeros)
            self.assertEqual(k1,known|single.orbit(split))
            self.assertTrue(z0.isdisjoint(k0) and z1.isdisjoint(k1))
        with self.assertRaises(AssertionError):single.assumptions(7,[7],[])
        with self.assertRaises(AssertionError):single.assumptions(15,[23],[23])

    def test_saved_witness_replay_and_tampering(self):
        files=list((base.HERE/'generated/pdual_refinement_v1').glob('*.json'))
        self.assertGreaterEqual(len(files),8)
        for path in files:
            record=base.read(path)
            if 'witness' not in record:continue
            with self.subTest(file=path.name):
                checker=multi.verify if 'initial_pivot' in record['witness'] else single.verify
                for name,digest in record['source_sha256'].items():self.assertEqual(base.sha(base.ROOT/name),digest)
                self.assertEqual(checker(record),record['witness']['rank'])
                forged=copy.deepcopy(record);forged['witness']['rank']+=1
                with self.assertRaises(AssertionError):checker(forged)
                forged=copy.deepcopy(record);forged['zero_indices']=[]
                with self.assertRaises(AssertionError):checker(forged)

    def test_multi_pivot_rank_against_all_length_five_binary_words(self):
        # Exponents modulo 255 acting on the subgroup of order five. Enumerate
        # every nonzero binary support, then all valid one-step append choices.
        from affine_wambach import gf_pow
        root=gf_pow(2,51)
        checked=0
        for word in range(1,32):
            support=[j for j in range(5) if (word>>j)&1]
            fourier=[]
            for e in range(255):
                value=0
                for j in support:value ^= gf_pow(root,e*j)
                fourier.append(value)
            zeros={e for e,v in enumerate(fourier) if not v}
            known=set(range(255))-zeros
            initial=min(known)
            pivots=sorted({min(single.orbit(e)) for e in known})[:3]
            for k in range(8):
                for shift in range(255):
                    moved=(initial*(1<<k)+shift)%255
                    if moved not in zeros:continue
                    for pivot in pivots:
                        witness=dict(initial_pivot=initial,steps=[(k,shift,pivot)],rank=2,
                                     independent_set=sorted([moved,pivot]))
                        self.assertEqual(multi.check(zeros,known,witness),2)
                        self.assertGreaterEqual(len(support),2)
                        checked+=1
        self.assertGreater(checked,1000)


if __name__=='__main__':unittest.main()
