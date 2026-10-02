"""Independent completed-receipt checks for the new GL15/drop-row10 code."""
import json
import os
from pathlib import Path
import unittest

import test_disjoint_pair_receipt as audited
import screen_paired_restriction as restriction
from flint import ctx

RECEIPT=Path(os.environ.get('SPIN_PAIRED15_RECEIPT',
    str(Path(__file__).resolve().parent.parent/'proof/paired-s15-drop10-whole-p256.json')))


class Paired15ReceiptTests(audited.CompleteReceiptTests):
    @classmethod
    def setUpClass(cls):
        ctx.prec=256
        cls.record=json.loads(RECEIPT.read_bytes())
        cls.wrapper,cls.map_record=restriction.prepare((10,))

    def test_complete_exact_geometry(self):
        record=self.record
        expected=dict(K=65536,N=131072,threshold=13107,
            target_minimum_distance=13108,target_margin_bits=40,
            groups=512,group_dimension=128,regions=64,packet_bits=4,
            physical_t=64,state_bits=15,physical_steps=2048,
            physical_steps_per_region=32,precision=256,
            zero_initial_state=True,final_flush=False,fresh_replay=True,
            whole_code_certificate=True,all_occupancies_covered=True,target_met=True,
            state_continuity='retained_across_every_step_and_region',
            inner_randomness='independent uniform GL15 matrix for each physical t64 step')
        for name,value in expected.items():self.assertEqual(record[name],value,name)
        self.assertEqual(set(record['occupancy_uppers']),set(map(str,range(1,513))))
        self.assertEqual(set(record['tail_choices']),set(map(str,range(3,513))))
        for first,last,tilts in record['tail_recipes']:
            for q in range(first,last+1):self.assertIn(record['tail_choices'][str(q)],tilts)

    def test_all_saved_source_hashes(self):
        super().test_all_saved_source_hashes()
        for name in ('reproduce_s15_drop10.py','screen_paired_restriction.py'):
            self.assertTrue(any(Path(p).name==name for p in self.record['source_sha256']),name)

    def test_receipt_map_is_fresh_declared_map(self):
        self.assertEqual(self.record['map_record'],json.loads(json.dumps(self.map_record)))
        self.assertEqual(self.record['map_record']['deleted_parent_rows'],[10])
        self.assertEqual(self.wrapper['bits'],15)
        self.assertEqual(self.wrapper['physical_steps'],2)
        self.assertEqual(self.wrapper['physical_windows'],16)
        self.assertEqual(self.wrapper['macro_windows'],32)

    # The inherited exact-integer union and fresh direct-physical q=1 checks
    # use this class's wrapper/receipt and do not import any replay endpoint.


if __name__=='__main__':
    unittest.main()
