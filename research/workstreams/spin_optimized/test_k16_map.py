"""Exact selected-map checks; optionally bind constants to an assembled proof."""
import argparse
from collections import Counter
import hashlib
import json
from pathlib import Path
import re
import unittest
from k16_map import A, B, header


class SelectedMapTest(unittest.TestCase):
    def test_expansion_spectrum(self):
        spectrum = Counter(sum((c & state).bit_count() & 1 for c in A) for state in range(1, 4096))
        self.assertEqual(spectrum, {24:286, 28:880, 32:1743, 36:912, 40:274})

    def test_feedback_full_rank(self):
        self.assertEqual(len(B), 64)
        for state in range(1, 4096):
            self.assertTrue(any((c & state).bit_count() & 1 for c in B))

    def test_generated_circuits_and_map(self):
        text = header()  # Symbolically checks all three straight-line circuits.
        for name, expected in [('columns', A), ('feedbackColumns', B)]:
            actual = re.search(r'\b'+name+r'\{([^}]+)\};', text)
            self.assertEqual(list(map(int, actual[1].split(','))), expected)


def bind(proof_path, header_path=None):
    proof = json.loads(proof_path.read_text())
    instance = proof['instance']
    inner = instance['inner']
    assert proof['full_distance_proved'] and proof['margin_bits'] > 40
    assert (instance['message_bits'], instance['output_bits'], instance['cutoff']) == (65536, 131072, 13107)
    assert (inner['t'], inner['s']) == (64,12)
    assert inner['transvection_rounds'] in (1,2)
    if inner['transvection_rounds']==2:
        assert proof['status']=='VERIFIED_K16_TWO_ROUND_FULL_DISTANCE_BOUND' and proof['margin_bits']>49
    assert inner['expansion_columns'] == A and inner['feedback_columns'] == B
    if header_path:
        expected = header().encode()
        assert header_path.read_bytes() == expected
        print('Generated header SHA256:', hashlib.sha256(expected).hexdigest())
    print('Fixed A/B match assembled full certificate:', proof['margin_bits'], 'bits')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--proof', type=Path)
    parser.add_argument('--header', type=Path)
    args = parser.parse_args()
    suite = unittest.defaultTestLoader.loadTestsFromTestCase(SelectedMapTest)
    if not unittest.TextTestRunner(verbosity=2).run(suite).wasSuccessful():
        raise SystemExit(1)
    if args.proof:
        bind(args.proof, args.header)
