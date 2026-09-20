"""Read-only verification of saved dependencies and low-occupation aggregates."""
from __future__ import annotations

import itertools
import json
import math
import sys
from fractions import Fraction
from pathlib import Path

sys.dont_write_bytecode = True
import numpy as np
from audit_bch_q1_full_arb import decode
from run_higher_endpoint_preflight import sha

ROOT = Path(__file__).resolve().parents[1]
GEN = ROOT / 'generated'


def main():
    checked = set()

    def read(path):
        value = json.loads(path.read_text())
        if path not in checked:
            checked.add(path)
            for name, expected in value.get('source_sha256', {}).items():
                aliases = {'diagnostic': GEN / 'random_inner_oa15_majorant_diagnostic.json',
                           'johnson_w52': GEN / 'johnson_n256_w52_d38.json',
                           'johnson_w54': GEN / 'johnson_n256_w54_d38.json'} if path.name == 'random_inner_threshold_outward.json' else {}
                dependency = aliases.get(name, ROOT / name)
                if path.name == 'bch256_q2_envelope_diagnostic.json':
                    assert Path(name).name == name
                    dependency = GEN / 'bch256_q2_envelope_diagnostic_sources' / name
                assert sha(dependency) == expected, name
                if dependency.suffix == '.json':
                    read(dependency)
        return value

    joint = read(GEN / 'bch256_q123_joint_evidence.json')
    q1 = read(GEN / 'bch256_q1_full_arb_evidence_audit.json')
    q2 = read(GEN / 'bch256_q2_positive_outward.json')
    q3 = read(GEN / 'bch256_q3_capped_outward.json')
    transfer = read(GEN / 'bch256_q1_full_arb_transfer.json')
    weights = []
    first = Fraction(0)
    for term in q1['terms']:
        subtotal = term['cap'] * sum((decode(transfer['coefficient_rows'][str(w)]['coefficient_upper'])
                                     for w in term['weights']), Fraction(0))
        assert subtotal == decode(term['upper'])
        first += subtotal
        weights.extend(term['weights'])
    assert len(weights) == len(set(weights)) == len(transfer['coefficient_rows'])
    values = [first]
    for occupation, receipt, cache_name, field, cap_field in (
        (2, q2, 'bch256_q2_positive_outward.npz', 'conditional_pair_upper', 'spectrum_envelope'),
        (3, q3, 'bch256_q3_capped_outward.npz', 'conditional_upper', 'weight_class_caps')):
        with np.load(GEN / cache_name) as cache:
            conditional = cache[field]
            caps = {int(w): cap for w, cap in receipt[cap_field].items() if int(w) and cap}
            subtotal = sum((math.prod(caps[w] for w in ws) * Fraction.from_float(conditional[ws])
                            for ws in itertools.product(sorted(caps), repeat=occupation)), Fraction(0))
        values.append(math.comb(8192, occupation) * subtotal)
    assert values == [decode(q1['full_Q1_upper']), decode(q2['Q2_upper']), decode(q3['Q3_upper'])]
    assert values == [decode(joint['occupation_uppers'][str(q)]) for q in (1, 2, 3)]
    total = sum(values, Fraction(0))
    assert total == decode(joint['Q1_through_Q3_upper'])
    residual = Fraction(1, 1 << 40) - total
    assert residual == decode(joint['residual_for_Q_ge_4'])
    assert total < Fraction(1, 1 << 40) and residual >= Fraction(1, 1 << 45)
    assert joint['simple_sufficient_tail_allowance'] == '2^-45'
    assert joint['joint_shell_claim_accepted'] and joint['all_three_shell_tests_passed']
    assert decode(joint['joint_shell_claim_false_accept_upper_under_ideal_IID']) == Fraction(1, 1 << 40)
    assert not joint['full_SPIN_40_bit_theorem_established']
    print(json.dumps(dict(status='saved low-occupation evidence verified',
                          JSON_dependencies_checked=len(checked),
                          Q1_Q2_Q3_independently_reaggregated=True,
                          residual_at_least_2_to_minus_45=True,
                          scope='Does not rerun the long sampling replays or prove Q>=4'), indent=2))


if __name__ == '__main__':
    main()
