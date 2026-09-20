"""Combine completed Q1 shell evidence with deterministic Q2/Q3 certificates."""
from __future__ import annotations

import json
import itertools
import sys
from fractions import Fraction
from pathlib import Path

sys.dont_write_bytecode=True
from audit_bch_q1_full_arb import decode,encode
from run_higher_endpoint_preflight import sha,write_new

ROOT=Path(__file__).resolve().parents[1]
GEN=ROOT/'generated'


def checked_receipt(path):
    payload=json.loads(path.read_text())
    for name,expected in payload['source_sha256'].items():
        assert sha(ROOT/name)==expected,f'Certificate dependency changed: {name}'
    return payload


def main():
    output=GEN/'bch256_q123_joint_evidence.json'
    assert not output.exists()
    paths=[GEN/'bch256_q1_full_arb_evidence_audit.json',GEN/'bch256_q2_positive_outward.json',
           GEN/'bch256_q3_capped_outward.json']
    q1,q2,q3=[checked_receipt(p) for p in paths]
    audit_paths=[GEN/'bch256_low_occupation_crosschecks.json',GEN/'bch256_deterministic_caps_reaudit.json']
    crosschecks,cap_audit=[checked_receipt(p) for p in audit_paths]
    assert crosschecks['both_symmetrized_aggregates_below_2_to_minus_60']
    assert set(cap_audit['checks'])=={str(w) for w in range(38,56,2)}
    assert q1['parameters']['distance']==q2['parameters']['distance_cutoff']==q3['parameters']['distance_cutoff']==209716
    assert q1['parameters']['outer_rows']==q2['parameters']['outer_rows']==q3['parameters']['outer_rows']==8192
    assert q1['parameters']['memory_bits']==q2['parameters']['memory_bits']==q3['parameters']['memory_bits']==22
    values=[decode(q1['full_Q1_upper']),decode(q2['Q2_upper']),decode(q3['Q3_upper'])]
    total=sum(values,Fraction(0))
    target=Fraction(1,1<<40)
    residual=target-total
    alphas=[Fraction(1,1<<40),Fraction(1,1<<41),Fraction(1,1<<41)]
    # Intersection-union soundness: accepting the joint claim requires every
    # test to pass. For each fixed nonempty set of false caps, acceptance is
    # contained in every false cap's false-acceptance event.
    null_cases=[subset for n in (1,2,3) for subset in itertools.combinations(range(3),n)]
    joint_error=max(min(alphas[i] for i in subset) for subset in null_cases)
    assert joint_error==Fraction(1,1<<40)
    all_caps_accepted=set(q1['accepted_shell_caps'])=={'38','40','42'}
    receipt=dict(classification='Exact rational aggregate of certified Q1/Q2/Q3 bounds; Q1 shell evidence is conditional statistical acceptance',
        accepted_shell_caps=q1['accepted_shell_caps'],
        occupation_uppers={str(j+1):encode(x) for j,x in enumerate(values)},
        Q1_through_Q3_upper=encode(total),margin_bits_diagnostic=-encode(total)['log2_diagnostic'],
        exact_Q1_through_Q3_below_target=total<target,
        residual_for_Q_ge_4=encode(residual),residual_fraction_of_target=float(residual/target),
        simple_sufficient_tail_allowance='2^-45' if residual>=Fraction(1,1<<45) else None,
        full_SPIN_40_bit_theorem_established=False,
        next_required_inequality='sum_{q=4}^{8192} F_q <= residual_for_Q_ge_4',
        randomness_qualification=q1['randomness_qualification'],
        statistical_familywise_false_accept_error_under_ideal_IID='at most 2^-39 for the three fixed-budget shell tests, separate from SPIN setup failure',
        all_three_shell_tests_passed=all_caps_accepted,
        joint_shell_claim_accepted=all_caps_accepted,
        joint_shell_claim_false_accept_upper_under_ideal_IID=encode(joint_error),
        joint_soundness_argument='The joint claim is accepted only when all three predeclared tests pass. Any fixed violating code has at least one false cap, so joint acceptance is contained in that test\'s false-acceptance event. The largest individual error bound is 2^-40. No independence between the three tests is required for this implication. This is not a posterior probability or a deterministic spectrum theorem.',
        scope='Fixed BCH-derived C=[256,128,38], RandomStepConv-M22, 2^20 message bits and 2^21 output bits; not an RM2Sub result',
        additional_checks='Independent zero-weight degeneracy checks, symmetrized reaggregation, exact LP dual-witness and cap-mapping rechecks',
        remaining_obligations=['Bound and sum Q>=4 within the residual','Replace the RandomStepConv transfer for RM2Sub application','A deterministic spectrum theorem is not supplied by statistical shell acceptance'],
        source_sha256={str(p.relative_to(ROOT)):sha(p) for p in paths+audit_paths+[Path(__file__)]})
    write_new(output,receipt)
    print(json.dumps({k:receipt[k] for k in ('accepted_shell_caps','margin_bits_diagnostic',
        'exact_Q1_through_Q3_below_target','residual_fraction_of_target','full_SPIN_40_bit_theorem_established')},indent=2))


if __name__=='__main__':
    main()
