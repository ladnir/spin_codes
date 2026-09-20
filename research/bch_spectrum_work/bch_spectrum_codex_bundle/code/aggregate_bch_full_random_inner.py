"""Combine all occupations into a conditional RandomStepConv distance bound."""
from __future__ import annotations
import json
import sys
from fractions import Fraction
from pathlib import Path
sys.dont_write_bytecode=True
sys.set_int_max_str_digits(0)
from audit_bch_q1_full_arb import decode,encode
from run_higher_endpoint_preflight import sha,write_new
from verify_bch_q123_evidence import main as verify_low_occupations
ROOT=Path(__file__).resolve().parents[1]
GEN=ROOT/'generated'


def checked(path):
    value=json.loads(path.read_text())
    for name,expected in value['source_sha256'].items():
        assert sha(ROOT/name)==expected,name
    return value


def main():
    output=GEN/'bch256_full_random_inner_conditional.json'
    assert not output.exists()
    verify_low_occupations()
    paths=[GEN/'bch256_q123_joint_evidence.json',GEN/'bch256_occupation_tail_outward.json',
           GEN/'bch256_tail_envelope_checks.json',GEN/'bch256_tail_transfer_crosschecks.json']
    low,tail,envelopes,transfer=[checked(p) for p in paths]
    assert [r['q'] for r in tail['rows']]==list(range(4,8193))
    for row in tail['rows']:
        value=decode(row['log2_upper'])
        assert row['dyadic_upper_exponent']==-(-value.numerator//value.denominator)
    independent_tail=sum((Fraction(1,1<<-r['dyadic_upper_exponent']) for r in tail['rows']),Fraction(0))
    assert independent_tail==decode(tail['tail_upper'])<Fraction(1,1<<83)
    assert tail['exact_tail_at_most_2_to_minus_45'] and transfer['all_8189_occupations_checked']
    assert envelopes['published_Table7_k71_matches_every_anchor_entry']
    assert envelopes['generator_containment_L71_in_Q123_verified']
    assert envelopes['anchor_dual_minimum_distance']==16
    assert low['joint_shell_claim_accepted']
    total=decode(low['Q1_through_Q3_upper'])+independent_tail
    target=Fraction(1,1<<40)
    assert total<target
    receipt=dict(classification='All-occupation first-moment distance bound, conditional on the three fixed BCH shell caps',
        parameters=tail['parameters'],accepted_shell_caps=low['accepted_shell_caps'],
        full_first_moment_upper=encode(total),full_margin_bits_diagnostic=-encode(total)['log2_diagnostic'],
        exact_full_bound_below_2_to_minus_40=True,
        tail_upper=encode(independent_tail),tail_exactly_below_2_to_minus_83=True,
        all_occupations_1_through_8192_covered=True,
        conditional_random_inner_distance_bound_established=True,
        bad_event='There exists a nonzero 2^20-bit message with encoded output weight at most 209716',
        probability_space='Independent uniform row-coordinate permutations and region permutations; ideal RandomStepConv-M22 maps sampled independently by position and shared by all messages',
        inner_interface='Initial state zero. At each position an independent uniform binary linear map sends the 22-bit state and one input bit to the next 22-bit state and one output bit.',
        method='First moment and Markov inequality; no independence between different messages and no variance estimate required',
        statistical_joint_false_accept_upper_under_ideal_IID=low['joint_shell_claim_false_accept_upper_under_ideal_IID'],
        randomness_qualification=low['randomness_qualification'],
        unconditional_deterministic_BCH_spectrum_theorem=False,RM2Sub_transfer_established=False,
        remaining_obligations=['Replace statistical shell acceptance by deterministic caps for an unconditional fixed-code theorem',
            'Transfer to RM2Sub or another concrete inner separately',
            'Any computational implementation of the ideal setup needs its own randomness/model justification'],
        source_sha256={str(p.relative_to(ROOT)):sha(p) for p in paths+[Path(__file__),ROOT/'code/verify_bch_q123_evidence.py']})
    write_new(output,receipt)
    print(json.dumps({k:receipt[k] for k in ('full_margin_bits_diagnostic','exact_full_bound_below_2_to_minus_40',
        'tail_exactly_below_2_to_minus_83','all_occupations_1_through_8192_covered',
        'conditional_random_inner_distance_bound_established')},indent=2))


if __name__=='__main__':
    main()
