"""Read-only dependency and exact-aggregation check of the full conditional bound."""
from __future__ import annotations
import json
import sys
from fractions import Fraction
from pathlib import Path
sys.dont_write_bytecode=True
sys.set_int_max_str_digits(0)
from audit_bch_q1_full_arb import decode
from run_higher_endpoint_preflight import sha
from verify_bch_q123_evidence import main as verify_low
ROOT=Path(__file__).resolve().parents[1]
GEN=ROOT/'generated'


def main():
    cache={}

    def checked(path):
        if path in cache:
            return cache[path]
        value=json.loads(path.read_text())
        cache[path]=value
        aliases={'diagnostic':GEN/'random_inner_oa15_majorant_diagnostic.json',
                 'johnson_w52':GEN/'johnson_n256_w52_d38.json',
                 'johnson_w54':GEN/'johnson_n256_w54_d38.json'} if path.name=='random_inner_threshold_outward.json' else {}
        for name,expected in value.get('source_sha256',{}).items():
            dependency=aliases.get(name,ROOT/name)
            if path.name=='bch256_q2_envelope_diagnostic.json':
                assert Path(name).name==name
                dependency=GEN/'bch256_q2_envelope_diagnostic_sources'/name
            assert sha(dependency)==expected,name
            if dependency.suffix=='.json':
                checked(dependency)
        return value

    full=checked(GEN/'bch256_full_random_inner_conditional.json')
    low=checked(GEN/'bch256_q123_joint_evidence.json')
    tail=checked(GEN/'bch256_occupation_tail_outward.json')
    transfer=checked(GEN/'bch256_tail_transfer_crosschecks.json')
    envelope=checked(GEN/'bch256_tail_envelope_checks.json')
    verify_low()
    rows=tail['rows']
    assert [row['q'] for row in rows]==list(range(4,8193))
    for row in rows:
        bound=decode(row['log2_upper'])
        assert row['dyadic_upper_exponent']==-(-bound.numerator//bound.denominator)<=-58
    minimum=min(row['dyadic_upper_exponent'] for row in rows)
    numerator=sum(1<<(row['dyadic_upper_exponent']-minimum) for row in rows)
    tail_total=Fraction(numerator,1<<-minimum)
    assert tail_total==decode(tail['tail_upper'])==decode(full['tail_upper'])<Fraction(1,1<<83)
    total=decode(low['Q1_through_Q3_upper'])+tail_total
    assert total==decode(full['full_first_moment_upper'])<Fraction(1,1<<40)
    assert full['parameters']==tail['parameters']
    assert full['parameters']==dict(outer_rows=8192,outer_length=256,message_bits=1<<20,
        output_bits=1<<21,distance_cutoff=209716,memory_bits=22)
    assert full['accepted_shell_caps']==low['accepted_shell_caps']
    assert full['conditional_random_inner_distance_bound_established']
    assert not full['unconditional_deterministic_BCH_spectrum_theorem']
    assert not full['RM2Sub_transfer_established']
    assert transfer['unnormalized_region_count_checks']==83
    assert transfer['direct_full_length_power_checks']==8106
    assert transfer['all_8189_occupations_checked']
    assert envelope['exact_prefix_comparisons_checked']==5397
    assert envelope['generator_containment_L71_in_Q123_verified']
    assert envelope['anchor_dual_minimum_distance']==16
    print(json.dumps(dict(status='full conditional random-inner distance bound verified',
        JSON_dependencies_checked=len(cache),all_occupations_covered=True,
        exact_full_bound_below_2_to_minus_40=True,
        scope='Source hashes and exact aggregates; independent transfer and polynomial checks are retained in linked receipts. Shell caps remain statistical, and RM2Sub is not covered.'),indent=2))


if __name__=='__main__':
    main()
