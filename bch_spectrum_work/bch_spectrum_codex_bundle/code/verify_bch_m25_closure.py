"""Read-only M25 verification, without relying on any statistical acceptance.

Checks exact BCH row derivation and dual certificates, Johnson witnesses,
saved higher-occupation arithmetic, all source hashes, and the full sum.
Run certify_bch_m25_closure.py --verify to recompute the Q1 transfer itself.
"""
from __future__ import annotations
import itertools
import json
import math
import sys
from fractions import Fraction
from pathlib import Path
import numpy as np
sys.dont_write_bytecode=True
sys.set_int_max_str_digits(0)
from audit_bch_q1_full_arb import decode
from audit_bch_closure_envelope import audit
from certify_bch_m25_closure import deterministic_caps
from export_constant_weight_delsarte import attach_qsopt_certificate
from run_higher_endpoint_preflight import sha

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
    full=checked(GEN/'bch256_m25_unconditional_distance.json')
    independent=checked(GEN/'bch256_m25_independent_transfer_check.json')
    envelope=checked(GEN/'bch256_closure_deterministic_envelope.json')
    assert audit()==envelope
    assert not envelope['orbit_lower_bounds_used']
    assert not envelope['statistical_caps_used']
    assert not envelope['lattice_rounding_used']
    assert independent['all_92_shells_independently_checked']
    assert independent['shells_checked']==92
    for w in (52,54):
        path=GEN/f'johnson_n256_w{w}_d38.json'
        stored=checked(path)
        assert (stored['n'],stored['weight'],stored['minimum_distance'])==(256,w,38)
        assert attach_qsopt_certificate(stored.copy(),256,w,19,path.with_suffix('.sol'))==stored
    caps=deterministic_caps()
    q1=Fraction(0)
    assert set(map(int,full['coefficient_rows']))=={w for w,c in caps.items() if c}
    for name,row in full['coefficient_rows'].items():
        w=int(name)
        assert row['weight']==w and row['cap']==caps[w]
        assert decode(row['s'])>0 and decode(row['coefficient_upper'])>0
        q1+=row['cap']*decode(row['coefficient_upper'])
    assert q1==decode(full['Q1_upper'])
    higher=[]
    r=Fraction((1<<20)+1,1<<20)
    for q,field,cap_field in ((2,'conditional_pair_upper','spectrum_envelope'),
                              (3,'conditional_upper','weight_class_caps')):
        stem='bch256_q2_positive_outward' if q==2 else 'bch256_q3_capped_outward'
        receipt=checked(GEN/(stem+'.json'))
        params=receipt['parameters']
        assert (params['memory_bits'],params['distance_cutoff'],params['occupation'],params['outer_rows'],
                params['outer_length'],params['message_bits'],params['output_bits'])==(22,209716,q,8192,256,1<<20,1<<21)
        old_caps={int(w):cap for w,cap in receipt[cap_field].items() if int(w) and cap}
        if q==2:
            assert set(old_caps)=={w for w,c in caps.items() if c}
            assert all(caps[w]<=r*old_caps[w] for w in old_caps)
        else:
            assert set(old_caps)==set(range(38,72,2))
            assert old_caps[70]==(1<<128)-1
            assert all(caps[w]<=r*old_caps[w] for w in range(38,70,2))
        with np.load(GEN/(stem+'.npz')) as data:
            array=data[field]
            value=math.comb(8192,q)*sum((math.prod(old_caps[w] for w in ws)*Fraction.from_float(array[ws])
                for ws in itertools.product(sorted(old_caps),repeat=q)),Fraction(0))
        assert value==decode(receipt[f'Q{q}_upper'])==decode(full[f'old_M22_Q{q}_upper'])
        higher.append(value)
    tail=checked(GEN/'bch256_occupation_tail_outward.json')
    assert tail['parameters']==dict(outer_rows=8192,outer_length=256,message_bits=1<<20,
        output_bits=1<<21,distance_cutoff=209716,memory_bits=22)
    assert [row['q'] for row in tail['rows']]==list(range(4,8193))
    for row in tail['rows']:
        value=decode(row['log2_upper'])
        assert row['dyadic_upper_exponent']==-(-value.numerator//value.denominator)<=-58
    exponent=min(row['dyadic_upper_exponent'] for row in tail['rows'])
    tail_sum=Fraction(sum(1<<(row['dyadic_upper_exponent']-exponent) for row in tail['rows']),1<<-exponent)
    assert tail_sum==decode(tail['tail_upper'])==decode(full['old_M22_tail_upper'])<Fraction(1,1<<83)
    assert checked(GEN/'bch256_tail_transfer_crosschecks.json')['all_8189_occupations_checked']
    assert checked(GEN/'bch256_tail_envelope_checks.json')['exact_prefix_comparisons_checked']==5397
    inflation=1/(1-8192*(r-1))
    assert inflation==Fraction(128,127)==decode(full['higher_occupation_inflation'])
    total=q1+inflation*(sum(higher,Fraction(0))+tail_sum)
    assert total==decode(full['full_first_moment_upper'])<Fraction(1,1<<41)
    assert full['parameters']==dict(message_bits=1<<20,outer_rows=8192,outer_length=256,
        output_bits=1<<21,distance_cutoff=209716,memory_bits=25)
    assert full['exact_full_bound_below_2_to_minus_40'] and full['all_occupations_covered']
    assert not full['statistical_shell_assumptions_used']
    assert not full['empirical_orbit_lower_bounds_used']
    assert not full['original_M22_target_closed']
    assert not full['paper_random_convolution_transfer_established']
    assert not full['RM2Sub_transfer_established']
    print(json.dumps(dict(status='M25 deterministic-spectrum closure verified',
        JSON_dependencies_checked=len(cache),BCH_constraint_rows_rederived=396,
        statistical_acceptance_used=False,all_occupations_covered=True,
        exact_full_bound_below_2_to_minus_41=True,original_M22_target_closed=False,
        scope='Published BCH anchors, exact LP witnesses, numerical certificates, and the explicit marginal-memory and envelope-inflation lemmas. No statistical test acceptance is used.'),indent=2))


if __name__=='__main__':
    main()
