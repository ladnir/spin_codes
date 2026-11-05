"""Fresh directed selected high-arity bound for scheduled Hlambda.

This computes one new regional bank per replay. It reuses reviewed numerical
operators, not an old outer receipt. Full independent nonzero D1 and fixed
invertible D2 Q Hlambda P give a unique original preimage; on total RS arity
A>=247 the expected raw counting density is at most255^(-8*247).
Only r=18 is covered, but every indexed source support and matching profile
at that r is included. This is not a whole-code certificate.
"""
from fractions import Fraction as F
from hashlib import sha256
import json
from pathlib import Path
from time import perf_counter

from flint import arb,ctx

import hd_two_sum_gf16_high_arity_certificate as arithmetic
from hd_lambda_outer import matrix
from hd_outer_pair_layers import field_multiply
from hd_pair_scheduled_outer import column_permutations

DIRECTORY=Path(__file__).resolve().parent
R,A0,MARK,ALPHA=18,247,F(691,700),F(1,5)
PINS={
    'hd_two_sum_gf16_high_arity_certificate.py':'9c9c2e4d1bfc95f2952d6bc743b9047489128d5cc26a43f5186ea721b1da4027',
    'hd_lambda_outer.py':'f0c063843e295a160a1186c6c084a07dd1c2a3745005ed68057fde59913c99e7',
    'hd_outer_pair_layers.py':'017c5ba20747ae8e25acf02393e9702a03b086038e1814ac2424cdfe49d66f60',
    'hd_pair_scheduled_outer.py':'64dc50bf18171d371ae16f4484087e9d42d89cc7c9b55a7576be55b1251e4741',
}
SCOPE=dict(construction='scheduled-lambda-bc-g16-rs16x8-h4-swapped-pairs-k20-rate-half',
    field_order=256,field_polynomial='0x11b',lambda_aes_byte=0xbc,
    input_schedule='P_j(i)=inv_GF16(i xor j); inv(0)=0; polynomial0x13',
    output_schedule='Q_j(r)=r xor j',column_map='D2_j Q_j Hlambda P_j D1_j',
    D1='full independent uniform nonzero coordinate multipliers',
    D2='full independent uniform nonzero coordinate multipliers, independent of D1',
    crossed_outer_layer=False,source_dimension=1048576,output_length=2097152,
    forbidden_weight=734003,outer_groups=8192,constituent_length=16,constituent_dimension=8,
    active_groups=R,total_arity_minimum=A0,total_arity_maximum=16*R,
    mark=str(MARK),alpha=str(ALPHA),packet_mantissa_bits=60,
    density_per_fixed_source='256^(256*r)/255^(8*A0)',
    regional_amplitudes='uniform field including zero',
    alternating_matchings=True,pair_swaps=True,independent_regional_pair_permutations=True,
    ordered_initial_zero_state=True,regions=256,inner_state_dimension=4,
    source_count_method='profile_count_upper',rowset_count_factor=1,cohort_count_factor=1,
    source_support_union_outside_alpha=True,projective_normalizations=1,
    all_selected_indexed_source_supports=True,all_selected_matching_profiles=True,
    old_outer_certificates_transferred=False,saved_receipts_used=False,
    low_total_arity_complement_covered=False,other_active_group_weights_covered=False,
    whole_code_certificate=False,directed=True)


def checked_dependencies():
    for name,want in PINS.items():
        if sha256((DIRECTORY/name).read_bytes()).hexdigest()!=want:
            raise ValueError('Reviewed dependency changed: '+name)
    arithmetic.checked_dependencies()
    h=matrix(0xbc)
    # Exact AES-field identity, independent of the fast mixer circuit.
    for i in range(16):
        for j in range(16):
            entry=0
            for k in range(16):entry^=field_multiply(h[i][k],h[k][j],modulus=0x11b)
            if entry!=int(i==j):raise ArithmeticError('Hlambda is not an involution')
    for j in range(16):
        before,after=column_permutations('input-inverse-output-xor',j)
        if set(before)!=set(range(16)) or set(after)!=set(range(16)):
            raise ArithmeticError('Nonbijective fixed schedule')
    return dict(PINS,**arithmetic.HELPER_PINS)


def receipt(precision=256,*,max_seconds=60,progress=None):
    if type(precision) is not int or precision not in (256,384):
        raise ValueError('Use fresh256/384-bit replay')
    if type(max_seconds) not in (int,float) or not 0<max_seconds<=60:
        raise ValueError('Positive bounded replay budget required')
    pins=checked_dependencies();started=perf_counter()
    def check_time():
        if perf_counter()-started>=max_seconds:
            raise TimeoutError('No complete high-arity receipt produced before budget')
    with ctx.workprec(precision):
        if progress:progress(dict(stage='fresh-regional-bank',precision_bits=precision))
        regions=arithmetic.uniform.regional_directed(MARK,ALPHA,R,60)
        check_time()
        terms=arithmetic.source_profile_terms(R,regions,check_time=check_time)
        prefix=arithmetic.high_prefix(R,A0,MARK,ALPHA)
        value=prefix+arithmetic.log2(sum((term[3] for term in terms),arb(0)))
        ceiling=arithmetic.uniform.rational_ceiling(value)
        rows=[dict(even_full=a,odd_full=b,indexed_source_count_upper=count,
                   failure_log2_upper=str(arithmetic.uniform.rational_ceiling(
                       prefix+arithmetic.log2(contribution))))
              for a,b,count,contribution in terms]
        check_time()
    if checked_dependencies()!=pins:raise ValueError('Dependencies changed during replay')
    return dict(SCOPE,precision_bits=precision,fresh_regional_banks=1,
        profiles=rows,profile_count=len(rows),failure_log2_upper=str(ceiling),
        failure_log2_enclosure=str(value),selected_class_reaches_40_bits=ceiling<-40,
        reviewed_source_sha256=pins,source_sha256=sha256(Path(__file__).read_bytes()).hexdigest(),
        elapsed_seconds=perf_counter()-started)


def checked_receipt(data,precision):
    expected=dict(SCOPE,precision_bits=precision,fresh_regional_banks=1,profile_count=99,
        reviewed_source_sha256=checked_dependencies(),
        source_sha256=sha256(Path(__file__).read_bytes()).hexdigest())
    arithmetic.require_fields(data,expected)
    expected_counts=arithmetic.profile_counts(R);observed={}
    for row in data.get('profiles',[]):
        a,b=row.get('even_full'),row.get('odd_full')
        if type(a) is not int or type(b) is not int or (a,b) not in expected_counts or (a,b) in observed:
            raise ValueError('Unexpected or repeated matching profile')
        arithmetic.require_fields(row,dict(indexed_source_count_upper=expected_counts[a,b]))
        if type(row.get('failure_log2_upper')) is not str:raise ValueError('Rational ceiling required')
        observed[a,b]=F(row['failure_log2_upper'])
    if set(observed)!=set(expected_counts) or type(data.get('failure_log2_upper')) is not str:
        raise ValueError('Incomplete selected class')
    return observed,F(data['failure_log2_upper'])


def verify(directory=DIRECTORY):
    checked=[];inputs=[]
    for precision in (256,384):
        path=Path(directory)/f'hd_lambda_high_arity_receipt{precision}.json'
        raw=path.read_bytes();data=json.loads(raw)
        checked.append(checked_receipt(data,precision))
        inputs.append(dict(file=path.name,sha256=sha256(raw).hexdigest(),precision_bits=precision))
    if checked[0]!=checked[1]:raise ValueError('Independent replay ceilings differ')
    profiles,raw=checked[0]
    rounded=arithmetic.rounded_union(profiles.values(),precision=512)
    if rounded<raw:raise ArithmeticError('Rounded profile union below raw enclosure')
    return dict(SCOPE,receipt_precisions=[256,384],sum_precision=512,
        all99_profile_ceilings_agree=True,raw_union_log2_upper=str(raw),
        rounded_profile_union_log2_upper=str(rounded),selected_class_reaches_40_bits=rounded<-40,
        reviewed_source_sha256=checked_dependencies(),inputs=inputs,
        trust='observed fresh directed replays; hashes alone do not prove arithmetic')


if __name__=='__main__':
    import argparse
    parser=argparse.ArgumentParser(description=__doc__)
    modes=parser.add_mutually_exclusive_group(required=True)
    modes.add_argument('--precision',type=int,choices=(256,384));modes.add_argument('--verify',action='store_true')
    parser.add_argument('--directory',type=Path,default=DIRECTORY)
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args()
    if args.output.exists():parser.error('Refuse to overwrite evidence')
    result=(verify(args.directory) if args.verify else receipt(args.precision,
        progress=lambda row:print(json.dumps(row),flush=True)))
    with args.output.open('x',encoding='utf-8') as handle:
        json.dump(result,handle,indent=2,allow_nan=False);handle.write('\n')
    print(json.dumps({key:value for key,value in result.items() if key not in
        ('profiles','reviewed_source_sha256')}),flush=True)
