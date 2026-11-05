"""Directed homogeneous-A12 shear bounds at several active-group counts.

The exact inventory fixes each original twelve-row set across all eight
column pairs. Summing its fourth-moment roots gives the raw group weight
B, including every labelled rowset and the RS comparison factor once.
For r independent active groups the outer factor is B**r/255.

For each requested r, compute ALL indexed matching-profile contributions
from a fresh uniform regional bank. The source union stays outside
theta=4/15. The bound is
 theta*(r*log2(B)-log2(255)) - alpha*W*log2(z) + log2(source_union).
Here alpha=1/5, z=691/700, and W=734003. No fitted mark or profile subset
is selected. Different r are disjoint classes, combined by a positive sum.

This extends only the homogeneous-A12 class of the SHEAR construction.
It does not transfer a certificate from another outer, cover mixed
arities, or assert a whole-code guarantee. No native operations.
"""
from fractions import Fraction as F
from hashlib import sha256
import json
from pathlib import Path
from time import perf_counter

from flint import arb, ctx

import hd_lambda_shear_moment_all_a12 as inventory

DIRECTORY = Path(__file__).resolve().parent
IDEAL, A = inventory.ideal, inventory.ideal.A
PINS = {'hd_lambda_shear_moment_all_a12.py':
    '3544a89d48894cf6c82ec478bfac4d96f3840df9f54d49b948e1c7d2f999651c'}
DEFAULT_GROUPS = tuple(range(16,21))


def dependencies():
    pins = inventory.dependencies()
    for name,digest in PINS.items():
        if sha256((DIRECTORY/name).read_bytes()).hexdigest() != digest:
            raise ValueError('Frozen A12 controller changed: '+name)
        if name in pins and pins[name] != digest:
            raise ValueError('Conflicting provenance: '+name)
        pins[name] = digest
    return dict(pins,**{Path(__file__).name:sha256(Path(__file__).read_bytes()).hexdigest()})


def checked_groups(groups):
    groups = tuple(groups)
    if (not groups or len(groups)>32 or tuple(sorted(set(groups))) != groups
            or any(type(r) is not int or not 1<=r<=32 for r in groups)):
        raise ValueError('Distinct increasing active-group counts in1..32 required')
    return groups


def log_bound(raw_group, r, uniform_source_log):
    checked_groups((r,))
    if not raw_group.is_finite() or not raw_group > 0 or not uniform_source_log.is_finite():
        raise ValueError('Positive finite raw weight and finite uniform log required')
    theta = A.uniform.rational(F(4,15))
    threshold = -A.uniform.rational(F(1,5))*734003*A.log2(A.uniform.rational(F(691,700)))
    return theta*(r*A.log2(raw_group)-A.log2(arb(255)))+threshold+uniform_source_log


def run(path,digest,*,groups=DEFAULT_GROUPS,precision=256,max_seconds=120,progress=None):
    groups = checked_groups(groups)
    if type(precision) is not int or precision not in (256,384):
        raise ValueError('Use independent256/384-bit directed replays')
    pins = dependencies()
    started,check = inventory._timer(max_seconds)
    tails,input_pin = inventory.checked_inventory(path,digest,check_time=check)
    records = []
    with ctx.workprec(precision):
        raw,details = inventory.raw_group_upper(tails)
        if progress:
            progress(dict(stage='fresh-regional-bank',maximum_active_groups=max(groups)))
        regions = A.uniform.regional_directed(F(691,700),F(1,5),max(groups),60)
        check()
        union = arb(0)
        for r in groups:
            check()
            terms = A.source_profile_terms(r,regions,check_time=check)
            reference = A.log2(sum((row[3] for row in terms),arb(0)))
            value = log_bound(raw,r,reference)
            upper = A.uniform.rational_ceiling(value)
            lower = -A.uniform.rational_ceiling(-value)
            union += arb(2)**A.uniform.rational(upper)
            row = dict(active_groups=r,homogeneous_arity=12,total_arity=12*r,
                matching_profiles=len(terms),all_matching_profiles=True,
                source_count_upper=sum(t[2] for t in terms),
                uniform_source_log2_enclosure=str(reference),
                failure_upper_log2_enclosure=str(value),
                failure_upper_log2_lower=str(lower),failure_upper_log2_upper=str(upper),
                selected_class_40_bits=bool(upper<=-40))
            records.append(row)
            if progress:
                progress(dict(row,elapsed_seconds=perf_counter()-started))
        union_upper = A.uniform.rational_ceiling(A.log2(union))
    check()
    if dependencies() != pins or sha256(Path(input_pin['file']).read_bytes()).hexdigest() != digest:
        raise ValueError('Frozen inputs changed during directed evaluation')
    return dict(complete=True,construction='scheduled-lambda-bc-staggered-column-shear-final-D2',
        requested_active_groups=list(groups),homogeneous_arity=12,entries=records,
        source_dimension=1048576,output_length=2097152,forbidden_weight=734003,
        uniform_alpha='1/5',outer_theta='4/15',mark='691/700',moment_order=4,
        all_labelled_original_rowsets=True,original_A_fixed_across_pairs=True,
        pair_count=8,minimum_active_pairs=5,exact_original_row_floors_retained=False,
        restricted_density_not_renormalized=True,source_count_outside_theta=True,
        rowset_sum_inside_theta=True,additional_rowset_binomial_factor=False,
        projective_normalizations=1,additional_beta_factor=False,
        full_independent_D1_comparison=True,permitted_gauges_per_pair=1,
        final_independent_D2_after_shear=True,arithmetic_enclosures_directed=True,
        precision_bits=precision,packet_mantissa_bits=60,fresh_regional_banks=1,
        regional_maximum=max(groups),mixed_arities_covered=False,
        unlisted_group_counts_covered=False,whole_code_certificate=False,
        failure_probability_upper_bound=True,selected_union_log2_upper=str(union_upper),
        selected_union_40_bits=bool(union_upper<=-40),raw_group_details=details,
        exact_inventory=input_pin,reviewed_source_sha256=pins,native_operations=0,
        elapsed_seconds=perf_counter()-started)


if __name__=='__main__':
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--inventory',type=Path,required=True)
    parser.add_argument('--inventory-sha256',required=True)
    parser.add_argument('--groups',type=int,nargs='+',default=DEFAULT_GROUPS)
    parser.add_argument('--precision',type=int,choices=(256,384),default=256)
    parser.add_argument('--max-seconds',type=float,default=120)
    parser.add_argument('--output',type=Path,required=True)
    args = parser.parse_args()
    if args.output.exists():
        parser.error('Refuse to overwrite evidence')
    result = run(args.inventory,args.inventory_sha256,groups=args.groups,precision=args.precision,
        max_seconds=args.max_seconds,progress=lambda row:print(json.dumps(row),flush=True))
    with args.output.open('x',encoding='utf-8') as handle:
        json.dump(result,handle,indent=2,allow_nan=False)
        handle.write('\n')
    print(json.dumps({k:v for k,v in result.items() if k not in
        ('entries','reviewed_source_sha256','raw_group_details')}),flush=True)
