"""Arb replay of the K18 two-step routing-conditioned fractional bound.

Imports authenticated outward physical endpoints and the retained complete
certificate. Recomputes new global endpoints at256-bit precision, replaces
only improved occupancy bounds, and sums all2048 dyadics exactly.
"""
from __future__ import annotations

import argparse
from fractions import Fraction as F
from math import comb, log2
from pathlib import Path
import hashlib
import json
import sys
from time import monotonic

from flint import arb, arb_mat, ctx

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[2]
sys.path.insert(0,str(HERE.parent/'k16_design'))
import packet_regional_power as placement

MAP='e2467b9d5de8403ccd8263f0f55c73df128a8a2dd4c4cc405bc3a8e057fc3581'
OLD_SHA='20cd124704ff83d7ac0e098f5b80ce97e59aabf48f547bd118b1b3b957d5b60b'
OLD=ROOT/'tmp/rs16-s22-k262144-whole-fresh-p256.json'


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def authenticate(pins):
    for path,expected in pins.items():
        if digest(path)!=expected:
            raise ArithmeticError(f'pin mismatch: {path}')


def aq(x):
    x=F(x)
    return arb(x.numerator)/x.denominator


def fraction(pair):
    n,e=pair
    if type(n)!=int or type(e)!=int or n<0:
        raise ValueError('nonnegative integer dyadic pair required')
    return F(n)*F(2)**e


def upper(x):
    if not x.is_finite() or not x>=0:
        raise ArithmeticError('finite nonnegative Arb enclosure required')
    return x.upper()


def endpoint(x):
    return list(map(int,upper(x).man_exp()))


def encode_exact(x):
    """An exact sum of dyadics is dyadic; do not round the final union."""
    x=F(x)
    den=x.denominator
    if den&(den-1):
        raise ArithmeticError('non-dyadic sum')
    return dict(numerator_hex=hex(x.numerator),exponent=-(den.bit_length()-1))


def bits(x):
    x=F(x)
    return log2(x.denominator)-log2(x.numerator)


def load_local(path):
    record=json.loads(Path(path).read_text(encoding='utf-8'))
    authenticate(record['source_sha256'])
    if (record.get('schema')!='k18-certified22-outward-physical-export-1'
            or record['map_record']['map_sha256']!=MAP
            or record['precision']<192):
        raise ValueError('authenticated exact certified22 physical export required')
    raw=record['matrices']
    if len(raw)!=17 or len(raw[0])!=11 or any(len(matrix)!=11 or any(len(row)!=11 for row in matrix) for matrix in raw):
        raise ValueError('complete physical0..16 / eleven-coordinate family required')
    matrices=[arb_mat([[aq(fraction(x)) for x in row] for row in matrix]) for matrix in raw]
    if matrices[0][0,0]!=1 or any(matrices[0][0,j]!=0 for j in range(1,11)):
        raise ValueError('zero input and zero state must remain zero')
    return matrices,F(record['tilt']),record


def thin(physical):
    size=physical[0].nrows()
    result=[]
    for j in range(17):
        value=arb_mat(size,size)
        for k in range(j+1):
            value+=aq(F(comb(j,k)*15**k,16**j))*physical[k]
        result.append(placement.rounded(value))
    return result


def fractional_blocks(physical,alpha):
    """No normalization needed: Arb supports arbitrarily small exponents."""
    alpha=F(alpha)
    if not 0<alpha<=1:
        raise ValueError('alpha in(0,1] required')
    exponent=aq(alpha)
    potential=thin(physical)
    size=potential[0].nrows()
    result=[]
    for total in range(33):
        value=arb_mat(size,size)
        for a in range(max(0,total-16),min(16,total)+1):
            product=placement.rounded(potential[a]*potential[total-a])
            powered=arb_mat([[arb(0) if product[i,j]==0 else upper(product[i,j]**exponent)
                             for j in range(size)] for i in range(size)])
            value+=aq(F(comb(16,a)*comb(16,total-a),comb(32,total)))*powered
        result.append(placement.rounded(value))
    return result


def scalar_bound(regional,q,tilt,alpha):
    beta=F(1<<256,65535**8)
    power=regional**64
    moment=sum((power[0,j] for j in range(power.ncols())),arb(0))
    prefactor=(aq(alpha)*(q*aq(beta).log()+aq(tilt)*52428)).exp()
    return upper(arb(comb(2048,q))*prefactor*moment)


def run(local,output,minimum=62,maximum=143,alpha=F(4,5),precision=256):
    if output.exists() or not 3<=minimum<=maximum<=384 or precision<256:
        raise ValueError('fresh output, qsubset3..384, precision>=256 required')
    start=monotonic()
    ctx.prec=precision
    physical,tilt,source=load_local(local)
    if digest(OLD)!=OLD_SHA:
        raise ArithmeticError('retained whole receipt changed')
    old=json.loads(OLD.read_text(encoding='utf-8'))
    authenticate(old['source_sha256'])
    if (old['map_record']['map_sha256']!=MAP or old['K']!=262144 or old['threshold']!=52428
            or not old['all_occupancies_covered'] or not old['whole_code_certificate']
            or set(old['occupancy_uppers'])!={str(q) for q in range(1,2049)}):
        raise ValueError('complete retained K18 certificate required')
    pins=dict(old['source_sha256'])
    for path,value in source['source_sha256'].items():
        if path in pins and pins[path]!=value:
            raise ArithmeticError('conflicting pins')
        pins[path]=value
    for path in (Path(local).resolve(),OLD,Path(__file__).resolve(),Path(placement.__file__).resolve(),HERE/'SLACK_REVIEW.md'):
        pins[str(path)]=digest(path)
    blocks=fractional_blocks(physical,alpha)
    print(f'Outward G2 operators ready in {monotonic()-start:.3f}s',flush=True)
    regionals=placement.placement_power(blocks,epochs=64,windows=32,maximum_groups=maximum)
    print(f'Outward regional placement through q{maximum} ready in {monotonic()-start:.3f}s',flush=True)
    candidates={str(q):endpoint(scalar_bound(regionals[q],q,tilt,alpha)) for q in range(minimum,maximum+1)}
    chosen={}
    improvements=[]
    for q,retained in old['occupancy_uppers'].items():
        candidate=candidates.get(q)
        if candidate is not None and fraction(candidate)<fraction(retained):
            chosen[q]=candidate
            improvements.append(int(q))
        else:
            chosen[q]=retained
    total=sum(map(fraction,chosen.values()),F(0))
    tail=sum((fraction(value) for q,value in chosen.items() if int(q)>=3),F(0))
    target=total<F(1,1<<60)
    authenticate(pins)
    result=dict(schema='k18-routing-conditioned-fractional-whole-1',
        whole_code_certificate=True,all_occupancies_covered=True,
        certificate_kind='authenticated retained endpoints plus fresh outward conditional bound',
        K=262144,N=524288,distance='1/10',threshold=52428,minimum_distance=52429,
        map_sha256=MAP,outer_group_dimension=128,packet_bits=4,state_bits=22,
        regions=64,slots_per_region=2048,physical_step_bits=64,
        physical_steps_per_macro=2,macro_slots=32,macros_per_region=64,
        continuous_state=True,zero_initial_state=True,final_flush=False,
        alpha=str(alpha),tilt=str(tilt),precision=precision,
        exact_geometry=True,thinning_before_fractional_power=True,
        chronological_product_before_power=True,placement_weights_outside_power=True,
        group_union_outside_power=True,local_receipt=str(Path(local).resolve()),
        retained_receipt=str(OLD),retained_receipt_sha256=OLD_SHA,
        local_census_recomputed_by_this_replay=False,
        new_occupancy_bounds=candidates,selected_occupancy_bounds=chosen,
        improved_occupancies=sorted(improvements),exact_union=encode_exact(total),
        exact_tail_union=encode_exact(tail),margin_bits_display=bits(total),tail_margin_bits_display=bits(tail),
        target_bits=60,target_met=target,source_sha256=pins,
        elapsed_seconds=monotonic()-start,
        scope='Whole-code setup-failure bound under the retained independent ideal ensemble; '
              'not a guarantee for every seeded instance. No construction or encoder change.')
    output.parent.mkdir(parents=True,exist_ok=True)
    with output.open('x',encoding='utf-8') as handle:
        json.dump(result,handle,indent=2)
        handle.write('\n')
    print(json.dumps({k:result[k] for k in ('target_met','margin_bits_display','tail_margin_bits_display','elapsed_seconds')}),flush=True)
    return result


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--local',type=Path,required=True)
    parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--min-q',type=int,default=62)
    parser.add_argument('--max-q',type=int,default=143)
    parser.add_argument('--alpha',type=F,default=F(4,5))
    args=parser.parse_args()
    run(args.local,args.output,args.min_q,args.max_q,args.alpha)
