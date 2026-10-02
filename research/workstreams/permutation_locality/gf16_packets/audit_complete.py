"""Exact scope and sum audit of a completed K=2^20 GF16 replay receipt.

This checks receipt consistency, not the underlying numerical inequalities.
Run assemble.py to regenerate those inequalities from the actual construction.
An optional independently regenerated sparse prefix must match exactly.
"""
import argparse
from fractions import Fraction as Q
import hashlib
import json
from pathlib import Path


def dyadic(pair):
    if (not isinstance(pair,(list,tuple)) or len(pair)!=2
            or any(type(x) is not int for x in pair) or pair[0]<0):
        raise ValueError('nonnegative integer dyadic endpoint required')
    return Q(pair[0])*Q(2)**pair[1]


def audit(report,dense_raw,distance,bits,prefix=None):
    """Check the requested strict distance and margin using exact arithmetic."""
    distance=Q(distance)
    if not 0<distance<1 or type(bits) is not int or bits<40:
        raise ValueError('distance in (0,1) and integer margin >=40 required')
    dense=json.loads(dense_raw)
    if (report.get('schema')!='gf16-complete-replay-1'
            or dense.get('schema')!='gf16-packet-scalar-cover-1'
            or report.get('dense_sha256')!=hashlib.sha256(dense_raw).hexdigest()):
        raise ValueError('receipt schema or dense input hash mismatch')
    if (type(report.get('precision')) is not int or report['precision']<128
            or type(report.get('updates')) is not int or report['updates'] not in (2,3,4)
            or type(report.get('threshold')) is not int or not 0<=report['threshold']<1<<21):
        raise ValueError('invalid replay claim')
    for key in ('updates','threshold'):
        if type(dense.get(key)) is not int or report[key]!=dense[key]:
            raise ValueError('dense and full claim mismatch')
    if type(dense.get('minimum_groups')) is not int or not 1<=dense['minimum_groups']<=2048:
        raise ValueError('invalid dense occupancy boundary')
    expected=dict(output_length=1<<21,minimum_distance=report['threshold']+1,
                  dense_through=2048,sparse_through=dense['minimum_groups']-1)
    if (any(type(report.get(k)) is not int or report[k]!=v for k,v in expected.items())
            or dense.get('unresolved')!={} or dense.get('screen_only',False)):
        raise ValueError('incomplete or inconsistent occupancy scope')
    leaves=dense.get('leaves')
    if (not isinstance(leaves,dict) or not leaves
            or any(set(p)-{'0','1'} for p in leaves)
            or any(not isinstance(v,dict) or not isinstance(v.get('witness'),dict) for v in leaves.values())):
        raise ValueError('binary paths with witnesses required')
    if (any(p[:i] in leaves for p in leaves for i in range(len(p)))
            or sum((Q(1,2)**len(p) for p in leaves),Q(0))!=1):
        raise ValueError('dense partition overlaps or leaves a gap')
    if Q(report['minimum_distance'],report['output_length'])<=distance:
        raise ValueError('requested strict distance not achieved')
    dense_upper=dyadic(report['dense_upper']);sparse_upper=dyadic(report['sparse_upper'])
    total=dyadic(report['total_upper'])
    if (dense_upper<=0 or (expected['sparse_through']>0 and sparse_upper<=0)
            or not dense_upper+sparse_upper<=total<Q(2)**-bits):
        raise ValueError('invalid outward sum or insufficient margin')
    if prefix is not None:
        parameters=prefix.get('parameters',{})
        if (prefix.get('schema')!='gf16-fresh-sparse-prefix-1'
                or prefix.get('covered_occupancies')!=[1,report['sparse_through']]
                or any(parameters.get(k)!=report.get(k) for k in
                       ('updates','threshold','precision','sparse_inner','exact_feedback',
                        'single_group_exact','sparse_joint_return_through','sparse_lazy_density_through'))
                or dyadic(prefix['upper'])!=sparse_upper):
            raise ValueError('independent sparse prefix scope or endpoint differs')
    return len(leaves)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('receipt',type=Path)
    parser.add_argument('dense',type=Path)
    parser.add_argument('--distance',required=True,help='Strict lower bound, as an exact fraction or decimal')
    parser.add_argument('--bits',type=int,default=40)
    parser.add_argument('--prefix',type=Path,help='Optional independent same-precision sparse prefix receipt')
    args=parser.parse_args()
    report=json.loads(args.receipt.read_bytes())
    prefix=None if args.prefix is None else json.loads(args.prefix.read_bytes())
    count=audit(report,args.dense.read_bytes(),args.distance,args.bits,prefix)
    print(f'PASS: matching full scope, {count} disjoint exhaustive dense intervals, '
          f'distance > {Q(args.distance)}, and total upper < 2^-{args.bits}.')
    if prefix is not None:print('Independent sparse prefix endpoint matches exactly.')
    print('Receipt consistency only; use assemble.py for fresh numerical proof replay.')


if __name__=='__main__':main()
