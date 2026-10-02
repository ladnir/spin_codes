"""Compose authenticated count duals to improve their exact box correction.

A dual bounds c.x by a constant plus a residual dot x. The ordinary check
uses 0<=x_i<=1. If another checked dual bounds x_i, substitute that dual
for a positive residual coefficient. The resulting multipliers prove the
same objective against the original constraints; no saved cap is assumed.
"""
import argparse
import hashlib
import json
from fractions import Fraction as Q
from pathlib import Path
from math import log2

import shared_support
from constraint_counts import Constraints, dual_residual, verify_dual
from count_refinements import refine_pair
from dual_moments import dual_shell_caps


def compose(inequalities,equalities,objective,y,z,known,variables,rounds=4):
    """Return improved exact multipliers using known coordinate duals.

    Every coordinate dual is verified anew. All substitutions in a round
    use the prior residual, so each source coefficient remains nonnegative.
    Adding valid duals introduces no extra premise or circular reasoning.
    """
    if type(rounds) is not int or not 0<=rounds<=8:
        raise ValueError('zero through eight composition rounds required')
    upper,correction=verify_dual(inequalities,equalities,objective,y,z,variables)
    if upper<0:raise ArithmeticError('negative initial bound')
    sources={}
    for i,(yy,zz) in known.items():
        if type(i) is not int or not 0<=i<variables:
            raise ValueError('valid coordinate required')
        bound,_=verify_dual(inequalities,equalities,{i:Q(1)},yy,zz,variables)
        if bound<0:raise ArithmeticError('negative coordinate bound')
        if bound<1:sources[i]=(list(map(Q,yy)),list(map(Q,zz)),bound)
    y,z=list(map(Q,y)),list(map(Q,z))
    trace=[dict(round=0,upper=str(upper),correction=str(correction))]
    for iteration in range(1,rounds+1):
        _,residual=dual_residual(inequalities,equalities,objective,y,z,variables)
        selected=[(i,residual[i],yy,zz,bound) for i,(yy,zz,bound) in sources.items()
                  if residual[i]>0]
        if not selected:break
        yy=y[:];zz=z[:]
        predicted=upper-sum((r*(1-bound) for _,r,_,_,bound in selected),Q(0))
        for _,r,sy,sz,_ in selected:
            for j,value in enumerate(sy):
                if value:yy[j]+=r*value
            for j,value in enumerate(sz):
                if value:zz[j]+=r*value
        candidate,remaining=verify_dual(inequalities,equalities,objective,yy,zz,variables)
        if candidate<0:raise ArithmeticError('composed dual contradicts the premises')
        if candidate>predicted:
            raise ArithmeticError('residual substitution did not obey its exact bound')
        trace.append(dict(round=iteration,upper=str(candidate),correction=str(remaining),
                          substituted=[i for i,_,_,_,_ in selected]))
        if candidate>=upper:break
        y,z,upper,correction=yy,zz,candidate,remaining
    return y,z,upper,trace


def expand(witness,system):
    system.verify_witness(witness)
    y=[Q(0)]*len(system.inequalities);z=[Q(0)]*len(system.equalities)
    for entries,target in ((witness['y'],y),(witness['z'],z)):
        for i,value in entries:target[i]=Q(value)
    return y,z


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--count-witnesses',nargs='+',type=Path,required=True)
    parser.add_argument('--supports',type=int,nargs='+',default=[115,116])
    parser.add_argument('--last',type=int,default=192)
    parser.add_argument('--rounds',type=int,default=4)
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args()
    if args.output.exists():parser.error('refusing to overwrite a receipt')
    if not 128<=args.last<=256 or any(not 38<=u<=args.last for u in args.supports):
        parser.error('valid support range required')
    contents=[path.read_bytes() for path in args.count_witnesses]
    records=[json.loads(data) for data in contents]
    provenance=[dict(path=str(path),sha256=hashlib.sha256(data).hexdigest())
                for path,data in zip(args.count_witnesses,contents)]
    if any(r['schema']!='shared-gf16-constraint-counts-1' or r['last']!=args.last
           or r['complement_symmetry'] is not False for r in records):
        parser.error('matching ordinary joint-count models required')
    _,ranks=shared_support.shared_counts(True)
    spectrum=shared_support.authenticated_caps()
    primal,dual,dims,ddims=refine_pair(ranks,spectrum)
    system=Constraints(primal,dual,dims,ddims,128,args.last,
                       spectra=(spectrum,dual_shell_caps()),ones=(False,False))
    known={};initial={}
    for record in records:
        for row in record['rows']:
            if 'dual_witness' not in row:continue
            witness=row['dual_witness'];h,u=witness['rank'],witness['support']
            y,z=expand(witness,system)
            i=system.indices[0,h,u]
            upper,_=verify_dual(system.inequalities,system.equalities,{i:Q(1)},y,z,len(system.keys))
            if i not in initial or upper<initial[i][0]:
                known[i]=(y,z);initial[i]=(upper,witness)
    results=[]
    for u in args.supports:
        i=system.indices[0,4,u]
        if i not in known:raise ValueError('a starting witness is required for each support')
        y,z,upper,trace=compose(system.inequalities,system.equalities,{i:Q(1)},
                              *known[i],known,len(system.keys),args.rounds)
        old=initial[i][0];prior=system.caps[0][3][u]
        cap=prior*min(Q(1),upper).numerator//min(Q(1),upper).denominator
        witness=dict(rank=4,support=u,inequality_count=len(y),equality_count=len(z),
                     y=[[j,str(v)] for j,v in enumerate(y) if v],
                     z=[[j,str(v)] for j,v in enumerate(z) if v])
        if system.verify_witness(witness)!=cap:raise ArithmeticError('saved dual replay mismatch')
        gain=log2(old)-log2(upper) if upper else None
        print('EXACT RESIDUAL COMPOSITION',u,'extra gain bits',gain,'rounds',len(trace)-1,flush=True)
        results.append(dict(rank=4,support=u,prior=str(prior),cap=str(cap),
                            verified_ratio=str(upper),extra_gain_bits=gain,trace=trace,
                            dual_witness=witness,status='exact composed dual checked'))
        if any(hashlib.sha256(path.read_bytes()).hexdigest()!=p['sha256']
               for path,p in zip(args.count_witnesses,provenance)):
            raise RuntimeError('source witness changed during replay')
        args.output.parent.mkdir(parents=True,exist_ok=True)
        args.output.write_text(json.dumps(dict(schema='shared-gf16-constraint-counts-1',
            last=args.last,complement_symmetry=False,composition_rounds=args.rounds,
            source_witnesses=provenance,rows=results),indent=2)+'\n')


if __name__=='__main__':main()
