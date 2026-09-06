"""Parameterized finite-size certificates using the frozen positive recurrence."""
import argparse
from fractions import Fraction as F
from functools import lru_cache
import math
from pathlib import Path

import numpy as np
from scipy.optimize import minimize_scalar
from flint import arb,ctx
import bridge as base
import larger_state_maps as maps
import polynomial_regions as poly
import scaled_adaptive as scaled
import positive_line_hull as hull
import activation_bridge as matrix
import certify_refresh_q1 as refresh
from certify_k30_q1 import polynomial_region
import k30_sparse as caps_module
from occupation_three import BANDS
from screen_constant_row_split import moment
from certify_larger_state_range import costs_for
from close_larger_state_gap import dyadic_power,as_fraction
from audit_bch_q1_full_arb import rational
from migrate_legacy_workspace import restore_or_verify


def parameters(m):
    assert 20<=m<=30
    rows=1<<(m-7)
    return rows,256*rows//10


def sources():
    files=[Path(__file__),Path(poly.__file__),Path(scaled.__file__),Path(hull.__file__),
        Path(matrix.__file__),Path(refresh.__file__),Path(caps_module.__file__),
        base.HERE/'certify_k30_q1.py',base.HERE/'general_occupancy.py',
        base.HERE/'tightened_occupancy.py',base.HERE/'general_batch_certificate.py',
        base.HERE/'certify_larger_state_range.py',base.HERE/'close_larger_state_gap.py',
        base.HERE/'screen_constant_row_split.py',base.HERE/'occupation_three.py',
        base.HERE/'MIGRATION_MANIFEST.json']
    return {p.relative_to(base.ROOT).as_posix():base.sha(p) for p in files}


def choose(region,q,tilt,rows,cutoff,caps):
    selected=np.array([[float(v.log()) if v>0 else -math.inf for v in row]
                       for row in region[:q+1]]).reshape(-1,3,3)
    ps,gammas=[],[]
    for band in BANDS:
        if band==(256,):ps.append(1.);gammas.append(0.);continue
        w=np.array(band);v=np.array([math.log(caps[x])-math.log(math.comb(256,x)) for x in band])
        def objective(theta):
            p=1/(1+math.exp(-theta))
            g=float(np.max(v-w*math.log(p)-(256-w)*math.log1p(-p)))
            return moment(selected,q,p)+q*g
        fit=minimize_scalar(objective,bounds=(-6,14),method='bounded',options={'xatol':1e-6})
        assert fit.success
        p=1/(1+math.exp(-float(fit.x)));ps.append(p)
        gammas.append(float(np.max(v-w*math.log(p)-(256-w)*math.log1p(-p))))
    ps=np.array(ps);roots=np.exp(np.array(gammas)/256)
    keep=hull.indices(roots*(1-ps),roots*ps)
    mantissas,exponents=scaled.initial(region[:q+1],outward=False)
    value=scaled.terminal_log(mantissas,exponents,ps[keep],roots[keep])
    value+=cutoff*math.exp(tilt/10)+math.log(math.comb(rows,q))+q*math.log(len(BANDS))
    return dict(anchor=q,tilt=tilt,p=[base.encode(F.from_float(float(p))) for p in ps],margin_bits_diagnostic=-value/math.log(2))


def evaluate(region,witness,lower,upper,rows,cutoff,caps):
    ps=[base.decode(p) for p in witness['p']]
    costs=costs_for(BANDS,ps,caps)
    roots=[((arb(c.numerator)/c.denominator).log()/256).exp().upper() for c in costs]
    probs=[arb(p.numerator)/p.denominator for p in ps]
    left=np.nextafter(np.array([float((r*(1-p)).upper()) for r,p in zip(roots,probs)]),np.inf)
    right=np.nextafter(np.array([float((r*p).upper()) for r,p in zip(roots,probs)]),np.inf)
    keep=hull.indices(left,right);mantissas,exponents=scaled.initial(region[:upper+1])
    correction=(cutoff*(arb(witness['tilt'])/10).exp()).exp()
    results=[]
    for q,mat,exponent in scaled.matrices(mantissas,exponents,left[keep],right[keep]):
        if q<lower:continue
        value=tuple(arb(float(v)) for v in mat.flat)
        for _ in range(8):value=matrix.positive_mul(value,value)
        bound=math.comb(rows,q)*len(BANDS)**q*rational((sum(value[:3],arb(0))*correction*arb(2)**(256*exponent)).upper())
        results.append(dict(occupation=q,upper_power=dyadic_power(bound),
            margin_bits_diagnostic=math.log2(bound.denominator)-math.log2(bound.numerator)))
    return results


def run(output,m,first,last,verify=False):
    output=output.resolve();saved=base.read(output) if verify else None
    if saved:
        m=saved['message_exponent'];first,last=saved['occupancy_range']
        for name,digest in saved['source_sha256'].items():assert base.sha(base.ROOT/name)==digest
    else:assert not output.exists()
    rows,cutoff=parameters(m);assert 2<=first<=last<=min(rows,8192)
    restore_or_verify(base.ROOT);ctx.prec=512 if verify else 256
    t,s,spectrum,kernel=maps.load('t64_s20');caps=caps_module.caps()
    @lru_cache(maxsize=6)
    def region(tilt,degree):
        return poly.regions(t,s,spectrum,kernel,(-(arb(tilt)/10).exp()).exp(),degree,length=rows)
    best={};shards=[]
    if saved:
        for shard in saved['shards']:
            lo,hi=shard['occupancy_range'];witness=shard['witness']
            values=evaluate(region(witness['tilt'],hi),witness,lo,hi,rows,cutoff,caps)
            assert len(values)==len(shard['rows'])
            for actual,old in zip(values,shard['rows']):
                assert actual['occupation']==old['occupation'] and actual['upper_power']<=old['upper_power']
                if old['upper_power']<=-70:best[old['occupation']]=min(best.get(old['occupation'],0),old['upper_power'])
            print('replay',lo,hi,flush=True)
        assert sorted(best)==list(range(first,last+1))
        assert [best[q] for q in sorted(best)]==saved['upper_powers']
        assert sum((as_fraction(v) for v in best.values()),F(0))==base.decode(saved['range_upper'])
        base.write_new(output.with_name(output.stem+'_replay.json'),dict(status='FRONTIER_SPARSE_512_BIT_REPLAY_PASSED',
            producer_sha256=base.sha(output),message_exponent=m,occupancy_range=[first,last]))
        return
    while len(best)<last-first+1:
        anchor=next(q for q in range(first,last+1) if q not in best)
        hi=min(last,max(anchor,int(anchor*1.22)))
        prediction=round((-76+6*math.log2(anchor)-10*(m-20)*math.log(2))/5)*5
        trials=[choose(region(tilt,hi),anchor,tilt,rows,cutoff,caps)
                for tilt in sorted(set(prediction+d for d in (-5,0,5)))]
        witness=max(trials,key=lambda r:r['margin_bits_diagnostic'])
        values=evaluate(region(witness['tilt'],hi),witness,anchor,hi,rows,cutoff,caps)
        shard=dict(occupancy_range=[anchor,hi],witness=witness,rows=values)
        shards.append(shard)
        for row in values:
            if row['upper_power']<=-70:best[row['occupation']]=min(best.get(row['occupation'],0),row['upper_power'])
        print('certify',anchor,hi,'margin',round(values[0]['margin_bits_diagnostic'],4),'covered',len(best),flush=True)
        if anchor not in best:
            base.write_new(output.with_name(output.stem+'_failed.json'),dict(status='FAILED_ANCHOR_NOT_CERTIFIED',
                message_exponent=m,anchor=anchor,shards=shards,source_sha256=sources()))
            return
    total=sum((as_fraction(best[q]) for q in best),F(0))
    base.write_new(output,dict(status='FRONTIER_SPARSE_OUTWARD_CERTIFICATE',message_exponent=m,
        parameters=dict(outer_rows=rows,cutoff=cutoff,step_bits=t,state_bits=s),occupancy_range=[first,last],
        upper_powers=[best[q] for q in sorted(best)],range_upper=base.encode(total),shards=shards,
        source_sha256=sources(),full_distance_proved=False))


def q1(output,m,verify=False):
    output=output.resolve();saved=base.read(output) if verify else None
    if saved:
        m=saved['message_exponent']
        for name,digest in saved['source_sha256'].items():assert base.sha(base.ROOT/name)==digest
    else:assert not output.exists()
    rows,cutoff=parameters(m);restore_or_verify(base.ROOT);ctx.prec=512 if verify else 256
    t,s,spectrum,_=maps.load('t64_s20');best={w:F(rows) for w in base.WEIGHTS}
    for a in (264,295,332,376):
        lam=arb(a)/(100*rows);epoch=refresh.epoch(t,s,spectrum,(-lam).exp(),arb)
        reg=polynomial_region(*epoch,rows//t) if verify else refresh.region(*epoch,rows//t,arb)
        co=refresh.coefficients(*reg,256,arb);factor=rows*(cutoff*lam).exp()
        for w in best:best[w]=min(best[w],rational((co[w]*factor).upper()))
    bound,factor,rest=base.bch_bound(best)
    assert bound<F(1,1<<40)
    if saved:
        assert all(best[w]<=base.decode(saved['coefficient_upper'][str(w)]) for w in best)
        assert base.bch_bound({int(w):base.decode(v) for w,v in saved['coefficient_upper'].items()})==tuple(
            base.decode(saved[k]) for k in ('Q1_upper','factor','rest'))
        base.write_new(output.with_name(output.stem+'_replay.json'),dict(status='FRONTIER_Q1_512_BIT_REPLAY_PASSED',
            producer_sha256=base.sha(output),message_exponent=m,coefficients_checked=92))
    else:
        base.write_new(output,dict(status='FRONTIER_Q1_OUTWARD_CERTIFICATE',message_exponent=m,
            parameters=dict(outer_rows=rows,cutoff=cutoff,step_bits=t,state_bits=s),Q1_upper=base.encode(bound),
            factor=base.encode(factor),rest=base.encode(rest),coefficient_upper={str(w):base.encode(v) for w,v in best.items()},
            margin_bits=math.log2(bound.denominator)-math.log2(bound.numerator),source_sha256=sources(),full_distance_proved=False))
    print('Q1',m,'margin',math.log2(bound.denominator)-math.log2(bound.numerator),flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True)
    p.add_argument('--m',type=int,default=24);p.add_argument('--first',type=int,default=2)
    p.add_argument('--last',type=int,default=1024);p.add_argument('--q1',action='store_true');p.add_argument('--verify',action='store_true')
    a=p.parse_args()
    if a.q1:q1(a.output,a.m,a.verify)
    else:run(a.output,a.m,a.first,a.last,a.verify)
