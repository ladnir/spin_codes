"""Outward boxes for the iid-domination bound; no sampled gaps are accepted."""
import argparse
from fractions import Fraction as F
from functools import lru_cache
import math
from pathlib import Path
import time

import numpy as np
from flint import arb,ctx
import bridge as base
import larger_state_maps as maps
import tightened_occupancy as transfer
import activation_bridge as matrix
import k30_sparse as caps_module
import k30_dense_unrestricted_screen as discovery
from k30_dense_screen import upper_hull
from frontier_sparse import parameters
from audit_bch_q1_full_arb import rational
from migrate_legacy_workspace import restore_or_verify


def aa(x):
    return arb(x.numerator)/x.denominator if isinstance(x,F) else arb(x)


def entropy(x):
    if x==0 or x==1:return arb(0)
    p=aa(x);return -p*p.log()-(1-p)*(1-p).log()


def hull_value(hull,x):
    for (a,b),(c,d) in zip(hull,hull[1:]):
        if a<=x<=c:return b+(d-b)*(x-a)/(c-a)
    raise AssertionError('outside hull')


def hull_max(hull,lo,hi):
    return max([hull_value(hull,lo),hull_value(hull,hi)]+[v for x,v in hull if lo<=x<=hi])


def terminal_float(mat,power):
    mat=mat.copy();scale=0.
    for _ in range(power):
        mat=mat@mat;maximum=float(mat.max());mat/=maximum
        scale=2*scale+math.log(maximum)
    value=float(mat[0].sum())
    return math.log(value)+scale if value>0 else math.inf


def monotone_parts(coefficients):
    """Exact decomposition c_j = ascending_j - descending_j."""
    ascending=[coefficients[0]];descending=[F(0)]
    for previous,current in zip(coefficients,coefficients[1:]):
        delta=current-previous
        ascending.append(ascending[-1]+max(F(0),delta))
        descending.append(descending[-1]+max(F(0),-delta))
    return ascending,descending


def bernstein_weights(t,r):
    return [math.comb(t,j)*r**j*(1-r)**(t-j) for j in range(t+1)]


class Checker:
    def __init__(self,m,ps):
        self.rows,self.cutoff=parameters(m);self.power=m-5
        self.t,self.s,self.spectrum,self.kernel=maps.load('t64_s20')
        counts=caps_module.caps()
        assert len(ps)==len(caps_module.BANDS) and ps[-1]==1 and all(0<p<1 for p in ps[:-1])
        costs=[rational(aa(caps_module.density_cost(b,p,counts)).log().upper())
               for b,p in zip(caps_module.BANDS[:-1],ps[:-1])]+[F(0)]
        assert all(c>=0 for c in costs)
        self.hull=upper_hull(list(zip(ps,costs)));self.pmin=min(ps)
        self.ps=ps;self.calls=0
        self.epoch=lru_cache(maxsize=160)(self._epoch)
        self.count_cost=lru_cache(maxsize=4096)(self._count_cost)

    def _epoch(self,tilt):
        lam=(arb(tilt)/10).exp()
        co=transfer.epoch_matrices(self.t,self.s,self.spectrum,self.kernel,(-lam).exp(),arb,self.t)
        co=[tuple(v.upper() for v in row) for row in co]
        parts=[monotone_parts([rational(row[k]) for row in co]) for k in range(9)]
        ascending=[[aa(v) for v in a] for a,d in parts]
        descending=[[aa(v) for v in d] for a,d in parts]
        return lam,ascending,descending,np.array([[float(v) for v in row] for row in co]).reshape(-1,3,3)

    def _count_cost(self,lo,hi):
        r=F(1,2) if 2*lo<=self.rows<=2*hi else F(hi if 2*hi<self.rows else lo,self.rows)
        return self.rows*entropy(r)+hi*arb(len(caps_module.BANDS)).log()+256*arb(self.rows+1).log()

    def density_interval(self,lo,hi,vlo,vhi):
        plo=self.pmin+(1-self.pmin)*vlo;phi=self.pmin+(1-self.pmin)*vhi
        return plo,phi,F(lo,self.rows)*plo,F(hi,self.rows)*phi

    def witness(self,lo,hi,vlo,vhi):
        _,_,rlo,rhi=self.density_interval(lo,hi,vlo,vhi)
        r=float((rlo+rhi)/2);x=(lo+hi)/(2*self.rows)
        weights=np.array([math.comb(64,j)*r**j*(1-r)**(64-j) for j in range(65)])
        prediction=round(10*math.log(x));best=None
        for tilt in range(prediction-15,prediction+29,3):
            lam,co,derivative,ap=self.epoch(tilt)
            value=terminal_float(np.einsum('i,ijk->jk',weights,ap),self.power)+self.cutoff*float(lam)
            if best is None or value<best[0]:best=value,tilt
        return best[1]

    def bound(self,lo,hi,vlo,vhi,tilt):
        self.calls+=1
        plo,phi,rlo,rhi=self.density_interval(lo,hi,vlo,vhi)
        lam,ascending,descending,_=self.epoch(tilt)
        low=bernstein_weights(64,aa(rlo));high=bernstein_weights(64,aa(rhi))
        # Increasing coefficient sequences have increasing binomial means.
        # Bound their difference by the high ascending mean minus the low
        # descending mean. No global derivative floor is introduced.
        mat=tuple(max(arb(0),(sum((w*c for w,c in zip(high,a)),arb(0))-
                              sum((w*c for w,c in zip(low,d)),arb(0))).upper())
                  for a,d in zip(ascending,descending))
        for _ in range(self.power):mat=matrix.positive_mul(mat,mat)
        return (self.count_cost(lo,hi)+hi*aa(hull_max(self.hull,plo,phi))+
                self.cutoff*lam+sum(mat[:3],arb(0)).log()).upper()


def source_hashes():
    paths=[Path(__file__),Path(transfer.__file__),Path(transfer.original.__file__),Path(matrix.__file__),
        Path(caps_module.__file__),Path(discovery.__file__),base.HERE/'k30_dense_screen.py',
        base.HERE/'frontier_sparse.py',base.HERE/'general_batch_certificate.py',
        base.HERE/'occupation_three.py',base.HERE/'MIGRATION_MANIFEST.json']
    return {p.relative_to(base.ROOT).as_posix():base.sha(p) for p in paths}


def run(output,m,first,last,slope,verify=False):
    output=output.resolve();saved=base.read(output) if verify else None
    if saved:
        m=saved['message_exponent'];first,last=saved['occupancy_range']
        for name,digest in saved['source_sha256'].items():assert base.sha(base.ROOT/name)==digest
        ps=[base.decode(p) for p in saved['row_probabilities']]
    else:
        assert not output.exists()
        ps,_=discovery.row_witnesses(caps_module.caps(),slope)
    restore_or_verify(base.ROOT);ctx.prec=512 if verify else 256
    checker=Checker(m,ps);assert 2<=first<=last<=checker.rows
    threshold=-80*arb(2).log();nodes=[];leaves=0;started=time.monotonic();last_print=started
    stack=[(first,last,F(0),F(1),0)];cursor=0
    while stack:
        lo,hi,vlo,vhi,depth=stack.pop()
        if saved:
            assert cursor<len(saved['tree']);code=saved['tree'][cursor];cursor+=1
        else:
            tilt=checker.witness(lo,hi,vlo,vhi)
            bound=checker.bound(lo,hi,vlo,vhi,tilt)
            if bound<threshold:code=tilt
            else:
                assert depth<45,f'Unclosed box {lo,hi,vlo,vhi}, log bound {bound}'
                # Separate q and the normalized active density. q splits are
                # disjoint integer ranges; density splits share only an endpoint.
                q_width=(hi-lo)/max(1,lo)
                p_width=float((1-checker.pmin)*(vhi-vlo)/(checker.pmin+(1-checker.pmin)*vlo))
                code='q' if lo<hi and q_width>=p_width else 'p'
        nodes.append(code)
        if code=='q':
            assert lo<hi;mid=(lo+hi)//2
            stack.append((mid+1,hi,vlo,vhi,depth+1));stack.append((lo,mid,vlo,vhi,depth+1))
        elif code=='p':
            mid=(vlo+vhi)/2
            stack.append((lo,hi,mid,vhi,depth+1));stack.append((lo,hi,vlo,mid,depth+1))
        else:
            assert isinstance(code,int)
            if saved:assert checker.bound(lo,hi,vlo,vhi,code)<threshold
            leaves+=1
        if time.monotonic()-last_print>15:
            print('replay' if saved else 'certify','leaves',leaves,'pending',len(stack),'Q',lo,hi,flush=True)
            last_print=time.monotonic()
    if saved:
        assert cursor==len(saved['tree']) and leaves==saved['leaves']
        base.write_new(output.with_name(output.stem+'_replay.json'),dict(status='FRONTIER_DENSE_512_BIT_REPLAY_PASSED',
            producer_sha256=base.sha(output),message_exponent=m,occupancy_range=[first,last],leaves=leaves))
    else:
        base.write_new(output,dict(status='FRONTIER_DENSE_INTERVAL_CERTIFICATE',message_exponent=m,
            occupancy_range=[first,last],row_probabilities=[base.encode(p) for p in ps],
            per_occupancy_upper_power=-80,range_upper=base.encode(F(last-first+1,1<<80)),
            tree=nodes,leaves=leaves,source_sha256=source_hashes(),full_distance_proved=False))
    print('Complete dense interval',first,last,'leaves',leaves,'seconds',round(time.monotonic()-started,2),flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True)
    p.add_argument('--m',type=int,default=24);p.add_argument('--first',type=int,default=4096)
    p.add_argument('--last',type=int,default=131072);p.add_argument('--slope',type=float,default=.75)
    p.add_argument('--verify',action='store_true');a=p.parse_args();run(a.output,a.m,a.first,a.last,a.slope,a.verify)
