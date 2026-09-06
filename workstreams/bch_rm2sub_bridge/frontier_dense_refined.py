"""Finer witness search for the unchanged outward dense interval inequality."""
import argparse
from fractions import Fraction as F
from pathlib import Path
import time
import numpy as np
from flint import arb,ctx
import frontier_dense as core


class Checker(core.Checker):
    def witness(self,lo,hi,vlo,vhi):
        center=super().witness(lo,hi,vlo,vhi)
        _,_,rlo,rhi=self.density_interval(lo,hi,vlo,vhi)
        weights=np.array(core.bernstein_weights(64,float((rlo+rhi)/2)))
        def score(tilt):
            lam,_,_,co=self.epoch(tilt)
            return core.terminal_float(np.einsum('i,ijk->jk',weights,co),self.power)+self.cutoff*float(lam)
        return min(range(center-2,center+3),key=score)


def run(output,m,first,last,slope,verify=False):
    base=core.base;output=output.resolve();saved=base.read(output) if verify else None
    if saved:
        m=saved['message_exponent'];first,last=saved['occupancy_range']
        for name,digest in saved['source_sha256'].items():assert base.sha(base.ROOT/name)==digest
        ps=[base.decode(p) for p in saved['row_probabilities']]
    else:
        assert not output.exists()
        ps,_=core.discovery.row_witnesses(core.caps_module.caps(),slope)
    core.restore_or_verify(base.ROOT);ctx.prec=512 if verify else 256
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
        hashes=core.source_hashes();hashes[Path(__file__).relative_to(base.ROOT).as_posix()]=base.sha(Path(__file__))
        base.write_new(output,dict(status='FRONTIER_DENSE_INTERVAL_CERTIFICATE',message_exponent=m,
            occupancy_range=[first,last],row_probabilities=[base.encode(p) for p in ps],
            per_occupancy_upper_power=-80,range_upper=base.encode(F(last-first+1,1<<80)),
            tree=nodes,leaves=leaves,source_sha256=hashes,full_distance_proved=False))
    print('Complete refined interval',first,last,'leaves',leaves,'seconds',round(time.monotonic()-started,2),flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True)
    p.add_argument('--m',type=int,default=24);p.add_argument('--first',type=int,default=4096)
    p.add_argument('--last',type=int,default=131072);p.add_argument('--slope',type=float,default=.75)
    p.add_argument('--verify',action='store_true');a=p.parse_args();run(a.output,a.m,a.first,a.last,a.slope,a.verify)
