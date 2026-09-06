"""Retain a tighter dense budget, refining a frozen certificate's interval tree."""
import argparse
from fractions import Fraction as F
from pathlib import Path
import time
from flint import arb,ctx
import frontier_dense_density as density
import frontier_ledger as ledger


def run(output,source,bits=100,verify=False):
    core=density.core;base=core.base;output=output.resolve()
    saved=base.read(output) if verify else None
    if saved:
        source=base.ROOT/saved['seed_certificate'];bits=-saved['per_occupancy_upper_power']
        for name,digest in saved['source_sha256'].items():assert base.sha(base.ROOT/name)==digest
    else:assert not output.exists()
    assert 80<=bits<=256
    source=source.resolve();seed=base.read(source);replay=source.with_name(source.stem+'_replay.json')
    assert seed['status']=='FRONTIER_DENSE_INTERVAL_CERTIFICATE'
    assert seed['method']=='sqrt_density_domination' and seed['per_occupancy_upper_power']==-80
    assert base.read(replay)['producer_sha256']==base.sha(source)
    for name,digest in seed['source_sha256'].items():assert base.sha(base.ROOT/name)==digest
    m=seed['message_exponent'];first,last=seed['occupancy_range']
    assert ledger.validate_tree(seed['tree'],first,last)==seed['leaves']
    if saved:
        assert saved['message_exponent']==m and saved['occupancy_range']==[first,last]
        assert saved['row_probabilities']==seed['row_probabilities']
    core.restore_or_verify(base.ROOT);ctx.prec=512 if saved else 256
    checker=density.Checker(m,[base.decode(p) for p in seed['row_probabilities']])
    threshold=-bits*arb(2).log();nodes=[];leaves=0;cursor=0;started=time.monotonic();printed=started
    tree=saved['tree'] if saved else seed['tree']
    stack=[(first,last,F(0),F(1),0,True)]
    while stack:
        lo,hi,vlo,vhi,depth,in_tree=stack.pop()
        if in_tree:
            assert cursor<len(tree);code=tree[cursor];cursor+=1
        else:code=checker.witness(lo,hi,vlo,vhi)
        if type(code) is int:
            bound=checker.bound(lo,hi,vlo,vhi,code)
            if saved:assert bound<threshold
            elif not bound<threshold:
                tilt=checker.witness(lo,hi,vlo,vhi)
                if tilt!=code:bound=checker.bound(lo,hi,vlo,vhi,tilt)
                code=tilt
                if not bound<threshold:
                    assert depth<50,f'Unclosed box {lo,hi,vlo,vhi}'
                    qw=(hi-lo)/max(1,lo)
                    pw=float((1-checker.pmin)*(vhi-vlo)/(checker.pmin+(1-checker.pmin)*vlo))
                    code='q' if lo<hi and qw>=pw else 'p'
                    in_tree=False
        nodes.append(code)
        if code=='q':
            assert lo<hi;mid=(lo+hi)//2
            stack.extend(((mid+1,hi,vlo,vhi,depth+1,in_tree),(lo,mid,vlo,vhi,depth+1,in_tree)))
        elif code=='p':
            mid=(vlo+vhi)/2
            stack.extend(((lo,hi,mid,vhi,depth+1,in_tree),(lo,hi,vlo,mid,depth+1,in_tree)))
        else:assert type(code) is int;leaves+=1
        if time.monotonic()-printed>15:
            print('replay' if saved else 'refine','leaves',leaves,'Q',lo,hi,flush=True);printed=time.monotonic()
    assert cursor==len(tree) and ledger.validate_tree(nodes,first,last)==leaves
    total=base.encode(F(last-first+1,1<<bits))
    if saved:
        assert saved['leaves']==leaves and saved['range_upper']==total
        base.write_new(output.with_name(output.stem+'_replay.json'),dict(status='FRONTIER_DENSE_512_BIT_REPLAY_PASSED',
            producer_sha256=base.sha(output),message_exponent=m,occupancy_range=[first,last],leaves=leaves))
    else:
        hashes=dict(seed['source_sha256'])
        for p in (Path(__file__),Path(ledger.__file__),source,replay):hashes[p.relative_to(base.ROOT).as_posix()]=base.sha(p)
        base.write_new(output,dict(status='FRONTIER_DENSE_INTERVAL_CERTIFICATE',method='retained_sqrt_density_domination',
            seed_certificate=source.relative_to(base.ROOT).as_posix(),message_exponent=m,occupancy_range=[first,last],
            row_probabilities=seed['row_probabilities'],per_occupancy_upper_power=-bits,range_upper=total,
            tree=nodes,leaves=leaves,source_sha256=hashes,full_distance_proved=False))
    print('Completed retained dense bound',first,last,'bits',bits,'leaves',leaves,'seconds',round(time.monotonic()-started,2),flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True)
    p.add_argument('--source',type=Path);p.add_argument('--bits',type=int,default=100);p.add_argument('--verify',action='store_true')
    a=p.parse_args();run(a.output,a.source,a.bits,a.verify)
