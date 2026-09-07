"""Outward anchor checks for candidate maps using exact fixed-weight regions.

Reuses the frozen activation-aware three-state inequality and directed
recurrence, with explicit map/length parameters. An anchor is not full coverage.
"""
import argparse
from pathlib import Path
import time
from flint import arb, ctx

import inner_candidate_screen as screen
import frontier_sparse as sparse

base=screen.base


def run(name,m,q,tilt,output,verify=False):
    saved=base.read(output) if verify else None
    if saved:
        name,m,q,tilt=(saved[k] for k in ('configuration','message_exponent','occupation','tilt'))
        for path,digest in saved['source_sha256'].items(): assert base.sha(base.ROOT/path)==digest,path
    elif output.exists(): raise FileExistsError('Use a fresh output path')
    t,s,spectrum,kernel=screen.load(name)
    if m<7: raise ValueError('Complete outer rows required')
    rows=1<<(m-7);cutoff=256*rows//10
    if rows<t or rows%t or not 2<=q<=rows: raise ValueError('Invalid epoch or occupancy geometry')
    started=time.monotonic();ctx.prec=512 if verify else 256
    caps=screen.caps_module.caps()
    region=sparse.poly.regions(t,s,spectrum,kernel,(-(arb(tilt)/10).exp()).exp(),q,length=rows)
    print(name,'region ready','Q',q,'precision',ctx.prec,'seconds',round(time.monotonic()-started,2),flush=True)
    witness=saved['witness'] if saved else sparse.choose(region,q,tilt,rows,cutoff,caps)
    result=sparse.evaluate(region,witness,q,q,rows,cutoff,caps)
    assert len(result)==1 and result[0]['occupation']==q
    power=result[0]['upper_power']
    upper=sparse.as_fraction(power)
    print('anchor result',result[0],flush=True)
    if saved:
        assert upper<=base.decode(saved['upper'])
        receipt=dict(status='CANDIDATE_ANCHOR_512_BIT_REPLAY_PASSED',configuration=name,occupation=q,
                     producer_sha256=base.sha(output),full_distance_proved=False)
        base.write_new(output.with_name(output.stem+'_replay.json'),receipt)
    else:
        # The frozen evaluator retains up to 80 dyadic margin bits.
        hashes=screen.source_hashes();hashes.update(sparse.sources())
        hashes[Path(__file__).relative_to(base.ROOT).as_posix()]=base.sha(Path(__file__))
        data=dict(status='CANDIDATE_FIXED_WEIGHT_ANCHOR_OUTWARD',configuration=name,message_exponent=m,
            occupation=q,tilt=tilt,rows=rows,cutoff=cutoff,witness=witness,upper=base.encode(upper),
            retained_margin_bits=-power,source_sha256=hashes,seconds=time.monotonic()-started,
            full_distance_proved=False)
        base.write_new(output,data)
    print(name,'Q',q,'complete','seconds',round(time.monotonic()-started,2),flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--configuration',choices=(*screen.maps.NAMES,*base.CONFIGS))
    p.add_argument('--m',type=int,default=20);p.add_argument('--occupation',type=int,default=8192)
    p.add_argument('--tilt',type=int,default=6);p.add_argument('--output',type=Path,required=True)
    p.add_argument('--verify',action='store_true');a=p.parse_args()
    if not a.verify and a.configuration is None:p.error('--configuration required')
    run(a.configuration,a.m,a.occupation,a.tilt,a.output,a.verify)
