"""Bounded K16 parameter checks; partial results are never full certificates."""
import argparse
from collections import Counter
from fractions import Fraction as F
import json
import math
from pathlib import Path
import time
from flint import arb,ctx
import model


def run(args):
    # Import the historical screen explicitly: this script also has that name.
    old_screen = model.study.core.grid.ladder.screen
    output = args.output.resolve()
    saved = json.loads(output.read_text()) if args.verify else None
    if saved:
        model.base.authenticate(saved)
        assert saved['parameters'] == [args.t,args.s,16] and saved['mode'] == args.mode
    elif output.exists():
        raise FileExistsError(output)
    ctx.prec = 512 if saved else 256
    e = model.Engine(args.t,args.s)
    result = dict(parameters=[args.t,args.s,16],mode=args.mode,instance=e.identity(),
        status='OUTWARD_K16_PARTIAL_SCREEN',full_distance_proved=False,
        feedback_weights=dict(Counter(v.bit_count() for v in e.columns)),
        kernel_low=e.kernel[:6],precision_bits=ctx.prec)
    start = time.monotonic()
    if args.mode == 'q1':
        result['q1'] = old_screen.q1(e,saved['q1'] if saved else None)
        print('Q1',args.t,args.s,result['q1']['margin_bits'],flush=True)
    elif args.mode == 'sparse':
        result['sparse'] = []
        for i,q in enumerate(args.q):
            row = old_screen.sparse(e,q,saved['sparse'][i] if saved else None)
            result['sparse'].append(row)
            print('sparse',args.t,args.s,q,row['margin_bits'],flush=True)
    else:
        c = model.Checker(args.t,args.s)
        result['points'] = []
        for i,q in enumerate(args.q):
            v = F(args.coordinate)
            witness = saved['points'][i]['witness'] if saved else c.witness(q,q,v,v)
            value = c.bound(q,q,v,v,witness)
            bound = model.base.pack(value.exp())
            if saved:
                assert value.exp() <= model.base.unpack(saved['points'][i]['upper'])
                bound = saved['points'][i]['upper']
            row = dict(q=q,coordinate=str(v),witness=witness,upper=bound,
                margin_bits=float(-model.base.unpack(bound).log()/arb(2).log()))
            result['points'].append(row)
            print('dense point',args.t,args.s,q,str(v),row['margin_bits'],flush=True)
    result['elapsed_seconds'] = time.monotonic()-start
    result['source_sha256'] = model.base.sources()
    if saved:
        result = dict(status='K16_PARTIAL_SCREEN_512_BIT_REPLAY_PASSED',
            producer_sha256=model.base.base.sha(output),full_distance_proved=False)
        output = output.with_name(output.stem+'_replay.json')
    output.parent.mkdir(parents=True,exist_ok=True)
    model.base.base.write_new(output,result)


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--t',type=int,required=True)
    p.add_argument('--s',type=int,required=True)
    p.add_argument('--mode',choices=['q1','sparse','dense'],required=True)
    p.add_argument('--q',type=int,nargs='+',default=[218])
    p.add_argument('--coordinate',default='21/128')
    p.add_argument('--output',type=Path,required=True)
    p.add_argument('--verify',action='store_true')
    run(p.parse_args())
