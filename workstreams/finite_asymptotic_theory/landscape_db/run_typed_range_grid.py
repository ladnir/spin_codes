"""Sequential whole-grid typed dense union diagnostics."""
import argparse
import csv
import json
import math
from pathlib import Path

import numpy as np

import activation_occupation as general
import dense_occupation_ranges as dense
import batched_typed_ranges as batched
import typed_dense_boxes as typed
import random_spectrum_variance as variance
import read_grid_receipts as coverage
import run_complete_q1_grid as grid
import run_occupation_grid as common
import composition_occupation as sparse


def intervals(length,minimum,count):
    if not 1<=minimum<=length or count<1: raise ValueError('invalid dense partition')
    # Integer arithmetic makes the partition deterministic at large L.
    edges=sorted({minimum+(length+1-minimum)*j//count for j in range(count+1)})
    return [(a,b-1) for a,b in zip(edges,edges[1:])]


def run_batch(directory,rows,arguments,sources):
    if (directory/'manifest.json').exists():
        receipt=common.verify(directory,arguments)
        if receipt['parameter_keys']!=[list(grid.key(r)) for r in rows]:raise ValueError('typed coverage changed')
        return receipt['row_count']
    if (directory/'occupations.csv').exists():raise ValueError('final CSV without receipt requires inspection')
    directory.mkdir(parents=True,exist_ok=True)
    map_path=grid.pilot.ROOT/rows[0]['map_source']
    if any(r['map_tag']!=rows[0]['map_tag'] for r in rows) or rows[0]['map_tag']!='nested-'+grid.pilot.sha(map_path)[:16]:
        raise ValueError('mixed or changed maps')
    config=json.loads(map_path.read_text());t=config['step_bits'];s=config['state_bits']
    ac={w:n for w,n in enumerate(config['a_counts']) if w and n}
    dependencies=dict(sources)
    dependencies[map_path.relative_to(grid.pilot.ROOT).as_posix()]=grid.pilot.sha(map_path)
    contexts={};prepared=[];models={}
    for row in rows:
        b,d,e=int(row['block_bits']),int(row['dimension']),int(row['message_exponent'])
        L=(1<<e)//d
        if L<arguments['minimum_occupation']:raise ValueError('dense range exceeds outer row count')
        model=(row['series'],b,d)
        if model not in models:
            event=''
            if row['outer_model']=='random-ensemble-average':
                payload=variance.evidence(b,d,arguments['setup_failure_bits']);counts=payload['counts']
                path=directory.parent/f'random{b}_caps.json';encoded=json.dumps(payload,sort_keys=True,indent=2)+'\n'
                if path.exists() and path.read_text()!=encoded:raise ValueError('setup event changed')
                if not path.exists():path.write_text(encoded)
                event=grid.pilot.sha(path);dependencies[path.relative_to(grid.pilot.ROOT).as_posix()]=event
            else:counts=common.exact_counts(row)
            models[model]=(counts,event)
        counts,event=models[model]
        key=(b,L,tuple(sorted(counts.items())))
        if key not in contexts:
            contexts[key]=batched.TypedPartition(counts,b,L,minimum=arguments['minimum_occupation'],
                depth=arguments['partition_depth'],scales=arguments['probability_scales'],
                eligibility_shifts=arguments['eligibility_shifts'])
        prepared.append((row,event,contexts[key]))
    for ordinal,tilt in enumerate(arguments['log_surprisals']):
        lam=math.exp(tilt);epoch=general.epoch_logs(t,s,ac,config['kernel_counts'],lam,t)
        for context in contexts.values():context.observe(epoch,t,lam,tilt)
        print(f't{t} s{s}: typed tilt {ordinal+1}/{len(arguments["log_surprisals"])}',flush=True)
    fields=list(rows[0])+['occupation_min','occupation_max','coverage_kind','setup_event_id','setup_failure_bits','witness_shift','range_witness_source']
    output=[];fallbacks=0
    for ordinal,(row,event,context) in enumerate(prepared):
        table=context.result();raw=table['log_union_upper']
        trivial=(1<<int(row['message_exponent']))*math.log(2)
        bound=min(raw,trivial);fallback=trivial<raw;fallbacks+=int(fallback)
        table.update(parameter_key=list(grid.key(row)),setup_event_id=event,
                     selected_method='all_message_count' if fallback else 'typed_coefficient_cover',
                     typed_margin_bits=-raw/math.log(2),margin_bits=-bound/math.log(2))
        path=directory/f'row{ordinal}_types.json';encoded=json.dumps(table,indent=2)+'\n'
        if path.exists() and path.read_text()!=encoded:raise ValueError('typed witness changed on resume')
        if not path.exists():path.write_text(encoded)
        dependencies[path.relative_to(grid.pilot.ROOT).as_posix()]=grid.pilot.sha(path)
        dominant=max(table['boxes'],key=lambda x:x['log_union_upper'])
        result=dict(row,occupation_min=table['occupation_min'],occupation_max=table['occupation_max'],
                    coverage_kind='occupation_union',margin_bits=table['margin_bits'],dominant_weight='',
                    dominant_log_surprisal=dominant['log_surprisal'] if not fallback else '',
                    dominant_witness_at_grid_edge=int(not fallback and dominant['log_surprisal'] in (min(arguments['log_surprisals']),max(arguments['log_surprisals']))),
                    setup_event_id=event,setup_failure_bits=arguments['setup_failure_bits'] if event else '',witness_shift='',
                    range_witness_source=path.relative_to(grid.pilot.ROOT).as_posix(),
                    notes=f'Whole dense range evaluated with a deterministic typed simplex-box cover; selected {table["selected_method"]}; raw typed margin {table["typed_margin_bits"]:.12g} bits; binary64 diagnostic.')
        if event:
            result['outer_model']='random-simultaneous-spectrum-caps'
            result['notes']+=' Add the authenticated shared setup failure once in the final union.'
        output.append(result)
    pending=directory/'occupations.pending.csv'
    with pending.open('w',newline='') as handle:
        writer=csv.DictWriter(handle,fieldnames=fields);writer.writeheader();writer.writerows(output)
    receipt=dict(schema='activation-aware-batched-typed-range-v1',status='BINARY64_DIAGNOSTIC',arguments=arguments,
                 parameter_keys=[list(grid.key(r)) for r in rows],row_count=len(output),message_count_fallbacks=fallbacks,
                 source_sha256=dependencies,csv_sha256=grid.pilot.sha(pending),
                 limitations=['Full evaluated coverage does not establish a useful positive bound.',
                              'Finite type-box widths and coefficient conditioning can leave substantial slack.',
                              'Every selected counting fallback is explicit; no outward certificate.'])
    (directory/'manifest.pending.json').write_text(json.dumps(receipt,indent=2)+'\n')
    pending.replace(directory/'occupations.csv');(directory/'manifest.pending.json').replace(directory/'manifest.json')
    print(f't{t} s{s}: {len(output)} typed range rows, {fallbacks} counting fallbacks',flush=True)
    return len(output)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output-dir',type=Path,required=True)
    parser.add_argument('--steps',nargs='+',type=int,default=list(grid.STEPS))
    parser.add_argument('--states',nargs='+',type=int)
    parser.add_argument('--message-exponents',nargs='+',type=int,default=list(grid.EXPONENTS))
    parser.add_argument('--blocks',nargs='+',type=int)
    parser.add_argument('--minimum-occupation',type=int,default=65)
    parser.add_argument('--partition-depth',type=int,default=7)
    parser.add_argument('--probability-scales',nargs='+',type=float,default=[.5,.75,1.])
    parser.add_argument('--eligibility-shifts',nargs='+',type=float,default=[-.5,0.,.5,1.])
    parser.add_argument('--log-surprisals',nargs='+',type=float,default=[-8.,-7.,-6.,-5.,-4.,-3.5,-3.,-2.5,-2.,-1.5,-1.,-.5,0.,.25,.5,.75,1.])
    parser.add_argument('--setup-failure-bits',type=int,default=60)
    args=parser.parse_args()
    if (args.minimum_occupation<1 or not 0<=args.partition_depth<=16 or args.setup_failure_bits<0
            or not all(math.isfinite(x) and x>=0 for x in args.probability_scales)
            or not all(math.isfinite(x) and abs(x)<=20 for x in args.eligibility_shifts)
            or not all(math.isfinite(x) for x in args.log_surprisals)
            or not set(args.steps).issubset(grid.STEPS) or not set(args.message_exponents).issubset(grid.EXPONENTS)):
        raise ValueError('invalid dense grid arguments')
    if (grid.OUT/'run.lock').exists(): raise ValueError('Q1 numerical producer is active')
    planned,observations,_=coverage.snapshot()
    selected=[r for r in planned if r['native'] and r['step_bits'] in args.steps
              and r['message_exponent'] in args.message_exponents
              and (args.states is None or r['state_bits'] in args.states)
              and (args.blocks is None or r['block_bits'] in args.blocks)]
    if not selected or any(grid.key(r) not in observations for r in selected): raise ValueError('empty grid or missing Q1 maps')
    directory=args.output_dir.resolve();directory.relative_to(grid.HERE);directory.mkdir(parents=True,exist_ok=True)
    arguments={k:v for k,v in vars(args).items() if k!='output_dir'}
    _,_,sources,_,_=grid.inventory()
    for path in (Path(__file__),Path(dense.__file__),Path(batched.__file__),Path(typed.__file__),Path(variance.__file__),Path(general.__file__),Path(coverage.__file__),Path(common.__file__),Path(sparse.__file__)):
        sources[path.relative_to(grid.pilot.ROOT).as_posix()]=grid.pilot.sha(path)
    lock=grid.HERE/'occupation_grid.lock'
    with lock.open('x') as handle:handle.write('Sequential dense occupation producer is active.\n')
    try:
        for t in args.steps:
            for s in range(t.bit_length(),21):
                rows=[observations[grid.key(r)] for r in selected if (r['step_bits'],r['state_bits'])==(t,s)]
                if rows:run_batch(directory/f't{t}_s{s}',rows,arguments,sources)
    finally:lock.unlink()


if __name__=='__main__':main()
