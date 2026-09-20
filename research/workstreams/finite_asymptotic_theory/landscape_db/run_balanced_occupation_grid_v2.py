"""Sequential balanced-probability occupation grid with variance spectrum caps.

This producer records each evaluated integer Q separately. A requested
ceiling below the outer row count is always partial coverage.
"""
import argparse
import csv
import json
import math
from pathlib import Path

import numpy as np

import activation_occupation as transfer
import balanced_occupation as balanced
import balanced_occupation_v2 as batched
import random_spectrum_variance as variance
import read_grid_receipts as coverage
import run_complete_q1_grid as grid


def exact_counts(row):
    key,c=next((key,c) for key,c in grid.pilot.CONSTITUENTS.items() if c.name+' exact'==row['series'])
    if c.spectrum_path.exists():
        counts=grid.pilot.load_spectrum(c)
    else:
        payload=json.loads((grid.pilot.SMALL/'spectra'/f'{key}_weight_counts.json').read_text())
        counts={int(w):int(n) for w,n in payload['weight_counts'].items()}
    if sum(counts.values())!=1<<c.dimension or counts.get(0)!=1:
        raise ValueError('incomplete outer spectrum')
    return {w:n for w,n in counts.items() if w and n}


def verify(directory,arguments):
    receipt=json.loads((directory/'manifest.json').read_text())
    if receipt['arguments']!=arguments:
        raise ValueError('occupation batch arguments changed')
    if grid.pilot.sha(directory/'occupations.csv')!=receipt['csv_sha256']:
        raise ValueError('occupation CSV changed')
    for path,digest in receipt['source_sha256'].items():
        if grid.pilot.sha(grid.pilot.ROOT/path)!=digest:
            raise ValueError(f'occupation dependency changed: {path}')
    return receipt


def run_batch(directory,rows,arguments,sources):
    if (directory/'manifest.json').exists():
        receipt=verify(directory,arguments)
        if receipt['parameter_keys']!=[list(grid.key(r)) for r in rows]:
            raise ValueError('balanced batch coverage changed')
        return receipt['row_count']
    if (directory/'occupations.csv').exists():
        raise ValueError('final CSV without receipt requires inspection')
    directory.mkdir(parents=True,exist_ok=True)
    map_path=grid.pilot.ROOT/rows[0]['map_source']
    if any(r['map_tag']!=rows[0]['map_tag'] for r in rows) or rows[0]['map_tag']!='nested-'+grid.pilot.sha(map_path)[:16]:
        raise ValueError('mixed or changed maps')
    config=json.loads(map_path.read_text());t=config['step_bits'];s=config['state_bits']
    a_counts={w:n for w,n in enumerate(config['a_counts']) if w and n}
    dependencies=dict(sources)
    dependencies[map_path.relative_to(grid.pilot.ROOT).as_posix()]=grid.pilot.sha(map_path)
    fields=list(rows[0])+['occupation','setup_event_id','setup_failure_bits','witness_probability_scale']
    geometries={};models={};prepared=[]
    for row in rows:
        b,d,e=int(row['block_bits']),int(row['dimension']),int(row['message_exponent'])
        length=(1<<e)//d;maximum=min(arguments['maximum_occupation'],length)
        key=(row['series'],b,d);event=''
        if key not in models:
            if row['outer_model']=='random-ensemble-average':
                payload=variance.evidence(b,d,arguments['setup_failure_bits']);counts=payload['counts']
                path=directory.parent/f'random{b}_caps.json'
                encoded=json.dumps(payload,sort_keys=True,indent=2)+'\n'
                if path.exists() and path.read_text()!=encoded:raise ValueError('random event changed')
                if not path.exists():path.write_text(encoded)
                event=grid.pilot.sha(path)
                dependencies[path.relative_to(grid.pilot.ROOT).as_posix()]=event
            else:counts=exact_counts(row)
            models[key]=(batched.PreparedModel(counts,b,arguments['probability_scales'],arguments['bands']),event)
        model,event=models[key]
        item=dict(row=row,length=length,maximum=maximum,cutoff=b*length//10,model=model,event=event,
                  best=np.full(maximum,np.inf),witness=np.zeros((maximum,2)))
        prepared.append(item)
        geometries.setdefault(length,[]).append(item)
    for ordinal,tilt in enumerate(arguments['log_surprisals']):
        lam=math.exp(tilt)
        epoch=transfer.epoch_logs(t,s,a_counts,config['kernel_counts'],lam,min(t,arguments['maximum_occupation']))
        regions=batched.region_ladder(epoch,t,geometries,arguments['maximum_occupation'])
        for length,items in geometries.items():
            for item in items:
                values,scales=item['model'].bounds(regions[length],length,item['cutoff'],lam)
                improve=values<item['best']
                item['best'][improve]=values[improve]
                item['witness'][improve,0]=tilt;item['witness'][improve,1]=scales[improve]
        if ordinal%8==7:print(f't{t} s{s}: {ordinal+1}/{len(arguments["log_surprisals"])} tilts',flush=True)
    output=[]
    for item in prepared:
        row=item['row'];event=item['event']
        if not np.isfinite(item['best']).all():raise ArithmeticError('nonfinite balanced bound')
        for q in range(2,item['maximum']+1):
            tilt,scale=map(float,item['witness'][q-1])
            result=dict(row,occupation=q,margin_bits=-float(item['best'][q-1])/math.log(2),dominant_weight='',
                        dominant_log_surprisal=tilt,
                        dominant_witness_at_grid_edge=int(tilt in (min(arguments['log_surprisals']),max(arguments['log_surprisals']))),
                        witness_probability_scale=scale,setup_event_id=event,
                        setup_failure_bits=arguments['setup_failure_bits'] if event else '',
                        notes=f'Activation-aware adaptive shell bound; fixed shell logits scaled by {scale}; batched witnesses and shared region powers; nearest binary64 diagnostic.')
            if event:
                result['outer_model']='random-simultaneous-spectrum-caps'
                result['notes']+=' Add the authenticated shared spectrum-event failure once in the final union.'
            output.append(result)
    pending=directory/'occupations.pending.csv'
    with pending.open('w',newline='') as handle:
        writer=csv.DictWriter(handle,fieldnames=fields);writer.writeheader();writer.writerows(output)
    receipt=dict(schema='activation-aware-balanced-occupation-grid-v2',status='BINARY64_DIAGNOSTIC',
                 arguments=arguments,parameter_keys=[list(grid.key(r)) for r in rows],row_count=len(output),
                 source_sha256=dependencies,csv_sha256=grid.pilot.sha(pending),
                 limitations=['Only explicitly listed integer occupations are covered.',
                              'Conditional random bounds require their shared setup failure charge.',
                              'No outward arithmetic and no extrapolation certificate.'])
    (directory/'manifest.pending.json').write_text(json.dumps(receipt,indent=2)+'\n')
    pending.replace(directory/'occupations.csv');(directory/'manifest.pending.json').replace(directory/'manifest.json')
    return len(output)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output-dir',type=Path,required=True)
    parser.add_argument('--maximum-occupation',type=int,default=64)
    parser.add_argument('--steps',nargs='+',type=int,default=list(grid.STEPS))
    parser.add_argument('--states',nargs='+',type=int)
    parser.add_argument('--message-exponents',nargs='+',type=int,default=list(grid.EXPONENTS))
    parser.add_argument('--blocks',nargs='+',type=int)
    parser.add_argument('--log-surprisals',nargs='+',type=float,default=[j/4 for j in range(-32,5)])
    parser.add_argument('--probability-scales',nargs='+',type=float,default=[0.,.05,.1,.15,.2,.25,.3,.4,.5,.75,1.])
    parser.add_argument('--bands',type=int,default=0,help='0 means singleton shells')
    parser.add_argument('--setup-failure-bits',type=int,default=60)
    args=parser.parse_args()
    if (args.maximum_occupation<2 or args.bands<0
            or not all(math.isfinite(x) and x>=0 for x in args.probability_scales)
            or not all(math.isfinite(x) for x in args.log_surprisals) or args.setup_failure_bits<0
            or not set(args.steps).issubset(grid.STEPS)
            or not set(args.message_exponents).issubset(grid.EXPONENTS)):
        raise ValueError('invalid occupation grid arguments')
    if (grid.OUT/'run.lock').exists():
        raise ValueError('Q1 producer is active: numerical jobs must run sequentially')
    planned,observations,_=coverage.snapshot()
    selected=[r for r in planned if r['native'] and r['step_bits'] in args.steps
              and r['message_exponent'] in args.message_exponents
              and (args.states is None or r['state_bits'] in args.states)
              and (args.blocks is None or r['block_bits'] in args.blocks)]
    if not selected or any(grid.key(r) not in observations for r in selected):
        raise ValueError('selected parameter grid is empty or lacks completed Q1 maps')
    directory=args.output_dir.resolve()
    directory.relative_to(grid.HERE)  # Receipts stay inside this workstream.
    directory.mkdir(parents=True,exist_ok=True)
    arguments={k:v for k,v in vars(args).items() if k!='output_dir'}
    sources={}
    # Authenticate the complete Q1 sources as provenance for all selected maps.
    for name in grid.BASES:
        sources.update(grid.verify(grid.HERE/name)['source_sha256'])
    for path in (Path(__file__),Path(transfer.__file__),Path(balanced.__file__),Path(batched.__file__),Path(variance.__file__),Path(coverage.__file__),Path(grid.__file__)):
        sources[path.relative_to(grid.pilot.ROOT).as_posix()]=grid.pilot.sha(path)
    lock=grid.HERE/'occupation_grid.lock'
    with lock.open('x') as handle: handle.write('Sequential occupation producer is active.\n')
    try:
        for t in args.steps:
            for s in range(t.bit_length(),21):
                rows=[observations[grid.key(r)] for r in selected if (r['step_bits'],r['state_bits'])==(t,s)]
                if rows:
                    count=run_batch(directory/f't{t}_s{s}',rows,arguments,sources)
                    print(f't{t} s{s}: verified {count} occupation rows',flush=True)
    finally:
        lock.unlink()


if __name__=='__main__':
    main()
