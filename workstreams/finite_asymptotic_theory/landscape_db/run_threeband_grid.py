"""Sequential fixed-three-band occupation grid with exact-region composition boxes.

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
import composition_boxes as boxes
import composition_occupation as sparse
import typed_dense_boxes as typed
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
        if receipt['parameter_keys']!=[list(grid.key(r)) for r in rows]:raise ValueError('composition coverage changed')
        return receipt['row_count']
    if (directory/'occupations.csv').exists():raise ValueError('final CSV without receipt requires inspection')
    directory.mkdir(parents=True,exist_ok=True)
    map_path=grid.pilot.ROOT/rows[0]['map_source']
    if any(r['map_tag']!=rows[0]['map_tag'] for r in rows) or rows[0]['map_tag']!='nested-'+grid.pilot.sha(map_path)[:16]:
        raise ValueError('mixed or changed map')
    config=json.loads(map_path.read_text());t=config['step_bits'];s=config['state_bits']
    ac={w:n for w,n in enumerate(config['a_counts']) if w and n}
    dependencies=dict(sources);dependencies[map_path.relative_to(grid.pilot.ROOT).as_posix()]=grid.pilot.sha(map_path)
    fields=list(rows[0])+['occupation','setup_event_id','setup_failure_bits','witness_low_probability','range_witness_source']
    output=[]
    for ordinal,row in enumerate(rows):
        b,d,e=int(row['block_bits']),int(row['dimension']),int(row['message_exponent'])
        L=(1<<e)//d;maximum=min(L,arguments['maximum_occupation']);H=b*L//10;event=''
        if row['outer_model']=='random-ensemble-average':
            payload=variance.evidence(b,d,arguments['setup_failure_bits']);counts=payload['counts']
            path=directory.parent/f'random{b}_caps.json';encoded=json.dumps(payload,sort_keys=True,indent=2)+'\n'
            if path.exists() and path.read_text()!=encoded:raise ValueError('random event changed')
            if not path.exists():path.write_text(encoded)
            event=grid.pilot.sha(path);dependencies[path.relative_to(grid.pilot.ROOT).as_posix()]=event
        else:counts=exact_counts(row)
        groups=[[w for w in counts if 16*w<3*b],[w for w in counts if 3*b<=16*w<=13*b],[w for w in counts if 16*w>13*b]]
        bands=[g for g in groups if g]
        models=[]
        for p in arguments['low_probabilities']:
            probabilities=[1. if g==[b] else prob for g,prob in zip(groups,[p,.5,arguments['high_probability']]) if g]
            models.append(boxes.CompositionBoxes(counts,b,bands,probabilities,maximum))
        choose=np.array([math.log(math.comb(L,q)) for q in range(1,maximum+1)])
        best=np.full(maximum,np.inf);selected=np.zeros(maximum,dtype=int);witnesses=[]
        for tilt in arguments['log_surprisals']:
            lam=math.exp(tilt);epoch=transfer.epoch_logs(t,s,ac,config['kernel_counts'],lam,min(t,maximum))
            regions=transfer.region_logs(epoch,t,L,maximum)
            for index,model in enumerate(models):
                current=regions;matrices=[]
                for q in range(1,maximum+1):
                    current=model.envelope.apply(current[:-1],current[1:]);matrices.append(current[0])
                values=choose+np.arange(1,maximum+1)*math.log(len(bands))+H*lam+sparse.terminal_logs(np.array(matrices),b)
                label=dict(log_surprisal=tilt,low_probability=arguments['low_probabilities'][index],high_probability=arguments['high_probability'])
                improve=values<best;best[improve]=values[improve];selected[improve]=len(witnesses)
                witnesses.append((model,regions,lam,label))
        for q in range(2,maximum+1):
            root=float(best[q-1]);winner=witnesses[int(selected[q-1])][3];bound=root;path_text='';method='three_band_adaptive'
            if q>=arguments['refine_minimum'] and root>-arguments['target_bits']*math.log(2):
                local=[w for w in witnesses if abs(w[3]['log_surprisal']-winner['log_surprisal'])<=arguments['witness_window']+1e-12]
                result=boxes.search(local,q,L,H,maximum_nodes=arguments['maximum_nodes'],target_bits=arguments['target_bits'])
                result.update(parameter_key=list(grid.key(row)),root_log_upper=root,bands=bands,setup_event_id=event)
                if result['log_union_upper']<bound:bound=result['log_union_upper'];method='composition_box_cover'
                path=directory/f'row{ordinal}_q{q}_boxes.json';encoded=json.dumps(result,indent=2)+'\n'
                if path.exists() and path.read_text()!=encoded:raise ValueError('box witness changed on resume')
                if not path.exists():path.write_text(encoded)
                path_text=path.relative_to(grid.pilot.ROOT).as_posix();dependencies[path_text]=grid.pilot.sha(path)
                print(f't{t} s{s} {row["series"]} Q{q}: {-bound/math.log(2):.3f} bits, {result["nodes_evaluated"]} nodes',flush=True)
            trivial=float(choose[q-1])+q*math.log((1<<d)-1)
            if trivial<bound:bound=trivial;method='message_count'
            if not math.isfinite(bound):raise ArithmeticError('nonfinite composition bound')
            output.append(dict(row,occupation=q,margin_bits=-bound/math.log(2),dominant_weight='',
                dominant_log_surprisal=winner['log_surprisal'],dominant_witness_at_grid_edge=int(winner['log_surprisal'] in (min(arguments['log_surprisals']),max(arguments['log_surprisals']))),
                witness_low_probability=winner['low_probability'],range_witness_source=path_text,
                outer_model='random-simultaneous-spectrum-caps' if event else row['outer_model'],
                setup_event_id=event,setup_failure_bits=arguments['setup_failure_bits'] if event else '',
                notes=f'Exact region coefficients with fixed three-band measures; selected {method}; CSV tilt describes adaptive baseline, refined witnesses are in the authenticated box file; nearest binary64 diagnostic.'))
    pending=directory/'occupations.pending.csv'
    with pending.open('w',newline='') as handle:
        writer=csv.DictWriter(handle,fieldnames=fields);writer.writeheader();writer.writerows(output)
    receipt=dict(schema='activation-aware-three-band-composition-box-grid-v1',status='BINARY64_DIAGNOSTIC',arguments=arguments,
                 parameter_keys=[list(grid.key(r)) for r in rows],row_count=len(output),source_sha256=dependencies,csv_sha256=grid.pilot.sha(pending),
                 limitations=['Only explicitly evaluated occupations are covered.',
                              'A box witness stays fixed within its count box; random setup failure is charged once.',
                              'No outward rounding or extrapolation certificate.'])
    (directory/'manifest.pending.json').write_text(json.dumps(receipt,indent=2)+'\n')
    pending.replace(directory/'occupations.csv');(directory/'manifest.pending.json').replace(directory/'manifest.json')
    return len(output)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output-dir',type=Path,required=True)
    parser.add_argument('--maximum-occupation',type=int,default=256)
    parser.add_argument('--steps',nargs='+',type=int,default=list(grid.STEPS))
    parser.add_argument('--states',nargs='+',type=int)
    parser.add_argument('--message-exponents',nargs='+',type=int,default=list(grid.EXPONENTS))
    parser.add_argument('--blocks',nargs='+',type=int)
    parser.add_argument('--log-surprisals',nargs='+',type=float,default=[j/5 for j in range(-25,6)])
    parser.add_argument('--low-probabilities',nargs='+',type=float,default=[.2,.225,.25,.275,.3,.35])
    parser.add_argument('--high-probability',type=float,default=.8677722630069483)
    parser.add_argument('--maximum-nodes',type=int,default=31)
    parser.add_argument('--refine-minimum',type=int,default=33)
    parser.add_argument('--target-bits',type=float,default=80.)
    parser.add_argument('--witness-window',type=float,default=.8)
    parser.add_argument('--exact-only',action='store_true')
    parser.add_argument('--setup-failure-bits',type=int,default=60)
    args=parser.parse_args()
    if (args.maximum_occupation<2 or args.maximum_nodes<1 or args.refine_minimum<2 or args.witness_window<0
            or not math.isfinite(args.witness_window) or not math.isfinite(args.target_bits) or args.target_bits<0
            or not 0<args.high_probability<1
            or not all(0<x<1 for x in args.low_probabilities)
            or not all(math.isfinite(x) for x in args.log_surprisals) or args.setup_failure_bits<0
            or not set(args.steps).issubset(grid.STEPS)
            or not set(args.message_exponents).issubset(grid.EXPONENTS)):
        raise ValueError('invalid occupation grid arguments')
    if (grid.OUT/'run.lock').exists():
        raise ValueError('Q1 producer is active: numerical jobs must run sequentially')
    planned,observations,_=coverage.snapshot()
    selected=[r for r in planned if r['native'] and r['step_bits'] in args.steps
              and (not args.exact_only or r['outer_model']!='random-ensemble-average')
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
    for path in (Path(__file__),Path(transfer.__file__),Path(balanced.__file__),Path(batched.__file__),Path(boxes.__file__),Path(sparse.__file__),Path(typed.__file__),Path(variance.__file__),Path(coverage.__file__),Path(grid.__file__)):
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
