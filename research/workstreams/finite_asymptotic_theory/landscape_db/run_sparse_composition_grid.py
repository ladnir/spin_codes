"""Sparse occupation grid with composition-specific witness selection."""
import argparse
import csv
import json
import math
from pathlib import Path

import numpy as np

import activation_occupation as general
import composition_occupation as sparse
import read_grid_receipts as coverage
import run_complete_q1_grid as grid
import run_occupation_grid as common


def run_batch(directory,rows,arguments,sources):
    if (directory/'manifest.json').exists():
        receipt=common.verify(directory,arguments)
        if receipt['parameter_keys']!=[list(grid.key(r)) for r in rows]:
            raise ValueError('sparse batch coverage changed')
        return receipt['row_count']
    if (directory/'occupations.csv').exists():
        raise ValueError('final sparse CSV without receipt requires inspection')
    directory.mkdir(parents=True,exist_ok=True)
    t,s=int(rows[0]['step_bits']),int(rows[0]['state_bits'])
    map_path=grid.pilot.ROOT/rows[0]['map_source']
    if any(r['map_tag']!=rows[0]['map_tag'] for r in rows) or rows[0]['map_tag']!='nested-'+grid.pilot.sha(map_path)[:16]:
        raise ValueError('sparse map identity mismatch')
    config=json.loads(map_path.read_text())
    a_counts={w:n for w,n in enumerate(config['a_counts']) if w and n}
    dependencies=dict(sources)
    dependencies[map_path.relative_to(grid.pilot.ROOT).as_posix()]=grid.pilot.sha(map_path)
    fields=list(rows[0])+['occupation','setup_event_id','setup_failure_bits','witness_shift','composition_witness_source']
    output=[];epoch_cache={};region_cache={}
    for ordinal,row in enumerate(rows):
        b,d,e=int(row['block_bits']),int(row['dimension']),int(row['message_exponent'])
        L=(1<<e)//d; cutoff=b*L//10
        event=''
        if row['outer_model']=='random-ensemble-average':
            counts=general.random_spectrum_caps(b,d,arguments['setup_failure_bits'])
            payload=dict(block_bits=b,dimension=d,setup_failure_bits=arguments['setup_failure_bits'],
                         counts=counts,method='simultaneous Markov shell caps; one reused full-rank subspace')
            encoded=json.dumps(payload,sort_keys=True,indent=2)+'\n'
            event_path=directory.parent/f'random{b}_caps.json'
            if event_path.exists() and event_path.read_text()!=encoded: raise ValueError('setup event changed')
            if not event_path.exists(): event_path.write_text(encoded)
            event=grid.pilot.sha(event_path)
            dependencies[event_path.relative_to(grid.pilot.ROOT).as_posix()]=event
        else:
            counts=common.exact_counts(row)
        for q in arguments['occupations']:
            if q>L: continue
            engine=sparse.SparseComposition(counts,b,q,arguments['tail_bands'],arguments['singleton_prefix'])
            best=np.full(len(engine.indices),np.inf)
            witness=np.zeros((len(best),2))
            for tilt in arguments['log_surprisals']:
                lam=math.exp(tilt)
                if (q,tilt) not in epoch_cache:
                    epoch_cache[q,tilt]=general.epoch_logs(t,s,a_counts,config['kernel_counts'],lam,min(t,q))
                if (L,q,tilt) not in region_cache:
                    region_cache[L,q,tilt]=general.region_logs(epoch_cache[q,tilt],t,L,q)
                region=region_cache[L,q,tilt]
                for shift in arguments['shifts']:
                    values=engine.components(region,cutoff,lam,shift)
                    improved=values<best
                    best[improved]=values[improved]
                    witness[improved,0]=tilt;witness[improved,1]=shift
            margin=-engine.aggregate(best,L)/math.log(2)
            if not math.isfinite(margin): raise ArithmeticError('nonfinite sparse bound')
            dominant=int(np.argmax(best+engine.log_mult))
            tilt,shift=map(float,witness[dominant])
            selected={}
            for i,pair in enumerate(witness):
                selected.setdefault(tuple(map(float,pair)),[]).append(i)
            evidence=dict(parameter_key=list(grid.key(row)),occupation=q,bands=engine.bands,
                          composition_order='itertools.combinations_with_replacement(range(len(bands)),occupation)',
                          composition_count=len(best),
                          witnesses=[dict(log_surprisal=a,shift=b,composition_indices=ids) for (a,b),ids in sorted(selected.items())],
                          dominant_composition=list(map(int,engine.indices[dominant])),margin_bits=margin,
                          setup_event_id=event,arithmetic='nearest binary64 diagnostic')
            witness_path=directory/f'row{ordinal}_q{q}_witness.json'
            encoded=json.dumps(evidence,indent=2)+'\n'
            if witness_path.exists() and witness_path.read_text()!=encoded:
                raise ValueError('sparse witness changed during resume')
            if not witness_path.exists(): witness_path.write_text(encoded)
            dependencies[witness_path.relative_to(grid.pilot.ROOT).as_posix()]=grid.pilot.sha(witness_path)
            result=dict(row,occupation=q,margin_bits=margin,dominant_weight='',
                        dominant_log_surprisal=tilt,
                        dominant_witness_at_grid_edge=int(tilt in (min(arguments['log_surprisals']),max(arguments['log_surprisals']))),
                        setup_event_id=event,setup_failure_bits=arguments['setup_failure_bits'] if event else '',
                        witness_shift=shift,composition_witness_source=witness_path.relative_to(grid.pilot.ROOT).as_posix(),
                        notes='Sparse band compositions retained across all regions; pointwise witness selection per composition; binary64 diagnostic.')
            if event:
                result['outer_model']='random-simultaneous-spectrum-caps'
                result['notes']+=' Conditional on one simultaneous spectrum event; add its failure budget once in the union.'
            output.append(result)
            print(f't{t} s{s} e{e} {row["series"]} Q{q}: {margin:.6f} bits, {len(best)} compositions',flush=True)
    pending=directory/'occupations.pending.csv'
    with pending.open('w',newline='') as handle:
        writer=csv.DictWriter(handle,fieldnames=fields)
        writer.writeheader();writer.writerows(output)
    receipt=dict(schema='activation-aware-sparse-composition-grid-v1',status='BINARY64_DIAGNOSTIC',
                  arguments=arguments,parameter_keys=[list(grid.key(r)) for r in rows],row_count=len(output),
                  source_sha256=dependencies,csv_sha256=grid.pilot.sha(pending),
                  limitations=['Only explicitly listed integer occupations are evaluated.',
                               'Auxiliary Bernoulli parameters pay the exact band density costs.',
                               'Conditional random bounds require their shared setup-failure term.'])
    (directory/'manifest.pending.json').write_text(json.dumps(receipt,indent=2)+'\n')
    pending.replace(directory/'occupations.csv')
    (directory/'manifest.pending.json').replace(directory/'manifest.json')
    return len(output)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output-dir',type=Path,required=True)
    parser.add_argument('--steps',nargs='+',type=int,default=list(grid.STEPS))
    parser.add_argument('--states',nargs='+',type=int)
    parser.add_argument('--message-exponents',nargs='+',type=int,default=list(grid.EXPONENTS))
    parser.add_argument('--blocks',nargs='+',type=int)
    parser.add_argument('--occupations',nargs='+',type=int,default=[3,4])
    parser.add_argument('--tail-bands',type=int,default=8)
    parser.add_argument('--singleton-prefix',type=int,default=4)
    parser.add_argument('--log-surprisals',nargs='+',type=float,default=[-8.,-7.,-6.,-5.5,-5.,-4.5,-4.,-3.,-2.,-1.,0.])
    parser.add_argument('--shifts',nargs='+',type=float,default=[0.,.5,1.,2.,4.])
    parser.add_argument('--setup-failure-bits',type=int,default=60)
    args=parser.parse_args()
    if (not args.occupations or min(args.occupations)<2 or max(args.occupations)>8
            or not set(args.steps).issubset(grid.STEPS)
            or not set(args.message_exponents).issubset(grid.EXPONENTS)):
        raise ValueError('invalid sparse grid arguments')
    if (grid.OUT/'run.lock').exists(): raise ValueError('Q1 numerical producer is active')
    planned,observations,_=coverage.snapshot()
    selected=[r for r in planned if r['native'] and r['step_bits'] in args.steps
              and r['message_exponent'] in args.message_exponents
              and (args.states is None or r['state_bits'] in args.states)
              and (args.blocks is None or r['block_bits'] in args.blocks)]
    if not selected or any(grid.key(r) not in observations for r in selected):
        raise ValueError('empty grid or missing completed Q1 maps')
    directory=args.output_dir.resolve();directory.relative_to(grid.HERE)
    directory.mkdir(parents=True,exist_ok=True)
    arguments={k:v for k,v in vars(args).items() if k!='output_dir'}
    _,_,sources,_,_=grid.inventory()
    for path in (Path(__file__),Path(sparse.__file__),Path(general.__file__),Path(coverage.__file__),Path(common.__file__)):
        sources[path.relative_to(grid.pilot.ROOT).as_posix()]=grid.pilot.sha(path)
    lock=grid.HERE/'occupation_grid.lock'
    with lock.open('x') as handle:handle.write('Sequential sparse occupation producer is active.\n')
    try:
        for t in args.steps:
            for s in range(t.bit_length(),21):
                rows=[observations[grid.key(r)] for r in selected if (r['step_bits'],r['state_bits'])==(t,s)]
                if rows:run_batch(directory/f't{t}_s{s}',rows,arguments,sources)
    finally:
        lock.unlink()


if __name__=='__main__':main()
