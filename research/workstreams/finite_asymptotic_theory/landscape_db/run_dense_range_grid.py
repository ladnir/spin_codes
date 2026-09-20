"""Sequential whole-grid dense interval union diagnostics."""
import argparse
import csv
import json
import math
from pathlib import Path

import numpy as np

import activation_occupation as general
import dense_occupation_ranges as dense
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
        if receipt['parameter_keys']!=[list(grid.key(r)) for r in rows]: raise ValueError('dense coverage changed')
        return receipt['row_count']
    if (directory/'occupations.csv').exists(): raise ValueError('final dense CSV without receipt requires inspection')
    directory.mkdir(parents=True,exist_ok=True)
    map_path=grid.pilot.ROOT/rows[0]['map_source']
    if any(r['map_tag']!=rows[0]['map_tag'] for r in rows) or rows[0]['map_tag']!='nested-'+grid.pilot.sha(map_path)[:16]:
        raise ValueError('dense map identity mismatch')
    config=json.loads(map_path.read_text());t=config['step_bits'];s=config['state_bits']
    a_counts={w:n for w,n in enumerate(config['a_counts']) if w and n}
    dependencies=dict(sources)
    dependencies[map_path.relative_to(grid.pilot.ROOT).as_posix()]=grid.pilot.sha(map_path)
    fields=list(rows[0])+['occupation_min','occupation_max','coverage_kind','setup_event_id','setup_failure_bits','witness_shift','range_witness_source']
    pairs=[(p,eta) for p in arguments['probabilities'] for eta in arguments['log_fugacities']]
    ps=np.array([x[0] for x in pairs]);etas=np.array([x[1] for x in pairs])
    epochs={};moments={};output=[]
    for ordinal,row in enumerate(rows):
        b,d,e=int(row['block_bits']),int(row['dimension']),int(row['message_exponent'])
        L=(1<<e)//d;cutoff=b*L//10
        if L<arguments['minimum_occupation']: continue
        event=''
        if row['outer_model']=='random-ensemble-average':
            counts=general.random_spectrum_caps(b,d,arguments['setup_failure_bits'])
            payload=dict(block_bits=b,dimension=d,setup_failure_bits=arguments['setup_failure_bits'],
                         counts=counts,method='simultaneous Markov shell caps; one reused full-rank subspace')
            event_path=directory.parent/f'random{b}_caps.json'
            encoded=json.dumps(payload,sort_keys=True,indent=2)+'\n'
            if event_path.exists() and event_path.read_text()!=encoded: raise ValueError('setup event changed')
            if not event_path.exists(): event_path.write_text(encoded)
            event=grid.pilot.sha(event_path)
            dependencies[event_path.relative_to(grid.pilot.ROOT).as_posix()]=event
        else: counts=common.exact_counts(row)
        partition=intervals(L,arguments['minimum_occupation'],arguments['interval_count'])
        left=np.array([a for a,b in partition]);right=np.array([b for a,b in partition])
        width=np.log(right-left+1)
        mass=(1<<d)-1
        mode=(L+1)*mass//(mass+1)
        peaks=np.minimum(right,np.maximum(left,mode))
        best=width+dense.log_choose(L,peaks)+peaks*math.log(mass)
        witnesses=[None]*len(partition)
        gammas={p:dense.density_cost(counts,b,p) for p in arguments['probabilities']}
        gamma=np.array([gammas[p] for p in ps])
        for tilt in arguments['log_surprisals']:
            lam=math.exp(tilt)
            if tilt not in epochs:
                epochs[tilt]=general.epoch_logs(t,s,a_counts,config['kernel_counts'],lam,t)
            if (e,tilt) not in moments:
                moments[e,tilt]=dense.serialized_moments(epochs[tilt],t,b*L,ps,etas)
            J=moments[e,tilt]
            lo=dense.point_logs(left[None,:],b,L,gamma[:,None],etas[:,None],J[:,None],cutoff,lam)
            hi=dense.point_logs(right[None,:],b,L,gamma[:,None],etas[:,None],J[:,None],cutoff,lam)
            values=np.maximum(lo,hi)+width[None,:]
            selected=np.argmin(values,axis=0)
            for interval,index in enumerate(selected):
                value=float(values[index,interval])
                if value<best[interval]:
                    best[interval]=value
                    witnesses[interval]=dict(log_surprisal=tilt,probability=float(ps[index]),log_fugacity=float(etas[index]))
        if not np.isfinite(best).all(): raise ArithmeticError('nonfinite dense interval bound')
        union=float(np.logaddexp.reduce(best));dominant=int(np.argmax(best))
        witness=witnesses[dominant]
        table=dict(parameter_key=list(grid.key(row)),occupation_min=int(left[0]),occupation_max=int(right[-1]),
                   intervals=[dict(lo=int(a),hi=int(b),log_union_upper=float(v),
                                   method='cauchy_endpoint' if w else 'message_count',witness=w)
                              for (a,b),v,w in zip(partition,best,witnesses)],
                   margin_bits=-union/math.log(2),setup_event_id=event,
                   arithmetic='nearest binary64 diagnostic',all_integers_covered=True)
        path=directory/f'row{ordinal}_ranges.json'
        encoded=json.dumps(table,indent=2)+'\n'
        if path.exists() and path.read_text()!=encoded: raise ValueError('dense witness changed during resume')
        if not path.exists():path.write_text(encoded)
        dependencies[path.relative_to(grid.pilot.ROOT).as_posix()]=grid.pilot.sha(path)
        result=dict(row,occupation_min=int(left[0]),occupation_max=int(right[-1]),coverage_kind='occupation_union',
                    margin_bits=table['margin_bits'],dominant_weight='',
                    dominant_log_surprisal=witness['log_surprisal'] if witness else '',
                    dominant_witness_at_grid_edge=int(witness is not None and witness['log_surprisal'] in (min(arguments['log_surprisals']),max(arguments['log_surprisals']))),
                    witness_shift='',setup_event_id=event,setup_failure_bits=arguments['setup_failure_bits'] if event else '',
                    range_witness_source=path.relative_to(grid.pilot.ROOT).as_posix(),
                    notes='Dense interval union; Cauchy coefficient bound and convex endpoint envelopes cover every integer; includes recorded counting fallbacks; binary64 diagnostic.')
        if event:
            result['outer_model']='random-simultaneous-spectrum-caps'
            result['notes']+=' Add the shared spectrum-event failure budget once in the final union.'
        output.append(result)
        print(f't{t} s{s} e{e} {row["series"]} Q{left[0]}..{right[-1]}: {table["margin_bits"]:.3f} bits',flush=True)
    pending=directory/'occupations.pending.csv'
    with pending.open('w',newline='') as handle:
        writer=csv.DictWriter(handle,fieldnames=fields);writer.writeheader();writer.writerows(output)
    receipt=dict(schema='activation-aware-dense-range-grid-v1',status='BINARY64_DIAGNOSTIC',arguments=arguments,
                  parameter_keys=[list(grid.key(r)) for r in rows],row_count=len(output),
                  source_sha256=dependencies,csv_sha256=grid.pilot.sha(pending),
                  limitations=['One whole-spectrum Bernoulli majorant can be loose.',
                               'Every integer in the recorded ranges is covered by a stated interval inequality.',
                               'Counting fallbacks are distinguished from Cauchy witnesses; no outward certificate.'])
    (directory/'manifest.pending.json').write_text(json.dumps(receipt,indent=2)+'\n')
    pending.replace(directory/'occupations.csv');(directory/'manifest.pending.json').replace(directory/'manifest.json')
    return len(output)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output-dir',type=Path,required=True)
    parser.add_argument('--steps',nargs='+',type=int,default=list(grid.STEPS))
    parser.add_argument('--states',nargs='+',type=int)
    parser.add_argument('--message-exponents',nargs='+',type=int,default=list(grid.EXPONENTS))
    parser.add_argument('--blocks',nargs='+',type=int)
    parser.add_argument('--minimum-occupation',type=int,default=17)
    parser.add_argument('--interval-count',type=int,default=64)
    parser.add_argument('--probabilities',nargs='+',type=float,default=[.1,.25,.5,.75,.9])
    parser.add_argument('--log-fugacities',nargs='+',type=float,default=[-12.,-9.,-7.,-5.,-3.,-2.,-1.,0.,1.,2.,3.,5.,7.,9.,12.])
    parser.add_argument('--log-surprisals',nargs='+',type=float,default=[-7.,-5.,-3.,-2.,-1.,0.])
    parser.add_argument('--setup-failure-bits',type=int,default=60)
    args=parser.parse_args()
    if (args.minimum_occupation<1 or args.interval_count<1 or args.setup_failure_bits<0
            or not all(0<p<1 for p in args.probabilities)
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
    for path in (Path(__file__),Path(dense.__file__),Path(general.__file__),Path(coverage.__file__),Path(common.__file__),Path(sparse.__file__)):
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
