"""Compare Q2 at s=12 and s=18 across three epoch lengths, with exact outers."""
import csv
import json
import math
from pathlib import Path

import numpy as np

import activation_q2 as q2
import run_activation_pilot as pilot
from refine_activation_grid import counts_for

HERE=Path(__file__).resolve().parent
OUT=HERE/'activation_pilot_q2_v1'


def main():
    if (OUT/'q2.csv').exists():
        raise ValueError('write-once output already exists')
    parents=[HERE/'activation_pilot_v1',HERE/'activation_pilot_small_state_v1']
    sources={}
    for parent in parents:
        manifest=json.loads((parent/'manifest.json').read_text())
        for path,digest in manifest['source_sha256'].items():
            if pilot.sha(pilot.ROOT/path)!=digest:
                raise ValueError(f'changed dependency: {path}')
        if pilot.sha(parent/'q1.csv')!=manifest['csv_sha256']:
            raise ValueError('changed parent CSV')
        sources.update(manifest['source_sha256'])
        for path in (parent/'manifest.json',parent/'q1.csv'):
            sources[path.relative_to(pilot.ROOT).as_posix()]=pilot.sha(path)
    keys=['ebch8','ebch32','xbch64','ebch128','rm37','rm49']
    families=[]
    for key in keys:
        c=pilot.CONSTITUENTS[key]
        row=dict(series=c.name+' exact',outer_model='fixed',block_bits=c.block_bits,dimension=c.dimension)
        families.append((row,counts_for(row)))
    kernel=q2.PairKernel()
    tilts=np.array([-8.,-7.,-6.,-5.5,-5.,-4.5,-4.,-3.,-2.,-1.,0.])
    lam=np.exp(tilts)
    fields=['series','outer_model','block_bits','dimension','message_exponent','message_bits','output_bits',
            'outer_rows','bad_weight','occupation','step_bits','state_bits','epochs_per_region','margin_bits',
            'dominant_weight','dominant_second_weight','dominant_log_surprisal','dominant_witness_at_grid_edge',
            'a_minimum_distance','kernel_minimum_distance','kernel_weight_four','map_tag','map_source','notes']
    OUT.mkdir(parents=True,exist_ok=True)
    completed=0
    with (OUT/'q2.csv').open('w',newline='') as handle:
        writer=csv.DictWriter(handle,fieldnames=fields)
        writer.writeheader()
        for t in (64,128,256):
            for s in (12,18):
                parent=parents[1] if s==12 else parents[0]
                path=parent/'maps'/f't{t}_s{s}.json'
                config=json.loads(path.read_text())
                generators=[int(g,16) for g in config['generator_words_hex']]
                columns=[sum(((g>>j)&1)<<i for i,g in enumerate(generators)) for j in range(t)]
                if 0 in columns or len(set(columns))!=t:
                    raise ValueError('Q2 requires distinct nonzero B columns')
                spectrum={w:n for w,n in enumerate(config['a_counts']) if w and n}
                for template,counts in families:
                    b,k=template['block_bits'],template['dimension']
                    e=16; L=(1<<e)//k; cutoff=(b*L)//10
                    if (1<<e)%k or L%t:
                        raise ValueError('non-native length')
                    region=q2.region_logs(t,s,spectrum,lam,L)
                    best=np.full((b+1,b+1),np.inf)
                    witness=np.zeros((b+1,b+1),dtype=int)
                    print(f'Q2 t{t} s{s} {template["series"]}',flush=True)
                    for i in range(len(tilts)):
                        value=np.minimum(0.,kernel.coefficients(region[i],b)+cutoff*lam[i])
                        improved=value<best
                        witness[improved]=i
                        np.minimum(best,value,out=best)
                    # Ordered local messages in two chosen row positions. The
                    # same realized spectrum is multiplied, never its expectation.
                    weights=np.array(sorted(counts),dtype=int)
                    logs=np.array([counts[w] for w in weights])
                    terms=best[np.ix_(weights,weights)]+logs[:,None]+logs[None,:]+math.log(math.comb(L,2))
                    margin=-float(np.logaddexp.reduce(terms.ravel()))/math.log(2)
                    a_index,b_index=np.unravel_index(np.argmax(terms),terms.shape)
                    a,bw=int(weights[a_index]),int(weights[b_index])
                    if not math.isfinite(margin):
                        raise ArithmeticError('nonfinite Q2 margin')
                    index=witness[a,bw]
                    row=dict(template,message_exponent=e,message_bits=1<<e,output_bits=b*L,outer_rows=L,
                             bad_weight=cutoff,occupation=2,step_bits=t,state_bits=s,epochs_per_region=L//t,
                             margin_bits=margin,dominant_weight=a,dominant_second_weight=bw,
                             dominant_log_surprisal=tilts[index],dominant_witness_at_grid_edge=int(index in (0,len(tilts)-1)),
                             a_minimum_distance=min(spectrum),kernel_minimum_distance=next(w for w,n in enumerate(config['kernel_counts']) if w and n),
                             kernel_weight_four=config['kernel_counts'][4],map_tag='nested-'+pilot.sha(path)[:16],
                             map_source=path.relative_to(pilot.ROOT).as_posix(),
                             notes='Activation-aware exact-support Q2; pointwise witness minimum before the realized-spectrum pair sum; binary64 only.')
                    writer.writerow(row); handle.flush(); completed+=1
                    print(f'Q2 margin {margin:.6f}, weights {a},{bw}, tilt {tilts[index]}, rows {completed}',flush=True)
    for path in (Path(__file__),HERE/'activation_q2.py',HERE/'activation_q2_kernel.cpp',HERE/'build_activation_q2.ps1'):
        sources[path.relative_to(pilot.ROOT).as_posix()]=pilot.sha(path)
    manifest=dict(schema='activation-aware-q2-pilot-v1',status='BINARY64_DIAGNOSTIC',occupation=2,
                  row_count=completed,source_sha256=sources,csv_sha256=pilot.sha(OUT/'q2.csv'),
                  log_surprisal_grid=list(tilts),native_binary_sha256=pilot.sha(HERE/'native_build/activation_q2_kernel.dll'),
                  compiler_flags='/O2 /EHsc /std:c++17 /fp:precise /LD; link /Brepro',
                  limitations=['Q2 only, exact fixed spectra; no random-outer Q2 reference.',
                               'All arithmetic nearest binary64; no full-distance certificate.',
                               'Kernel spectrum at weights >=3 is not exercised by Q2.',
                               'Coarse common witness grid is not a proved optimizer.'])
    (OUT/'manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')


if __name__=='__main__':
    main()
