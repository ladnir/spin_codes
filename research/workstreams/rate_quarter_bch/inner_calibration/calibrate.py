"""Controlled nested RM2Sub calibration for the fixed [128,32] outer.

No encoder is modified. Q1 screens alone are not complete margin claims.
"""
import argparse
from fractions import Fraction
import json
import math
from pathlib import Path
import random
import sys

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parent))
import evaluate_fixed_inner as fixed
import smaller_outer
sys.path.insert(0,str(fixed.ROOT/'workstreams/finite_asymptotic_theory/small_k_replay'))
import generate_rm2sub_calibration_constituent as algebra


def selected(t,s):
    directory = fixed.BRIDGE/('generated/larger_state_inputs_v1' if (t,s) in ((64,20),(128,19)) else 'inputs')
    path = directory/f't{t}_s{s}_selection.json'
    manifest = directory/'manifest.json'
    fixed.outer.require(fixed.sha(path)==json.loads(manifest.read_text())['snapshot_sha256'][path.name],'changed selected map')
    return json.loads(path.read_text()),[path,manifest]


def audit(t,s,masks):
    m = t.bit_length()-1
    monomials = [(i,j) for i in range(m) for j in range(i+1,m)]
    columns = algebra.coordinate_columns(m,masks,monomials)
    rows = algebra.generator_words(columns,s)
    fixed.outer.require(len(masks)==s-m-1 and algebra.rank(rows)==s,'invalid rank')
    fixed.outer.require(len(set(columns))==t and 0 not in columns,'invalid columns')
    fixed.outer.require(not any((x&y).bit_count()%2 for x in rows for y in rows),'BA != 0')
    spectrum = algebra.enumerate_spectrum(rows,s,t)
    kernel = fixed.outer.macwilliams(spectrum,s)
    fixed.outer.require(kernel[4]==algebra.kernel_weight_four_count(columns),'kernel-four check mismatch')
    # Cheap implementation work indicators, not a timing model. The actual
    # compiled emitter and quadratic circuit can share additional work.
    groups = (s+3)//4
    lookups = sum(sum(bool((c>>(4*g))&15) for g in range(groups)) for c in columns)
    zeta_xors = sum(d for d in [1<<i for i in range(m)] for start in range(0,t,2*d)
                    if (start//(2*d)).bit_count()<=2)
    return dict(step_bits=t,state_bits=s,quadratic_masks_hex=[hex(x) for x in masks],
        columns=columns,generator_rows_hex=[hex(x) for x in rows],a_counts=spectrum,kernel_counts=kernel,
        a_distance=next(w for w in range(1,t+1) if spectrum[w]),
        kernel_distance=next(w for w in range(1,t+1) if kernel[w]),kernel_weight_four=kernel[4],
        cost_indicators=dict(epochs_at_k20=(4*(1<<20))//t,
            field_row_words_at_k20=(4*(1<<20))*s//t,
            emission_lookups_per_output=lookups/t,
            field_lookups_per_output=s*groups/t,zeta_xors_per_output=zeta_xors/t))


def build_maps():
    directory = HERE/'maps'
    directory.mkdir(parents=True,exist_ok=True)
    results,sources = [],[]
    for t,base_s in ((64,20),(128,19),(256,14)):
        source,paths = selected(t,base_s); sources+=paths
        masks = [int(x,16) for x in source['selected']['quadratic_masks_hex']]
        m = t.bit_length()-1
        rng = random.Random(20260911+t)
        while len(masks)<20-m-1:
            candidate = rng.randrange(1,1<<(m*(m-1)//2))
            if algebra.rank(masks+[candidate])>len(masks): masks.append(candidate)
        for s in range(14,21):
            tag = f't{t}_s{s}_nested'
            record = dict(tag=tag,origin=f'prefix/extension of selected t{t}_s{base_s}',
                **audit(t,s,masks[:s-m-1]))
            if s==base_s:
                fixed.outer.require(record['columns']==[int(x,16) for x in source['selected']['B_columns_hex']],'baseline identity changed')
                record['existing_selected_map']=True
            path = directory/f'{tag}.json'
            path.write_text(json.dumps(record,indent=2)+'\n',encoding='utf-8',newline='\n')
            results.append((path,record))
    for t,s in ((64,16),(128,15)):
        source,paths = selected(t,s); sources+=paths
        record = dict(tag=f't{t}_s{s}_selected',origin='separately selected existing control',existing_selected_map=True,
            **audit(t,s,[int(x,16) for x in source['selected']['quadratic_masks_hex']]))
        fixed.outer.require(record['columns']==[int(x,16) for x in source['selected']['B_columns_hex']],'control identity mismatch')
        path = directory/f"{record['tag']}.json"
        path.write_text(json.dumps(record,indent=2)+'\n',encoding='utf-8',newline='\n')
        results.append((path,record))
    return results,sources


def q1_screen(record,counts):
    t,s = record['step_bits'],record['state_bits']
    L,N = (1<<20)//32,4*(1<<20)
    a = {w:n for w,n in enumerate(record['a_counts']) if w and n}
    zs = np.arange(-720,1)/40
    lam = np.exp(zs)
    moments = fixed.q1.coefficient_logs(*fixed.q1.region_logs(*fixed.q1.epoch_logs(t,s,a,lam),L//t),128)
    rows = []
    for delta,target in ((Fraction(33,200),40.),(Fraction(19,100),30.)):
        values = np.minimum(0.,moments+(N*delta.numerator//delta.denominator)*lam[:,None])
        indices = np.argmin(values,axis=0)
        best = values[indices,np.arange(129)]
        weights = sorted(counts)
        terms = [math.log(L)+math.log(counts[w])+best[w] for w in weights]
        dominant = weights[int(np.argmax(terms))]
        margin = -float(np.logaddexp.reduce(terms))/math.log(2)
        rows.append(dict(distance=str(delta),target_bits=target,q1_margin_bits=margin,
            q1_screen_pass=margin>=target,dominant_weight=dominant,log_surprisal=float(zs[indices[dominant]]),
            witness_at_edge=bool(indices[dominant] in (0,len(zs)-1))))
    return rows


def main():
    argparse.ArgumentParser(description=__doc__).parse_args()
    smaller_outer.construction()
    counts = {w:n for w,n in enumerate(smaller_outer.spectrum()) if w and n}
    maps,sources = build_maps()
    results = []
    for path,record in maps:
        rows = q1_screen(record,counts)
        results.append(dict(tag=record['tag'],step_bits=record['step_bits'],state_bits=record['state_bits'],
            map_source=path.relative_to(fixed.ROOT).as_posix(),map_sha256=fixed.sha(path),
            a_distance=record['a_distance'],kernel_distance=record['kernel_distance'],kernel_weight_four=record['kernel_weight_four'],
            cost_indicators=record['cost_indicators'],results=rows))
        print(record['tag'],[round(r['q1_margin_bits'],4) for r in rows],f"kernel d={record['kernel_distance']}",flush=True)
    sources += [Path(__file__),Path(algebra.__file__),Path(fixed.__file__),Path(fixed.q1.__file__),
        Path(fixed.outer.__file__),Path(smaller_outer.__file__),fixed.HERE/'sources/EBCH128_29.wd',fixed.HERE/'sources/EBCH128_36.wd']
    payload = dict(status='Q1_ONLY_BINARY64_CALIBRATION',outer=[128,32,32],message_exponent=20,
        output_bits=4*(1<<20),state_grid=list(range(14,21)),step_grid=[64,128,256],extension_seed_base=20260911,
        results=results,source_sha256={p.relative_to(fixed.ROOT).as_posix():fixed.sha(p) for p in sources})
    (HERE/'Q1_SCREEN.json').write_text(json.dumps(payload,indent=2)+'\n',encoding='utf-8',newline='\n')


if __name__=='__main__': main()
