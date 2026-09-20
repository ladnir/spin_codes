"""Check sparse occupations and transport the retained dense cover to new inners.

All results use nearest binary64, not outward-certified arithmetic.
Dense box geometry depends on the outer, not on the inner. Every transported
witness is reevaluated with the candidate's exact constituent spectra.
"""
import argparse
import copy
from fractions import Fraction
import json
import math
from pathlib import Path

import numpy as np
import calibrate as calibration

fixed = calibration.fixed
import activation_q2 as pair
import typed_dense_boxes as typed
from verify_fixed_inner_results import check_coverage

L, B, N = 32768, 128, 4194304
DISTANCES = (Fraction(33,200), Fraction(19,100))


def sparse(record, counts, maximum):
    t,s = record['step_bits'],record['state_bits']
    assert L%t == 0 and maximum >= 4
    a = {w:n for w,n in enumerate(record['a_counts']) if w and n}
    kernel = record['kernel_counts']
    cutoffs = [N*d.numerator//d.denominator for d in DISTANCES]
    best = np.full((2,maximum),np.inf)
    for i,row in enumerate(calibration.q1_screen(record,counts)):
        best[i,0] = -row['q1_margin_bits']*math.log(2)
    weights = np.array(sorted(counts),dtype=int)
    log_counts = np.array([math.log(counts[w]) for w in weights])
    pair_best = np.full((2,len(weights),len(weights)),np.inf)
    engine = pair.PairKernel()
    zs = np.arange(-44,-7)/4
    regions = pair.region_logs(t,s,a,np.exp(zs),L)
    for j,z in enumerate(zs):
        moment = engine.coefficients(regions[j],B)[np.ix_(weights,weights)]
        for i,cutoff in enumerate(cutoffs):
            np.minimum(pair_best[i],np.minimum(0.,moment+cutoff*math.exp(z)),out=pair_best[i])
    for i in range(2):
        terms = pair_best[i]+log_counts[:,None]+log_counts[None,:]
        best[i,1] = math.log(math.comb(L,2))+float(np.logaddexp.reduce(terms.ravel()))
    print(record['tag'],'Q1/Q2',(-best[:,:2]/math.log(2)).tolist(),flush=True)
    models = {q:fixed.composition.SparseComposition(counts,B,q) for q in (3,4)}
    components = {q:np.full((2,len(model.indices)),np.inf) for q,model in models.items()}
    grid = np.arange(-112,1,dtype=float)/8
    for j,z in enumerate(grid):
        lam = math.exp(z)
        # Longer tilts only improved the small-occupation part of the original
        # search. Retain its finer small-Q search and its Q<=128 extension.
        active_maximum = min(maximum,64) if z < -7 else maximum
        epoch = fixed.general.epoch_logs(t,s,a,kernel,lam,min(t,active_maximum))
        regions = fixed.general.region_logs(epoch,t,L,active_maximum)
        for i,cutoff in enumerate(cutoffs):
            for shift in (-.5,0.,.5,1.):
                values = fixed.general.occupation_bounds(regions,counts,B,L,cutoff,lam,shift,8)
                np.minimum(best[i,2:active_maximum],values[2:],out=best[i,2:active_maximum])
                for q,model in models.items():
                    np.minimum(components[q][i],model.components(regions,cutoff,lam,shift),out=components[q][i])
        if j%8 == 0:
            print(record['tag'],f'sparse tilts {j+1}/{len(grid)}',flush=True)
    for q,model in models.items():
        for i in range(2):
            best[i,q-1] = min(best[i,q-1],model.aggregate(components[q][i],L))
    return [dict(distance_target=str(d),bad_weight=cutoffs[i],covered_occupations=[1,maximum],
                 occupation_margins_bits=(-best[i]/math.log(2)).tolist(),
                 sparse_union_margin_bits=-float(np.logaddexp.reduce(best[i]))/math.log(2),
                 worst_occupation=int(np.argmax(best[i]))+1)
            for i,d in enumerate(DISTANCES)]


def transport_dense(record, counts, original, minimum):
    t,s = record['step_bits'],record['state_bits']
    a = {w:n for w,n in enumerate(record['a_counts']) if w and n}
    kernel = record['kernel_counts']
    result = copy.deepcopy(original['dense'])
    assert minimum >= result['occupation_min']
    boxes,epochs = [],{}
    for box in result['selected_boxes']:
        upper = np.array(box['upper']); lower = np.array(box['lower'])
        upper[0] = min(upper[0],L-minimum)
        corners = typed.vertices(lower,upper,L)
        if not len(corners): continue
        box['lower'],box['upper'] = corners.min(axis=0).tolist(),corners.max(axis=0).tolist()
        w = box['witness']; w.pop('maximum_vertex',None)
        z = w['log_surprisal']; lam = math.exp(z)
        if z not in epochs:
            epochs[z] = fixed.general.epoch_logs(t,s,a,kernel,lam,t)
        bank = result['probability_banks'][w['probability_bank']]
        p,costs,proposal = np.array(bank['probabilities']),np.array(bank['log_density_costs']),np.array(w['proposal'])
        assert np.all(proposal>0) and abs(sum(proposal)-1)<1e-12
        for g,band in enumerate(result['bands'],1):
            if band == [B]:
                assert p[g] == 1 and costs[g] == 0
            else:
                cost = max(math.log(counts[v])-math.log(math.comb(B,v))-v*math.log(p[g])-(B-v)*math.log1p(-p[g]) for v in band)
                assert abs(cost-costs[g])<1e-10
        moment,_ = typed.moment_and_density(epochs[z],t,N,float(proposal@p))
        value = max(typed.point_logs(corners,L,B,costs,proposal,moment,original['bad_weight'],lam))
        value += typed.lattice_log_count(box['lower'],box['upper'])
        box['own_log_bound'] = float(value)
        boxes.append(box)
    check_coverage(boxes,L,minimum)
    result['selected_boxes'] = boxes
    result['occupation_min'] = minimum
    result['log_union_upper'] = float(np.logaddexp.reduce([b['own_log_bound'] for b in boxes]))
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--tags',nargs='+',default=['t256_s18_nested'])
    parser.add_argument('--maximum',type=int,default=128)
    parser.add_argument('--dense-only',action='store_true')
    args = parser.parse_args()
    source = fixed.HERE/'SMALLER_MARGIN_CLOSED.json'
    retained = json.loads(source.read_text())
    for name,digest in retained['source_sha256'].items():
        fixed.outer.require(fixed.sha(fixed.ROOT/name)==digest,'changed retained dependency: '+name)
    screen_path = calibration.HERE/'Q1_SCREEN.json'
    screen = json.loads(screen_path.read_text())
    for name,digest in screen['source_sha256'].items():
        fixed.outer.require(fixed.sha(fixed.ROOT/name)==digest,'changed screen dependency: '+name)
    counts = {w:n for w,n in enumerate(calibration.smaller_outer.spectrum()) if w and n}
    for tag in args.tags:
        map_path = calibration.HERE/'maps'/f'{tag}.json'
        entry = next(r for r in screen['results'] if r['tag']==tag)
        assert fixed.sha(map_path)==entry['map_sha256']
        record = json.loads(map_path.read_text())
        rows = [] if args.dense_only else sparse(record,counts,args.maximum)
        for i,original in enumerate(retained['results']):
            minimum = max(args.maximum+1,original['dense']['occupation_min'])
            dense = transport_dense(record,counts,original,minimum)
            if args.dense_only: rows.append(dict(distance_target=original['distance_target']))
            row = rows[i]
            assert row['distance_target']==original['distance_target']
            row.update(dense=dense,dense_union_margin_bits=-dense['log_union_upper']/math.log(2),exact_integer_dense_coverage_checked=True)
            if not args.dense_only:
                assert minimum==args.maximum+1
                row['combined_margin_bits'] = -float(np.logaddexp(dense['log_union_upper'],-row['sparse_union_margin_bits']*math.log(2)))/math.log(2)
            print(tag,row['distance_target'],{k:v for k,v in row.items() if k.endswith('margin_bits')},flush=True)
        paths = [Path(__file__),Path(calibration.__file__),map_path,screen_path,source,
                 Path(fixed.general.__file__),Path(fixed.composition.__file__),Path(pair.__file__),
                 fixed.LANDSCAPE/'activation_q2_kernel.cpp',Path(typed.__file__),Path(typed.dense.__file__),
                 fixed.HERE/'verify_fixed_inner_results.py']
        payload = dict(status='DENSE_ONLY_BINARY64' if args.dense_only else 'FULL_OCCUPATION_BINARY64_NOT_OUTWARD_CERTIFICATE',
                       tag=tag,outer=[128,32,32],message_exponent=20,results=rows,
                       source_sha256={p.relative_to(fixed.ROOT).as_posix():fixed.sha(p) for p in paths})
        output = calibration.HERE/f"{tag}_{'dense' if args.dense_only else 'full'}.json"
        output.write_text(json.dumps(payload,indent=2)+'\n',encoding='utf-8',newline='\n')


if __name__=='__main__': main()
