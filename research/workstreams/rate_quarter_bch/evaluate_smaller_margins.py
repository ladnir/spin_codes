"""Test the proposed 16.5%/40-bit and 19%/30-bit smaller-outer points.

The fixed inner is unchanged. Exact outer spectra and exact pair supports
feed activation-aware transfer bounds, evaluated in nearest binary64.
"""
import argparse
from fractions import Fraction
import json
import math
from pathlib import Path

import numpy as np
import evaluate_fixed_inner as fixed
import smaller_outer
import activation_q2 as pair
import typed_dense_boxes as typed
from verify_fixed_inner_results import check_coverage


def sparse_results(counts, a, kernel, exponent, distances, maximum):
    B, L = 128, (1 << exponent)//32
    fixed.outer.require((1 << exponent)%32 == 0 and L%128 == 0, 'incomplete geometry')
    N = B*L
    cutoffs = [N*d.numerator//d.denominator for d in distances]
    zs = np.arange(-720,1)/40
    lam = np.exp(zs)
    moments = fixed.q1.coefficient_logs(*fixed.q1.region_logs(*fixed.q1.epoch_logs(128,19,a,lam),L//128),B)
    weights = np.array(sorted(counts),dtype=int)
    log_counts = np.array([math.log(counts[w]) for w in weights])
    best = np.full((len(distances),maximum),np.inf)
    for i, cutoff in enumerate(cutoffs):
        values = np.min(np.minimum(0.,moments+cutoff*lam[:,None]),axis=0)
        best[i,0] = math.log(L)+float(np.logaddexp.reduce(log_counts+values[weights]))
    pairs = pair.PairKernel()
    pair_best = np.full((len(distances),len(weights),len(weights)),np.inf)
    pair_zs = np.arange(-44,-7)/4
    regions = pair.region_logs(128,19,a,np.exp(pair_zs),L)
    for j,z in enumerate(pair_zs):
        moment = pairs.coefficients(regions[j],B)[np.ix_(weights,weights)]
        for i,cutoff in enumerate(cutoffs):
            np.minimum(pair_best[i],np.minimum(0.,moment+cutoff*math.exp(z)),out=pair_best[i])
    for i in range(len(distances)):
        terms = pair_best[i]+log_counts[:,None]+log_counts[None,:]
        best[i,1] = math.log(math.comb(L,2))+float(np.logaddexp.reduce(terms.ravel()))
    print('Exact-support Q1/Q2 margins:',(-best[:,:2]/math.log(2)).tolist(),flush=True)
    compositions = {q:fixed.composition.SparseComposition(counts,B,q) for q in range(3,min(4,maximum)+1)}
    component_best = {q:np.full((len(distances),len(model.indices)),np.inf) for q,model in compositions.items()}
    sparse_zs = np.arange(-112,1,dtype=float)/8
    shifts = [-.5,0.,.5,1.]
    for z in sparse_zs:
        lam = math.exp(z)
        epoch = fixed.general.epoch_logs(128,19,a,kernel,lam,min(128,maximum))
        region = fixed.general.region_logs(epoch,128,L,maximum)
        for i,cutoff in enumerate(cutoffs):
            for shift in shifts:
                values = fixed.general.occupation_bounds(region,counts,B,L,cutoff,lam,shift,8)
                np.minimum(best[i,2:],values[2:],out=best[i,2:])
                for q,model in compositions.items():
                    values = model.components(region,cutoff,lam,shift)
                    np.minimum(component_best[q][i],values,out=component_best[q][i])
    for q,model in compositions.items():
        for i in range(len(distances)):
            best[i,q-1] = min(best[i,q-1],model.aggregate(component_best[q][i],L))
    rows = []
    for i,delta in enumerate(distances):
        rows.append(dict(message_exponent=exponent,message_bits=1<<exponent,output_bits=N,outer_rows=L,
            distance_target=str(delta),bad_weight=cutoffs[i],covered_occupations=[1,maximum],
            occupation_margins_bits=(-best[i]/math.log(2)).tolist(),
            sparse_union_margin_bits=-float(np.logaddexp.reduce(best[i]))/math.log(2),
            worst_occupation=int(np.argmax(best[i]))+1))
    return rows, dict(q1_log_surprisals=zs.tolist(),q2_log_surprisals=pair_zs.tolist(),
                     sparse_log_surprisals=sparse_zs.tolist(),sparse_shifts=shifts)


def replay_dense(result, counts, a, kernel, row):
    B, L, N = 128,row['outer_rows'],row['output_bits']
    fixed.outer.require(result['occupation_max']==L,'dense endpoint mismatch')
    check_coverage(result['selected_boxes'],L,result['occupation_min'])
    logs, epochs = [], {}
    for box in result['selected_boxes']:
        w = box['witness']; z = w['log_surprisal']; lam = math.exp(z)
        if z not in epochs:
            epochs[z] = fixed.general.epoch_logs(128,19,a,kernel,lam,128)
        bank = result['probability_banks'][w['probability_bank']]
        p,costs,proposal = np.array(bank['probabilities']),np.array(bank['log_density_costs']),np.array(w['proposal'])
        for g,band in enumerate(result['bands'],1):
            if band == [B]:
                fixed.outer.require(p[g]==1 and costs[g]==0,'invalid all-one atom')
                continue
            cost = max(math.log(counts[v])-math.log(math.comb(B,v))-v*math.log(p[g])-(B-v)*math.log1p(-p[g]) for v in band)
            fixed.outer.require(abs(cost-costs[g])<1e-10,'invalid density cost')
        moment,_ = typed.moment_and_density(epochs[z],128,N,float(proposal@p))
        corners = typed.vertices(box['lower'],box['upper'],L)
        value = max(typed.point_logs(corners,L,B,costs,proposal,moment,row['bad_weight'],lam))
        value += typed.lattice_log_count(box['lower'],box['upper'])
        fixed.outer.require(abs(value-box['own_log_bound'])<2e-6,'dense replay mismatch')
        logs.append(value)
    total = float(np.logaddexp.reduce(logs))
    fixed.outer.require(abs(total-result['log_union_upper'])<2e-6,'dense union mismatch')
    return total


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--exponent',type=int,default=20)
    parser.add_argument('--distances',type=Fraction,nargs='+',default=[Fraction(33,200),Fraction(19,100)])
    parser.add_argument('--maximum-sparse',type=int,default=64)
    parser.add_argument('--dense-nodes',type=int,default=511)
    parser.add_argument('--sparse-only',action='store_true')
    parser.add_argument('--output',type=Path,default=fixed.HERE/'SMALLER_MARGIN_SCREEN.json')
    args = parser.parse_args()
    fixed.outer.require(args.maximum_sparse>=4 and all(0<d<Fraction(1,2) for d in args.distances),'invalid grid')
    construction = smaller_outer.construction()
    values = smaller_outer.spectrum()
    counts = {w:n for w,n in enumerate(values) if w and n}
    a,kernel,sources = fixed.load_inner()
    rows,grid = sparse_results(counts,a,kernel,args.exponent,args.distances,args.maximum_sparse)
    for row in rows:
        print({key:row[key] for key in ('distance_target','sparse_union_margin_bits','worst_occupation')},flush=True)
    if not args.sparse_only:
        evaluator = typed.TypedDense(counts,128,128,19,a,kernel,rows[0]['outer_rows'],
            np.arange(-12.,1.01,.5),probability_scales=(1.,),bands=[sorted(w for w in counts if w!=128),[128]])
        for row in rows:
            evaluator.cutoff = row['bad_weight']
            result = evaluator.search(args.maximum_sparse+1,args.dense_nodes,target_bits=60)
            total = replay_dense(result,counts,a,kernel,row)
            row['dense'] = result
            row['dense_union_margin_bits'] = -total/math.log(2)
            row['combined_margin_bits'] = -float(np.logaddexp(total,-row['sparse_union_margin_bits']*math.log(2)))/math.log(2)
            row['exact_integer_coverage_checked'] = row['dense_witness_replay_passed'] = True
            print({key:row[key] for key in ('distance_target','dense_union_margin_bits','combined_margin_bits')},flush=True)
    sources += [Path(__file__),Path(smaller_outer.__file__),Path(fixed.__file__),Path(fixed.outer.__file__),
        Path(pair.__file__),fixed.LANDSCAPE/'activation_q2_kernel.cpp',fixed.LANDSCAPE/'build_activation_q2.ps1',
        Path(fixed.q1.__file__),Path(fixed.general.__file__),Path(fixed.composition.__file__),Path(fixed.maps.__file__),
        Path(typed.__file__),Path(typed.dense.__file__),fixed.HERE/'verify_fixed_inner_results.py',
        fixed.HERE/'sources/EBCH128_29.wd',fixed.HERE/'sources/EBCH128_36.wd']
    payload = dict(status='BINARY64_DIAGNOSTIC_NOT_OUTWARD_CERTIFICATE',inner='existing optimized t128_s19',
        inner_reoptimized=False,construction=construction,spectrum={str(w):str(n) for w,n in counts.items()},
        arithmetic='nearest binary64',grid=grid,results=rows,
        arguments={k:str(v) if isinstance(v,Path) else [str(x) for x in v] if k=='distances' else v for k,v in vars(args).items()},
        source_sha256={p.relative_to(fixed.ROOT).as_posix():fixed.sha(p) for p in sources})
    args.output.write_text(json.dumps(payload,indent=2,default=lambda x:x.tolist())+'\n',encoding='utf-8',newline='\n')


if __name__ == '__main__':
    main()
