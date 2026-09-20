"""Sparse IMT compositions with reference logits scaled toward one half."""
import argparse
import math
from pathlib import Path

import numpy as np
import parameter_sparse_compositions as prior

full = prior.full


def bank(problem,scale):
    roots,active,inactive = [],[],[]
    for band in problem.bands:
        if band == [problem.block]:
            roots.append(math.log(problem.counts[problem.block])/problem.block)
            active.append(0.); inactive.append(-math.inf)
            continue
        p = (min(band)+max(band))/(2*problem.block)
        eta = scale*(math.log(p)-math.log1p(-p))
        lp,ln = -float(np.logaddexp(0.,-eta)),-float(np.logaddexp(0.,eta))
        roots.append(max(math.log(problem.counts[w])-math.log(math.comb(problem.block,w))
                         -w*lp-(problem.block-w)*ln for w in band)/problem.block)
        active.append(lp); inactive.append(ln)
    return (prior.composition.probability_logs(problem.indices,np.array(active),np.array(inactive)),
            problem.block*np.sum(np.array(roots)[problem.indices],axis=1))


def run(output,maximum=4):
    assert not output.exists()
    b,t,s,exponent = 128,64,20,20
    record,bs,kernel,caps,low = full.prepare(t,s)
    counts,_ = full.grid.outer(b)
    rows = (1 << exponent)//(b//2)
    cutoff = b*rows//10
    problems = [prior.composition.SparseComposition(counts,b,q,tail_bands=6,singleton_prefix=3)
                for q in range(2,maximum+1)]
    scales = (0.,0.5,1.)
    banks = [[bank(p,scale) for scale in scales] for p in problems]
    best = [np.full(len(p.indices),np.inf) for p in problems]
    witnesses = [[None]*len(p.indices) for p in problems]
    lo,hi = math.log(2/rows)-1,math.log(maximum/rows)+3
    tilts = np.arange(math.floor(10*lo),math.ceil(10*hi)+1,dtype=float)/10
    for i,tilt in enumerate(tilts):
        lam = math.exp(tilt)
        epochs = full.g.epochs(record['spectrum'],kernel,record['feedback_columns'],caps,low,lam,maximum=min(t,maximum))
        regions = full.g.regions(epochs,t,rows,maximum)
        for k,problem in enumerate(problems):
            for scale,(probabilities,cost) in zip(scales,banks[k]):
                matrices = np.full((len(probabilities),regions.shape[1],regions.shape[2]),-np.inf)
                for j in range(problem.occupation+1):
                    np.logaddexp(matrices,probabilities[:,j,None,None]+regions[j],out=matrices)
                values = np.minimum(problem.trivial,prior.terminal(matrices,b)+cost+cutoff*lam)
                for j in np.flatnonzero(values < best[k]):
                    witnesses[k][int(j)] = dict(tilt=float(tilt),scale=scale)
                best[k] = np.minimum(best[k],values)
        if i % 10 == 0:
            print('balanced IMT compositions',i+1,'/',len(tilts),
                [-p.aggregate(v,rows)/math.log(2) for p,v in zip(problems,best)],flush=True)
    results = [dict(q=p.occupation,bands=p.bands,indices=p.indices.tolist(),
        component_log_bounds=v.tolist(),witnesses=w,log_bound=p.aggregate(v,rows),
        margin_bits=-p.aggregate(v,rows)/math.log(2)) for p,v,w in zip(problems,best,witnesses)]
    sources = full.grid.ladder.candidate.sources()
    _,path = full.grid.outer(b)
    sources[path.relative_to(full.grid.ROOT).as_posix()] = full.grid.maps.sha(path)
    full.grid.ladder.model.base.write_new(output,dict(status='BINARY64_IMT_BALANCED_COMPOSITIONS',
        geometry=[b,t,s,exponent],inner=record,results=results,
        full_distance_proved=False,source_sha256=sources))
    print('balanced sparse margins',[r['margin_bits'] for r in results],flush=True)


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output',type=Path,required=True)
    p.add_argument('--maximum',type=int,choices=(2,3,4),default=4)
    a = p.parse_args()
    run(a.output.resolve(),a.maximum)
