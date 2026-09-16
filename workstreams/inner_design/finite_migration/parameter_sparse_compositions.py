"""Keep small-Q band compositions fixed across the whole IMT computation."""
import argparse
import math
from pathlib import Path

import numpy as np
import parameter_full as full
import composition_occupation as composition


def terminal(matrices,length):
    current = np.full_like(matrices,-np.inf)
    n = matrices.shape[1]
    current[:,np.arange(n),np.arange(n)] = 0
    power = matrices
    while length:
        if length & 1:
            current = full.grid.wm.product(current,power)
        length >>= 1
        if length:
            power = full.grid.wm.product(power,power)
    return np.logaddexp.reduce(current[:,0],axis=1)


def components(problem,regions,cutoff,lam,shift):
    if shift not in problem.auxiliary_cache:
        roots,active,inactive = composition.general.density_roots(problem.counts,problem.block,problem.bands,shift)
        problem.auxiliary_cache[shift] = (
            composition.probability_logs(problem.indices,active,inactive),
            problem.block*np.sum(roots[problem.indices],axis=1))
    probabilities,cost = problem.auxiliary_cache[shift]
    matrices = np.full((len(probabilities),regions.shape[1],regions.shape[2]),-np.inf)
    for j in range(problem.occupation+1):
        np.logaddexp(matrices,probabilities[:,j,None,None]+regions[j],out=matrices)
    return np.minimum(problem.trivial,terminal(matrices,problem.block)+cost+cutoff*lam)


def run(output,b,t,s,exponent,maximum):
    assert not output.exists()
    record,bs,kernel,caps,low = full.prepare(t,s)
    counts,_ = full.grid.outer(b)
    rows = (1 << exponent)//(b//2)
    cutoff = b*rows//10
    problems = [composition.SparseComposition(counts,b,q,tail_bands=6,singleton_prefix=3)
                for q in range(2,maximum+1)]
    best = [np.full(len(p.indices),np.inf) for p in problems]
    witnesses = [[None]*len(p.indices) for p in problems]
    lo,hi = math.log(2/rows)-1,math.log(maximum/rows)+3
    tilts = np.arange(math.floor(10*lo),math.ceil(10*hi)+1,dtype=float)/10
    for i,tilt in enumerate(tilts):
        lam = math.exp(tilt)
        epochs = full.g.epochs(record['spectrum'],kernel,record['feedback_columns'],caps,low,lam,maximum=min(t,maximum))
        regions = full.g.regions(epochs,t,rows,maximum)
        for shift in (-0.5,0.,0.5):
            for k,problem in enumerate(problems):
                values = components(problem,regions,cutoff,lam,shift)
                improved = values < best[k]
                for j in np.flatnonzero(improved):
                    witnesses[k][int(j)] = dict(tilt=float(tilt),shift=shift)
                best[k] = np.minimum(best[k],values)
        if i % 10 == 0:
            print('IMT fixed compositions',i+1,'/',len(tilts),
                  [-p.aggregate(v,rows)/math.log(2) for p,v in zip(problems,best)],flush=True)
    results = [dict(q=p.occupation,bands=p.bands,indices=p.indices.tolist(),
        component_log_bounds=v.tolist(),witnesses=w,log_bound=p.aggregate(v,rows),
        margin_bits=-p.aggregate(v,rows)/math.log(2)) for p,v,w in zip(problems,best,witnesses)]
    sources = full.grid.ladder.candidate.sources()
    _,path = full.grid.outer(b)
    sources[path.relative_to(full.grid.ROOT).as_posix()] = full.grid.maps.sha(path)
    full.grid.ladder.model.base.write_new(output,dict(status='BINARY64_IMT_SPARSE_COMPOSITIONS',
        geometry=[b,t,s,exponent],inner=record,results=results,
        full_distance_proved=False,source_sha256=sources))
    print('completed sparse margins',[r['margin_bits'] for r in results],flush=True)


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output',type=Path,required=True)
    p.add_argument('--maximum',type=int,choices=(2,3,4),default=4)
    a = p.parse_args()
    run(a.output.resolve(),128,64,20,20,a.maximum)
