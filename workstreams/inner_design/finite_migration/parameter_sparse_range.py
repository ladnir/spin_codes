"""Numerical IMT bounds for every occupancy in a moderate sparse range."""
import argparse
import math
from pathlib import Path

import numpy as np
import parameter_full as full
import positive_regions


def run(output,b,t,s,exponent,first,last):
    assert not output.exists() and 2 <= first <= last
    record,bs,kernel,caps,low = full.prepare(t,s)
    counts,_ = full.grid.outer(b)
    rows = (1 << exponent)//(b//2)
    assert last < rows and rows % t == 0
    cutoff = b*rows//10
    best = np.full(last-first+1,np.inf)
    witnesses = [None]*len(best)
    banks = []
    for scale in (0.,0.5,1.):
        bands,roots,active,inactive = full.balanced.density_roots(counts,b,scale)
        banks.append((scale,len(bands),full.balanced.Envelope(roots,active,inactive)))
    lo,hi = math.log(first/rows)-1,math.log(last/rows)+3
    tilts = np.arange(math.floor(8*lo),math.ceil(8*hi)+1,dtype=float)/8
    for i,tilt in enumerate(tilts):
        lam = math.exp(tilt)
        epoch = full.g.epochs(record['spectrum'],kernel,record['feedback_columns'],caps,low,lam,maximum=min(t,last))
        regions = positive_regions.regions(epoch,t,rows,last)
        for scale,band_count,envelope in banks:
            current = regions.copy()
            for q in range(1,last+1):
                current = envelope.apply(current[:-1],current[1:])
                if q < first:
                    continue
                bound = (math.log(math.comb(rows,q))+q*math.log(band_count)+cutoff*lam
                         +full.g.terminal(current[0],b))
                if bound < best[q-first]:
                    best[q-first] = bound
                    witnesses[q-first] = dict(tilt=float(tilt),scale=scale)
        if i % 8 == 0:
            print('sparse range',i+1,'/',len(tilts),'worst margin',-float(max(best))/math.log(2),flush=True)
    sources = full.grid.ladder.candidate.sources()
    _,path = full.grid.outer(b)
    sources[path.relative_to(full.grid.ROOT).as_posix()] = full.grid.maps.sha(path)
    result = dict(status='BINARY64_IMT_SPARSE_RANGE',geometry=[b,t,s,exponent],inner=record,
        first=first,last=last,log_bounds=best.tolist(),witnesses=witnesses,
        margin_bits=-float(np.logaddexp.reduce(best))/math.log(2),
        full_distance_proved=False,source_sha256=sources)
    full.grid.ladder.model.base.write_new(output,result)
    print('sparse range done',result['margin_bits'],flush=True)


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output',type=Path,required=True)
    p.add_argument('--b',type=int,choices=(64,128),default=128)
    p.add_argument('--t',type=int,choices=(64,128,256),default=64)
    p.add_argument('--s',type=int,default=20)
    p.add_argument('--m',type=int,default=20)
    p.add_argument('--first',type=int,default=5)
    p.add_argument('--last',type=int,default=64)
    a = p.parse_args()
    run(a.output.resolve(),a.b,a.t,a.s,a.m,a.first,a.last)
