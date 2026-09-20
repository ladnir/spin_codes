"""Select numerical witnesses for a separate outward quarter-rate verifier.

Q2 uses band compositions, not the native pair-support implementation.
Q>=3 uses the simpler adaptive bound. No numerical margin is trusted by
the outward verifier; this file only chooses tilts and counting measures.
"""
from fractions import Fraction as F
import json
import math
from pathlib import Path

import numpy as np
import evaluate_fixed_inner as fixed
import smaller_outer


def main():
    source = fixed.HERE/'SMALLER_MARGIN_CLOSED.json'
    old = json.loads(source.read_text())
    for name,digest in old['source_sha256'].items():
        assert fixed.sha(fixed.ROOT/name)==digest,name
    a,kernel,sources = fixed.load_inner()
    counts = {w:n for w,n in enumerate(smaller_outer.spectrum()) if w and n}
    L,N,B = 32768,4194304,128
    zs = np.arange(-720,1)/40
    lam = np.exp(zs)
    moments = fixed.q1.coefficient_logs(*fixed.q1.region_logs(*fixed.q1.epoch_logs(128,19,a,lam),L//128),B)
    pair = fixed.composition.SparseComposition(counts,B,2)
    rows = []
    best = [np.full(r['covered_occupations'][1],np.inf) for r in old['results']]
    pair_best = np.full((2,len(pair.indices)),np.inf)
    pair_witness = [[None]*len(pair.indices) for _ in range(2)]
    for row in old['results']:
        indices = np.argmin(moments+row['bad_weight']*lam[:,None],axis=0)
        rows.append(dict(distance_target=row['distance_target'],bad_weight=row['bad_weight'],
            maximum_sparse=row['covered_occupations'][1],
            q1=[dict(weight=w,log_tilt=str(F(int(indices[w])-720,40))) for w in counts],
            adaptive=[None]*row['covered_occupations'][1],dense=row['dense']))
    shifts = [-.5,0.,.5,1.]
    for index,numerator in enumerate(range(-112,1)):
        z = numerator/8; lam = math.exp(z)
        maximum = 64 if z < -7 else 128
        region = fixed.general.region_logs(fixed.general.epoch_logs(128,19,a,kernel,lam,maximum),128,L,maximum)
        for i,row in enumerate(rows):
            limit = min(maximum,row['maximum_sparse'])
            for shift in shifts:
                values = fixed.general.occupation_bounds(region[:limit+1],counts,B,L,row['bad_weight'],lam,shift,8)
                improved = values<best[i][:limit]
                best[i][:limit] = np.minimum(best[i][:limit],values)
                for q in np.flatnonzero(improved):
                    row['adaptive'][int(q)] = dict(log_tilt=str(F(numerator,8)),shift=str(F(shift)))
                values = pair.components(region,row['bad_weight'],lam,shift)
                improved = values<pair_best[i]
                np.minimum(pair_best[i],values,out=pair_best[i])
                for c in np.flatnonzero(improved):
                    pair_witness[i][int(c)] = dict(log_tilt=str(F(numerator,8)),shift=str(F(shift)))
        if index%8==0: print('sparse witness search',index+1,'/113',flush=True)
    for i,row in enumerate(rows):
        row['adaptive'] = [dict(occupation=q+1,**w) for q,w in enumerate(row['adaptive']) if q>=2]
        row['q2'] = [dict(band_indices=list(map(int,c)),**w) for c,w in zip(pair.indices,pair_witness[i])]
        row['q2_margin_diagnostic'] = -pair.aggregate(pair_best[i],L)/math.log(2)
        row['adaptive_margin_diagnostic'] = -float(np.logaddexp.reduce(best[i][2:]))/math.log(2)
        print(row['distance_target'],row['q2_margin_diagnostic'],row['adaptive_margin_diagnostic'],flush=True)
    sources += [Path(__file__),source,Path(smaller_outer.__file__),Path(fixed.__file__),Path(fixed.q1.__file__),
                Path(fixed.general.__file__),Path(fixed.composition.__file__)]
    payload = dict(status='WITNESS_SELECTION_ONLY_NOT_CERTIFIED',outer=[128,32,32],message_exponent=20,
                   q2_bands=pair.bands,adaptive_bands=fixed.general.bands_for(counts,B,8),results=rows,
                   source_sha256={p.relative_to(fixed.ROOT).as_posix():fixed.sha(p) for p in sources})
    (fixed.HERE/'SMALLER_OUTWARD_WITNESSES.json').write_text(json.dumps(payload,indent=2)+'\n',encoding='utf-8',newline='\n')


if __name__=='__main__': main()
