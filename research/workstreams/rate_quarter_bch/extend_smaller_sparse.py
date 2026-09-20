"""Close the small-occupation boundary without changing the outer or inner."""
import json
import math
from pathlib import Path

import numpy as np
import evaluate_smaller_margins as base
from refine_smaller_dense import Refiner


def main():
    source = base.fixed.HERE/'SMALLER_MARGIN_REFINED.json'
    output = base.fixed.HERE/'SMALLER_MARGIN_CLOSED.json'
    payload = json.loads(source.read_text())
    for name,digest in payload['source_sha256'].items():
        base.fixed.outer.require(base.fixed.sha(base.fixed.ROOT/name)==digest,'changed input: '+name)
    a,kernel,_ = base.fixed.load_inner()
    counts = {w:n for w,n in enumerate(base.smaller_outer.spectrum()) if w and n}
    row = next(r for r in payload['results'] if r['distance_target']=='19/100')
    maximum = 128
    length = row['outer_rows']
    previous = row['covered_occupations'][1]
    best = np.full(maximum,np.inf)
    best[:previous] = -np.array(row['occupation_margins_bits'])*math.log(2)
    zs = np.arange(-56,1,dtype=float)/8
    shifts = [-.5,0.,.5,1.]
    witnesses = np.zeros((maximum,2))
    for index,z in enumerate(zs):
        lam = math.exp(z)
        epoch = base.fixed.general.epoch_logs(128,19,a,kernel,lam,128)
        regions = base.fixed.general.region_logs(epoch,128,length,maximum)
        for shift in shifts:
            values = base.fixed.general.occupation_bounds(regions,counts,128,length,row['bad_weight'],lam,shift,8)
            improve = values<best
            best[improve] = values[improve]
            witnesses[improve] = [z,shift]
        if index%8==0: print(f'extended sparse tilts: {index+1}/{len(zs)}',flush=True)
    row['occupation_margins_bits'] = (-best/math.log(2)).tolist()
    row['sparse_union_margin_bits'] = -float(np.logaddexp.reduce(best))/math.log(2)
    row['worst_occupation'] = int(np.argmax(best))+1
    row['covered_occupations'] = [1,maximum]
    row['sparse_extension'] = dict(previous_maximum=previous,log_surprisals=zs.tolist(),shifts=shifts,
                                 extra_witnesses=witnesses[previous:].tolist())
    # Intersect each disjoint dense box with the new occupation domain.
    # Its old fixed witness remains valid on that smaller box.
    model = Refiner(row,a,kernel)
    clipped = []
    for box in row['dense']['selected_boxes']:
        lower,upper = np.array(box['lower']),np.array(box['upper'])
        upper[0] = min(upper[0],length-maximum-1)
        corners = base.typed.vertices(lower,upper,length)
        if not len(corners): continue
        witness = dict(box['witness'])
        witness.pop('maximum_vertex',None)  # old search hint need not lie in the clipped box
        new = dict(lower=corners.min(axis=0).tolist(),upper=corners.max(axis=0).tolist(),witness=witness)
        new['own_log_bound'] = model.direct(new,new['witness'])[0]
        clipped.append(new)
    row['dense']['selected_boxes'] = clipped
    row['dense']['occupation_min'] = maximum+1
    row['dense']['log_union_upper'] = float(np.logaddexp.reduce([b['own_log_bound'] for b in clipped]))
    dense = base.replay_dense(row['dense'],counts,a,kernel,row)
    row['dense_union_margin_bits'] = -dense/math.log(2)
    row['combined_margin_bits'] = -float(np.logaddexp(dense,-row['sparse_union_margin_bits']*math.log(2)))/math.log(2)
    row['dense_clipped_to_new_sparse_boundary'] = True
    payload['source_sha256'][Path(__file__).relative_to(base.fixed.ROOT).as_posix()] = base.fixed.sha(Path(__file__))
    payload['source_sha256'][source.relative_to(base.fixed.ROOT).as_posix()] = base.fixed.sha(source)
    output.write_text(json.dumps(payload,indent=2)+'\n',encoding='utf-8',newline='\n')
    print('FINAL',[(r['distance_target'],r['sparse_union_margin_bits'],r['dense_union_margin_bits'],r['combined_margin_bits']) for r in payload['results']],flush=True)


if __name__ == '__main__':
    main()
