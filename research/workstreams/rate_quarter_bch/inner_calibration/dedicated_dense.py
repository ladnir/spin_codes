"""Fresh, candidate-specific dense search with retained numerical witnesses.

Uses the existing three-category convex counting bound and adaptive partition.
No old box or witness seeds the fresh search. All arithmetic is nearest binary64.
"""
import argparse
import copy
import json
import math
from pathlib import Path
import time

import numpy as np
import check_candidates as check


class LoggedDense(check.typed.TypedDense):
    def __init__(self,*args,**kwargs):
        self.evaluations = 0
        self.started = time.monotonic()
        super().__init__(*args,**kwargs)

    def evaluate(self,lower,upper):
        result = super().evaluate(lower,upper)
        self.evaluations += 1
        if self.evaluations%16 == 1:
            print(dict(nodes=self.evaluations,elapsed_seconds=round(time.monotonic()-self.started,1),
                       current_box_margin_bits=-result['own_log_bound']/math.log(2) if result else None),flush=True)
        return result


def validate_receipt(path):
    payload = json.loads(path.read_text())
    for name,digest in payload['source_sha256'].items():
        check.fixed.outer.require(check.fixed.sha(check.fixed.ROOT/name)==digest,'changed dependency: '+name)
    return payload


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--tag',default='t256_s18_nested')
    parser.add_argument('--nodes',type=int,default=1023)
    parser.add_argument('--distances',nargs='+',default=['33/200','19/100'])
    args = parser.parse_args()
    folder = check.calibration.HERE
    source = folder/f'{args.tag}_full.json'
    retained = validate_receipt(source)
    map_path = folder/'maps'/f'{args.tag}.json'
    assert map_path.relative_to(check.fixed.ROOT).as_posix() in retained['source_sha256']
    record = json.loads(map_path.read_text())
    t,s = record['step_bits'],record['state_bits']
    a = {w:n for w,n in enumerate(record['a_counts']) if w and n}
    counts = {w:n for w,n in enumerate(check.calibration.smaller_outer.spectrum()) if w and n}
    tilts = np.arange(-32,5,dtype=float)/4
    print('Preparing fresh model',args.tag,'tilts',tilts.tolist(),flush=True)
    model = LoggedDense(counts,check.B,t,s,a,record['kernel_counts'],check.L,tilts,
                       probability_scales=(1.,),bands=[sorted(w for w in counts if w!=check.B),[check.B]])
    paths = [Path(__file__),source,map_path,Path(check.__file__),Path(check.typed.__file__),
             Path(check.typed.dense.__file__),Path(check.fixed.general.__file__),Path(check.fixed.composition.__file__),
             Path(check.fixed.q1.__file__),check.fixed.HERE/'verify_fixed_inner_results.py']
    payload = dict(status='FRESH_DENSE_BINARY64_NOT_OUTWARD_CERTIFICATE',tag=args.tag,outer=[128,32,32],
                   message_exponent=20,arithmetic='nearest binary64',
                   arguments=vars(args),grid=dict(log_surprisals=tilts.tolist(),probability_scales=[1.],dense_target_bits=60),
                   results=[],source_sha256={p.relative_to(check.fixed.ROOT).as_posix():check.fixed.sha(p) for p in paths})
    output = folder/f'{args.tag}_dedicated_dense.json'
    for distance in args.distances:
        original = next(r for r in retained['results'] if r['distance_target']==distance)
        row = {k:copy.deepcopy(v) for k,v in original.items() if k not in ('dense','dense_union_margin_bits','combined_margin_bits')}
        model.cutoff = original['bad_weight']
        model.evaluations = 0; model.started = time.monotonic()
        print('START',distance,flush=True)
        dense = model.search(original['covered_occupations'][1]+1,args.nodes,target_bits=60)
        row['dense'] = dense
        replay = check.transport_dense(record,counts,row,dense['occupation_min'])
        # Replay tightens redundant coordinate bounds to simplex extrema.
        # The represented integer set is unchanged, but the lattice-count
        # penalty can decrease. This is a valid improvement, not a mismatch.
        assert replay['log_union_upper']<=dense['log_union_upper']+2e-6
        for b,c in zip(replay['selected_boxes'],dense['selected_boxes']):
            np.testing.assert_array_equal(check.typed.vertices(b['lower'],b['upper'],check.L),
                                          check.typed.vertices(c['lower'],c['upper'],check.L))
            penalty_change = (check.typed.lattice_log_count(b['lower'],b['upper'])
                              -check.typed.lattice_log_count(c['lower'],c['upper']))
            assert abs(b['own_log_bound']-c['own_log_bound']-penalty_change)<2e-6
        row['search_log_union_upper_before_coordinate_tightening'] = dense['log_union_upper']
        row['dense'] = dense = replay
        row['dense_union_margin_bits'] = -dense['log_union_upper']/math.log(2)
        row['combined_margin_bits'] = -float(np.logaddexp(dense['log_union_upper'],-row['sparse_union_margin_bits']*math.log(2)))/math.log(2)
        row['dense_witness_replay_passed'] = True
        row['fresh_search_elapsed_seconds'] = time.monotonic()-model.started
        row['old_transported_dense_margin_bits'] = original['dense_union_margin_bits']
        payload['results'].append(row)
        output.write_text(json.dumps(payload,indent=2,default=lambda x:x.tolist())+'\n',encoding='utf-8',newline='\n')
        print('FINISHED',distance,'dense margin',row['dense_union_margin_bits'],'combined',row['combined_margin_bits'],flush=True)


if __name__=='__main__': main()
