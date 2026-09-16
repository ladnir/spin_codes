"""Setup-independent distance ceiling from a restricted zero-state subspace.

Assumes a binary linear [B,D] outer on L independent rows, followed only
by coordinate permutations before an inner with X_i -> X_i+A q_i,
q_{i+1}=alpha_i q_i+B_i X_i, q_0=0. Each B_i has at most s output bits.
Restrict to Q outer rows and impose B_i X_i=0 in all E=N/t epochs.
The resulting subspace has dimension at least r=DQ-sE. Its support is
at most BQ, and averaging nonzero words bounds d_min by
floor(BQ * 2^(r-1)/(2^r-1)). No setup probability or spectrum is needed.
"""
import json
from fractions import Fraction
from pathlib import Path

import check_candidates as check


def ceiling(block,dimension,rows,t,s):
    assert 0<dimension<=block and rows>0 and t>0 and s>0
    n = block*rows
    assert n%t==0
    epochs = n//t
    constraints = s*epochs
    q = constraints//dimension+1
    if q>rows:
        return None
    rank_lower = dimension*q-constraints
    support = block*q
    bound = support*(1<<(rank_lower-1))//((1<<rank_lower)-1)
    return dict(outer_block=block,outer_dimension=dimension,outer_rows=rows,
                output_bits=n,step_bits=t,state_bits=s,epochs=epochs,
                maximum_syndrome_constraints=constraints,restricted_outer_rows=q,
                zero_state_subspace_dimension_lower=rank_lower,
                support_upper=support,minimum_distance_upper=bound,
                relative_distance_upper=str(Fraction(bound,n)))


def main():
    rows = []
    for s in range(14,29):
        row = ceiling(128,32,check.L,256,s)
        row['targets_impossible_for_every_setup'] = [str(d) for d in check.DISTANCES
             if row['minimum_distance_upper']<=check.N*d.numerator//d.denominator]
        rows.append(row)
    baseline = ceiling(128,32,check.L,128,19)
    payload = dict(status='DETERMINISTIC_DISTANCE_UPPER_BOUND_EXACT_INTEGER',
                   construction_scope=__doc__,results=rows,baseline=baseline,
                   limitations=['A ceiling above a target does not show that target is achievable.',
                                'Changing the pre-inner transform beyond coordinate permutations requires a new support argument.'],
                   source_sha256={p.relative_to(check.fixed.ROOT).as_posix():check.fixed.sha(p)
                                  for p in (Path(__file__),Path(check.__file__))})
    output = check.calibration.HERE/'ZERO_STATE_RANK_CEILING.json'
    output.write_text(json.dumps(payload,indent=2)+'\n',encoding='utf-8',newline='\n')
    for row in rows:
        print(row['state_bits'],row['minimum_distance_upper'],float(Fraction(row['relative_distance_upper'])),
              row['targets_impossible_for_every_setup'],flush=True)


if __name__=='__main__': main()
