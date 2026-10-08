"""Proof gates for lower-work GF256 constructions at 80--84 percent rate.

This reuses the audited full-Cauchy punctured local transfer unchanged.
Only the outer length/dimension and packet geometry vary. Results are
optimized floating upper bounds, not distance certificates or evidence
that a negative bound describes an actual low-distance code.
"""
from __future__ import annotations

import argparse
from dataclasses import dataclass
from fractions import Fraction as F
import json

from gf8_punctured_screen import profile,search


@dataclass(frozen=True)
class Candidate:
    name:str
    length:int
    dimension:int
    width:int
    outputs:int
    state:int

    @property
    def rate(self):
        return F(self.dimension*self.width,self.length*self.outputs)

    def geometry(self):
        if not F(4,5)<=self.rate<=F(21,25):
            raise ValueError('candidate is outside the 80--84 percent target')
        g=profile(self.width,self.outputs,self.length,self.dimension)
        if self.width+self.state>128:
            raise ValueError('chosen Cauchy row/column node sets would overlap')
        return g


CANDIDATES=(
    Candidate('control128-80-w8-s8-o6',128,80,8,6,8),
    Candidate('outer64-40-w8-s8-o6',64,40,8,6,8),
    Candidate('outer64-39-w8-s8-o6',64,39,8,6,8),
    Candidate('outer64-45-w8-s8-o7',64,45,8,7,8),
    Candidate('outer128-80-w4-s4-o3',128,80,4,3,4),
    Candidate('outer128-78-w4-s4-o3',128,78,4,3,4),
    Candidate('outer64-40-w4-s4-o3',64,40,4,3,4),
    Candidate('outer64-39-w4-s4-o3',64,39,4,3,4),
    Candidate('outer128-80-w8-s4-o6',128,80,8,6,4),
    Candidate('outer128-78-w8-s4-o6',128,78,8,6,4),
    Candidate('outer64-45-w8-s4-o7',64,45,8,7,4),
)


def screen(candidate,value,sparse=False,iterations=70,start=None):
    candidate.geometry()
    result=search(value,sparse,candidate.width,candidate.outputs,candidate.state,
                  candidate.length,candidate.dimension,iterations,start)
    result.update(candidate=candidate.name,nominal_rate_exact=str(candidate.rate),
                  requested_rate=F(2**20,result['output_n']).__str__(),
                  output_distance='1/10',proof_target_bits=40,
                  local_transfer='audited full-Cauchy positive two-state zero-flat cover',
                  setup='independent uniform regional permutations and nonzero GF256 diagonals',
                  native_claim=False)
    return result


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--candidate',nargs='+',choices=[c.name for c in CANDIDATES])
    parser.add_argument('--sparse',nargs='*',type=int,default=[64,128])
    parser.add_argument('--density',nargs='*',type=float,default=[.03,.10])
    parser.add_argument('--iterations',type=int,default=70)
    args=parser.parse_args()
    for c in CANDIDATES:
        if args.candidate and c.name not in args.candidate:
            continue
        for sparse,values in ((True,args.sparse),(False,args.density)):
            previous=None
            for value in values:
                result=screen(c,value,sparse,args.iterations,previous)
                previous=result['log_parameters']
                print(json.dumps(result),flush=True)
