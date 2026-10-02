"""Exact local GF16 feedback census; not a distance certificate."""
from collections import Counter
from fractions import Fraction
from itertools import combinations
from pathlib import Path
import sys

sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from group_moment import maps


def rank(columns):
    basis={}
    for x in columns:
        while x:
            bit=x.bit_length()-1
            if bit not in basis:
                basis[bit]=x
                break
            x^=basis[bit]
    return len(basis)


def census(columns,max_packets=4):
    if not columns or len(columns)%4:raise ValueError('complete four-bit packets required')
    spans=[]
    for start in range(0,len(columns),4):
        values=[0]
        for c in columns[start:start+4]:values += [v^c for v in values]
        spans.append(Counter(values[1:]))
    # A uniform nonzero input has each of these 15 preimages equally likely.
    zero_counts=Counter(v[0] for v in spans)
    pair_zero=Counter()
    pair_max=Counter()
    for a,b in combinations(spans,2):
        values=Counter()
        for x,nx in a.items():
            for y,ny in b.items():values[x^y]+=nx*ny
        pair_zero[values[0]]+=1
        pair_max[max(values.values())]+=1
    ranks={j:Counter(rank([c for w in windows for c in columns[4*w:4*w+4]])
                     for windows in combinations(range(len(spans)),j))
           for j in range(1,min(max_packets,len(spans))+1)}
    atom_caps={j:str(Fraction(2**(4*j-min(counts)),15**j)) for j,counts in ranks.items()}
    return dict(packet_count=len(spans),distinct_feedback_images=dict(Counter(len(v) for v in spans)),
                zero_preimages_out_of_15=dict(zero_counts),
                pair_zero_preimages_out_of_225=dict(pair_zero),
                pair_largest_atom_out_of_225=dict(pair_max),
                packet_set_rank_histograms={j:dict(v) for j,v in ranks.items()},
                worst_case_atom_caps=atom_caps)


if __name__=='__main__':
    _,columns,_=maps()
    print(census(columns))
