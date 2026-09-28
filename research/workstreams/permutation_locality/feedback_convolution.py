"""Exact independent-window feedback convolutions for a future state envelope.

These counts assume two independent window/mask choices. They are not yet
integrated with SPIN's without-replacement placement or tilted moments.
"""
from collections import Counter
from fractions import Fraction
from itertools import combinations_with_replacement

from two_column_moment import census


def counts(atoms):
    result={}
    for a,b in combinations_with_replacement(sorted(atoms),2):
        convolution=Counter()
        for x,cx in atoms[a].items():
            for y,cy in atoms[b].items():
                convolution[x^y]+=cx*cy
        denominator=sum(atoms[a].values())*sum(atoms[b].values())
        assert sum(convolution.values())==denominator
        assert convolution[0]==sum(count*atoms[b][x] for x,count in atoms[a].items())
        peak=max(count for x,count in convolution.items() if x)
        result[a,b]=(Fraction(convolution[0],denominator),Fraction(peak,denominator))
    return result


def main():
    data=census(1)
    atoms={shape:Counter({x:c for x,c in values.items() if c}) for shape,values in data[3].items()}
    single=max(Fraction(max(values.values()),sum(values.values())) for values in atoms.values())
    result=counts(atoms)
    print('Independent single feedback maximum atom:',single,flush=True)
    for shapes,(zero,nonzero) in result.items():
        print('weights',shapes,'pair zero',zero,'pair maximum nonzero atom',nonzero,flush=True)
    print('Worst two-input nonzero atom:',max(p for _,p in result.values()),flush=True)
    print('Counts only: no without-replacement or tilted-state certificate.')


if __name__=='__main__':
    main()
