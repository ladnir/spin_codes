"""Small exact checks for compressed high-occupancy shape envelopes."""
from itertools import product
from math import comb
from collections import Counter
from fractions import Fraction
from flint import ctx,arb,arb_mat

from occupancy_memory import shape_classes,coarse_epoch,prepare,epoch_operators,F,Z,C,TERMINAL
from occupancy_model import local_data
from group_moment import maps


def restricted_pairs(data):
    _,columns,_=maps()
    syndromes=[]
    for j in range(32):
        row=[0]*16
        for mask in range(1,16):
            bit=mask&-mask
            row[mask]=row[mask^bit]^columns[4*j+bit.bit_length()-1]
        syndromes.append(row)
    peak=max(Fraction(row[2],row[0]) for row in data[1].values())
    tests=0
    for available in ((0,1),(0,31),(0,7,31),tuple(range(0,32,2)),tuple(range(16,32))):
        for a,b in product(range(1,5),repeat=2):
            counts=Counter(syndromes[j][ma]^syndromes[k][mb]
                           for j in available for k in available if j!=k
                           for ma in range(1,16) if ma.bit_count()==a
                           for mb in range(1,16) if mb.bit_count()==b)
            denominator=len(available)*(len(available)-1)*comb(4,a)*comb(4,b)
            assert counts[0]==0 and sum(counts.values())==denominator
            assert Fraction(max(counts.values()),denominator)<=peak*32*31/(len(available)*(len(available)-1))
            tests+=1
    print('Restricted-window pair tests passed:',tests,'global peak',peak,flush=True)


def main():
    for q in range(1,9):
        exact={(sum(weights),max(comb(4,w) for w in weights)) for weights in product(range(1,5),repeat=q)}
        assert shape_classes(q)==exact
        full={(sum(weights),max(comb(4,w) for w in weights),weights.count(4)) for weights in product(range(1,5),repeat=q)}
        assert shape_classes(q,True)==full
    print('Compressed shape classes match every weight tuple at q=1..8',flush=True)
    prepared=prepare(local_data(4))
    restricted_pairs(prepared[0][1])
    ctx.prec=192
    spectrum=prepared[0][0][0]
    terminal=arb_mat([[int(x)] for x in TERMINAL])
    for q in (3,4,5,8,16,32):
        t=coarse_epoch(spectrum,q,'0')
        mass=t*terminal
        assert all(mass[i,0].upper()>=terminal[i,0] for i in range(9))
    for tilt in ('.0016','.008'):
        detailed=epoch_operators(prepared,tilt,detailed=True)
        weighted=epoch_operators(prepared,tilt,detailed=True,input_penalty='1.01',full_penalty='.5')
        for q in range(1,5):
            for weights,t in detailed[q].items():
                scale=arb('1.01')**sum(weights)*arb('.5')**weights.count(4)
                for i in range(9):
                    for j in range(9):
                        assert abs(weighted[q][weights][i,j]-t[i,j]*scale)<arb(2)**-165
        for q in (3,4):
            coarse=coarse_epoch(spectrum,q,tilt)
            for weights,t in detailed[q].items():
                for i in range(9):
                    for j in range(9):
                        # These two entries additionally exploit the fresh
                        # source's 1/32 density cap, unlike the older bound.
                        if i==F and j in (Z,C):
                            continue
                        assert coarse[i,j]+arb(2)**-170>=t[i,j]
    print('Tail mass checks and lower-occupancy envelope comparisons passed',flush=True)


if __name__=='__main__':
    main()
