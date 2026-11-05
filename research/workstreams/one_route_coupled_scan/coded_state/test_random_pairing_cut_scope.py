"""Exact target-set sensitivity and joint-cardinality relaxation controls.

Standard library only, serial <=30s, no output files. The joint cut
barrier concerns a proof relaxation, not the actual code's distance.
"""
from fractions import Fraction as F
from itertools import product
import shutil
import tempfile
import time
import math


def main():
    assert shutil.disk_usage(tempfile.gettempdir()).free>=20*(1<<30),'disk guard'
    began=time.monotonic();sensitivity=0;mixtures=0
    counts=(0,1,2,3,4,4);M=len(counts);R=4;parity=2
    for a in range(1<<M):
        ga=parity*a.bit_count()-sum(c for j,c in enumerate(counts) if a>>j&1)
        for b in range(1<<M):
            gb=parity*b.bit_count()-sum(c for j,c in enumerate(counts) if b>>j&1)
            removed=(a&~b).bit_count();added=(b&~a).bit_count()
            assert gb-ga <= (R-parity)*removed+parity*added
            sensitivity+=1
    assert -47942+128*374==-70
    assert 128*655-83887==-47
    assert 964689-128*7537==-47
    for M in range(2,17):
        for s in range(1,M):
            for parity in (1,2,128):
                empty_rows=min(2*s,M);full_rows=2*s-empty_rows
                assert F(empty_rows+full_rows,2)==s
                # Empty state: cost=-parity*|T|. Full state: cost=+parity*|T|.
                mean_cost=F(-parity*empty_rows+parity*full_rows,2)
                assert mean_cost==-parity*min(s,M-s)
                mixtures+=1
    # Small exact distinction: hard half-full layers versus dual averages.
    # Two regions of four rows, one four-coordinate update per region.
    route=((0,1,2,3),(2,3,0,1));n=4;M=4;parity=1
    all_records=[];hard={s:None for s in range(1,M)}
    for H1,H2 in product(range(16),repeat=2):
        live=0
        for pair in (3,12):
            if H1&pair:live|=pair
        visible=(H1,live|H2);counts=[0]*M
        for step,mask in enumerate(visible):
            for j in range(n):counts[route[step][j]]+=(mask>>j)&1
        deaths=(H1&~H2).bit_count()
        for T in range(16):
            s=T.bit_count()
            cost=deaths+sum(counts[i]-parity for i in range(M) if T>>i&1)
            all_records.append((cost,H1.bit_count(),H2.bit_count(),s))
            if 0<s<M and H1.bit_count()==2 and H2.bit_count()==2:
                hard[s]=cost if hard[s] is None else min(hard[s],cost)
    dual_checks=0
    for s in range(1,M):
        barrier=-parity*min(s,M-s)
        assert hard[s]>barrier
        for l1,l2,gamma in product((-2,-1,0,1,2),repeat=3):
            dual=min(c+l1*(a-2)+l2*(b-2)+gamma*(t-s) for c,a,b,t in all_records)
            assert dual<=barrier
            dual_checks+=1
    assert time.monotonic()-began<25
    total=math.comb(8192,1974)
    ball=sum(math.comb(1974,j)*math.comb(6218,j) for j in range(188))
    cover_lower_power=total.bit_length()-1-ball.bit_length()
    print(dict(target_sensitivity_checks=sensitivity,mixture_identities=mixtures,
               dual_price_checks=dual_checks,toy_hard_minima=hard,
               equal_size_ball_swaps=187,cover_centers_lower_power_of_two=cover_lower_power,
               seconds=time.monotonic()-began,construction_obstruction=False))


if __name__=='__main__':main()
