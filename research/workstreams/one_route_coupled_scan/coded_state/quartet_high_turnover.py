"""Exact all-row-count union for covered inactive positions plus deaths.

Replays the scalar union from the established pin-and-shortening count.
No floating point, generated files, external modules, or whole-code claim.
"""
from fractions import Fraction as F
from itertools import product
from pathlib import Path
import shutil


def certificate():
    q=1<<32;N=1<<21;M=8192;W=(9*N)//20
    assert 1/(1-(N//4)*F(280,q))==F(1024,989)<2
    # Coordinate support/selected-zero automaton: [[1,1],[1,2]].
    potential=(F(1),F(13,8));rate=F(21,8)
    assert potential[0]+potential[1]<=rate*potential[0]
    assert potential[0]+2*potential[1]<=rate*potential[1]
    assert min(potential)>=1
    # Independent short enumeration of the support/selected-zero automaton.
    vector=(1,0)
    for steps in range(1,9):
        vector=(vector[0]+vector[1],vector[0]+2*vector[1])
        exact=0
        for suffix in product((0,1),repeat=steps):
            path=(0,)+suffix
            exact+=2**sum(u*v for u,v in zip(path,path[1:]))
        assert exact==sum(vector) and exact<=rate**steps
    # All 512 coordinate histories, initially zero, with free terminal:
    # count <= rate^(512*4096) <= 2^history_bits.
    assert N%1024==0
    assert 21**1024 < 2**(1426+3*1024)
    history_bits=1426*(N//1024)
    slack=(history_bits+M+64)//32
    assert 32*slack==history_bits+M+64
    log2_union=1+M+history_bits-32*slack
    assert log2_union==-63
    all_deaths_rows=(W+slack+127)//128
    threshold=lambda n:max(0,W-128*n+slack)
    assert threshold(all_deaths_rows)==0 and threshold(all_deaths_rows-1)>0
    free_gap=256*2100-(threshold(6092)-1)
    assert free_gap==282137 and N-free_gap==1815015
    print(dict(certified_failure_less_than='2^-63',all_row_counts=True,
               field_symbol_slack=slack,death_threshold_at_inactive_2100=threshold(6092),
               bounded_score='covered_inactive_coordinates + support_deaths',
               automaton_enumeration_checks=8,free_minus_deaths_floor_at_inactive_2100=free_gap,
               death_threshold_is_a_sufficient_corollary=True,
               every_low_output_history_covered_when_active_rows_at_least=all_deaths_rows,
               free_terminal=True,whole_code_distance=False,
               arithmetic='exact rational and integer',files_written=0))
    return slack


if __name__=='__main__':
    assert shutil.disk_usage(Path(__file__).resolve().anchor).free>=20*(1<<30)
    certificate()
