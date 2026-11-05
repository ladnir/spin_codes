"""Exhaustive all-pair-type temporal tests; no distance/performance claim."""
from itertools import product
from fractions import Fraction as F
import time
from layered_fixed_rank import path_nullity
from layered_component_transfer import STATES,coordinate_step
from layered_general_components import ALL_TERMS,ALL_LAST,rational_operator


def run():
    started=time.monotonic();width=4;cases=0;by_type=[0]*4;total=F(0)
    for s1,s2 in product(range(16),repeat=2):
        full=[];fixed_base={j:0 for j in range(width)};counts=[0,0,0]
        for t,mask in enumerate((s1,s2),1):
            bit=t%2
            for j in range(width):
                if j>>bit&1:continue
                k=j^(1<<bit);a=(mask>>j)&1;c=(mask>>k)&1
                by_type[2*a+c]+=1
                if a and c:full.append((t,j,k));counts[2]+=1
                elif a:fixed_base[t*width+j]=1
                elif c:counts[1]+=1
                else:fixed_base[t*width+j]=fixed_base[t*width+k]=0;counts[0]+=1
        inactive={t*width+j for t,mask in enumerate((s1,s2)) for j in range(width) if not mask>>j&1}
        common=[j for j in range(width) if (s1&s2)>>j&1]
        births=(s2&~s1).bit_count()
        forced_internal=0
        for j in range(width):
            if not j>>1&1 and not (s1>>j&1) and not (s1>>(j^2)&1):
                forced_internal+=int(not (s2>>j&1))+int(not (s2>>(j^2)&1))
        assert forced_internal>=max(0,2*counts[0]-births-width)
        for emask in range(1<<len(full)):
            fixed=fixed_base.copy()
            for k,(t,i,j) in enumerate(full):
                if emask>>k&1:fixed[t*width+i]=0;fixed[t*width+j]=1
            for zmask in range(1<<len(common)):
                ze={width+j for k,j in enumerate(common) if zmask>>k&1}
                kappa=path_nullity(width,3,inactive,ze,fixed);scanned=0
                for j in range(width):
                    src=STATES.index(((s1>>j)&1,0));sn=(s2>>j)&1
                    zero=(width+j in ze) or not ((s1|s2)>>j&1)
                    items=[v for v in coordinate_step(src,fixed.get(width+j,-1))
                           if STATES[v[0]][0]==sn and v[3]==int(not zero)]
                    assert len(items)==1
                    dst,gain,_,_=items[0]
                    scanned+=gain+STATES[dst][1]+(2*width+j not in fixed)
                assert scanned==kappa
                weight=s1.bit_count()+(s1|s2).bit_count()-len(ze)
                # q=2^32 and theta=31/32; nontrivial markers for all types.
                total+=F(2)**kappa*F(2)**(-31*len(ze))*F(2,3)**counts[0]*F(3,4)**counts[1]*F(2)**counts[2]*F(1,3)**emask.bit_count()*F(2,3)**births*F(1,2)**weight
                cases+=1
    toy=[F(2,3),F(3,4),F(2),F(1,3),F(2,3),F(1,2),F(2),F(3),F(3),F(4)]
    # Positive z/r select all source rows; numerical profile counts are unused
    # by rational_operator, which constructs only the local positive operator.
    h,rows,ends=rational_operator((1,1,1,1,0,0,0,0),toy,31,depth=1)
    initial=((1+h[0][1]*toy[5])*(1+h[1][1]*toy[5]))**2
    bound=initial*max(rows[0].values())*max(rows[1].values())*max(ends[0].values())**2
    assert total<=bound
    # Every row/branch has the required deterministic fixed/free coordinates.
    assert len({v[0] for v in ALL_LAST})==9
    assert all(v[3]==0 or v[2]==0 for v in ALL_LAST)  # not both00 and01
    shortening_checks=0
    for R in range(129,257):
        for g in range(R+1):
            assert R*max(0,g-128)<=(R-128)*g
            shortening_checks+=1
    assert time.monotonic()-started<20
    return dict(exhaustive_histories=cases,pair_type_occurrences=by_type,
                local_terms=len(ALL_TERMS),terminal_terms=len(ALL_LAST),
                shortening_checks=shortening_checks,
                exact_component_checks=True,exact_partition_bound=True)


if __name__=='__main__':print(run())
