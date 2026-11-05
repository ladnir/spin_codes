"""Exact small checks for the confined-state projection argument; <=30s."""
from collections import Counter
from fractions import Fraction as F
from itertools import product,combinations
import math
import time
from diagnose_uniform_pairs import field
from diagnose_random_pairing import crossing_law,all_matchings,matchings


def tests():
    began=time.monotonic();q=8;mul=field(q,11);rs_cases=0
    words=[]
    for coefficients in product(range(q),repeat=3):
        word=[]
        for x in range(6):
            value=0
            for c in reversed(coefficients):value=mul[value][x]^c
            word.append(value)
        words.append(word)
    for m in range(7):
        for outside in combinations(range(6),m):
            exact=sum((F(1,q**sum(word[i]!=0 for i in outside)) for word in words),F(0))
            upper=1+2**m*F(q)**(3-m)
            assert exact<=upper
            if m<3:assert exact==F(q)**(3-m)*(2-F(1,q))**m
            rs_cases+=1
    subsets=list(combinations(range(4),2));na_cases=0
    for cutoff in range(4):
        for tilt in (F(5,4),F(3,2),F(2)):
            exact=F(0)
            for regions in product(subsets,repeat=3):
                counts=[sum(i in region for region in regions) for i in range(4)]
                exact+=tilt**sum(max(0,b-cutoff) for b in counts)
            exact/=len(subsets)**3
            marginal=sum((F(math.comb(3,b),8)*tilt**max(0,b-cutoff) for b in range(4)),F(0))
            assert exact<=marginal**4;na_cases+=1
    delta=2
    good=[(a,b,c,d) for a,b,c,d in product(range(1,q),repeat=4)
          if mul[a][d]^mul[b][c]==1 and a!=delta and d!=delta
          and mul[a^delta][d^delta]^mul[b][c]!=0]
    # At most q matrices map a fixed nonzero input to a specified pair.
    density=F(q**3,len(good));pair_cases=0;moving_cases=0
    fixtures=list(product(range(2),repeat=4))+[(2,)*4,(3,)*4]
    for w in range(1,5):
        H=(1<<w)-1;state=tuple(int(i<w) for i in range(4));law=crossing_law(4,w)
        counters=[]
        for matching in all_matchings(tuple(range(4))):
            pair_images=[]
            for i,j in matching:
                if not (state[i] or state[j]):pair_images.append(Counter({(0,0):len(good)}))
                else:
                    pair_images.append(Counter((mul[a][state[i]]^mul[b][state[j]],
                                                mul[c][state[i]]^mul[d][state[j]]) for a,b,c,d in good))
            census=Counter()
            for images in product(*(list(c) for c in pair_images)):
                image=[0]*4;weight=1
                for pair,values,counter in zip(matching,images,pair_images):
                    image[pair[0]],image[pair[1]]=values;weight*=counter[values]
                census[tuple(image)]+=weight
            counters.append(census)
        for X in fixtures:
            outside=sum(X[i]!=0 for i in range(w,4))
            endpoint=law.get(outside,F(0))/math.comb(4-w,outside)
            zero_counts=Counter();joint=Counter()
            for census in counters:
                for image,multiplicity in census.items():
                    nxt=tuple(x^a for x,a in zip(X,image))
                    next_mask=sum((v!=0)<<i for i,v in enumerate(nxt))
                    zeros=sum((nxt[i]==mul[delta][state[i]])<<i for i in range(w))
                    joint[next_mask,zeros]+=multiplicity
                    if next_mask==H:zero_counts[zeros]+=multiplicity
            for selected in range(1<<w):
                count=sum(c for zeros,c in zero_counts.items() if zeros&selected==selected)
                exact=F(count,matchings(4)*len(good)**2)
                upper=endpoint*density**2/F(q**(outside+selected.bit_count()))
                assert exact<=upper,(w,X,selected,exact,upper)
                pair_cases+=1
            for next_mask in range(16):
                births=(next_mask&~H).bit_count();deaths=(H&~next_mask).bit_count()
                common=H&next_mask
                outside_union=sum(X[i]!=0 and not ((H|next_mask)>>i&1) for i in range(4))
                endpoint=sum((math.comb(births,j)*law.get(outside_union+j,F(0))
                              /math.comb(4-w,outside_union+j) for j in range(births+1)
                              if outside_union+j<=4-w),F(0))
                worst=max(prob/math.comb(4-w,k) for k,prob in law.items())
                assert endpoint<=min(F(1),2**births*worst)
                for selected in range(16):
                    if selected&~common:continue
                    count=sum(c for (mask,zeros),c in joint.items()
                              if mask==next_mask and zeros&selected==selected)
                    exact=F(count,matchings(4)*len(good)**2)
                    upper=endpoint*density**2/F(q**(outside_union+deaths+selected.bit_count()))
                    assert exact<=upper,(w,X,next_mask,selected,exact,upper)
                    moving_cases+=1
            if time.monotonic()-began>25:raise TimeoutError('projection test budget')
    print(dict(rs_projection_cases=rs_cases,regional_mgf_cases=na_cases,
               accepted_pair_constraints=pair_cases,moving_support_constraints=moving_cases,
               seconds=time.monotonic()-began))


if __name__=='__main__':tests()
