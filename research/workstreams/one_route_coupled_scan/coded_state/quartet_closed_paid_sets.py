"""Exact toy check: charge only paid inputs covered by the canonical support.

Width4, four updates, two regions. All12544 (inactive,paid) pairs and all
65536 allowed masks. No code-distance claim, source only, no generated files.
"""
from collections import defaultdict
import itertools
import math


def canonical(allowed):
    masks=[(allowed>>(4*t))&15 for t in range(4)]
    free=[masks[t] if masks[t+1]==15 else 0 for t in range(3)]+[masks[-1]]
    selected=[0]*4
    for j in range(4):
        t=0
        while t<4:
            if not (free[t]>>j)&1:t+=1;continue
            end=t
            while end<4 and (free[end]>>j)&1:end+=1
            if end==4 or 2*(end-t)-5>0:
                for u in range(t,end):selected[u]|=1<<j
            t=end
    cover=[selected[t]|(15 if t and selected[t-1] else 0) for t in range(4)]
    assert all((c&~a)==0 for c,a in zip(cover,masks))
    score=2*sum(s.bit_count() for s in selected)
    score-=5*sum((selected[t]&~selected[t+1]).bit_count() for t in range(3))
    return tuple(selected),sum(c<<(4*t) for t,c in enumerate(cover)),score


def tests():
    full=(1<<16)-1;direct=defaultdict(int);closure_checks=0;closed=0
    inactive_regions=[sum(1<<j for j in positions) for positions in itertools.combinations(range(8),2)]
    for low,high in itertools.product(inactive_regions,repeat=2):
        I=low|(high<<8);bits=[j for j in range(16) if (I>>j)&1]
        for choice in range(16):
            paid=sum(1<<j for k,j in enumerate(bits) if (choice>>k)&1)
            U=(full^I)|paid;support,cover,H=canonical(U);reduced=I&cover
            support2,cover2,H2=canonical((full^I)|reduced)
            assert support2==support and H2==H and I&cover2==reduced
            assert reduced&~paid==0;closure_checks+=1
            if paid==reduced:
                direct[paid.bit_count(),H]+=1;closed+=1
    grouped=defaultdict(int)
    for U in range(1<<16):
        support,cover,H=canonical(U);multiplicity=1;total=0
        for shift in (0,8):
            a=((U>>shift)&255).bit_count();c=a-6;v=((cover>>shift)&255).bit_count()
            if not 0<=c<=2 or c>v:multiplicity=0;break
            multiplicity*=math.comb(v,c);total+=c
        if multiplicity:grouped[total,H]+=multiplicity
    assert dict(grouped)==dict(direct)
    print(dict(exact_closure_pairs=closure_checks,closed_pairs=closed,allowed_masks=1<<16,
               exact_regional_closed_count_coefficients=len(grouped),local_only=True),flush=True)
    return dict(grouped)


if __name__=='__main__':tests()
