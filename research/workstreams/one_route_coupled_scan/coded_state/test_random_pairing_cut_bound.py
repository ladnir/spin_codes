"""Exhaustive small checks of the cut formulation and half-full lower bound."""
import shutil
import tempfile
import time
import numpy as np
from random_pairing_cut_bound import cut_bound


def main():
    assert shutil.disk_usage(tempfile.gettempdir()).free>=20*(1<<30),'disk guard'
    began=time.monotonic();rng=np.random.default_rng(20261027);checks=0
    n=4
    for L in (1,2,3):
        for fixture in range(3):
            target=rng.random((L,n))<.4
            partners=np.empty((L,n),dtype=np.int32)
            for t in range(L):
                order=rng.permutation(n);partners[t,order]=order[np.arange(n)^1]
            candidates=[];symmetric_candidates=[];half_full=[]
            for mask in range(1<<(L*n)):
                states=[(mask>>(n*t))&15 for t in range(L)]
                value=sum(bool(states[0]>>i&1) and target[0,i] for i in range(n))
                for t in range(1,L):
                    P=states[t-1];H=states[t]
                    live=P|sum(((P>>int(partners[t,i]))&1)<<i for i in range(n))
                    value+=sum(bool((live|H)>>i&1) and target[t,i] for i in range(n))
                    value+=(P&~H).bit_count()
                candidates.append((int(value),mask.bit_count()))
                symmetric_candidates.append((16*int(value)+8*(states[-1].bit_count()-states[0].bit_count()),mask.bit_count()))
                if all(H.bit_count()==n//2 for H in states):half_full.append(int(value))
            for mu in (0,3,8,16):
                result=cut_bound(target,partners,mu)
                exact=min(16*c-mu*a for c,a in candidates)+mu*L*n
                assert result['flow']==exact
                assert result['lower_scaled']<=16*min(half_full)
                checks+=1
                symmetric=cut_bound(target,partners,mu,symmetric=True)
                exact=min(c-mu*a for c,a in symmetric_candidates)+mu*L*n
                assert symmetric['flow']==exact
                assert symmetric['lower_scaled']<=16*min(half_full)
                checks+=1
    assert time.monotonic()-began<25
    print(dict(exhaustive_cut_checks=checks,seconds=time.monotonic()-began))


if __name__=='__main__':main()
