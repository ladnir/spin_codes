"""Exact two-orbit feedback collisions, not an output-distance certificate."""
from fractions import Fraction
from math import log2
import numpy as np
from scipy.sparse import csr_matrix
from rank_two_moment import census


def main():
    _,allowed,_,atoms,_=census()
    shapes=sorted(allowed)
    rows=[];columns=[];counts=[]
    for i,s in enumerate(shapes):
        for syndrome,count in atoms[s].items():
            if syndrome:
                rows.append(i);columns.append(syndrome);counts.append(count)
                assert syndrome.bit_count()%2==sum(s)%2
    frequencies=csr_matrix((np.array(counts,dtype=np.int64),(rows,columns)),shape=(len(shapes),1<<19))
    pairs=(frequencies@frequencies.T).toarray()
    assert np.array_equal(pairs,pairs.T)
    for i,s in enumerate(shapes):
        assert pairs[i,i]==sum(c*c for q,c in atoms[s].items() if q)
    zero=same_parity_zero=0
    for i,a in enumerate(shapes):
        for j,b in enumerate(shapes):
            if not pairs[i,j]:
                zero+=1
                same_parity_zero+=(sum(a)-sum(b))%2==0
            if (sum(a)-sum(b))%2:
                assert pairs[i,j]==0
    print('Ordered orbit pairs:',len(shapes)**2,'zero collision:',zero,
          'same-parity zero collision:',same_parity_zero,flush=True)
    selected=((1,1,1,1),(2,2,2,2),(3,3,3,3),(4,4,4,4),(1,2,2,3),(0,2,3,3))
    print('source,target,exact_probability,coarse_probability,improvement_bits',flush=True)
    for a in selected:
        i=shapes.index(a)
        for b in selected:
            j=shapes.index(b)
            na=8*len(allowed[a]);nb=8*len(allowed[b])
            exact=Fraction(int(pairs[i,j]),na*nb)
            coarse=Fraction((na-atoms[a][0])*max(c for q,c in atoms[b].items() if q),na*nb)
            assert exact<=coarse
            print(a,b,str(exact),str(coarse),log2(coarse/exact) if exact else 'infinite',flush=True)
    print('Exact identities for independent local orbit draws only. A recurrence must also account for refreshes, accumulated states, and output weights.')


if __name__=='__main__':
    main()
