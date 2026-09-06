"""Discovery only: iid domination of shuffled Poisson-binomial region inputs.

For L independent Bernoulli bits with mean rL, each count probability is
at most (L+1) times Binomial(L,r). Shuffling makes this a pointwise density
bound on the full region. Pay this factor for each of the 256 regions.
The upper concave hull of row density costs covers all type proportions.
Sampling that hull and occupancy is NOT interval certification.
"""
import argparse
from fractions import Fraction as F
import math
from pathlib import Path

import numpy as np
from scipy.special import gammaln
import bridge as base
import larger_state_maps as maps
import k30_kernel as general
import k30_sparse as sparse
from certify_k30_q1 import ROWS,CUTOFF
from migrate_legacy_workspace import restore_or_verify


def upper_hull(points):
    hull=[]
    for point in sorted(points):
        while len(hull)>=2:
            a,b=hull[-2:]
            if (b[1]-a[1])*(point[0]-b[0])>(point[1]-b[1])*(b[0]-a[0]):break
            hull.pop()
        hull.append(point)
    return hull


def terminal_log(matrix):
    matrix=matrix.copy();scale=0.
    for _ in range(25):
        matrix=matrix@matrix
        factor=float(matrix.max())
        if factor<=0:return float('nan')
        matrix/=factor;scale=2*scale+math.log(factor)
    value=float(matrix[0].sum())
    return math.log(value)+scale if value>0 else float('nan')


def run(output):
    assert not output.exists()
    restore_or_verify(base.ROOT)
    t,s,spectrum,kernel=maps.load('t64_s20');counts=sparse.caps()
    ps=[F(min(b)+max(b),512) for b in sparse.BANDS[:-1]]+[F(1)]
    costs=[math.log(sparse.density_cost(b,p,counts)) for b,p in zip(sparse.BANDS[:-1],ps[:-1])]+[0.]
    hull=upper_hull(list(zip(map(float,ps),costs)))
    grid=sorted(set(np.linspace(float(ps[0]),1,97).tolist()+[x for x,y in hull]))
    hull_cost=np.interp(grid,[x for x,y in hull],[y for x,y in hull])
    comb=np.array([math.comb(t,j) for j in range(t+1)],dtype=float)
    rows=[]
    for exponent in range(10,24):
        q=1<<exponent;x=q/ROWS
        rates=x*np.array(grid)
        weights=comb[None,:]*rates[:,None]**np.arange(t+1)[None,:]*(1-rates[:,None])**np.arange(t,-1,-1)[None,:]
        best=np.full(len(grid),np.inf);witness=np.zeros(len(grid))
        for lam in np.geomspace(max(1e-7,x/100),min(6,100*x),37):
            epoch=np.array(general.epochs(t,s,spectrum,kernel,math.exp(-float(lam)),t,float)).reshape(t+1,4,4)
            matrices=np.einsum('ij,jkl->ikl',weights,epoch)
            value=np.array([terminal_log(m) for m in matrices])+lam*CUTOFF
            mask=value<best;best[mask]=value[mask];witness[mask]=lam
        combinatorial=gammaln(ROWS+1)-gammaln(q+1)-gammaln(ROWS-q+1)+q*math.log(len(ps))+256*math.log(ROWS+1)
        margins=-(combinatorial+q*hull_cost+best)/math.log(2)
        index=int(np.argmin(margins))
        row=dict(occupation=q,margin_bits_diagnostic=float(margins[index]),
            worst_sample_mean_active_density=grid[index],tilt_witness=float(witness[index]))
        rows.append(row);print(row,flush=True)
    base.write_new(output,dict(status='DISCOVERY_ONLY_IID_DOMINATION_GRID',full_distance_proved=False,
        iid_region_density_factor=ROWS+1,row_probabilities=[base.encode(p) for p in ps],
        upper_hull=hull,rows=rows,warning='No coverage between sampled occupancies or type proportions. Binary64 only.',
        source_sha256={p.relative_to(base.ROOT).as_posix():base.sha(p) for p in
            (Path(__file__),Path(general.__file__),Path(sparse.__file__))}))


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,required=True)
    run(parser.parse_args().output)
