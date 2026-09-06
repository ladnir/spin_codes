"""Discovery-only variant with row probabilities tilted against activation loss."""
import argparse
from fractions import Fraction as F
import math
from pathlib import Path

import numpy as np
from scipy.optimize import minimize_scalar
from scipy.special import gammaln
import bridge as base
import larger_state_maps as maps
import k30_kernel as general
import k30_sparse as sparse
import k30_dense_screen as dense
from certify_k30_q1 import ROWS,CUTOFF
from migrate_legacy_workspace import restore_or_verify


def run(output,slope):
    assert not output.exists()
    restore_or_verify(base.ROOT)
    t,s,spectrum,kernel=maps.load('t64_s20');counts=sparse.caps()
    ps=[]
    for band in sparse.BANDS[:-1]:
        def objective(p):
            return max(math.log(counts[w])-math.log(math.comb(256,w))-w*math.log(p)-(256-w)*math.log1p(-p)
                       for w in band)-256*slope*p
        fit=minimize_scalar(objective,bounds=(.001,.999),method='bounded')
        assert fit.success
        ps.append(F.from_float(float(fit.x)))
    ps.append(F(1))
    costs=[math.log(sparse.density_cost(b,p,counts)) for b,p in zip(sparse.BANDS[:-1],ps[:-1])]+[0.]
    hull=dense.upper_hull(list(zip(map(float,ps),costs)))
    grid=sorted(set(np.linspace(min(map(float,ps)),1,97).tolist()+[x for x,y in hull]))
    hull_cost=np.interp(grid,[x for x,y in hull],[y for x,y in hull])
    comb=np.array([math.comb(t,j) for j in range(t+1)],dtype=float)
    rows=[]
    for exponent in range(10,24):
        q=1<<exponent;x=q/ROWS;rates=x*np.array(grid)
        weights=comb[None,:]*rates[:,None]**np.arange(t+1)[None,:]*(1-rates[:,None])**np.arange(t,-1,-1)[None,:]
        best=np.full(len(grid),np.inf);witness=np.zeros(len(grid))
        for lam in np.geomspace(max(1e-7,x/100),min(6,100*x),65):
            epoch=np.array(general.epochs(t,s,spectrum,kernel,math.exp(-float(lam)),t,float)).reshape(t+1,4,4)
            matrices=np.einsum('ij,jkl->ikl',weights,epoch)
            value=np.array([-256*ROWS*lam if r==1 else dense.terminal_log(m) for r,m in zip(rates,matrices)])+lam*CUTOFF
            mask=value<best;best[mask]=value[mask];witness[mask]=lam
        combinatorial=gammaln(ROWS+1)-gammaln(q+1)-gammaln(ROWS-q+1)+q*math.log(len(ps))+256*math.log(ROWS+1)
        margins=-(combinatorial+q*hull_cost+best)/math.log(2)
        index=int(np.argmin(margins))
        row=dict(occupation=q,margin_bits_diagnostic=float(margins[index]),
            worst_sample_mean_active_density=grid[index],tilt_witness=float(witness[index]))
        rows.append(row);print(row,flush=True)
    base.write_new(output,dict(status='DISCOVERY_ONLY_TILTED_IID_DOMINATION_GRID',full_distance_proved=False,
        slope=slope,iid_region_density_factor=ROWS+1,row_probabilities=[base.encode(p) for p in ps],
        upper_hull=hull,rows=rows,warning='No coverage between sampled occupancies or type proportions. Binary64 only.',
        source_sha256={p.relative_to(base.ROOT).as_posix():base.sha(p) for p in
            (Path(__file__),Path(dense.__file__),Path(general.__file__),Path(sparse.__file__))}))


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--slope',type=float,default=1.)
    args=parser.parse_args();run(args.output,args.slope)
