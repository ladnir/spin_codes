"""Discovery-only dense comparison using the existing three-state envelope.

Keeping the unrestricted nonzero-state emission moment avoids paying the
separate kernel/nonkernel upper bounds in k30_kernel. The old transfer's
refresh loss is small relative to dense-occupancy exponents.
"""
import argparse
from fractions import Fraction as F
import math
from pathlib import Path

import numpy as np
from scipy.optimize import minimize_scalar
from scipy.special import gammaln
import bridge as base
import larger_state_maps as maps
import tightened_occupancy as transfer
import k30_sparse as sparse
import k30_dense_screen as dense
from certify_k30_q1 import ROWS,CUTOFF
from migrate_legacy_workspace import restore_or_verify


def row_witnesses(counts,slope):
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
    return ps,dense.upper_hull(list(zip(map(float,ps),costs)))


def run(output):
    assert not output.exists()
    restore_or_verify(base.ROOT)
    t,s,spectrum,kernel=maps.load('t64_s20');counts=sparse.caps()
    witnesses=[]
    for slope in (0.,.25,.5,.75,1.,1.25):
        ps,hull=row_witnesses(counts,slope)
        grid=sorted(set(np.linspace(min(map(float,ps)),1,97).tolist()+[x for x,y in hull]))
        cost=np.interp(grid,[x for x,y in hull],[y for x,y in hull])
        witnesses.append((slope,ps,hull,np.array(grid),cost))
    comb=np.array([math.comb(t,j) for j in range(t+1)],dtype=float)
    rows=[]
    for exponent in range(10,24):
        q=1<<exponent;x=q/ROWS
        tilts=np.geomspace(max(1e-7,x/100),min(6,100*x),65)
        epochs=[np.array(transfer.epoch_matrices(t,s,spectrum,kernel,math.exp(-float(lam)),float,t)).reshape(t+1,3,3)
                for lam in tilts]
        combinatorial=gammaln(ROWS+1)-gammaln(q+1)-gammaln(ROWS-q+1)+q*math.log(len(sparse.BANDS))+256*math.log(ROWS+1)
        selected=None;diagnostics=[]
        for slope,ps,hull,grid,hull_cost in witnesses:
            rates=x*grid
            weights=comb[None,:]*rates[:,None]**np.arange(t+1)[None,:]*(1-rates[:,None])**np.arange(t,-1,-1)[None,:]
            best=np.full(len(grid),np.inf);tilt_witness=np.zeros(len(grid))
            for lam,epoch in zip(tilts,epochs):
                matrices=np.einsum('ij,jkl->ikl',weights,epoch)
                value=np.array([-256*ROWS*lam if r==1 else dense.terminal_log(m) for r,m in zip(rates,matrices)])+lam*CUTOFF
                mask=value<best;best[mask]=value[mask];tilt_witness[mask]=lam
            margins=-(combinatorial+q*hull_cost+best)/math.log(2)
            assert np.isfinite(margins).all()
            index=int(np.argmin(margins));margin=float(margins[index])
            diagnostics.append(dict(slope=slope,margin_bits_diagnostic=margin))
            if selected is None or margin>selected['margin_bits_diagnostic']:
                selected=dict(occupation=q,margin_bits_diagnostic=margin,slope=slope,
                    worst_sample_mean_active_density=float(grid[index]),tilt_witness=float(tilt_witness[index]))
        selected['slope_comparisons']=diagnostics
        rows.append(selected);print({k:v for k,v in selected.items() if k!='slope_comparisons'},flush=True)
    base.write_new(output,dict(status='DISCOVERY_ONLY_UNRESTRICTED_IID_GRID',full_distance_proved=False,
        iid_region_density_factor=ROWS+1,rows=rows,
        witnesses=[dict(slope=c,row_probabilities=[base.encode(p) for p in ps],upper_hull=hull)
                   for c,ps,hull,grid,cost in witnesses],
        warning='No coverage between sampled occupancies or type proportions. Binary64 only.',
        source_sha256={p.relative_to(base.ROOT).as_posix():base.sha(p) for p in
            (Path(__file__),Path(dense.__file__),Path(transfer.__file__),Path(transfer.original.__file__),Path(sparse.__file__))}))


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,required=True)
    run(parser.parse_args().output)
