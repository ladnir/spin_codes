"""Diagnostic common potentials for avoiding exponential support-box grids.

If R_r h <= a b^r h and h >= eta*terminal, then the full support union
at occupancy q is <= binom(2048,q) exp(lambda*T) a^256 S(b)^q / eta,
where S(b) is the outer support-CDF bound evaluated at decreasing b^u.
This script only searches binary64 witnesses; it does not certify them.
"""
import argparse
from math import comb,log
import numpy as np
from scipy.optimize import linprog
from scipy.special import logsumexp
from flint import ctx

from bch_joint_support import authenticated_caps,support_caps
from shortened_bound import dimension_caps
from basis_lattice import improve_caps
from occupancy_model import local_data
from occupancy_memory import prepare,build,TERMINAL
from occupancy_cdf import differences


def potential(region,a,b):
    n=len(TERMINAL)
    rows=[]
    for r,operator in enumerate(region):
        row=np.column_stack((operator-a*b**r*np.eye(n),np.zeros(n)))
        scale=np.max(np.abs(row),axis=1)
        nonzero=scale>0
        rows.extend(row[nonzero]/scale[nonzero,None])
    for i,entry in enumerate(TERMINAL):
        if entry:
            row=np.zeros(n+1)
            row[i]=-1
            row[-1]=1
            rows.append(row)
    objective=np.zeros(n+1)
    objective[-1]=-1
    fitted=linprog(objective,A_ub=rows,b_ub=np.zeros(len(rows)),
                   bounds=[(1,1)]+[(0,None)]*(n-1)+[(0,1)],method='highs')
    if not fitted.success or fitted.x[-1]<=0:
        return None
    h,eta=fitted.x[:-1],fitted.x[-1]
    # Reject large solver residuals. This is not a substitute for outward replay.
    if any(np.any(t@h>a*b**r*h+1e-8*np.maximum(1,t@h)) for r,t in enumerate(region)):
        return None
    return eta,h


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--groups',type=int,choices=(3,4),default=4)
    args=parser.parse_args()
    spectrum=authenticated_caps()
    caps=improve_caps(support_caps(spectrum,g=4,dimensions=dimension_caps()),spectrum)
    increments=differences([sum(row[u] for row in caps) for u in range(257)])
    prepared=prepare(local_data(args.groups))
    ctx.prec=192
    best=None
    for tilt in ('.00064','.001','.001375','.0016','.002','.0025'):
        _,region=build(prepared,tilt)
        for a in (1.,1.001,1.01,1.05,1.1,1.2):
            for b in (.01,.015,.02,.025,.03,.04,.05,.06,.075,.1,.15,.2,.3,.5,.75):
                fit=potential(region,a,b)
                if fit is None:
                    continue
                eta,h=fit
                outer=logsumexp([log(c)+u*log(b) for u,c in enumerate(increments) if c])
                score=(log(comb(2048,args.groups))+float(tilt)*209715+256*log(a)
                       +args.groups*outer-log(eta))/log(2)
                if best is None or score<best[0]:
                    best=score,tilt,a,b,eta,h
        print('Potential screen tilt',tilt,'best log2',None if best is None else best[0],flush=True)
    print('Best binary64 candidate:',best,flush=True)
    print('No certificate: common-potential inequalities still need outward verification.')


if __name__=='__main__':
    main()
