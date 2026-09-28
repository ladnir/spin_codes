"""Counterfactual diagnostics of where the three-group transfer bound loses.

The altered transfers are NOT bounds for the implemented encoder. They only
measure sensitivity to terms in the existing envelope. No certificate.
"""
from math import comb,log
import numpy as np
from flint import ctx

from bch_joint_support import authenticated_caps,support_caps
from shortened_bound import dimension_caps
from basis_lattice import improve_caps
from occupancy_model import local_data,epoch_operators
from occupancy_screen import optimize


def placement(operators):
    coefficients=[np.eye(7)]
    for _ in range(64):
        coefficients=[sum((coefficients[r-k]@operators[k])*comb(32,k)
                          for k in range(r+1) if r-k<len(coefficients))
                      for r in range(4)]
    return [value/comb(2048,r) for r,value in enumerate(coefficients)]


def main():
    spectrum=authenticated_caps()
    caps=improve_caps(support_caps(spectrum,g=4,dimensions=dimension_caps()),spectrum)
    data=local_data(3)
    ctx.prec=192
    states=(1<<19)-1
    levels=data[0][0]
    distribution=np.array([levels[w]/states for w in sorted(levels)])
    cases=('baseline','no triple zero feedback','no active cancellation','refresh every active epoch')
    best={(case,u):(0.,None) for case in cases for u in (76,80,83,84,96)}
    for tilt in ('.001','.001125','.00125','.001375','.0015','.0016'):
        exact=epoch_operators(data,tilt)
        operators=[np.array([[float(t[i,j]) for j in range(7)] for i in range(7)]) for t in exact]
        for case in cases:
            variant=[t.copy() for t in operators]
            if case=='no triple zero feedback':
                variant[3][0,0]=0.
            if case=='no active cancellation':
                for t in variant[1:]:
                    t[1:,0]=t[1:,1]/states
            if case=='refresh every active epoch':
                for t in variant[1:]:
                    t[:,2:]+=t[:,1,None]*distribution[None,:]
                    t[:,1]=0.
            region=placement(variant)
            for u in (76,80,83,84,96):
                value,_=optimize(region,(u,)*3,float(tilt))
                if value<best[case,u][0]:
                    best[case,u]=value,tilt
        print('Counterfactual point diagnostic tilt',tilt,flush=True)
    for u in (76,80,83,84,96):
        count=sum(row[u] for row in caps)
        offset=3*log(count)+log(comb(2048,3))
        print('support',u,'log2 point contributions:',
              [(case,round((best[case,u][0]+offset)/log(2),4),best[case,u][1]) for case in cases],flush=True)
    print('Altered transfers are counterfactual, not bounds for the current encoder.')


if __name__=='__main__':
    main()
