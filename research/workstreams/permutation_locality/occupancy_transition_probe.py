"""Counterfactual sensitivity of selected support points; NOT valid upper bounds.

Artificially halving transition coefficients isolates useful proof targets.
No modified coefficient is asserted to bound the actual encoder.
"""
import argparse
from math import comb,log
import numpy as np
from flint import arb,arb_mat,ctx
from bch_joint_support import authenticated_caps,support_caps
from shortened_bound import dimension_caps
from basis_lattice import improve_caps
from occupancy_allones import weighted_cdf
from occupancy_model import local_data,placement
from occupancy_memory import prepare,rounded,TERMINAL,Z,F,M,C,U
from occupancy_fresh_moment import fresh_census,refine
from occupancy_window_average import prepare_inputs,averages,refine as window_refine
from occupancy_multi_average import refine as multi_refine
from occupancy_screen import optimize


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--groups',type=int,default=13)
    parser.add_argument('--support',type=int,default=152)
    parser.add_argument('--tilt',default='.008')
    parser.add_argument('--penalty',default='.625')
    args=parser.parse_args()
    spectrum=authenticated_caps();dimensions=dimension_caps()
    caps=improve_caps(support_caps(spectrum,g=4,dimensions=dimensions),spectrum)
    counts=weighted_cdf(spectrum,caps,dimensions,args.penalty,prefix_flags=True)
    prepared=prepare(local_data(4));fresh=fresh_census(prepared)
    inputs=prepare_inputs();ctx.prec=192
    values=averages(inputs,args.tilt)
    ops=refine(prepared,fresh,args.tilt,args.groups,args.penalty)
    ops=window_refine(prepared,values,args.tilt,args.groups,args.penalty,ops)
    ops=multi_refine(ops,fresh,args.tilt,args.penalty)
    tests={'baseline':[], 'density cancellation':[(C,Z)],
           'density persistence':[(C,C)], 'fresh cancellation':[(F,Z)],
           'fresh to density':[(F,C)], 'zero return':[(Z,Z)],
           'mature lazy mass':[(M,M)],
           'uniform lazy density':[(U+i,C) for i in range(5)]}
    for name,entries in tests.items():
        altered=[arb_mat(t) for t in ops]
        for j,t in enumerate(altered):
            if not j:
                continue
            for i,k in entries:
                t[i,k]*=arb('.5')
        regions=placement(altered,rounding=rounded)
        arrays=[np.array([[float(t[i,j]) for j in range(9)] for i in range(9)]) for t in regions]
        value,ps=optimize(arrays,(args.support,)*args.groups,float(args.tilt),TERMINAL)
        score=(value+args.groups*log(counts[args.support])+log(comb(2048,args.groups)))/log(2)
        print('COUNTERFACTUAL except baseline:',name,'log2 score',score,flush=True)
    print('Sensitivity only. Halved entries are not proved and cannot be used for a certificate.')


if __name__=='__main__':
    main()
