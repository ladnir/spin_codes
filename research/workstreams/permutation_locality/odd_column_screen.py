"""Selected-point odd-column weighting experiment. Binary64, not a certificate."""
import argparse
from math import comb,log
import numpy as np
from flint import ctx,arb
from occupancy_sensitivity import (prepare,local_data,fresh_census,prepare_inputs,census,pair_census,
    zero_census,averages,fresh_refine,window_refine,multi_refine,collision_refine,
    zero_refine,lift,transform,authenticated_caps,dimension_caps,support_caps,improve_caps)
from shortening_moments import improve
from odd_column_caps import weighted_support,self_test
from occupancy_memory import epoch_operators
from mature_tail import self_test as inner_test
from round_sensitivity import score


def array(ops):
    return np.array([[[float(t[i,j]) for j in range(t.ncols())] for i in range(t.nrows())] for t in ops])


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--groups',type=int,default=64)
    parser.add_argument('--supports',type=int,nargs='+',default=[96,128,160])
    parser.add_argument('--tilt',default='.032')
    parser.add_argument('--penalty',default='.75')
    parser.add_argument('--etas',nargs='+',default=['1','1.1','1.25','1.5'])
    parser.add_argument('--full-feedback',type=int,choices=range(2,11),default=6)
    args=parser.parse_args();ctx.prec=192;self_test()
    spectrum=authenticated_caps();dimensions=dimension_caps()
    caps,_=improve(improve_caps(support_caps(spectrum,g=4,dimensions=dimensions),spectrum),dimensions)
    counts={(u,eta):weighted_support(spectrum,caps,dimensions,u,eta,args.penalty)
            for u in args.supports for eta in args.etas}
    prepared=prepare(local_data(4));fresh=fresh_census(prepared)
    inputs=prepare_inputs();tails=census(inputs,prepared)
    pairs=pair_census(inputs,prepared);zeros=zero_census()
    values=averages(inputs,args.tilt)
    from full_feedback_census import census as full_census
    from full_feedback_refinement import refine as full_refine
    full=full_census(args.full_feedback)
    unweighted=epoch_operators(prepared,args.tilt,detailed=True,full_penalty=args.penalty)
    for eta in args.etas:
        # Every detailed shape must acquire exactly its odd-column factor.
        detailed=epoch_operators(prepared,args.tilt,detailed=True,full_penalty=args.penalty,odd_penalty=eta)
        for j in range(1,5):
            for shape,t in detailed[j].items():
                expected=unweighted[j][shape]*(arb(eta)**sum(w%2 for w in shape))
                assert np.allclose(array([t]),array([expected]),rtol=1e-14,atol=0)
        kw={'odd_penalty':eta};local=min(args.groups,32)
        ops=fresh_refine(prepared,fresh,args.tilt,local,args.penalty,**kw)
        ops=window_refine(prepared,values,args.tilt,local,args.penalty,ops,**kw)
        ops=multi_refine(ops,fresh,args.tilt,args.penalty,**kw)
        ops=collision_refine(ops,prepared,args.tilt,args.penalty,**kw)
        ops=zero_refine(ops,zeros,args.tilt,args.penalty,**kw)
        ops=lift(ops,inputs[3],tails,args.tilt,args.penalty,64,pairs,**kw)
        ops=transform(ops,inputs[3],args.tilt,2)
        if eta=='1' and args.groups==64 and args.tilt=='.032' and args.penalty=='.75':
            with np.load('tmp/r2-q64-sensitivity.npz',allow_pickle=False) as saved:
                assert np.allclose(array(ops),saved['base'],rtol=1e-13,atol=0)
            print('Default local operators reproduce pre-odd-weight diagnostic cache',flush=True)
        ops=full_refine(ops,full,inputs[3],args.tilt,args.penalty,2,odd_penalty=eta)
        inner_test(inputs,ops,args.tilt,args.penalty,rounds=2,odd_penalty=eta)
        for u in args.supports:
            outer=args.groups*log(counts[u,eta])+log(comb(2048,args.groups))
            value,p=score(array(ops),args.groups,u,args.tilt,outer)
            print('ODD-WEIGHT DIAGNOSTIC eta/support/log2-bound/p',eta,u,value,p,
                  'outer log2/group',log(counts[u,eta])/log(2),flush=True)


if __name__=='__main__':main()
