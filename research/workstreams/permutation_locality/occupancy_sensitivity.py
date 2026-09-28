"""Locate slack in the two-update bound. Altered operators are NOT bounds.

Binary64 normalized placement is diagnostic only. It must not be used by
the outward certificate or interpreted as evidence of a bad codeword.
"""
import argparse
from math import comb,log
import numpy as np
from scipy.optimize import minimize_scalar
from scipy.special import gammaln
from flint import ctx,fmpq_mat
from occupancy_model import placement,local_data
from occupancy_memory import prepare,Z,F,M,C,U
from occupancy_fresh_moment import fresh_census,refine as fresh_refine
from occupancy_window_average import prepare_inputs,averages,refine as window_refine
from occupancy_multi_average import refine as multi_refine
from fresh_collision import refine as collision_refine
from zero_moment import census as zero_census,refine as zero_refine
from mature_tail import census,lift,TAIL_TERMINAL
from pair_tail import census as pair_census
from mixing_rounds import transform
from bch_joint_support import authenticated_caps,support_caps
from shortened_bound import dimension_caps
from basis_lattice import improve_caps
from occupancy_allones import weighted_cdf
from occupancy_screen import matrix_for_probabilities
from two_group_screen import log_power_moment,log_binomial_mass


def float_placement(operators,degree,epochs=64,windows=32):
    """Stable hypergeometric recurrence; no outward rounding is claimed."""
    size=operators.shape[1]
    current=np.eye(size)[None,:,:]
    for epoch in range(1,epochs+1):
        limit=min(degree,epoch*windows)
        nxt=np.zeros((limit+1,size,size))
        total=epoch*windows;previous=total-windows
        for k in range(min(windows,len(operators)-1,limit)+1):
            rs=np.arange(k,min(limit,k+len(current)-1)+1)
            old=rs-k
            # All arguments here are nonnegative and feasible.
            log_weight=(log(comb(windows,k))+gammaln(previous+1)-gammaln(old+1)
                        -gammaln(previous-old+1)-gammaln(total+1)
                        +gammaln(rs+1)+gammaln(total-rs+1))
            nxt[rs]+=np.exp(log_weight)[:,None,None]*(current[old]@operators[k])
        current=nxt
    assert np.isfinite(current).all() and (current>=0).all()
    return current


def self_test():
    ops=[fmpq_mat([[1,1],[0,1]]),fmpq_mat([[1,0],[1,1]]),fmpq_mat([[2,1],[0,1]])]
    floats=np.array([[[float(t[i,j]) for j in range(2)] for i in range(2)] for t in ops])
    for epochs in (1,2,3,7):
        exact=placement(ops,epochs,2,fmpq_mat,lambda x:x,maximum_groups=2*epochs)
        actual=float_placement(floats,2*epochs,epochs,2)
        expected=np.array([[[float(t[i,j]) for j in range(2)] for i in range(2)] for t in exact])
        assert np.allclose(actual,expected,rtol=1e-12,atol=1e-12)
    print('Diagnostic normalized placement matches exact rational examples',flush=True)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--groups',type=int,default=64)
    parser.add_argument('--support',type=int,default=128)
    parser.add_argument('--tilt',default='.032')
    parser.add_argument('--penalty',default='.75')
    parser.add_argument('--snapshot',help='Optional diagnostic-only NumPy cache path')
    parser.add_argument('--load-snapshot',help='Reuse a matching diagnostic-only cache')
    args=parser.parse_args();self_test()
    if args.load_snapshot:
        with np.load(args.load_snapshot,allow_pickle=False) as saved:
            assert tuple(saved['parameters'])==(str(args.groups),str(args.support),args.tilt,args.penalty)
            run_tests(saved['base'],args,float(saved['outer']))
        return
    spectrum=authenticated_caps();dimensions=dimension_caps()
    caps=improve_caps(support_caps(spectrum,g=4,dimensions=dimensions),spectrum)
    counts=weighted_cdf(spectrum,caps,dimensions,args.penalty,prefix_flags=True)
    prepared=prepare(local_data(4));fresh=fresh_census(prepared)
    inputs=prepare_inputs();tails=census(inputs,prepared)
    pairs=pair_census(inputs,prepared);zeros=zero_census();ctx.prec=192
    values=averages(inputs,args.tilt)
    local=min(32,args.groups)
    ops=fresh_refine(prepared,fresh,args.tilt,local,args.penalty)
    ops=window_refine(prepared,values,args.tilt,local,args.penalty,ops)
    ops=multi_refine(ops,fresh,args.tilt,args.penalty)
    ops=collision_refine(ops,prepared,args.tilt,args.penalty)
    ops=zero_refine(ops,zeros,args.tilt,args.penalty)
    ops=lift(ops,inputs[3],tails,args.tilt,args.penalty,64,pairs)
    ops=transform(ops,inputs[3],args.tilt,2)
    n=len(TAIL_TERMINAL)
    base=np.array([[[float(t[i,j]) for j in range(n)] for i in range(n)] for t in ops])
    outer=args.groups*log(counts[args.support])+log(comb(2048,args.groups))
    if args.snapshot:
        np.savez_compressed(args.snapshot,base=base,outer=outer,
                            parameters=np.array([str(args.groups),str(args.support),args.tilt,args.penalty]))
    run_tests(base,args,outer)


def run_tests(base,args,outer):
    n=len(TAIL_TERMINAL)
    tests={'baseline':[], 'density cancellation':[(C,Z)],'fresh cancellation':[(F,Z)],
           'direct zero return':[(Z,Z)],
           'uniform cancellation':[(U+i,Z) for i in range(5)],
           'mature refresh cancellation':[(i,Z) for i in (M,9,10)],
           'non-mature returns':[(i,Z) for i in (Z,F,C)]+[(U+i,Z) for i in range(5)],
           'all returns to zero':[(i,Z) for i in range(n)],
           'density persistence':[(C,C)],'fresh density':[(F,C)],
           'uniform density':[(U+i,C) for i in range(5)],
           'mature lazy mass':[(i,M) for i in (M,9,10)]}
    for name,entries in tests.items():
        for factor in ([1.] if name=='baseline' else [.5,0.]):
            changed=base.copy()
            for i,j in entries:changed[1:,i,j]*=factor
            region=float_placement(changed,args.groups)
            def objective(z):
                p=1/(1+np.exp(-z))
                return (log_power_moment(matrix_for_probabilities(region,[p]*args.groups),256,TAIL_TERMINAL)
                        +float(args.tilt)*209715-args.groups*log_binomial_mass(256,args.support,p))
            result=minimize_scalar(objective,bounds=(-8.,16.),method='bounded')
            score=(result.fun+outer)/log(2)
            print('DIAGNOSTIC NOT A CERTIFICATE:',name,'factor',factor,'log2 score',score,
                  'p',1/(1+np.exp(-result.x)),flush=True)


if __name__=='__main__':main()
