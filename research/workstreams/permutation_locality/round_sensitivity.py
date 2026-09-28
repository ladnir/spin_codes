"""Diagnostic only: change the update count in a cached r=2 envelope.

No implementation or goal is changed. Timings of the two-update code do
not apply to these alternative distributions. All scores are binary64.
"""
import argparse
from math import log
import numpy as np
from scipy.optimize import minimize_scalar
from occupancy_memory import Z,F,M,U
from mature_tail import TAIL_TERMINAL
from group_moment import maps
from occupancy_sensitivity import float_placement
from occupancy_screen import matrix_for_probabilities
from two_group_screen import log_power_moment,log_binomial_mass


def score(operators,groups,support,tilt,outer):
    region=float_placement(operators,groups)
    def objective(z):
        p=1/(1+np.exp(-z))
        return (log_power_moment(matrix_for_probabilities(region,[p]*groups),256,TAIL_TERMINAL)
                +float(tilt)*209715-groups*log_binomial_mass(256,support,p))
    fit=minimize_scalar(objective,bounds=(-8.,16.),method='bounded')
    return (fit.fun+outer)/log(2),1/(1+np.exp(-fit.x))


def rescale(base,spectrum,tilt,rounds):
    # r=0 denotes the limiting uniform refresh model, not zero updates.
    probability=2.**(-rounds) if rounds else 0.
    alpha=probability/.25;beta=(1-probability)/.75
    result=base.copy();levels=sorted(spectrum)
    for occupancy,old in enumerate(base):
        for i in range(1,old.shape[0]):
            result[occupancy,i,:]=old[i,:]*alpha
            result[occupancy,i,U:U+5]=old[i,U:U+5]*beta
            if occupancy:
                if i==M or i in (9,10):
                    result[occupancy,i,Z]=beta*old[i,Z]
                elif i==F or U<=i<U+5:
                    refresh=min(old[i,U+k]/spectrum[v] for k,v in enumerate(levels))
                    result[occupancy,i,Z]=alpha*old[i,Z]+(beta-alpha)*refresh
            elif U<=i<U+5:
                v=levels[i-U]
                result[occupancy,i,i]=beta*old[i,i]-(beta-alpha)*np.exp(-float(tilt)*v)/4
    assert (result>=-1e-15).all()
    return np.maximum(result,0.)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--snapshot',default='tmp/r2-q64-sensitivity.npz')
    parser.add_argument('--full-feedback',type=int,choices=range(2,11))
    parser.add_argument('--ablate',action='store_true',help='Invalid-bound diagnostics deleting selected zero returns')
    args=parser.parse_args()
    with np.load(args.snapshot,allow_pickle=False) as saved:
        base=saved['base'];outer=float(saved['outer'])
        groups,support,tilt,penalty=map(str,saved['parameters'])
    groups=int(groups);support=int(support);spectrum=maps()[2]
    data=None
    if args.full_feedback:
        from flint import ctx
        from full_feedback_census import census
        ctx.prec=192
        data=census(args.full_feedback)
    assert np.array_equal(rescale(base,spectrum,tilt,2),base)
    for rounds in (2,3,4,6,0):
        ops=rescale(base,spectrum,tilt,rounds)
        if data:
            from flint import arb,arb_mat
            from full_feedback_refinement import refine
            matrices=[arb_mat([[arb(float(x)) for x in row] for row in t]) for t in ops]
            changed=refine(matrices,data,spectrum,tilt,penalty,rounds or None)
            revised=np.array([[[float(t[i,j]) for j in range(11)] for i in range(11)] for t in changed])
            assert (revised<=ops+1e-15).all()
            ops=revised
        value,p=score(ops,groups,support,tilt,outer)
        print('DIAGNOSTIC ALTERNATIVE ENSEMBLE:',rounds if rounds else 'uniform refresh limit',
              'log2 score',value,'p',p,
              'outer improvement bits/group to reach -40',max(0.,value+40)/groups,flush=True)
        if args.ablate and rounds in (2,0):
            for name,rows in [('zero source',[Z]),('nonzero sources',list(range(1,11))),
                              ('uniform sources',list(range(U,U+5))),('all sources',list(range(11)))]:
                modified=ops.copy();modified[1:,rows,Z]=0
                value,p=score(modified,groups,support,tilt,outer)
                print('ABLATION INVALID AS BOUND:',rounds if rounds else 'uniform refresh limit',
                      name,'log2 score',value,'p',p,flush=True)


if __name__=='__main__':main()
