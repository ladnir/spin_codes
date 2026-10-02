"""Selected row-weight classes using actual reference packet probabilities.

Reuses the independent-row Bernoulli domination lemma. This binary64 search
does not cover mixed row classes or certify a distance. It tests whether
forgetting packet-weight frequencies is a significant remaining loss.
"""
import argparse
from math import comb, factorial, log

import numpy as np

import attack_cache
import shape_potential as shape
import shape_grid
import mass_density_screen as screen
from row_counts import row_gamma_function
from row_verify import packet_law


def unpenalized(families, rho=.9):
    """Remove the proof-only all-four weight, exactly through shape cutoff."""
    if not 0 < rho <= 1:
        raise ValueError('penalty must be in (0,1]')
    result=[]
    for j,family in enumerate(families):
        if len(family)>1:
            catalog=sorted(shape.expected_shapes(j))
            if len(catalog)!=len(family):
                raise ValueError('complete shape catalog required')
            counts=np.array([s.count(4) for s in catalog])
        else:
            # Every shape has at most j full packets. This is conservative.
            counts=np.array([j])
        result.append(family/rho**counts[:,None,None])
    return result


def average(families, theta):
    theta=np.asarray(theta,dtype=float)
    if theta.shape!=(4,) or (theta<0).any() or not np.isclose(theta.sum(),1):
        raise ValueError('nonzero-packet weight distribution required')
    result=[]
    for j,family in enumerate(families):
        if len(family)==1:
            result.append(family[0]);continue
        catalog=sorted(shape.expected_shapes(j))
        if len(catalog)!=len(family):raise ValueError('complete shape catalog required')
        weights=[]
        for s in catalog:
            value=factorial(j)
            for b,p in enumerate(theta,1):
                n=s.count(b)
                value*=p**n/factorial(n)
            weights.append(value)
        assert np.isclose(sum(weights),1)
        result.append(np.einsum('s,sij->ij',weights,family))
    return np.array(result)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--tilts',nargs='+',default=['.072','.088'])
    parser.add_argument('--rounds',type=int,nargs='+',default=[2,3,4])
    parser.add_argument('--groups',type=int,nargs='+',default=[96,128])
    parser.add_argument('--active-rows',type=int,nargs='+',default=[2,4])
    parser.add_argument('--weights',type=int,nargs='+',default=[64,80,96,112,128])
    parser.add_argument('--probabilities',type=float,nargs='+',default=[.3,.4,.5,.6,.7,.8,.9])
    parser.add_argument('--directory',default='tmp/four-bit-attack-cache')
    args=parser.parse_args()
    if (any(not 2<=r<=32 for r in args.rounds) or any(not 1<=q<=2048 for q in args.groups)
            or any(not 1<=r<=4 for r in args.active_rows) or any(not 1<=w<256 for w in args.weights)
            or any(not 0<p<1 for p in args.probabilities)):
        parser.error('invalid row-class or witness parameters')
    print('SELECTED ROW-CLASS BINARY64 SCREEN, NOT A CERTIFICATE',flush=True)
    data=attack_cache.build(args.tilts,directory=args.directory)
    caps=screen.baseline.authenticated_caps()
    gamma={w:row_gamma_function(caps,w,w,exclude_zero=True) for w in args.weights}
    best={}
    for tilt,(base,feedback,details) in data.items():
        shaped=shape.shape_matrices(base,feedback,details,tilt,'.9',2)
        families=shape.as_families(base,shaped)
        for rounds in args.rounds:
            mixed=unpenalized(shape_grid.retarget(families,details['spectrum'],tilt,rounds))
            for r in args.active_rows:
                for p in args.probabilities:
                    active,theta=packet_law(p,r)
                    local=average(mixed,theta)
                    region=screen.float_placement(local,max(args.groups))
                    for q in args.groups:
                        matrix=screen.baseline.matrix_for_probabilities(region,[active]*q)
                        inner=screen.baseline.log_power_moment(matrix,256,screen.baseline.TAIL_TERMINAL)
                        for w in args.weights:
                            value=(inner+float(tilt)*193986+q*r*gamma[w](p)
                                   +log(comb(2048,q))+q*log(comb(4,r)))/log(2)
                            key=(q,rounds,r,w)
                            if key not in best or value<best[key][0]:best[key]=(value,tilt,p)
            print('ROW SCREEN progress tilt/R',tilt,rounds,flush=True)
    for key,result in sorted(best.items()):
        print('ROW BEST q/R/active_rows/weight',*key,'log2/tilt/p',*result,flush=True)
    print('These row classes overlap union-support classes; their bounds cannot simply be added to those covers.')


if __name__=='__main__':main()
