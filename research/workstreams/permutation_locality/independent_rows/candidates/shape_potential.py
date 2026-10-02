"""Retain packet shapes until applying a continuation potential.

Research-only binary64 diagnostics. The shape matrices are rebuilt using
the existing outward local inequalities; global floating-point evaluation
is NOT a certificate. No encoder or production verifier is changed.
"""
import argparse
from contextlib import redirect_stdout
from fractions import Fraction as Q
from io import StringIO
from itertools import combinations_with_replacement
from math import comb, log

import numpy as np
from scipy.special import gammaln

import mass_density_screen as screen
from occupancy_memory import Z, F, C, U
from group_rank_one_verify import up


def expected_shapes(j):
    return set(combinations_with_replacement(range(1,5),j))


def shape_matrices(base, feedback, details, tilt, penalty, rounds, maximum=8):
    """Intersect existing all-shape and same-shape component bounds.

The input weight tilt is one. Every occupancy up to maximum must contain
all shapes, including the all-one shape. Higher occupancies retain base.
No mass/density column replacement or entrywise min of coupled alternative
decompositions is performed.
"""
    if type(maximum) is not int or not 1<=maximum<=8:
        raise ValueError('complete shape range must be in 1..8')
    if len(base)!=33 or any(t.nrows()!=11 or t.ncols()!=11 for t in base):
        raise ValueError('all 33 eleven-coordinate occupancy operators required')
    if not 0<Q(penalty)<=1 or Q(tilt)<=0 or type(rounds) is not int or rounds<1:
        raise ValueError('positive tilt, penalty in (0,1], positive update count required')
    if details.get('parameters')!=(Q(tilt),Q(penalty),rounds):
        raise ValueError('shape records must match the tilt, penalty and update count')
    moments=details['moments']
    spectrum=details['spectrum']
    result={}
    with redirect_stdout(StringIO()):
        for j in range(1,maximum+1):
            records={}
            for shape in sorted(expected_shapes(j)):
                counts=tuple(shape.count(b) for b in range(1,5))
                if counts not in moments:
                    raise ValueError('missing shape output moment')
                if j<=6 and (shape not in feedback or shape not in details['density']):
                    raise ValueError('missing shape feedback/density record')
                if j<=3 and shape not in details['cancellations'][0]:
                    raise ValueError('missing shape joint cancellation record')
                candidate=list(base)
                candidate[j]=base[j]*1
                if j==1:
                    # A zero state emits the input packet itself and then
                    # enters its fresh feedback distribution. This entry is
                    # independent of the number of nonzero-state updates.
                    b=shape[0]
                    activation=((-screen.baseline.arb(tilt)*b).exp()
                                *screen.baseline.arb(penalty)**int(b==4))
                    candidate[j][Z,F]=min(candidate[j][Z,F],up(activation))
                candidate=screen.window_histogram.refine(candidate,details['histograms'],
                           {counts:moments[counts]},penalty,rounds,'1')
                if j<=6:
                    candidate=screen.full_feedback_refinement.refine(
                        candidate,{shape:feedback[shape]},spectrum,tilt,penalty,rounds,'1')
                    record=details['density'][shape]
                    scale=screen.baseline.arb(penalty)**shape.count(4)
                    for source,value in zip((C,*range(U,U+5)),
                           (record['density'],*(record['uniform'][v] for v in sorted(spectrum)))):
                        candidate[j][source,C]=min(candidate[j][source,C],up(value*scale))
                if j<=3:
                    joint=({shape:details['cancellations'][0][shape]},spectrum)
                    candidate=screen.cancellation_joint.refine(candidate,joint,tilt,penalty,rounds,'1')
                matrix=candidate[j]
                if any(matrix[a,b]<0 or matrix[a,b]>base[j][a,b]
                       for a in range(11) for b in range(11)):
                    raise ValueError('shape refinement must be nonnegative and no larger than base')
                records[shape]=matrix
            result[j]=records
    return result


def as_families(base, shaped):
    """Check completeness before dropping shape labels for numerical work."""
    if any(type(j) is not int or not 1<=j<=8 for j in shaped):
        raise ValueError('shape-specific occupancy keys must be in 1..8')
    arrays=screen.as_array(base)
    result=[]
    for j,matrix in enumerate(arrays):
        if j in shaped:
            if set(shaped[j])!=expected_shapes(j):
                raise ValueError('incomplete shape family')
            family=screen.as_array([shaped[j][s] for s in sorted(shaped[j])])
        else:
            family=matrix[None,:,:]
        if not np.isfinite(family).all() or (family<0).any() or (family>matrix).any():
            raise ValueError('invalid diagnostic family')
        result.append(family)
    return result


def placement_weights(degree,epochs=64,windows=32):
    """Normalized, without-replacement placement weights, binary64 only."""
    if any(type(x) is not int for x in (degree,epochs,windows)) or not 0<=degree<=epochs*windows or min(epochs,windows)<1:
        raise ValueError('valid integer placement geometry required')
    result=[]
    for epoch in range(1,epochs+1):
        limit=min(degree,epoch*windows)
        previous=(epoch-1)*windows
        records=[]
        for k in range(min(windows,limit)+1):
            rs=np.arange(k,min(limit,k+previous)+1)
            old=rs-k
            value=(log(comb(windows,k))+gammaln(previous+1)-gammaln(old+1)
                   -gammaln(previous-old+1)-gammaln(epoch*windows+1)
                   +gammaln(rs+1)+gammaln(epoch*windows-rs+1))
            records.append((k,rs,old,np.exp(value)))
        result.append((limit,records))
    return result


def region_action(families,potential,weights):
    """Backward continuation bound; maximize a whole row after applying v.

Shape choices may depend on remaining occupancy and envelope coordinate.
This is an upper relaxation, not an assertion about actual adaptive setup.
"""
    current=np.asarray(potential,dtype=float)[None,:]
    for limit,records in weights:
        nxt=np.zeros((limit+1,current.shape[1]))
        for k,rs,old,prob in records:
            values=np.matmul(families[k],current[old].T).max(axis=0).T
            nxt[rs]+=prob[:,None]*values
        current=nxt
    return current


def common_potential(families,weights,masses,initial,regions=256,iterations=12,terminal=None):
    """Collatz upper proposals, including the finite terminal factor.

Every iterate is usable if its inequalities are later replayed outward;
convergence is an optimization aid, not a validity assumption.
"""
    terminal=screen.baseline.TAIL_TERMINAL if terminal is None else np.asarray(terminal,dtype=float)
    if any(type(n) is not int or n<1 for n in (regions,iterations)):
        raise ValueError('positive integer region and iteration counts required')
    v=np.asarray(initial,dtype=float).copy()
    if v.shape!=terminal.shape or not np.isfinite(v).all() or (v<=0).any():
        raise ValueError('strictly positive potential matching the terminal required')
    best=None
    for iteration in range(iterations):
        v/=v[0]
        value=masses@region_action(families,v,weights)
        if not np.isfinite(value).all() or (value<=0).any():
            raise ValueError('nonpositive or nonfinite continuation; binary64 screen unavailable')
        ratios=value/v
        lam=float(ratios.max())
        constant=float((terminal/v).max())
        proposal=regions*log(lam)+log(constant)
        if best is None or proposal<best['log_moment']:
            best={'log_moment':proposal,'lambda':lam,'terminal_factor':constant,
                  'potential':v.tolist(),'iteration':iteration,
                  'spread':float(ratios.max()/ratios.min()-1)}
        if ratios.max()/ratios.min()-1<1e-10:
            break
        v=value
    return best


def matrix_potential(matrix):
    # Positive terminal perturbation supplies nonzero auxiliary coefficients.
    v=np.ones(matrix.shape[0])
    for _ in range(128):
        v=matrix@v
        v/=v[0]
    return v


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--rounds',type=int,default=2)
    parser.add_argument('--tilt',default='.056')
    parser.add_argument('--penalty',default='.9')
    parser.add_argument('--groups',type=int,nargs='+',default=[80,96,128])
    parser.add_argument('--support',type=int,default=200)
    parser.add_argument('--cutoff',type=int,default=193986)
    parser.add_argument('--maximum',type=int,default=8)
    parser.add_argument('--iterations',type=int,default=12)
    args=parser.parse_args()
    if (not 1<=args.rounds<=32 or Q(args.tilt)<=0 or not 0<Q(args.penalty)<=1
        or any(not 1<=q<=2048 for q in args.groups) or not 38<=args.support<256
        or not 0<=args.cutoff<(1<<21) or not 1<=args.maximum<=8 or args.iterations<1):
        parser.error('invalid shape-potential screen parameters')
    print('SHAPE POTENTIAL: selected homogeneous support events; binary64, NOT a certificate.',flush=True)
    grid,feedback,details=screen.epoch_grid([args.tilt],args.penalty,rounds=args.rounds,details=True)
    base=grid[args.tilt]
    shaped=shape_matrices(base,feedback,details[args.tilt],args.tilt,args.penalty,args.rounds,args.maximum)
    families=as_families(base,shaped)
    collapsed=np.array([a.max(axis=0) for a in families])
    print('Complete retained shape counts:',{j:len(x) for j,x in shaped.items()},flush=True)
    baseline_region=screen.float_placement(screen.as_array(base),max(args.groups))
    refined_region=screen.float_placement(collapsed,max(args.groups))
    caps=screen.baseline.authenticated_caps()
    cdf=screen.baseline.integer_cdf(screen.baseline.weighted_cdf_upper(caps,1<<128,full_weight=1/Q(args.penalty)))
    shells=screen.baseline.weighted_union_shells(caps,full_weight=1/Q(args.penalty))
    count=min(Q(cdf[args.support]),shells[args.support])
    for q in args.groups:
        old,p=screen.score(baseline_region,q,args.support,count,args.tilt,cutoff=args.cutoff)
        masses=np.array([1.])
        for _ in range(q):masses=np.convolve(masses,[1-p,p])
        matrix=np.einsum('r,rij->ij',masses,refined_region[:q+1])
        refined=screen.baseline.log_power_moment(matrix,256,screen.baseline.TAIL_TERMINAL)
        outer=(float(args.tilt)*args.cutoff-q*screen.baseline.log_binomial_mass(256,args.support,p)
               +q*(log(count.numerator)-log(count.denominator))+log(comb(2048,q)))
        initial=matrix_potential(matrix)
        ratios=(matrix@initial)/initial
        linear_common=256*log(float(ratios.max()))+log(float((screen.baseline.TAIL_TERMINAL/initial).max()))
        print('SHAPE CONTROL q/u',q,args.support,'p',p,'original',old,
              'intersect-then-max',(refined+outer)/log(2),
              'linear-common',(linear_common+outer)/log(2),flush=True)
        common=common_potential(families,placement_weights(q),masses,initial,iterations=args.iterations)
        print('SHAPE COUPLED q/u',q,args.support,'log2 proposal',(common['log_moment']+outer)/log(2),
              'gain vs linear common bits',(linear_common-common['log_moment'])/log(2),
              'witness',common,flush=True)
    print('No complete support cover, outward global replay, or new distance certificate.',flush=True)


if __name__=='__main__':main()
