"""Floating rank-potential screen: striped MDS / prefix / global route / 2 states.

This is a research bound, not a directed or complete certificate. It counts
exact supports and their shortened source dimensions, not field messages.
The global route is uniform on ALL N positions; these numbers cannot be
imported for the faster regional routing distribution.
"""
import argparse
import json
import math

import numpy as np
from scipy.optimize import minimize, minimize_scalar
from scipy.special import logsumexp


def rank_transfers(y, g):
    matrices = []
    for active in (0,1):
        matrix = np.zeros((4,4))
        for old in range(4):
            for new in range(4):
                if not active and new != old:
                    continue
                if active and any(not(old & bit) and not(new & bit) for bit in (1,2)):
                    continue
                returns = (old & ~new).bit_count() if active else 0
                for live_output in (0,1):
                    if not live_output:
                        if not active and old in (1,2):
                            continue
                        if active and not old:
                            continue
                        if (active or old) and not new:
                            continue  # required local determinant excludes this
                    elif not active and not old:
                        continue
                    output_pivot = int(not live_output and bool(active or old))
                    matrix[old,new] += g**(returns+output_pivot)*y**live_output
        matrices.append(matrix)
    return tuple(matrices)


def row_log_transfer(b, k, logtau, g):
    """Coordinate-major one-prefix source-support transfer.

    Per prescribed MDS support w> b-k, its shortened dimension is w-b+k.
    An active-to-active prefix transition has weight1, a return has g.
    """
    scale = k*logtau/2
    source = np.array([-scale]+[-math.inf]*(b-k)+[
        (w-b+k)*logtau-scale for w in range(b-k+1,b+1)])
    weighted = np.full((b+1,b+1),-math.inf)
    weighted[0] = source
    for size in range(1,b+1):
        count = b-size+1
        weighted[size,:count] = np.logaddexp(weighted[size-1,:count],
                                             weighted[size-1,1:count+1])
    combinations = [np.array([math.log(math.comb(n,j)) for j in range(n+1)])
                    for n in range(b+1)]
    out = np.full((b+1,b+1),-math.inf)
    for incoming in range(b+1):
        births = np.arange(b-incoming+1)
        for returns in range(incoming+1):
            following = incoming-returns+births
            terms = (combinations[b-incoming]+combinations[incoming][returns]
                     +returns*math.log(g)+weighted[incoming-returns,births+returns])
            out[incoming,following] = np.logaddexp(out[incoming,following],terms)
    return out


def first_log_moment(matrix, logy, b, k, logtau, rows, steps=72):
    matrix = matrix+np.arange(b+1)[None,:]*logy
    vector = np.zeros(b+1)
    for _ in range(steps):
        following = logsumexp(matrix+vector[None,:],axis=1)
        vector = following-float(following.max())
    upper = float((logsumexp(matrix+vector[None,:],axis=1)-vector).max())
    return rows*(upper+k*logtau/2)+b*math.log(2)+float(vector.max()-vector.min())


def post_log_moment(logx, logy, n, g):
    t0,t1 = rank_transfers(math.exp(logy),g)
    matrix = t0+math.exp(logx)*t1
    vector = np.ones(4)
    for _ in range(96):
        following = matrix @ vector
        vector = following/float(following.max())
    upper = float(np.max((matrix @ vector)/vector))
    return n*math.log(upper)+math.log(vector[0]/float(vector.min()))


def double_prefix_rank_transfers(y, g):
    # Exact combinatorial support transitions; g charges each closed return.
    def options(old, active):
        if not active:
            return ((old,1.),)
        if not old:
            return ((1,1.),)
        return ((0,g),(1,1.))
    result = []
    for active in (0,1):
        matrix = np.zeros((4,4))
        for old1 in (0,1):
            for old2 in (0,1):
                for new1,p1 in options(old1,active):
                    for new2,p2 in options(old2,new1):
                        matrix[2*old1+old2,2*new1+new2] += p1*p2*y**new2
        result.append(matrix)
    return tuple(result)


def product_first_bound(marks,b,k,logtau,rows,g):
    """Product-vector envelope for the two-prefix row transfer.

    Majorize the MDS row-support polynomial by the sum of the zero support
    and tau^(-b+k) times all supports weighted by tau per active position.
    For a product vector v^tensor b, the ratio is a product plus a product.
    It is convex in the histogram of entering coordinate states, so its
    maximum is bounded by a homogeneous histogram. This avoids claiming
    independent states or a second random permutation.
    """
    t0,t1 = double_prefix_rank_transfers(math.exp(marks[0]),g)
    logs = np.array([0.,*marks[1:]])
    v = np.exp(logs)
    a = (t0 @ v)/v
    c = ((t0+math.exp(logtau)*t1) @ v)/v
    upper = float(np.logaddexp(b*np.log(a),-(b-k)*logtau+b*np.log(c)).max())
    return rows*upper+b*(math.log(4)+float(logs.max()-logs.min()))


def middle_class(alpha, tau, b=128, q=65536, k_log=20, distance=.4, pre=1, sharp=False):
    k = b//2
    rows = 2**k_log//k
    n = b*rows
    weight = round(alpha*n)
    cutoff = math.ceil(distance*n)-1
    if not 0 < weight <= n or tau < 1:
        raise ValueError('invalid weight or rank tilt')
    logtau = math.log(tau)
    g = (max(1/(q-1),1/tau) if sharp
         else 1/(q-1)+(1-1/(q-1))/tau)
    if pre == 1:
        row = row_log_transfer(b,k,logtau,g)
        first = lambda mark: first_log_moment(row,mark,b,k,logtau,rows)-weight*mark
        fit = minimize_scalar(first,bounds=(-logtau-6,0),method='bounded',
                              options={'maxiter':40})
    elif pre == 2:
        first = lambda marks: product_first_bound(marks,b,k,logtau,rows,g)-weight*marks[0]
        first_fits = [minimize(first, start, method='Powell',
                              bounds=((-logtau-6,0),(-20,20),(-20,20),(-20,20)),
                              options={'maxiter':24,'ftol':1e-9})
                      for start in ((-1.,0.,0.,0.),(-logtau/2,-1.,-1.,-1.))]
        fit = min(first_fits,key=lambda r:r.fun)
    else:
        raise ValueError('pre order must be one or two')
    logchoose = math.lgamma(n+1)-math.lgamma(weight+1)-math.lgamma(n-weight+1)
    second = lambda mark: (post_log_moment(*mark,n,g)-weight*mark[0]
                           -cutoff*mark[1]-logchoose)
    fits = [minimize(second,mark,method='L-BFGS-B',bounds=((-25,25),(-25,0)),
                     options={'maxiter':50,'ftol':1e-12})
            for mark in ((0.,-1.),(logtau/2,-logtau/2),(logtau,-logtau))]
    best = min(fits,key=lambda r:r.fun)
    # Condition the independent-coefficient experiment on every local map
    # being invertible. This factor is paid once, not at each zero equation.
    conditional = -n*math.log1p(-(q-2)/(q-1)**2)
    bound = fit.fun+best.fun-logtau+conditional
    return dict(outer=[b,k],q=q,n=n,weight=weight,tau=tau,pre=pre,sharp_live_pivot=sharp,
                first_mark=np.asarray(fit.x).tolist(),post_marks=best.x.tolist(),
                margin_bits=-float(bound)/math.log(2),
                conditioning_cost_bits=conditional/math.log(2),
                first_converged=bool(fit.success),post_converged=bool(best.success),
                rigorous=False,whole_code=False,route='uniform global, not regional')


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--outer',type=int,default=128)
    p.add_argument('--pre',type=int,choices=(1,2),default=1)
    p.add_argument('--sharp',action='store_true')
    p.add_argument('--tau',type=float,nargs='+',default=[4,16,64])
    p.add_argument('--density',type=float,nargs='+',default=[.2,.4,.6,.8])
    a = p.parse_args()
    for tau in a.tau:
        for alpha in a.density:
            print(json.dumps(middle_class(alpha,tau,a.outer,pre=a.pre,sharp=a.sharp)),flush=True)
