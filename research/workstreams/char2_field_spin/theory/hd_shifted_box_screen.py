"""FLOAT neighborhood bounds for the alternating-pair construction.

Uses exact support-count envelopes, directed scalar normalizers, and FLOAT
positive aliases for the two regional coefficient matrices. This supplies
coverage proposals only, never a directed whole-code certificate.
"""
from fractions import Fraction
from math import log2

import numpy as np
from scipy.optimize import minimize

from hd_randomized_outer_state_profile import multiply,power,total
from hd_shifted_box_normalizer import histogram_ranges,scalar_envelope
from hd_shifted_general_count import ProfileBox,box_count_upper
from hd_shifted_general_screen import general_coefficients
from hd_shifted_general_sweep import CHECKPOINTS
from hd_shifted_good_pair_screen import (good_step_coefficients,
    good_pair_coefficients,low_impurity_frequency,paired_frequency)


def support_count_markers(box,pairs=4096):
    """FLOAT proposes positive rational tilts; exact count evaluation follows."""
    midpoint = np.array([(lo+hi)/2 for lo,hi in box.ranges])
    def objective(theta):
        x,y,p,q = theta
        matrix = tuple(tuple(v*x+np.logaddexp2(0.,y+u*p+v*q) for v in range(2))
                       for u in range(2))
        trace = power(matrix,pairs)
        return float(total((trace[0][0],trace[1][1]))-midpoint@theta)
    result = minimize(objective,np.zeros(4),method='BFGS',options={'maxiter':150,'gtol':1e-6})
    theta = np.clip(result.x,-32.,32.)
    exact = tuple(Fraction(float(2**x)) for x in theta)
    return tuple(value.limit_denominator(10**8) or value for value in exact)


def coefficient_envelope(z,alpha,odds,ranges,grid=64):
    """Unnormalized regional alias envelope for all target triples in ranges.

    Truncate above the largest possible degree, then bound every residue
    that the neighborhood can request. Histogram and tilt denominators
    are deliberately NOT divided out; scalar_envelope handles them jointly.
    """
    if grid<64 or grid&(grid-1):
        raise ValueError('Expected a power-of-two grid at least 64')
    target = tuple(hi for lo,hi in ranges)
    dirty_coefficients = general_coefficients(z,alpha,odds)
    dirty_coefficients[:,0,0,:] = 0
    dirty_coefficients[:,1,0,:] = 0
    dirty_coefficients[:,0,1,:] = 0
    for axis,maximum in enumerate(target,1):
        if maximum<32:
            index = [slice(None)]*4
            index[axis] = slice(maximum+1,None)
            dirty_coefficients[tuple(index)] = 0
    dirty = np.fft.fftn(dirty_coefficients,s=(grid,grid,grid),axes=(1,2,3))
    good = low_impurity_frequency(good_step_coefficients(z,alpha),odds,grid,target)
    composed = low_impurity_frequency(good_pair_coefficients(z,alpha),odds,grid,target)
    values = paired_frequency(dirty,good,composed)
    del dirty,good,composed,dirty_coefficients
    a,b,c,d = values[:,0,0,0].real
    eigenvalue = (a+d+np.sqrt((a-d)**2+4*b*c))/2
    vector = np.array([b,eigenvalue-a] if eigenvalue-a>=eigenvalue-d else [eigenvalue-d,c])
    if not eigenvalue>0 or not np.all(vector>0):
        raise ArithmeticError('Degenerate positive similarity scaling')
    vector /= np.max(vector)
    values *= (vector[None,:]/vector[:,None]/eigenvalue).ravel()[:,None,None,None]
    for _ in range(6):
        a,b,c,d = values
        trace,bc = a+d,b*c
        values = np.stack((a*a+bc,b*trace,c*trace,d*d+bc))
    aliases = np.fft.ifftn(values,axes=(1,2,3))
    residues = [np.arange(grid) if hi-lo+1>=grid else np.unique(np.arange(lo,hi+1)%grid)
                for lo,hi in ranges]
    requested = aliases[np.ix_(np.arange(4),*residues)]
    if np.max(abs(requested.imag))>1e-10 or np.min(requested.real)<-1e-10:
        raise ArithmeticError('FLOAT neighborhood alias diagnostic failed')
    cap = np.max(requested.real,axis=(1,2,3))
    if np.min(cap)<1e-10:
        raise ArithmeticError('FLOAT neighborhood coefficient below cancellation floor')
    logs = np.log2(cap+1e-12).reshape((2,2))
    logs += 64*log2(eigenvalue)+np.log2(vector[:,None]/vector[None,:])
    return logs,dict(residue_counts=[len(row) for row in residues],
                     minimum_normalized_entry_cap=float(np.min(cap)),
                     maximum_imaginary=float(np.max(abs(requested.imag))),directed=False)


def bound(box,z,alpha,odds,grid=64):
    exact = tuple(tuple(Fraction(float(2**x)) for x in row) for row in odds)
    markers = tuple(tuple(value.limit_denominator(10**9) or value for value in row) for row in exact)
    # The FLOAT operator uses exactly the rationalized marker proposals used
    # by the directed scalar calculation, rather than a nearby second tilt.
    odds = tuple(tuple(log2(marker) for marker in row) for row in markers)
    histograms = histogram_ranges(box)
    first,d0 = coefficient_envelope(z,alpha,odds[0],histograms[0][1:],grid)
    second,d1 = coefficient_envelope(z,alpha,odds[1],histograms[1][1:],grid)
    scalar = scalar_envelope(4096,box,Fraction(z),Fraction(alpha),markers)
    support_markers = support_count_markers(box)
    count = box_count_upper(4096,box,support_markers)
    if count==0:
        raise ValueError('Empty support box')
    run = power(multiply(tuple(map(tuple,first)),tuple(map(tuple,second))),128)
    value = log2(count)+float(scalar)+total(run[0])
    return dict(box=box.ranges,mark=z,alpha=alpha,grid=grid,
                marker_log_odds=odds,support_markers=list(map(str,support_markers)),
                log2_support_count_upper=log2(count),log2_scalar_upper=str(scalar),
                log2_upper_screen=float(value),regional_diagnostics=[d0,d1],
                directed=False,whole_code_certificate=False)


if __name__=='__main__':
    import argparse,json
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--radius',type=int,nargs='+',default=(0,1,2))
    parser.add_argument('--checkpoint',choices=tuple(CHECKPOINTS),default='asymmetric')
    parser.add_argument('--grid',type=int,default=64)
    args = parser.parse_args()
    profile,z,alpha,odds = CHECKPOINTS[args.checkpoint]
    for radius in args.radius:
        if radius<0:
            raise ValueError('Radius must be nonnegative')
        box = ProfileBox(*((max(0,n-radius),min(4096,n+radius)) for n in profile))
        print(json.dumps(bound(box,z,alpha,odds,args.grid)),flush=True)
