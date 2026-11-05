"""Directed selective-pair receipts for the alternating adjacent-pair route.

Good-good state paths are composed before alpha. Other cases retain the
positive GD+DG+DD envelope. Exact packet caps and interval radix-2 FFTs
enclose a class contribution; this is never a whole-code certificate.
"""
from fractions import Fraction as F
from math import comb,isfinite,sqrt

import numpy as np
from flint import arb,ctx

from hd_grouped_packet_screen import OCC
from hd_shifted_alias_certificate import (LIMIT,packet_table,rational,
    binary64_enclosure,interval_fft,interval_matrix_power,selected_inverse,
    _add,_mul,_cadd,_cmul,_slice,_map,_finite,_multinomial)
from hd_shifted_general_certificate import general_coefficients
from hd_shifted_general_count import exact_joint_count,matching_histograms
from hd_shifted_general_sweep import CHECKPOINTS
from hd_shifted_good_pair_screen import good_geometry
from hd_sumrank_paired_certificate import scalar_power


def good_coefficients(z,alpha,bits=60):
    """Arb coefficients C[k,j] in X^k y^j, X=x10+x01, before tilts."""
    table = packet_table(F(z),bits)
    impurity,degree,counts,indices = good_geometry()
    groups = [(int(k),int(j),int(n),tuple(rational(v) for v in table[tuple(OCC[index])]))
              for k,j,n,index in zip(impurity,degree,counts,indices)]
    exponent = rational(alpha)
    unit_exponent = F(alpha)==1
    result = {}
    for k,j,n,left in groups:
        for ell,h,m,right in groups:
            key = k+ell,j+h
            if key not in result:
                result[key] = [arb(0) for _ in range(4)]
            for r in range(2):
                for s in range(2):
                    value = left[2*r]*right[s]+left[2*r+1]*right[2+s]
                    if not value.is_zero():
                        result[key][2*r+s] += n*m*(value if unit_exponent else value**exponent)
    return result


def components(single,paired,markers,target):
    """Sparse tilted G,D,C coefficients, truncated only above target degrees."""
    markers = tuple(map(F,markers))
    if len(markers)!=3 or min(markers)<=0 or len(target)!=3 or min(target)<0:
        raise ValueError('Invalid markers or target')
    powers = [[rational(marker)**j for j in range(65)] for marker in markers]
    good,dirty,composed = {},{},{}
    for a in range(33):
        for b in range(33-a):
            for c in range(33-a-b):
                if any(n>maximum for n,maximum in zip((a,b,c),target)):
                    continue
                index = a*33**2+b*33+c
                values = tuple(single[entry][index] for entry in range(4))
                if all(value.is_zero() for value in values):
                    continue
                marker = powers[0][a]*powers[1][b]*powers[2][c]
                destination = good if a+b<=1 else dirty
                destination[a,b,c] = tuple(value*marker for value in values)
    for (k,c),values in paired.items():
        for a in range(k+1):
            b = k-a
            if any(n>maximum for n,maximum in zip((a,b,c),target)):
                continue
            marker = comb(k,a)*powers[0][a]*powers[1][b]*powers[2][c]
            composed[a,b,c] = tuple(value*marker for value in values)
    return good,dirty,composed


def matrix_product(left,right):
    entries = [_cadd(_cmul(_slice(left,2*r),_slice(right,s)),
                     _cmul(_slice(left,2*r+1),_slice(right,2+s)))
               for r in range(2) for s in range(2)]
    return tuple(tuple(np.stack([entry[part][end] for entry in entries])
                       for end in range(2)) for part in range(2))


def interval_polynomial(coefficients,factors,grid,dimensions=3):
    """Convert and fold positive coefficients modulo the FFT grid."""
    lower = np.zeros((4,)+(grid,)*dimensions)
    upper = np.zeros_like(lower)
    for degree,values in coefficients.items():
        if len(degree)!=dimensions:
            raise ValueError('Polynomial dimension mismatch')
        index = tuple(int(n)%grid for n in degree)
        for entry,value in enumerate(values):
            if value.lower()<0:
                raise ArithmeticError('Expected a known nonnegative coefficient')
            if not value.is_zero():
                selected = (entry,)+index
                bounds = binary64_enclosure(value*factors[entry])
                lower[selected],upper[selected] = _add((lower[selected],upper[selected]),bounds)
    return ((lower,upper),(np.zeros_like(lower),np.zeros_like(lower)))


def inverse_alias_box(powered,target_ranges):
    """Enclose the entrywise maximum of positive aliases in a target box."""
    shape = powered[0][0].shape
    if len(shape)!=4 or shape[0]!=4 or len(set(shape[1:]))!=1:
        raise ValueError('Expected four entries on a cubic FFT grid')
    grid = shape[1]
    if len(target_ranges)!=3 or any(not isinstance(lo,int) or not isinstance(hi,int)
                                   or not 0<=lo<=hi for lo,hi in target_ranges):
        raise ValueError('Expected three nonnegative target intervals')
    # The inverse FFT is unnormalized. Its factor is an exact binary power,
    # but enclose the multiplication to retain subnormals.
    aliases = interval_fft(powered,(1,2,3),inverse=True)
    aliases = tuple(_mul(part,(1./grid**3,1./grid**3)) for part in aliases)
    residues = [np.arange(grid) if hi-lo+1>=grid else np.unique(np.arange(lo,hi+1)%grid)
                for lo,hi in target_ranges]
    requested = _slice(aliases,np.ix_(np.arange(4),*residues))
    if np.any(requested[1][0]>0) or np.any(requested[1][1]<0):
        raise ArithmeticError('A neighborhood imaginary coefficient excludes zero')
    # max of upper endpoints majorizes every requested positive alias.
    real = tuple(np.max(endpoint,axis=(1,2,3)) for endpoint in requested[0])
    imag = (np.min(requested[1][0],axis=(1,2,3)),
            np.max(requested[1][1],axis=(1,2,3)))
    return real,imag


def regional_power(single,paired,markers,target,grid=64,progress=None):
    """Normalized interval transform of the regional polynomial.

    The returned transform can be inverted once and reused for every
    target whose coordinates do not exceed the supplied truncation cap.
    A cap of (64,64,64) omits no single-step or composed coefficient.
    """
    if grid<64 or grid&(grid-1):
        raise ValueError('Expected a power-of-two grid at least64')
    good,dirty,composed = components(single,paired,markers,target)
    return component_power(good,dirty,composed,grid,progress)


def component_power(good,dirty,composed,grid=64,progress=None,dimensions=3):
    """Shared directed transform for already tilted positive components."""
    if grid<64 or grid&(grid-1) or dimensions not in (2,3):
        raise ValueError('Invalid regional grid or dimension')
    sums = lambda coeff:tuple(sum((value[i] for value in coeff.values()),arb(0)) for i in range(4))
    g,d,c = map(sums,(good,dirty,composed))
    product = lambda a,b:tuple(a[2*r]*b[s]+a[2*r+1]*b[2+s] for r in range(2) for s in range(2))
    positive = tuple(sum(parts,arb(0)) for parts in zip(c,product(g,d),product(d,g),product(d,d)))
    a,b,c,d = map(float,positive)
    eigenvalue = (a+d+sqrt((a-d)**2+4*b*c))/2
    proposal = [b,eigenvalue-a] if eigenvalue-a>=eigenvalue-d else [eigenvalue-d,c]
    if not all(isfinite(value) and value>0 for value in proposal):
        raise ArithmeticError('Invalid positive similarity proposal')
    vector = [arb(value/max(proposal)) for value in proposal]
    row_cap = max(((positive[2*r]*vector[0]+positive[2*r+1]*vector[1])/vector[r]).upper()
                  for r in range(2))
    # rho is an exact binary64 witness. Its square is exact at precision>=128.
    rho = arb(float(np.nextafter(float(row_cap.sqrt().upper()),np.inf)))
    if not rho*rho>=row_cap:
        raise ArithmeticError('Shared two-step scaling does not dominate row sums')
    factor = [vector[s]/vector[r]/rho for r in range(2) for s in range(2)]
    paired_factor = [value/rho for value in factor]
    if progress:
        progress(dict(stage='good-pair-interval-fft',grid=grid))
    with np.errstate(over='raise',invalid='raise',divide='raise',under='ignore'):
        transformed = []
        for coeff,scale in ((good,factor),(dirty,factor),(composed,paired_factor)):
            transformed.append(interval_fft(interval_polynomial(coeff,scale,grid,dimensions),
                                            tuple(range(1,dimensions+1))))
        g,d,c = transformed
        values = _cadd(c,_cadd(matrix_product(g,d),
                       _cadd(matrix_product(d,g),matrix_product(d,d))))
        del g,d,c,transformed
        _finite(values)
        powered = interval_matrix_power(values,64)
        del values
    return powered,vector,rho


def regional_numerator(single,paired,markers,target,grid=64,progress=None,target_ranges=None):
    if target_ranges is not None and (len(target_ranges)!=3
            or any(not isinstance(lo,int) or not isinstance(hi,int) or not 0<=lo<=hi
                   for lo,hi in target_ranges)
            or tuple(hi for lo,hi in target_ranges)!=tuple(target)):
        raise ValueError('Target must be the upper endpoint of each target range')
    powered,vector,rho = regional_power(single,paired,markers,target,grid,progress)
    with np.errstate(over='raise',invalid='raise',divide='raise',under='ignore'):
        if target_ranges is None:
            selected = selected_inverse(powered,target)
        else:
            selected = inverse_alias_box(powered,target_ranges)
    if np.any(selected[1][0]>0) or np.any(selected[1][1]<0):
        raise ArithmeticError('Imaginary coefficient excludes zero')
    values = tuple(arb(max(0.,selected[0][1][i]))*rho**128*vector[i//2]/vector[i%2]
                   for i in range(4))
    return values,dict(grid=grid,composition_steps=2,paired_blocks=64,
        good_maximum_singletons=1,similarity_vector=list(map(str,vector)),rho=str(rho),
        normalized_alias_real_lower=selected[0][0].tolist(),
        normalized_alias_real_upper=selected[0][1].tolist(),
        normalized_alias_imag_lower=selected[1][0].tolist(),
        normalized_alias_imag_upper=selected[1][1].tolist(),
        target_ranges=target_ranges,directed=True)


def receipt(profile,z,alpha,markers,*,grid=64,precision=256,bits=60,threshold=-40,progress=None):
    profile = tuple(profile)
    z,alpha = F(z),F(alpha)
    markers = tuple(tuple(map(F,row)) for row in markers)
    if not 0<z<1 or not 0<alpha<=1 or precision<128 or len(markers)!=2:
        raise ValueError('Invalid directed witness')
    count = exact_joint_count(4096,*profile)
    if count<=0 or not sum(profile[:2]):
        raise ValueError('Empty or zero-message profile')
    histograms = matching_histograms(4096,*profile)
    with ctx.workprec(precision):
        if progress:
            progress(dict(stage='unmarked-single-coefficients'))
        single = general_coefficients(z,alpha,(F(1),)*3,bits=bits,progress=progress)
        if progress:
            progress(dict(stage='ordered-good-good-coefficients',geometry_pairs=425**2))
        paired = good_coefficients(z,alpha,bits)
        regions,diagnostics = [],[]
        for parity,(histogram,row) in enumerate(zip(histograms,markers)):
            if progress:
                progress(dict(stage='regional-operator',parity=parity))
            region,diagnostic = regional_numerator(single,paired,row,histogram[1:],grid,progress)
            normalizer = arb(_multinomial(histogram))
            for marker,degree in zip(row,histogram[1:]):
                normalizer *= rational(marker)**degree
            regions.append(tuple(value/normalizer for value in region))
            diagnostics.append(diagnostic)
        first,second = regions
        paired_region = tuple(first[2*r]*second[s]+first[2*r+1]*second[2+s]
                              for r in range(2) for s in range(2))
        product = scalar_power(paired_region,128)
        log2 = lambda value:value.log()/arb(2).log()
        exponent = rational(alpha)*(sum(profile[:2])*(2048-128*log2(arb(255)))
                                    -log2(arb(255))-LIMIT*log2(rational(z)))
        bound = log2(arb(count))+exponent+log2(product[0]+product[1])
        ceiling = F(int((bound.upper()*10**12).upper().ceil().fmpz()),10**12)
        return dict(profile=profile,mark=str(z),alpha=str(alpha),
                    markers=[list(map(str,row)) for row in markers],grid=grid,
                    precision=precision,packet_mantissa_bits=bits,
                    log2_upper_enclosure=str(bound),log2_upper_rational=str(ceiling),
                    accepted=bool(bound<threshold),threshold=threshold,
                    regional_diagnostics=diagnostics,directed=True,whole_code_certificate=False)


if __name__=='__main__':
    import argparse,json
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--checkpoint',choices=tuple(CHECKPOINTS),default='asymmetric')
    parser.add_argument('--grid',type=int,default=64)
    parser.add_argument('--precision',type=int,default=256)
    args = parser.parse_args()
    profile,z,alpha,odds = CHECKPOINTS[args.checkpoint]
    markers = tuple(tuple(F(float(2**v)).limit_denominator(10**9) or F(float(2**v)) for v in row) for row in odds)
    print(json.dumps(receipt(profile,F(z),F(alpha),markers,grid=args.grid,precision=args.precision,
        progress=lambda item:print(json.dumps(item),flush=True))),flush=True)
