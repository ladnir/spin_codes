"""Feedback class counts and atom bounds without full Walsh inversion.

For an input-shape feedback count n and its unnormalized Walsh transform
n_hat, n(t) <= sum_a |n_hat(a)| / S, where S is the state-space size.
For any expansion class A, sum_{t in A} n(t) equals
sum_a n_hat(a) * 1_A_hat(a) / S. Characters with identical packet-window
histograms have the same n_hat for every shape, so both sums can use the
existing compressed character polynomials.

The returned peak is an upper bound on ALL atoms, not an exact nonzero
peak. Zero and expansion-class counts are exact. Polynomial construction
keeps its existing signed-int64 guards; every final dot product below uses
unbounded Python integers. This candidate does not change a verifier.
"""
import argparse
from fractions import Fraction as Q

import numpy as np

from mass_density_screen import baseline
from feedback_character_census import character_polynomials
from feedback_density import class_transforms
from group_moment import maps
from mass_density import conditioned_peaks


def summaries(rows, denominators, multiplicities, inverse, expansion):
    rows=np.asarray(rows)
    inverse=np.asarray(inverse)
    expansion=np.asarray(expansion)
    multiplicities=np.asarray(multiplicities)
    if (rows.ndim!=2 or rows.dtype.kind!='i' or inverse.ndim!=1 or inverse.dtype.kind not in 'iu'
            or multiplicities.ndim!=1 or multiplicities.dtype.kind not in 'iu'
            or rows.shape[1]!=len(multiplicities) or rows.shape[0]!=len(denominators)
            or len(inverse)!=len(expansion) or not len(inverse)
            or np.any(inverse<0) or np.any(inverse>=len(multiplicities))):
        raise ValueError('consistent exact integer character data required')
    size=len(inverse)
    if size*size >= 1<<63:
        raise ValueError('class-transform compression exceeds its signed-int64 guard')
    if not np.array_equal(np.bincount(inverse.astype(np.int64),minlength=len(multiplicities)),multiplicities):
        raise ValueError('character histogram multiplicities do not match')
    levels,sizes,transforms=class_transforms(expansion)
    compressed={}
    for v in levels:
        values=np.zeros(len(multiplicities),dtype=np.int64)
        # At most S summands, each of magnitude at most S. No float cast.
        np.add.at(values,inverse,transforms[v])
        compressed[v]=tuple(map(int,values))
    multiples=tuple(map(int,multiplicities))
    expected=-multiplicities.astype(np.int64).copy()
    expected[int(inverse[0])]+=size
    assert np.array_equal(sum((np.array(x,dtype=np.int64) for x in compressed.values()),
                              np.zeros(len(multiplicities),dtype=np.int64)),expected)
    result=[]
    for row,denominator in zip(rows,denominators):
        denominator=int(denominator)
        coefficients=tuple(map(int,row))
        if (denominator<1 or coefficients[int(inverse[0])]!=denominator
                or max(map(abs,coefficients))>denominator):
            raise ValueError('invalid shape denominator or trivial character')
        zero_numerator=sum(a*b for a,b in zip(coefficients,multiples))
        if zero_numerator%size or not 0<=zero_numerator//size<=denominator:
            raise ValueError('invalid zero-feedback count')
        zero=zero_numerator//size
        absolute=sum(abs(a)*b for a,b in zip(coefficients,multiples))
        peak=min(denominator,(absolute+size-1)//size)
        classes={}
        for v,indicator in compressed.items():
            numerator=sum(a*b for a,b in zip(coefficients,indicator))
            if numerator%size or not 0<=numerator//size<=denominator:
                raise ValueError('invalid expansion-class feedback count')
            classes[v]=numerator//size
        assert zero+sum(classes.values())==denominator
        assert zero<=peak<=denominator
        result.append((zero,peak,denominator,classes))
    return result


def compare_known(result,known,maximum):
    checks=0
    for shape,(zero,peak,den,classes) in known.items():
        if not shape or len(shape)>maximum:
            continue
        z,p,d,c=result[shape]
        assert z*den==zero*d and all(c[v]*den==n*d for v,n in classes.items())
        assert p*den>=peak*d
        checks+=1
    return checks


def combined_records(spectral,known,through=6):
    """Retain exact short-census peaks and intersect longer atom bounds.

    Both inputs must describe the same maps and packet distribution. The
    spectral records supply exact zero/class counts; the conditioned short
    census supplies another bound on the largest atom, including zero.
    The integer ceiling is taken only after multiplying by that shape's
    exact counting denominator.
    """
    maximum=max(map(len,spectral))
    compare_known(spectral,known,min(through,maximum))
    peaks=conditioned_peaks(known,min(through,maximum),maximum)
    result={}
    for shape,(zero,peak,den,classes) in spectral.items():
        if len(shape)<=through:
            result[shape]=known[shape]
            continue
        scaled=peaks[shape]*den
        ceiling=(scaled.numerator+scaled.denominator-1)//scaled.denominator
        bounded=min(peak,ceiling)
        assert zero<=bounded<=den
        result[shape]=(zero,bounded,den,classes)
    return result


def build(maximum=10, known=None):
    if not isinstance(maximum,int) or isinstance(maximum,bool) or not 1<=maximum<=10:
        raise ValueError('this bounded polynomial construction supports one through ten windows')
    shapes,rows,denominators,multiplicities,inverse=character_polynomials(maximum)
    images,_,_=maps()
    expansion=np.array([image.bit_count() for image in images],dtype=np.int64)
    records=summaries(rows,[denominators[s] for s in shapes],multiplicities,inverse,expansion)
    result={tuple(b for b in range(1,5) for _ in range(shape[b-1])):record
            for shape,record in zip(shapes,records) if sum(shape)}
    if known is not None:
        checks=compare_known(result,known,maximum)
        print('COMPRESSED SPECTRAL FEEDBACK: independent census comparisons',checks,flush=True)
    for j in range(1,maximum+1):
        worst=max((Q(peak,den),shape) for shape,(_,peak,den,_) in result.items() if len(shape)==j)
        print('SPECTRAL ATOM UPPER',j,worst,flush=True)
    print('Exact zero and expansion-class counts; peak fields are upper bounds, not exact peaks.',flush=True)
    return result


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--maximum',type=int,default=10)
    args=parser.parse_args()
    from full_feedback_census import census
    build(args.maximum,census(min(2,args.maximum)))
