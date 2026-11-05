"""Cached directed regional bounds and disjoint local support coverage.

The interval transform depends on the witness, not the support profile.
Invert it once, then reuse its positive alias coefficients across boxes.
The driver certifies only its declared root box, never the entire code.
"""
from dataclasses import dataclass
from fractions import Fraction as F
from math import comb

import numpy as np
from flint import arb,ctx

from hd_shifted_alias_certificate import interval_fft,_mul,rational
from hd_shifted_general_certificate import general_coefficients
from hd_shifted_general_count import (ProfileBox,transfer_value,joint_count_upper,
                                      joint_is_feasible)
from hd_shifted_good_pair_certificate import good_coefficients,regional_power
from hd_shifted_box_normalizer import histogram_ranges,scalar_envelope
from hd_shifted_box_certificate import DEFAULT_COUNT_MARKERS
from hd_sumrank_paired_certificate import scalar_power


def rational_ceiling(value):
    return F(int((value.upper()*10**12).upper().ceil().fmpz()),10**12)


@dataclass(frozen=True)
class RegionalAliases:
    """Upper endpoints of a complete, untruncated regional alias table."""
    upper: np.ndarray
    factors: tuple

    @classmethod
    def from_power(cls,powered,vector,rho):
        grid = powered[0][0].shape[1]
        dimensions = len(powered[0][0].shape)-1
        with np.errstate(over='raise',invalid='raise',divide='raise',under='ignore'):
            inverse = interval_fft(powered,tuple(range(1,dimensions+1)),inverse=True)
            inverse = tuple(_mul(part,(1./grid**dimensions,1./grid**dimensions)) for part in inverse)
        if np.any(inverse[1][0]>0) or np.any(inverse[1][1]<0):
            raise ArithmeticError('Imaginary alias excludes zero')
        upper = np.maximum(0.,inverse[0][1])
        if not np.all(np.isfinite(upper)):
            raise ArithmeticError('Nonfinite alias bound')
        upper.setflags(write=False)
        factors = tuple(rho**128*vector[r]/vector[s] for r in range(2) for s in range(2))
        return cls(upper,factors)

    def numerator(self,ranges):
        dimensions = self.upper.ndim-1
        if len(ranges)!=dimensions or any(not isinstance(lo,int) or not isinstance(hi,int)
                                or not 0<=lo<=hi for lo,hi in ranges):
            raise ValueError('Nonnegative target ranges matching the atlas dimension required')
        grid = self.upper.shape[1]
        residues = [np.arange(grid) if hi-lo+1>=grid else np.unique(np.arange(lo,hi+1)%grid)
                    for lo,hi in ranges]
        upper = np.max(self.upper[np.ix_(np.arange(4),*residues)],axis=tuple(range(1,dimensions+1)))
        return tuple(arb(float(value))*factor for value,factor in zip(upper,self.factors))


class SupportCounts:
    """Cache the exact generating value and both parity-binomial sums."""
    def __init__(self,pairs,markers):
        self.pairs = pairs
        self.markers = tuple(map(F,markers))
        self.value = transfer_value(pairs,self.markers)
        self.prefix = [0]
        for n in range(pairs+1):
            self.prefix.append(self.prefix[-1]+comb(pairs,n))

    def upper(self,box):
        if any(not isinstance(n,int) or not 0<=n<=self.pairs
               for interval in box.ranges for n in interval) or any(lo>hi for lo,hi in box.ranges):
            raise ValueError('Invalid support box')
        monomial = F(1)
        for marker,(lo,hi) in zip(self.markers,box.ranges):
            monomial *= marker**(lo if marker>=1 else hi)
        value = self.value/monomial
        upper = value.numerator//value.denominator
        for_parity = [self.prefix[hi+1]-self.prefix[lo] for lo,hi in (box.even,box.odd)]
        upper = min(upper,for_parity[0]*for_parity[1])
        if all(lo==hi for lo,hi in box.ranges):
            upper = min(upper,joint_count_upper(self.pairs,*(lo for lo,hi in box.ranges)))
        return upper


class DirectedAtlas:
    """One fixed rational witness, two untruncated regional alias tables."""
    def __init__(self,z,alpha,markers,*,grid=64,precision=256,bits=60,
                 count_markers=DEFAULT_COUNT_MARKERS,progress=None):
        self.z,self.alpha = F(z),F(alpha)
        self.markers = tuple(tuple(map(F,row)) for row in markers)
        if (not 0<self.z<1 or not 0<self.alpha<=1 or precision<128
                or len(self.markers)!=2 or any(len(row)!=3 or min(row)<=0 for row in self.markers)):
            raise ValueError('Invalid directed witness')
        self.precision,self.grid,self.bits = precision,grid,bits
        self.counts = SupportCounts(4096,count_markers)
        with ctx.workprec(precision):
            if progress:
                progress(dict(stage='atlas-single-coefficients'))
            single = general_coefficients(self.z,self.alpha,(F(1),)*3,bits=bits)
            if progress:
                progress(dict(stage='atlas-good-pair-coefficients'))
            paired = good_coefficients(self.z,self.alpha,bits)
            self.regions = []
            for parity,row in enumerate(self.markers):
                if progress:
                    progress(dict(stage='atlas-regional-transform',parity=parity))
                # No term is removed at cap 64: single-step degree <=32,
                # and the selective two-step term has degree <=64.
                power,vector,rho = regional_power(single,paired,row,(64,64,64),grid,progress)
                self.regions.append(RegionalAliases.from_power(power,vector,rho))
                del power

    def receipt(self,box):
        if box.even[0]+box.odd[0]==0:
            raise ValueError('Box must exclude the zero message')
        count = self.counts.upper(box)
        if count==0:
            raise ValueError('Empty support box')
        with ctx.workprec(self.precision):
            histograms = histogram_ranges(box)
            first,second = [region.numerator(ranges[1:]) for region,ranges in zip(self.regions,histograms)]
            paired = tuple(first[2*r]*second[s]+first[2*r+1]*second[2+s]
                           for r in range(2) for s in range(2))
            value = scalar_power(paired,128)
            scalar = scalar_envelope(4096,box,self.z,self.alpha,self.markers,precision=self.precision)
            bound = (arb(count)*(value[0]+value[1])).log()/arb(2).log()+scalar
            ceiling = rational_ceiling(bound)
        return dict(box=box.ranges,log2_upper_rational=str(ceiling),
                    log2_upper_enclosure=str(bound),support_count_upper=str(count),
                    directed=True,whole_code_certificate=False)

    def witness(self):
        return dict(z=str(self.z),alpha=str(self.alpha),markers=[list(map(str,row)) for row in self.markers],
                    count_markers=list(map(str,self.counts.markers)),grid=self.grid,
                    precision=self.precision,packet_mantissa_bits=self.bits,untruncated=True)


def certainly_empty(pairs,box):
    """Cheap necessary histogram conditions; false does not assert feasibility."""
    try:
        histogram_ranges(box,pairs)
    except ValueError:
        return True
    if all(lo==hi for lo,hi in box.ranges):
        return not joint_is_feasible(pairs,*(lo for lo,hi in box.ranges))
    return False


def cover(root,evaluate,*,pairs=4096,threshold=-80,max_evaluations=256):
    """Partition a declared box into accepted, empty, and unresolved leaves.

    Bisection gives disjoint children whose union is their parent. No
    sampled point can certify a box; evaluate must return a directed bound.
    The exact tree also accounts for unfinished work at the evaluation cap.
    """
    if (not isinstance(root,ProfileBox) or any(not isinstance(n,int) or not 0<=n<=pairs
            for interval in root.ranges for n in interval) or any(lo>hi for lo,hi in root.ranges)
            or root.even[0]+root.odd[0]==0 or max_evaluations<1):
        raise ValueError('Invalid nonzero root or evaluation budget')
    pending,accepted,empty,unresolved = [root],[],[],[]
    evaluations = 0
    while pending and evaluations<max_evaluations:
        box = pending.pop()
        if certainly_empty(pairs,box):
            empty.append(box.ranges)
            continue
        result = evaluate(box)
        evaluations += 1
        if result.get('directed') is not True or tuple(result.get('box',()))!=box.ranges:
            raise ValueError('Evaluator did not certify this box with directed arithmetic')
        upper = F(result['log2_upper_rational'])
        if upper<=threshold:
            accepted.append(result)
        elif all(lo==hi for lo,hi in box.ranges):
            unresolved.append(box.ranges)
        else:
            pending.extend(reversed(box.split()))
    unresolved.extend(box.ranges for box in pending)
    with ctx.workprec(256):
        total = sum((arb(2)**rational(F(item['log2_upper_rational'])) for item in accepted),arb(0))
        union = str(rational_ceiling(total.log()/arb(2).log())) if accepted else None
    return dict(root=root.ranges,accepted=accepted,empty=empty,unresolved=unresolved,
                evaluations=evaluations,complete_local_cover=not unresolved,
                accepted_union_log2_upper=union,directed=True,whole_code_certificate=False)


if __name__=='__main__':
    import argparse,json
    from hd_shifted_general_sweep import CHECKPOINTS
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--radius',type=int,default=4)
    parser.add_argument('--radii',type=int,nargs='+',help='Reuse one atlas for several local covers')
    parser.add_argument('--precision',type=int,default=256)
    parser.add_argument('--max-evaluations',type=int,default=256)
    parser.add_argument('--threshold',type=float,default=-80)
    parser.add_argument('--details',action='store_true',help='Print every accepted and unresolved leaf')
    args = parser.parse_args()
    radii = args.radii if args.radii is not None else [args.radius]
    if min(radii)<0:
        raise ValueError('Negative radius')
    profile,z,alpha,odds = CHECKPOINTS['asymmetric']
    markers = tuple(tuple(F(float(2**v)).limit_denominator(10**9) for v in row) for row in odds)
    atlas = DirectedAtlas(z,alpha,markers,precision=args.precision,
                          progress=lambda item:print(json.dumps(item),flush=True))
    for radius in radii:
        root = ProfileBox(*((max(0,n-radius),min(4096,n+radius)) for n in profile))
        result = cover(root,atlas.receipt,threshold=args.threshold,max_evaluations=args.max_evaluations)
        result['witness'] = atlas.witness()
        if not args.details:
            result['accepted_leaves'] = len(result.pop('accepted'))
            result['empty_leaves'] = len(result.pop('empty'))
            result['unresolved_leaves'] = len(result.pop('unresolved'))
        print(json.dumps(result),flush=True)
