"""Cached directed bounds on unsigned support boxes, including boundaries.

The box is intersected with the nonnegative-histogram domain. Convex
normalizers need at most 24 rational vertices, not a staircase of boxes.
All receipts cover every lane balance, but not the whole code by default.
"""
from fractions import Fraction as F
from functools import lru_cache
from itertools import product

from flint import arb,ctx

from hd_shifted_alias_certificate import rational,LIMIT
from hd_shifted_coverage import RegionalAliases,rational_ceiling
from hd_shifted_general_certificate import general_coefficients
from hd_shifted_good_pair_certificate import good_coefficients,component_power
from hd_shifted_unoriented_certificate import (collapse_coefficients,unsigned_components,
    positive_region,matrix_product)
from hd_shifted_pairs_count import profile_count_upper
from hd_sumrank_paired_certificate import scalar_power

CONSTRUCTION = 'toeplitz-gf256-h4-alternating-pairs-k20-rate-half'


def validate_box(ranges,pairs=4096):
    ranges = tuple(tuple(row) for row in ranges)
    if (len(ranges)!=3 or any(len(row)!=2 or any(not isinstance(n,int) for n in row)
                             or not 0<=row[0]<=row[1] for row in ranges)
            or not 1<=ranges[0][0]<=ranges[0][1]<=2*pairs
            or ranges[1][1]>pairs or ranges[2][1]>pairs):
        raise ValueError('Invalid nonzero unsigned support box')
    return ranges


def histogram_vertices(ranges,pairs=4096):
    """Vertices sufficient for a convex maximum on box intersect hist>=0."""
    (rl,rh),(al,ah),(bl,bh) = validate_box(ranges,pairs)
    lo,hi = max(rl,2*al,2*bl),min(rh,pairs+ah,pairs+bh)
    if lo>hi:
        return ()
    breaks = sorted({lo,hi}|{r for r in (pairs+al,2*ah,pairs+bl,2*bh) if lo<=r<=hi})
    result = set()
    for r in breaks:
        aa = (F(max(al,r-pairs)),min(F(ah),F(r,2)))
        bb = (F(max(bl,r-pairs)),min(F(bh),F(r,2)))
        for a,b in product(aa,bb):
            result.add((F(r),a,b))
    return tuple(sorted(result))


@lru_cache(16384)
def profile_count(ranges,pairs=4096):
    """Count feasible integer triples in this box, without lane balances."""
    (rl,rh),(al,ah),(bl,bh) = validate_box(ranges,pairs)
    result = 0
    for r in range(max(rl,2*al,2*bl),min(rh,pairs+ah,pairs+bh)+1):
        alo,ahi = max(al,r-pairs),min(ah,r//2)
        blo,bhi = max(bl,r-pairs),min(bh,r//2)
        if alo<=ahi and blo<=bhi:
            result += (ahi-alo+1)*(bhi-blo+1)
            # A nonconstant cycle cannot be full under both matchings.
            if r<2*pairs and r%2==0 and alo<=r//2<=ahi and blo<=r//2<=bhi:
                result -= 1
    return result


@lru_cache(8)
def unmarked(z,alpha,precision,bits):
    with ctx.workprec(precision):
        single = collapse_coefficients(general_coefficients(z,alpha,(F(1),)*3,bits=bits))
        paired = good_coefficients(z,alpha,bits)
    return single,paired


class CachedBound:
    """Fixed exact witness, with positive or positive-alias numerators."""
    construction = CONSTRUCTION

    @classmethod
    def unmarked_coefficients(cls,z,alpha,precision,bits):
        """Coefficient-provider hook; the original construction is unchanged."""
        return unmarked(z,alpha,precision,bits)

    def __init__(self,z,alpha,markers,*,precision=256,bits=60,grid=None,uniform_aliases=False,
                 run_markers=None,progress=None):
        self.z,self.alpha = F(z),F(alpha)
        self.markers = tuple(tuple(map(F,row)) for row in markers)
        if (not 0<self.z<1 or not 0<self.alpha<=1 or precision<128 or len(self.markers)!=2
                or any(len(row)!=2 or min(row)<=0 for row in self.markers)):
            raise ValueError('Invalid directed witness')
        if uniform_aliases and grid is None:
            raise ValueError('Uniform alias coefficients require an explicit grid')
        self.precision,self.bits,self.grid = precision,bits,grid
        self.uniform_aliases = bool(uniform_aliases)
        self.run_markers = None if run_markers is None else tuple(map(F,run_markers))
        if self.run_markers is not None and (len(self.run_markers)!=2 or min(self.run_markers)<=0):
            raise ValueError('Two positive run-count markers required')
        self.log_factorials = {}
        with ctx.workprec(precision):
            if progress:
                progress(dict(stage='cached-unmarked-coefficients'))
            single,paired = self.unmarked_coefficients(self.z,self.alpha,precision,bits)
            self.regions = []
            for parity,row in enumerate(self.markers):
                if progress:
                    progress(dict(stage='cached-region',parity=parity,grid=grid))
                if grid is None:
                    self.regions.append(positive_region(single,paired,row))
                else:
                    pieces = unsigned_components(single,paired,row)
                    power,vector,rho = component_power(*pieces,grid,dimensions=2)
                    atlas = RegionalAliases.from_power(power,vector,rho)
                    self.regions.append(atlas.numerator(((0,grid-1),(0,grid-1)))
                                        if self.uniform_aliases else atlas)
            log2 = lambda value:value.log()/arb(2).log()
            self.log_markers = tuple(tuple(log2(rational(m)) for m in row) for row in self.markers)
            self.slope = rational(self.alpha)*(2048-128*log2(arb(255)))
            self.constant = -rational(self.alpha)*(log2(arb(255))+LIMIT*log2(rational(self.z)))

    def log_factorial(self,n):
        n = F(n)
        if n not in self.log_factorials:
            self.log_factorials[n] = (rational(n)+1).lgamma()/arb(2).log()
        return self.log_factorials[n]

    def scalar(self,point,count_matching=None):
        """Profile scalar, with optional one-matching support count."""
        r,a,b = point
        histograms = tuple((4096-r+f,r-2*f,f) for f in (a,b))
        value = self.constant+self.slope*rational(r)
        for parity,(hist,row) in enumerate(zip(histograms,self.log_markers)):
            log_mult = self.log_factorial(4096)-sum((self.log_factorial(n) for n in hist),arb(0))
            value -= 128*(log_mult+sum((rational(n)*mark for n,mark in zip(hist[1:],row)),arb(0)))
            if parity==count_matching:
                value += log_mult+rational(hist[1])
        return value

    def receipt(self,ranges):
        ranges = validate_box(ranges)
        count = profile_count(ranges)
        if not count:
            raise ValueError('Empty unsigned support box')
        vertices = histogram_vertices(ranges)
        with ctx.workprec(self.precision):
            if self.grid is None or self.uniform_aliases:
                regions = self.regions
            else:
                regions = []
                for parity,region in enumerate(self.regions):
                    s = [r-2*p[parity+1] for p in vertices for r in (p[0],)]
                    f = [p[parity+1] for p in vertices]
                    # Round outward before requesting integer coefficient aliases.
                    floor = lambda v:v.numerator//v.denominator
                    ceil = lambda v:-((-v.numerator)//v.denominator)
                    target = ((floor(min(s)),ceil(max(s))),(floor(min(f)),ceil(max(f))))
                    regions.append(region.numerator(target))
            run = scalar_power(matrix_product(*regions),128)
            log2 = lambda value:value.log()/arb(2).log()
            operator = log2(run[0]+run[1])
            # Both one-matching counts are valid; minimize two complete envelopes.
            envelope = min(max(self.scalar(p,matching).upper() for p in vertices)
                           for matching in (0,1))
            bound = operator+envelope+log2(arb(count))
            count_method = 'matching-profile-envelope'
            if self.run_markers is not None:
                from hd_shifted_run_count import box_log2_upper,transfer_value
                support_log = box_log2_upper(4096,ranges,*self.run_markers,precision=self.precision)
                # This counts ALL indexed supports in the box, not profiles.
                # Do not multiply by the number of profile triples again.
                run_bound = operator+max(self.scalar(p).upper() for p in vertices)+support_log
                if run_bound.upper()<bound.upper():
                    bound = run_bound
                    count_method = 'whole-box-cyclic-runs'
                x,y = self.run_markers
                logx,logy = log2(rational(x)),log2(rational(y))
                # Count and scalar tilts must be combined before maximizing:
                # sum_U 2^S(U) <= G(x,y) max_U 2^S(U)/(x^r y^(r-a-b)).
                # Subtracting an affine form preserves convexity.
                tilted = max((self.scalar(p)-rational(p[0])*logx
                              -rational(p[0]-p[1]-p[2])*logy).upper() for p in vertices)
                weighted_bound = operator+log2(rational(transfer_value(4096,x,y)))+tilted
                if weighted_bound.upper()<bound.upper():
                    bound = weighted_bound
                    count_method = 'tilted-cyclic-run-sum'
            if all(lo==hi for lo,hi in ranges):
                point = tuple(lo for lo,hi in ranges)
                exact_upper = profile_count_upper(4096,*point)
                if exact_upper:
                    point_bound = operator+self.scalar(point)+log2(arb(exact_upper))
                    if point_bound.upper()<bound.upper():
                        bound = point_bound
                        count_method = 'joint-profile-count'
            ceiling = rational_ceiling(bound)
        return dict(ranges=ranges,integer_profile_triples=count,convex_vertices=len(vertices),
                    log2_upper_rational=str(ceiling),log2_upper_enclosure=str(bound),
                    support_count_method=count_method,
                    directed=True,covers_all_even_odd_splits=True,whole_code_certificate=False,
                    construction=self.construction,field_order=256,forbidden_weight=LIMIT)

    def witness(self):
        return dict(z=str(self.z),alpha=str(self.alpha),markers=[list(map(str,row)) for row in self.markers],
                    precision=self.precision,packet_mantissa_bits=self.bits,grid=self.grid,
                    uniform_aliases=self.uniform_aliases,
                    run_markers=None if self.run_markers is None else list(map(str,self.run_markers)),
                    coefficient_method='positive' if self.grid is None else 'positive-alias')
