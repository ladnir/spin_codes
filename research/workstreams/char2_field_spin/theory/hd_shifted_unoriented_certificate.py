"""Directed bounds that cover every lane balance at fixed (r,a,b).

For singleton count s and full-pair count f, divide each oriented local
coefficient C[u,s-u,f] by binomial(s,u), then take entrywise upper maxima.
This gives a two-variable coefficient envelope independent of lane counts.
It covers all indexed supports counted by profile_count_upper, not just
one even/odd split. It remains a partial, not whole-code, certificate.
"""
from fractions import Fraction as F
from itertools import product
from math import comb,factorial,prod

from flint import arb,ctx

from hd_shifted_alias_certificate import rational,LIMIT
from hd_shifted_general_certificate import general_coefficients
from hd_shifted_good_pair_certificate import good_coefficients,component_power
from hd_shifted_coverage import RegionalAliases,rational_ceiling
from hd_shifted_pairs_count import profile_count_upper
from hd_sumrank_paired_certificate import scalar_power


def collapse_coefficients(single):
    """Arb upper endpoints of max_u C[u,s-u,f]/binomial(s,u)."""
    result = {}
    for s in range(33):
        for f in range(33-s):
            candidates = [tuple(row[u*33**2+(s-u)*33+f]/comb(s,u) for row in single)
                          for u in range(s+1)]
            upper = tuple(max(values[entry].upper() for values in candidates) for entry in range(4))
            if any(not value.is_zero() for value in upper):
                result[s,f] = upper
    return result


def unsigned_components(single,paired,markers):
    markers = tuple(map(F,markers))
    if len(markers)!=2 or min(markers)<=0:
        raise ValueError('Positive singleton and full-pair markers required')
    powers = [[rational(marker)**degree for degree in range(65)] for marker in markers]
    good,dirty,composed = {},{},{}
    for (s,f),values in single.items():
        destination = good if s<=1 else dirty
        destination[s,f] = tuple(value*powers[0][s]*powers[1][f] for value in values)
    for (s,f),values in paired.items():
        composed[s,f] = tuple(value*powers[0][s]*powers[1][f] for value in values)
    return good,dirty,composed


def matrix_product(first,second):
    return tuple(first[2*i]*second[j]+first[2*i+1]*second[2+j]
                 for i in range(2) for j in range(2))


def positive_region(single,paired,markers):
    """The unnormalized regional matrix at positive real markers."""
    parts = unsigned_components(single,paired,markers)
    g,d,c = [tuple(sum((row[i] for row in values.values()),arb(0)) for i in range(4))
             for values in parts]
    two = tuple(sum(entries,arb(0)) for entries in
                zip(c,matrix_product(g,d),matrix_product(d,g),matrix_product(d,d)))
    return scalar_power(two,64)


def global_box_receipt(ranges,z,alpha,markers,*,precision=256,bits=60,progress=None):
    """Directed all-lane-balance union over an unsigned (r,a,b) box.

    Use the first matching's support count, then bound the resulting
    convex logarithmic expression by its eight corner values. The count
    of integer triples is a union factor, not a support multiplicity.
    """
    ranges = tuple(tuple(interval) for interval in ranges)
    z,alpha = F(z),F(alpha)
    markers = tuple(tuple(map(F,row)) for row in markers)
    if (len(ranges)!=3 or any(len(interval)!=2 or any(not isinstance(n,int) for n in interval)
                             or not 0<=interval[0]<=interval[1] for interval in ranges)
            or not 1<=ranges[0][0]<=ranges[0][1]<=8192 or ranges[1][1]>4096 or ranges[2][1]>4096
            or not 0<z<1 or not 0<alpha<=1 or precision<128 or len(markers)!=2
            or any(len(row)!=2 or min(row)<=0 for row in markers)):
        raise ValueError('Invalid unsigned box or rational witness')
    corners = list(product(*(tuple(dict.fromkeys(interval)) for interval in ranges)))
    histogram = lambda r,f:(4096-r+f,r-2*f,f)
    if any(min(histogram(r,f))<0 for r,a,b in corners for f in (a,b)):
        raise ValueError('Every histogram must be nonnegative throughout the box')
    with ctx.workprec(precision):
        if progress:
            progress(dict(stage='unsigned-global-coefficients'))
        single = collapse_coefficients(general_coefficients(z,alpha,(F(1),)*3,bits=bits))
        paired = good_coefficients(z,alpha,bits)
        regions = [positive_region(single,paired,row) for row in markers]
        run = scalar_power(matrix_product(*regions),128)
        log2 = lambda value:value.log()/arb(2).log()
        operator = log2(run[0]+run[1])
        log_factorials = {}
        def multinomial_log(hist):
            for n in (4096,)+hist:
                if n not in log_factorials:
                    log_factorials[n] = log2(arb(factorial(n)))
            return log_factorials[4096]-sum((log_factorials[n] for n in hist),arb(0))
        values = []
        for r,a,b in corners:
            histograms = (histogram(r,a),histogram(r,b))
            # N(r,a,b) <= Mult(H0)*2^(r-2a), including all lane balances.
            value = operator+multinomial_log(histograms[0])+histograms[0][1]
            value += rational(alpha)*(r*(2048-128*log2(arb(255)))
                                      -log2(arb(255))-LIMIT*log2(rational(z)))
            for hist,row in zip(histograms,markers):
                value -= 128*(multinomial_log(hist)+sum((n*log2(rational(mark))
                                                       for n,mark in zip(hist[1:],row)),arb(0)))
            values.append(value)
        count = prod(hi-lo+1 for lo,hi in ranges)
        upper = max(value.upper() for value in values)+log2(arb(count))
        ceiling = rational_ceiling(upper)
    return dict(ranges=ranges,integer_profile_triples=count,covers_all_even_odd_splits=True,
                log2_upper_rational=str(ceiling),log2_upper_enclosure=str(upper),
                corner_log2_upper=[str(rational_ceiling(v)) for v in values],
                directed=True,whole_code_certificate=False,
                witness=dict(z=str(z),alpha=str(alpha),markers=[list(map(str,row)) for row in markers],
                             precision=precision,packet_mantissa_bits=bits,coefficient_method='positive'))


class UnorientedAtlas:
    def __init__(self,z,alpha,markers,*,grid=128,precision=256,bits=60,progress=None):
        self.z,self.alpha = F(z),F(alpha)
        self.markers = tuple(tuple(map(F,row)) for row in markers)
        self.grid,self.precision,self.bits = grid,precision,bits
        if (not 0<self.z<1 or not 0<self.alpha<=1 or precision<128
                or len(self.markers)!=2 or any(len(row)!=2 or min(row)<=0 for row in self.markers)):
            raise ValueError('Invalid unoriented witness')
        with ctx.workprec(precision):
            if progress:
                progress(dict(stage='oriented-local-coefficients'))
            oriented = general_coefficients(self.z,self.alpha,(F(1),)*3,bits=bits)
            single = collapse_coefficients(oriented)
            del oriented
            if progress:
                progress(dict(stage='good-good-coefficients'))
            paired = good_coefficients(self.z,self.alpha,bits)
            self.regions = []
            for parity,row in enumerate(self.markers):
                if progress:
                    progress(dict(stage='unsigned-regional-transform',parity=parity,grid=grid))
                parts = unsigned_components(single,paired,row)
                powered,vector,rho = component_power(*parts,grid,progress,dimensions=2)
                self.regions.append(RegionalAliases.from_power(powered,vector,rho))

    def receipt(self,profile,pairs=4096):
        if pairs!=4096 or len(profile)!=3 or any(not isinstance(n,int) for n in profile):
            raise ValueError('Expected an integer K20 support profile (r,a,b)')
        r,a,b = profile
        count = profile_count_upper(pairs,r,a,b)
        histograms = tuple((pairs-r+f,r-2*f,f) for f in (a,b))
        if r<=0 or count<=0 or any(min(row)<0 for row in histograms):
            raise ValueError('Empty or zero-message profile')
        with ctx.workprec(self.precision):
            regions = []
            for atlas,histogram,markers in zip(self.regions,histograms,self.markers):
                values = atlas.numerator(tuple((n,n) for n in histogram[1:]))
                normalizer = arb(factorial(pairs))
                for n in histogram:
                    normalizer /= factorial(n)
                for n,marker in zip(histogram[1:],markers):
                    normalizer *= rational(marker)**n
                regions.append(tuple(value/normalizer for value in values))
            first,second = regions
            paired = tuple(first[2*i]*second[j]+first[2*i+1]*second[2+j]
                           for i in range(2) for j in range(2))
            value = scalar_power(paired,128)
            log2 = lambda value:value.log()/arb(2).log()
            prefix = rational(self.alpha)*(r*(2048-128*log2(arb(255)))
                                           -log2(arb(255))-LIMIT*log2(rational(self.z)))
            bound = log2(arb(count))+prefix+log2(value[0]+value[1])
            ceiling = rational_ceiling(bound)
            count_bits = str(log2(arb(count)))
        return dict(profile=profile,covers_all_even_odd_splits=True,
                    log2_upper_rational=str(ceiling),log2_upper_enclosure=str(bound),
                    log2_support_count_upper=count_bits,directed=True,whole_code_certificate=False,
                    witness=dict(z=str(self.z),alpha=str(self.alpha),
                                 markers=[list(map(str,row)) for row in self.markers],
                                 grid=self.grid,precision=self.precision,packet_mantissa_bits=self.bits))


if __name__=='__main__':
    import argparse,json
    from hd_shifted_general_sweep import CHECKPOINTS
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--grid',type=int,default=128)
    parser.add_argument('--precision',type=int,default=256)
    parser.add_argument('--coarse',choices=('middle','dense'),help='Replay a bulk positive-marker box')
    args = parser.parse_args()
    if args.coarse:
        if args.coarse=='middle':
            ranges = ((4032,4160),(1536,1664),(1536,1664))
            z,odds = .053587729568108665,(-.31797964871576334,.9361653072364948)
        else:
            ranges = ((6016,6272),(2372,2628),(2372,2628))
            z,odds = .012405036188127124,(1.8582139933645587,3.4036987243707566)
        alpha = .11372673628845778
        row = tuple(F(float(2**v)).limit_denominator(10**9) for v in odds)
        print(json.dumps(global_box_receipt(ranges,z,alpha,(row,row),precision=args.precision,
                         progress=lambda item:print(json.dumps(item),flush=True))),flush=True)
        raise SystemExit
    profile,z,alpha,odds = CHECKPOINTS['asymmetric']
    original = tuple(tuple(F(float(2**v)).limit_denominator(10**9) for v in row) for row in odds)
    markers = tuple((row[0]+row[1],row[2]) for row in original)
    atlas = UnorientedAtlas(z,alpha,markers,grid=args.grid,precision=args.precision,
                           progress=lambda item:print(json.dumps(item),flush=True))
    print(json.dumps(atlas.receipt((profile[0]+profile[1],profile[2],profile[3]))),flush=True)
