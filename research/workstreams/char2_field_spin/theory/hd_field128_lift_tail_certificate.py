"""Directed positive-marker, common-potential, box, and terminal bounds.

The setup is the genuine full-field lift. Receipts cover only their stated
weights or boxes. All input witnesses are exact rationals, all arithmetic
bounds use Arb, and no old GF256 numerical ceiling is imported.
"""

import ast
from fractions import Fraction as F
from hashlib import sha256
from itertools import product
from math import comb
from pathlib import Path

from flint import arb,arb_poly,ctx

import hd_field128_lift_certificate as base


def source_pins():
    directory = Path(__file__).resolve().parent
    pending,seen = [Path(__file__).name],{}
    while pending:
        name = pending.pop()
        if name in seen:
            continue
        path = directory/name
        raw = path.read_bytes()
        seen[name] = sha256(raw).hexdigest()
        for node in ast.walk(ast.parse(raw,filename=str(path))):
            modules = ([a.name for a in node.names] if isinstance(node,ast.Import)
                       else [node.module] if isinstance(node,ast.ImportFrom) and node.module else [])
            for module in modules:
                target = module.split('.')[0]+'.py'
                if (directory/target).is_file() and target not in seen:
                    pending.append(target)
    return dict(sorted(seen.items()))


def log2arb(value):
    return value.log()/arb(2).log()


def full_counts(pairs,weight):
    if type(pairs) is not int or type(weight) is not int or pairs<1 or not 0<=weight<=2*pairs:
        raise ValueError('Invalid integer matching size or weight')
    return range(max(0,weight-pairs),weight//2+1)


def one_pairing_count(pairs,weight,full):
    single = weight-2*full
    if full not in full_counts(pairs,weight):
        return 0
    return comb(pairs,full)*comb(pairs-full,single)*(1<<single)


def circular_run_counts(length,weight):
    result = [0]*(weight+1)
    if weight in (0,length):
        result[0] = 1
        return result
    result[1] = length
    for runs in range(1,min(weight,length-weight)):
        numerator = result[runs]*(weight-runs)*(length-weight-runs)
        result[runs+1],remainder = divmod(numerator,runs*(runs+1))
        if remainder:
            raise ArithmeticError('Nonintegral cycle count recurrence')
    return result


def contraction(entries,ratio):
    """Common state vector (1,ratio), with a rigorous contraction ceiling."""
    ratio = base.exact(ratio)
    if ratio<=0:
        raise ValueError('A positive exact state ratio is required')
    v = base.rational(ratio)
    a,b,c,d = entries
    return max(base.upper(a+b*v),base.upper(c/v+d))


def support_sum(pairs,weight,regions,ratio,*,method='marginal',repetitions=256,budget=None):
    """Every joint matching profile, using a common potential across them."""
    ratio = base.exact(ratio)
    counts = full_counts(pairs,weight)
    if (type(repetitions) is not int or repetitions<2 or repetitions%2
            or set(regions)!=set(counts) or ratio<=0):
        raise ValueError('Complete histogram range and positive even repetitions required')
    budget = budget or base.Budget()
    contractions = {f:contraction(regions[f],ratio) for f in counts}
    value = arb(0)
    if method=='marginal':
        for f in counts:
            budget.check()
            value += one_pairing_count(pairs,weight,f)*contractions[f]**repetitions
    elif method=='runs':
        coefficients = [arb(0)]*(weight//2+1)
        for f in counts:
            budget.check()
            coefficients[f] = base.upper(contractions[f]**(repetitions//2))
        convolution = arb_poly(coefficients)**2
        by_runs = circular_run_counts(2*pairs,weight)
        for degree in range(2*counts.start,2*(counts.stop-1)+1):
            budget.check()
            if by_runs[weight-degree]:
                value += by_runs[weight-degree]*base.upper(convolution[degree])
    else:
        raise ValueError('Support method must be marginal or runs')
    return base.upper(value*base.rational(max(F(1),1/ratio)))


def validate_box(ranges,pairs=4096):
    ranges = tuple(tuple(row) for row in ranges)
    if (len(ranges)!=3 or any(len(row)!=2 or any(type(n) is not int for n in row)
                            or not 0<=row[0]<=row[1] for row in ranges)
            or not 1<=ranges[0][0]<=ranges[0][1]<=2*pairs
            or ranges[1][1]>pairs or ranges[2][1]>pairs):
        raise ValueError('Invalid unsigned support box')
    return ranges


def histogram_vertices(ranges,pairs=4096):
    (rl,rh),(al,ah),(bl,bh) = validate_box(ranges,pairs)
    lo,hi = max(rl,2*al,2*bl),min(rh,pairs+ah,pairs+bh)
    if lo>hi:
        return ()
    breaks = sorted({lo,hi}|{r for r in (pairs+al,2*ah,pairs+bl,2*bh) if lo<=r<=hi})
    result = set()
    for r in breaks:
        aa = F(max(al,r-pairs)),min(F(ah),F(r,2))
        bb = F(max(bl,r-pairs)),min(F(bh),F(r,2))
        result.update((F(r),a,b) for a,b in product(aa,bb))
    return tuple(sorted(result))


def profile_count(ranges,pairs=4096):
    (rl,rh),(al,ah),(bl,bh) = validate_box(ranges,pairs)
    result = 0
    for r in range(max(rl,2*al,2*bl),min(rh,pairs+ah,pairs+bh)+1):
        alo,ahi = max(al,r-pairs),min(ah,r//2)
        blo,bhi = max(bl,r-pairs),min(bh,r//2)
        if alo<=ahi and blo<=bhi:
            result += (ahi-alo+1)*(bhi-blo+1)
            if r<2*pairs and r%2==0 and alo<=r//2<=ahi and blo<=r//2<=bhi:
                result -= 1
    return result


def evaluate_polynomial(coefficients,markers):
    x,y = map(base.rational,markers)
    result = [arb(0)]*4
    for (s,f),values in coefficients.items():
        monomial = x**s*y**f
        for i in range(4):
            result[i] += values[i]*monomial
    return tuple(base.upper(value) for value in result)


def positive_pair(single,composed,markers):
    """GG was composed before alpha; GD, DG, DD are positive products."""
    if composed is None:
        one = evaluate_polynomial(single,markers)
        return base.matrix_product(one,one)
    good = evaluate_polynomial({k:v for k,v in single.items() if k[0]<=1},markers)
    dirty = evaluate_polynomial({k:v for k,v in single.items() if k[0]>=2},markers)
    pieces = (evaluate_polynomial(composed,markers),base.matrix_product(good,dirty),
              base.matrix_product(dirty,good),base.matrix_product(dirty,dirty))
    return tuple(base.upper(sum((row[i] for row in pieces),arb(0))) for i in range(4))


class PositiveBank:
    """One exact witness and marker bank, reusable over every source weight."""

    def __init__(self,z,alpha,markers,*,q=base.FIELD_ORDER,distance=F(7,20),
                 precision=256,bits=80,outer_measure='invertible-toeplitz',max_seconds=120.):
        self.z,self.alpha,self.distance = map(base.exact,(z,alpha,distance))
        self.q,self.field_bits = q,base.checked_q(q)
        self.markers = tuple(dict.fromkeys(tuple(map(base.exact,row)) for row in markers))
        if (not 0<self.z<1 or not 0<self.alpha<=1 or not 0<self.distance<F(1,2)
                or not 1<=len(self.markers)<=64 or any(len(row)!=2 or min(row)<=0 for row in self.markers)
                or type(precision) is not int or not 128<=precision<=4096
                or type(bits) is not int or not 32<=bits<=512):
            raise ValueError('Invalid exact bank witness')
        self.precision,self.bits,self.outer_measure = precision,bits,outer_measure
        self.forbidden = int(self.distance*base.OUTPUT_LENGTH)
        self.pins = source_pins()
        budget = base.Budget(max_seconds)
        with ctx.workprec(precision):
            # A step has at most64 active positions; a good-good pair has128.
            # Retaining maximum128 here means no pair coefficient is omitted.
            single,composed = base.directed_parts(self.z,self.alpha,128,q,bits,budget,True)
            self.numerators = []
            for row in self.markers:
                budget.check()
                self.numerators.append(base.matrix_power(positive_pair(single,composed,row),64))
            self.density_log = log2arb(base.rational(base.outer_density(q,outer_measure)))
            self.slope = base.rational(self.alpha)*self.density_log
            self.constant = -base.rational(self.alpha)*(log2arb(arb(q-1))
                                                      +self.forbidden*log2arb(base.rational(self.z)))
            self.log_markers = tuple(tuple(log2arb(base.rational(x)) for x in row) for row in self.markers)
        budget.check()
        self.check_pins()

    def check_pins(self):
        if source_pins()!=self.pins:
            raise RuntimeError('A pinned transitive helper changed during the replay')

    def regional(self,weight):
        if type(weight) is not int or not 1<=weight<=base.GROUPS:
            raise ValueError('Nonzero source weight in1..8192 required')
        result = {}
        for full in full_counts(base.PAIRS,weight):
            single = weight-2*full
            count = comb(base.PAIRS,full)*comb(base.PAIRS-full,single)
            candidates = []
            for markers,numerator in zip(self.markers,self.numerators):
                denominator = count*base.rational(markers[0])**single*base.rational(markers[1])**full
                candidates.append(tuple(base.upper(entry/denominator) for entry in numerator))
            result[full] = tuple(min(row[i] for row in candidates) for i in range(4))
        return result

    def metadata(self):
        return dict(construction=base.CONSTRUCTION,field_order=self.q,field_bits=self.field_bits,
                    output_length=base.OUTPUT_LENGTH,source_dimension=base.DIMENSION,outer_groups=base.GROUPS,
                    forbidden_weight=self.forbidden,target_distance=str(self.distance),outer_measure=self.outer_measure,
                    pair_swaps=True,state_persists_across_regions=True,full_field_random_coefficients=True,
                    projective_divisor=str(self.q-1),directed=True,whole_code_certificate=False,
                    source_sha256=self.pins,pin_scope='complete recursive static local-import closure, including lazy imports',
                    witness=dict(z=str(self.z),alpha=str(self.alpha),markers=[list(map(str,row)) for row in self.markers],
                                 precision=self.precision,packet_mantissa_bits=self.bits),
                    coefficient_method='positive-rational-marker-evaluation-with-Arb')

    def weight_receipt(self,weight,ratio,*,method='marginal',max_seconds=60.):
        ratio = base.exact(ratio)
        budget = base.Budget(max_seconds)
        with ctx.workprec(self.precision):
            regions = self.regional(weight)
            support = support_sum(base.PAIRS,weight,regions,ratio,method=method,budget=budget)
            value = self.constant+self.slope*weight+log2arb(support)
            ceiling = base.rational_ceiling(value)
            enclosure = str(value)
        budget.check()
        self.check_pins()
        return dict(self.metadata(),active_groups=weight,regional_histograms=len(regions),
                    state_potential_ratio=str(ratio),support_method=method,
                    log2_upper_rational=str(ceiling),log2_upper_enclosure=enclosure,
                    covers_all_supports_at_weight=True,covers_all_even_odd_splits=True,
                    meets_40_bit_selected_weight_target=ceiling<=-40)

    def scalar(self,point,indices,count_matching=None):
        r,a,b = point
        value = self.constant+self.slope*base.rational(r)
        for parity,(full,index) in enumerate(zip((a,b),indices)):
            hist = base.PAIRS-r+full,r-2*full,full
            multinomial_log = (arb(base.PAIRS+1).lgamma()
                               -sum(((base.rational(n)+1).lgamma() for n in hist),arb(0)))/arb(2).log()
            value -= 128*(multinomial_log+sum((base.rational(n)*marker for n,marker in
                                              zip(hist[1:],self.log_markers[index])),arb(0)))
            if parity==count_matching:
                value += multinomial_log+base.rational(hist[1])
        return value

    def box_receipt(self,ranges,*,marker_indices=(0,0),max_seconds=60.):
        ranges = validate_box(ranges)
        count = profile_count(ranges)
        if (not count or len(marker_indices)!=2 or any(type(i) is not int or not 0<=i<len(self.markers)
                                                      for i in marker_indices)):
            raise ValueError('A nonempty box and two valid fixed marker indices are required')
        budget = base.Budget(max_seconds)
        vertices = histogram_vertices(ranges)
        with ctx.workprec(self.precision):
            matrices = tuple(self.numerators[i] for i in marker_indices)
            run = base.matrix_power(base.matrix_product(*matrices),128)
            envelopes = []
            for matching in (0,1):
                values = []
                for point in vertices:
                    budget.check()
                    values.append(self.scalar(point,marker_indices,matching).upper())
                envelopes.append(max(values))
            value = log2arb(run[0]+run[1])+min(envelopes)+log2arb(arb(count))
            ceiling = base.rational_ceiling(value)
            enclosure = str(value)
        self.check_pins()
        return dict(self.metadata(),ranges=ranges,integer_profile_triples=count,convex_vertices=len(vertices),
                    marker_indices=marker_indices,support_count_method='matching-profile-convex-envelope',
                    log2_upper_rational=str(ceiling),log2_upper_enclosure=enclosure,
                    covers_entire_clipped_box=True,covers_all_even_odd_splits=True)


def terminal_receipt(*,q=base.FIELD_ORDER,distance=F(7,20),first=None,failure_bits=40,precision=256):
    """Simultaneously cover all weights >= first via unique inverse images."""
    field_bits,distance = base.checked_q(q),base.exact(distance)
    if (not 0<distance<F(1,2) or type(failure_bits) is not int or failure_bits<0
            or type(precision) is not int or not 128<=precision<=4096):
        raise ValueError('Invalid terminal target')
    limit = int(distance*base.OUTPUT_LENGTH)
    if not limit:
        raise ValueError('Positive forbidden weight required')
    pins = source_pins()
    mark = F(limit,(base.OUTPUT_LENGTH-limit)*(q-1))
    with ctx.workprec(precision):
        z = base.rational(mark)
        logqm1 = log2arb(arb(q-1))
        ball = base.OUTPUT_LENGTH*log2arb(1+(q-1)*z)-limit*log2arb(z)
        loss = 96*field_bits+32*logqm1
        if first is None:
            first = max(1,int(((ball-logqm1+failure_bits)/loss).upper().ceil().fmpz()))
        if type(first) is not int or not 1<=first<=base.GROUPS:
            raise ValueError('Terminal first weight must lie in1..8192')
        value = ball-first*loss-logqm1
        ceiling = base.rational_ceiling(value)
        enclosure = str(value)
    if source_pins()!=pins:
        raise RuntimeError('A pinned helper changed during the replay')
    return dict(construction=base.CONSTRUCTION,field_order=q,field_bits=field_bits,
                weight_range=[first,base.GROUPS],output_length=base.OUTPUT_LENGTH,
                source_dimension=base.DIMENSION,forbidden_weight=limit,target_distance=str(distance),
                witness=dict(mark=str(mark),precision=precision),outer_measure='invertible-toeplitz',
                projective_divisor=str(q-1),log2_upper_rational=str(ceiling),log2_upper_enclosure=enclosure,
                covers_all_supports_in_range=True,covers_all_even_odd_splits=True,
                no_support_subset_union=True,method='inverse-image-terminal-range',
                meets_requested_range_target=ceiling<=-failure_bits,requested_failure_bits=failure_bits,
                source_sha256=pins,pin_scope='complete recursive static local-import closure, including lazy imports',
                directed=True,whole_code_certificate=False)


def main():
    import argparse
    import json
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--terminal',action='store_true')
    parser.add_argument('--first',type=int)
    parser.add_argument('--r',type=int)
    parser.add_argument('--box',type=int,nargs=6,metavar=('RLO','RHI','ALO','AHI','BLO','BHI'))
    parser.add_argument('--marker-indices',type=int,nargs=2,default=(0,0))
    parser.add_argument('--z',default='1/2')
    parser.add_argument('--alpha',default='1/128')
    parser.add_argument('--marker',nargs=2,action='append',metavar=('SINGLE','FULL'))
    parser.add_argument('--ratio',default='1/256')
    parser.add_argument('--method',choices=('marginal','runs'),default='marginal')
    parser.add_argument('--field-bits',type=int,default=128)
    parser.add_argument('--distance',default='7/20')
    parser.add_argument('--precision',type=int,default=256)
    parser.add_argument('--bits',type=int,default=80)
    parser.add_argument('--max-seconds',type=float,default=120.)
    parser.add_argument('--output',type=Path)
    args = parser.parse_args()
    if args.output and args.output.exists():
        parser.error('Output exists; choose an unused filename')
    if not 8<=args.field_bits<=128 or args.field_bits%8:
        parser.error('--field-bits must be a multiple of eight in8..128')
    if sum((args.terminal,args.r is not None,args.box is not None))!=1:
        parser.error('Select exactly one of --terminal, --r, --box')
    if args.terminal:
        result = terminal_receipt(q=1<<args.field_bits,distance=F(args.distance),first=args.first,precision=args.precision)
    else:
        budget = base.Budget(args.max_seconds)
        bank = PositiveBank(F(args.z),F(args.alpha),args.marker or (('1/16','1/16'),),
                            q=1<<args.field_bits,distance=F(args.distance),precision=args.precision,
                            bits=args.bits,max_seconds=args.max_seconds)
        remaining = args.max_seconds-(base.monotonic()-budget.started)
        budget.check()
        if args.box is not None:
            ranges = tuple(tuple(args.box[i:i+2]) for i in (0,2,4))
            result = bank.box_receipt(ranges,marker_indices=tuple(args.marker_indices),max_seconds=remaining)
        else:
            result = bank.weight_receipt(args.r,F(args.ratio),method=args.method,max_seconds=remaining)
    encoded = json.dumps(result,indent=2,allow_nan=False)
    if args.output:
        with args.output.open('x',encoding='utf-8') as stream:
            stream.write(encoded+'\n')
    print(encoded)


if __name__=='__main__':
    main()
