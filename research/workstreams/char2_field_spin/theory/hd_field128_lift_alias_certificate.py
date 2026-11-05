"""Directed residue-alias banks for genuine full-field paired routing.

Exact q-dependent packet caps feed the retained interval radix-2 DFT.
The binary64 transform uses outward-rounded intervals, never FLOAT FFT
coefficients or an empirical cushion. Arb handles coefficient preparation,
normalization, similarity witnesses, restoration, and probability bounds.
"""

import ast
from fractions import Fraction as F
from hashlib import sha256
from pathlib import Path

import numpy as np
from flint import arb,ctx

import hd_field128_lift_certificate as base
import hd_field128_lift_tail_certificate as tail
from hd_shifted_alias_certificate import interval_fft,interval_matrix_power,_cadd,_mul,_finite
from hd_shifted_good_pair_certificate import interval_polynomial,matrix_product as interval_matrix_product
from hd_shifted_coverage import RegionalAliases


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


def components(single,composed,markers):
    powers = tuple(tuple(base.rational(mark)**n for n in range(65)) for mark in markers)
    good,dirty,paired = {},{},{}
    for (s,f),values in single.items():
        destination = dirty if composed is None or s>=2 else good
        destination[s,f] = tuple(base.upper(x*powers[0][s]*powers[1][f]) for x in values)
    if composed is not None:
        for (s,f),values in composed.items():
            paired[s,f] = tuple(base.upper(x*powers[0][s]*powers[1][f]) for x in values)
    return good,dirty,paired


def similarity(positive):
    """Stable Arb Perron proposal, verified by directed row inequalities."""
    a,b,c,d = tuple(base.upper(x) for x in positive)
    if not b>0 or not c>0:
        raise ArithmeticError('Positive off-diagonal entries required for similarity')
    delta = a-d
    radical = (delta*delta+4*b*c).sqrt()
    ratio = 2*c/(delta+radical) if a>=d else (radical-delta)/(2*b)
    ratio = ratio.mid()
    if not ratio>0 or not ratio.is_finite():
        raise ArithmeticError('Invalid positive similarity ratio')
    vector = (arb(1),ratio)
    row_cap = max(base.upper(a+b*ratio),base.upper(c/ratio+d))
    # A small exact inflation makes the final inequality verifiable even
    # when squaring the rounded square-root endpoint incurs new rounding.
    rho = base.upper(row_cap.sqrt().upper()*(1+arb(2)**(-ctx.prec//2)))
    if not rho*rho>=row_cap:
        raise ArithmeticError('Similarity scale fails the directed row bound')
    return vector,rho


def component_alias(good,dirty,composed,grid,budget):
    if type(grid) is not int or grid<64 or grid>256 or grid&(grid-1):
        raise ValueError('Grid must be64,128,or256')
    sums = lambda entries:tuple(sum((row[i] for row in entries.values()),arb(0)) for i in range(4))
    g,d,c = map(sums,(good,dirty,composed))
    positive = tuple(sum(row,arb(0)) for row in zip(c,base.matrix_product(g,d),
                                                    base.matrix_product(d,g),base.matrix_product(d,d)))
    vector,rho = similarity(positive)
    single_factor = tuple(vector[s]/vector[r]/rho for r in range(2) for s in range(2))
    pair_factor = tuple(x/rho for x in single_factor)
    with np.errstate(over='raise',invalid='raise',divide='raise',under='ignore'):
        transformed = []
        for coefficients,factors in ((good,single_factor),(dirty,single_factor),(composed,pair_factor)):
            budget.check()
            transformed.append(interval_fft(interval_polynomial(coefficients,factors,grid,dimensions=2),(1,2)))
        g,d,c = transformed
        values = _cadd(c,_cadd(interval_matrix_product(g,d),
                      _cadd(interval_matrix_product(d,g),interval_matrix_product(d,d))))
        _finite(values)
        budget.check()
        powered = interval_matrix_power(values,64)
        budget.check()
        atlas = RegionalAliases.from_power(powered,vector,rho)
    # Store the positive-evaluation comparator as an independent fallback.
    # Taking its minimum with an interval alias remains a valid bound.
    fallback = base.matrix_power(positive,64)
    return atlas,fallback,dict(similarity_vector=list(map(str,vector)),rho=str(rho),
                               interval_arithmetic='outward-rounded-binary64-with-Arb-checked-conversions')


class AliasBank(tail.PositiveBank):
    """Directed normalized atlas minima, one common potential per weight."""

    def __init__(self,z,alpha,markers,*,q=base.FIELD_ORDER,distance=F(7,20),grid=128,
                 precision=256,bits=80,outer_measure='invertible-toeplitz',max_seconds=180.,progress=None):
        self.z,self.alpha,self.distance = map(base.exact,(z,alpha,distance))
        self.q,self.field_bits = q,base.checked_q(q)
        self.markers = tuple(dict.fromkeys(tuple(map(base.exact,row)) for row in markers))
        if (not 0<self.z<1 or not 0<self.alpha<=1 or not 0<self.distance<F(1,2)
                or not 1<=len(self.markers)<=128 or any(len(row)!=2 or min(row)<=0 for row in self.markers)
                or type(precision) is not int or not 128<=precision<=4096
                or type(bits) is not int or not 32<=bits<=512
                or type(grid) is not int or grid<64 or grid>256 or grid&(grid-1)):
            raise ValueError('Invalid exact alias-bank witness')
        self.precision,self.bits,self.outer_measure,self.grid = precision,bits,outer_measure,grid
        self.forbidden = int(self.distance*base.OUTPUT_LENGTH)
        self.pins = source_pins()
        budget = base.Budget(max_seconds)
        with ctx.workprec(precision):
            if progress:
                progress(dict(stage='exact-full-field-coefficients'))
            single,composed = base.directed_parts(self.z,self.alpha,128,q,bits,budget,True)
            self.atlases,self.numerators,self.diagnostics = [],[],[]
            for index,row in enumerate(self.markers):
                budget.check()
                if progress:
                    progress(dict(stage='directed-interval-alias',marker=index,grid=grid))
                atlas,fallback,diagnostic = component_alias(*components(single,composed,row),grid,budget)
                self.atlases.append(atlas)
                self.numerators.append(fallback)
                self.diagnostics.append(diagnostic)
            self.density_log = tail.log2arb(base.rational(base.outer_density(q,outer_measure)))
            self.slope = base.rational(self.alpha)*self.density_log
            self.constant = -base.rational(self.alpha)*(tail.log2arb(arb(q-1))
                           +self.forbidden*tail.log2arb(base.rational(self.z)))
            self.log_markers = tuple(tuple(tail.log2arb(base.rational(x)) for x in row) for row in self.markers)
        budget.check()
        self.check_pins()

    def check_pins(self):
        if source_pins()!=self.pins:
            raise RuntimeError('A pinned transitive alias helper changed during replay')

    def numerator(self,index,ranges):
        alias = self.atlases[index].numerator(ranges)
        return tuple(min(base.upper(x),base.upper(y)) for x,y in zip(alias,self.numerators[index]))

    def regional(self,weight):
        if type(weight) is not int or not 1<=weight<=base.GROUPS:
            raise ValueError('Nonzero source weight in1..8192 required')
        result = {}
        for full in tail.full_counts(base.PAIRS,weight):
            single = weight-2*full
            count = base.comb(base.PAIRS,full)*base.comb(base.PAIRS-full,single)
            candidates = []
            for index,markers in enumerate(self.markers):
                denominator = count*base.rational(markers[0])**single*base.rational(markers[1])**full
                numerator = self.numerator(index,((single,single),(full,full)))
                candidates.append(tuple(base.upper(entry/denominator) for entry in numerator))
            result[full] = tuple(min(row[i] for row in candidates) for i in range(4))
        return result

    def metadata(self):
        result = super().metadata()
        result['coefficient_method'] = 'positive-residue-aliases-directed-interval-radix2-DFT'
        result['witness']['grid'] = self.grid
        result['regional_diagnostics'] = self.diagnostics
        return result

    def box_receipt(self,ranges,*,marker_indices=(0,0),max_seconds=60.):
        ranges = tail.validate_box(ranges)
        count = tail.profile_count(ranges)
        if (not count or len(marker_indices)!=2 or any(type(i) is not int or not 0<=i<len(self.markers)
                                                      for i in marker_indices)):
            raise ValueError('A nonempty box and two fixed marker indices are required')
        budget = base.Budget(max_seconds)
        vertices = tail.histogram_vertices(ranges)
        with ctx.workprec(self.precision):
            numerators = []
            floor = lambda x:x.numerator//x.denominator
            ceil = lambda x:-((-x.numerator)//x.denominator)
            for parity,index in enumerate(marker_indices):
                singles = [p[0]-2*p[parity+1] for p in vertices]
                full = [p[parity+1] for p in vertices]
                targets = ((floor(min(singles)),ceil(max(singles))),(floor(min(full)),ceil(max(full))))
                numerators.append(self.numerator(index,targets))
            run = base.matrix_power(base.matrix_product(*numerators),128)
            envelopes = []
            for matching in (0,1):
                values = []
                for point in vertices:
                    budget.check()
                    values.append(self.scalar(point,marker_indices,matching).upper())
                envelopes.append(max(values))
            value = tail.log2arb(run[0]+run[1])+min(envelopes)+tail.log2arb(arb(count))
            ceiling,enclosure = base.rational_ceiling(value),str(value)
        budget.check()
        self.check_pins()
        return dict(self.metadata(),ranges=ranges,integer_profile_triples=count,convex_vertices=len(vertices),
                    marker_indices=marker_indices,support_count_method='matching-profile-convex-envelope',
                    log2_upper_rational=str(ceiling),log2_upper_enclosure=enclosure,
                    covers_entire_clipped_box=True,covers_all_even_odd_splits=True)


def main():
    import argparse
    import json
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--r',type=int,required=True)
    parser.add_argument('--z',required=True)
    parser.add_argument('--alpha',required=True)
    parser.add_argument('--marker',nargs=2,action='append',required=True)
    parser.add_argument('--ratio',required=True)
    parser.add_argument('--method',choices=('marginal','runs'),default='runs')
    parser.add_argument('--distance',default='7/20')
    parser.add_argument('--grid',type=int,default=128)
    parser.add_argument('--precision',type=int,default=256)
    parser.add_argument('--bits',type=int,default=80)
    parser.add_argument('--max-seconds',type=float,default=180.)
    parser.add_argument('--output',type=Path)
    args = parser.parse_args()
    if args.output and args.output.exists():
        parser.error('Output exists; choose an unused filename')
    budget = base.Budget(args.max_seconds)
    bank = AliasBank(F(args.z),F(args.alpha),args.marker,distance=F(args.distance),grid=args.grid,
                     precision=args.precision,bits=args.bits,max_seconds=args.max_seconds,
                     progress=lambda row:print(json.dumps(row),flush=True))
    remaining = args.max_seconds-(base.monotonic()-budget.started)
    budget.check()
    result = bank.weight_receipt(args.r,F(args.ratio),method=args.method,max_seconds=remaining)
    encoded = json.dumps(result,indent=2,allow_nan=False)
    if args.output:
        with args.output.open('x',encoding='utf-8') as stream:
            stream.write(encoded+'\n')
    print(encoded)


if __name__=='__main__':
    main()
