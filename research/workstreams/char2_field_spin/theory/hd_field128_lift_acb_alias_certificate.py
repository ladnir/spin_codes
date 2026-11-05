"""Arbitrary-precision directed residue aliases using native ACB DFTs.

This backend uses the same exact full-field packet coefficients and
positive alias argument as the binary64-interval backend. Native ACB
transforms retain arbitrarily small positive contributions and allow the
transform precision itself to increase. No FLOAT FFT output is imported.
"""

import ast
from fractions import Fraction as F
from hashlib import sha256
from pathlib import Path

from flint import acb,arb,ctx

import hd_field128_lift_certificate as base
import hd_field128_lift_tail_certificate as tail
import hd_field128_lift_alias_certificate as alias


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


def dft2(entries,grid,*,inverse=False,budget=None):
    """Two native one-dimensional DFT passes; inverse is normalized."""
    budget = budget or base.Budget()
    if len(entries)!=4 or any(len(row)!=grid*grid for row in entries):
        raise ValueError('Four square coefficient arrays are required')
    result = []
    for values in entries:
        budget.check()
        rows = [x for r in range(grid) for x in acb.dft(values[r*grid:(r+1)*grid],inverse=inverse)]
        following = [None]*(grid*grid)
        for c in range(grid):
            column = acb.dft(rows[c::grid],inverse=inverse)
            for r,value in enumerate(column):
                if not value.is_finite():
                    raise ArithmeticError('Nonfinite native ACB transform')
                following[r*grid+c] = value
        result.append(following)
    return tuple(result)


def polynomial(coefficients,factors,grid):
    result = [[acb(0) for _ in range(grid*grid)] for _ in range(4)]
    for (s,f),values in coefficients.items():
        index = (s%grid)*grid+(f%grid)
        for entry,value in enumerate(values):
            if value.lower()<0:
                raise ArithmeticError('Nonnegative coefficient caps are required')
            result[entry][index] += acb(value*factors[entry])
    return tuple(result)


def paired_power(good,dirty,composed,*,budget=None):
    """At each frequency compute (C+GD+DG+DD)^64 in ACB."""
    budget = budget or base.Budget()
    size = len(good[0])
    result = [[None]*size for _ in range(4)]
    for i in range(size):
        if not i%1024:
            budget.check()
        a,b,c,d = (row[i] for row in good)
        e,f,g,h = (row[i] for row in dirty)
        p,q,r,s = (row[i] for row in composed)
        p += a*e+b*g+e*a+f*c+e*e+f*g
        q += a*f+b*h+e*b+f*d+e*f+f*h
        r += c*e+d*g+g*a+h*c+g*e+h*g
        s += c*f+d*h+g*b+h*d+g*f+h*h
        for _ in range(6):
            trace,qr = p+s,q*r
            p,q,r,s = p*p+qr,q*trace,r*trace,s*s+qr
        if not all(x.is_finite() for x in (p,q,r,s)):
            raise ArithmeticError('Nonfinite ACB regional power')
        result[0][i],result[1][i],result[2][i],result[3][i] = p,q,r,s
    return tuple(result)


class RegionalAliases:
    """Directed positive upper endpoints, with no binary64 conversion."""

    def __init__(self,entries,grid,vector,rho):
        self.grid = grid
        self.upper = []
        for row in entries:
            following = []
            for value in row:
                if not value.imag.contains(0):
                    raise ArithmeticError('Imaginary native alias excludes zero')
                following.append(max(arb(0),value.real.upper()))
            self.upper.append(tuple(following))
        self.upper = tuple(self.upper)
        self.factors = tuple(rho**128*vector[r]/vector[s] for r in range(2) for s in range(2))

    def numerator(self,ranges):
        if len(ranges)!=2 or any(type(lo) is not int or type(hi) is not int or not 0<=lo<=hi
                                for lo,hi in ranges):
            raise ValueError('Two nonnegative integer target intervals are required')
        grid = self.grid
        residues = [range(grid) if hi-lo+1>=grid else tuple(sorted({n%grid for n in range(lo,hi+1)}))
                    for lo,hi in ranges]
        indices = tuple(s*grid+f for s in residues[0] for f in residues[1])
        return tuple(max(row[i] for i in indices)*factor for row,factor in zip(self.upper,self.factors))


def component_alias(good,dirty,composed,grid,budget):
    if type(grid) is not int or grid<64 or grid>256 or grid&(grid-1):
        raise ValueError('Grid must be64,128,or256')
    sums = lambda entries:tuple(sum((row[i] for row in entries.values()),arb(0)) for i in range(4))
    g,d,c = map(sums,(good,dirty,composed))
    positive = tuple(sum(row,arb(0)) for row in zip(c,base.matrix_product(g,d),
                                                    base.matrix_product(d,g),base.matrix_product(d,d)))
    vector,rho = alias.similarity(positive)
    single_factor = tuple(vector[s]/vector[r]/rho for r in range(2) for s in range(2))
    pair_factor = tuple(x/rho for x in single_factor)
    transformed = [dft2(polynomial(coefficients,factors,grid),grid,budget=budget)
                   for coefficients,factors in ((good,single_factor),(dirty,single_factor),(composed,pair_factor))]
    powered = paired_power(*transformed,budget=budget)
    del transformed
    inverse = dft2(powered,grid,inverse=True,budget=budget)
    del powered
    atlas = RegionalAliases(inverse,grid,vector,rho)
    fallback = base.matrix_power(positive,64)
    return atlas,fallback,dict(similarity_vector=list(map(str,vector)),rho=str(rho),
                               interval_arithmetic='native-ACB-DFT',transform_precision=ctx.prec)


class AliasBank(alias.AliasBank):
    """The same weight/box interface with arbitrary-precision aliases."""

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
            raise ValueError('Invalid exact native-ACB alias-bank witness')
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
                    progress(dict(stage='directed-native-ACB-alias',marker=index,grid=grid,precision=precision))
                atlas,fallback,diagnostic = component_alias(*alias.components(single,composed,row),grid,budget)
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
            raise RuntimeError('A pinned transitive native-ACB alias helper changed during replay')

    def metadata(self):
        result = super().metadata()
        result['coefficient_method'] = 'positive-residue-aliases-directed-native-ACB-DFT'
        return result

