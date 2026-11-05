"""Directed checks of selected histogram regions, not full coverage.

Decimal endpoints enclose every operation. Floating search supplies only
candidate centers and matrices, which are rechecked with interval arithmetic.
Standard library, serial, bounded, no output files.
"""
import math
import shutil
import sys
import time
from decimal import Decimal as D,getcontext
from fractions import Fraction as F
from functools import lru_cache
from itertools import product, combinations_with_replacement
from pathlib import Path
from random_pairing_flow_envelope import (
    TYPES,LAW,P,OUTPUT,LOSS,LOG_WEIGHT_NUMERATORS,dot,variational,moments,
    tables,probability_colors,local_data,switches,affine_ellipsoid,eigensystem,
)

getcontext().prec=65
VERIFY_DEADLINE=math.inf


class I:
    def __init__(self,x=0,hi=None):
        if isinstance(x,I):self.lo,self.hi=x.lo,x.hi
        elif isinstance(x,F):
            assert hi is None
            enclosure=I(x.numerator)/x.denominator
            self.lo,self.hi=enclosure.lo,enclosure.hi
        else:self.lo=D(x);self.hi=D(x if hi is None else hi)
    @staticmethod
    def rounded(lo,hi):return I(lo.next_minus(),hi.next_plus())
    def __add__(self,y):
        y=I(y);return I.rounded(self.lo+y.lo,self.hi+y.hi)
    __radd__=__add__
    def __neg__(self):return I(-self.hi,-self.lo)
    def __sub__(self,y):return self+-I(y)
    def __rsub__(self,y):return I(y)+-self
    def __mul__(self,y):
        y=I(y);v=[a*b for a in (self.lo,self.hi) for b in (y.lo,y.hi)]
        return I.rounded(min(v),max(v))
    __rmul__=__mul__
    def __truediv__(self,y):
        y=I(y);assert y.lo>0 or y.hi<0
        return self*I.rounded(1/y.hi,1/y.lo)
    def __rtruediv__(self,y):return I(y)/self
    def exp(self):return I.rounded(self.lo.exp(),self.hi.exp())
    def ln(self):
        assert self.lo>0;return I.rounded(self.lo.ln(),self.hi.ln())
    def sqrt(self):
        assert self.hi>=0
        return I(max(D(0),self.lo).sqrt().next_minus() if self.lo>0 else 0,self.hi.sqrt().next_plus())
    def sq(self):
        lo=0 if self.lo<=0<=self.hi else min(self.lo*self.lo,self.hi*self.hi)
        return I(max(D(0),D(lo).next_minus()),max(self.lo*self.lo,self.hi*self.hi).next_plus())
    def absupper(self):return max(abs(self.lo),abs(self.hi))
    def upper(self):return I(self.hi)
    def midfloat(self):return float((self.lo+self.hi)/2)


def length(v):return sum((x.sq() for x in v),I()).sqrt()


def arithmetic_tests():
    checks=0
    for a in range(-4,5):
        for b in range(1,5):
            x=I(a)/7;y=I(b)/11
            for actual,exact in ((x+y,F(a,7)+F(b,11)),(x-y,F(a,7)-F(b,11)),
                                 (x*y,F(a,7)*F(b,11)),(x/y,F(a,7)/F(b,11)),(x.sq(),F(a,7)**2)):
                assert F(actual.lo)<=exact<=F(actual.hi);checks+=1
    root=I(2).sqrt();assert root.lo*root.lo<=2<=root.hi*root.hi
    inverse=(I(3).ln()).exp();assert inverse.lo<=3<=inverse.hi
    for value in (F(1,3),F(-7,11),F(257,2)):
        enclosure=I(value);assert F(enclosure.lo)<=value<=F(enclosure.hi);checks+=1
    try:positive_logdet([[I(-1)]])
    except AssertionError:pass
    else:raise AssertionError('negative matrix accepted')
    return checks+3


def allowance_interval(log_weights=LOG_WEIGHT_NUMERATORS):
    u=[F(x,10**9) for x in log_weights];k=F(2,25);ordered=sorted(u)
    for size in range(1,6):
        threshold=(2*k+sum(ordered[:size]))/size
        if size==5 or threshold<=ordered[size]:break
    fractions=[max(F(0),(threshold-x)/(2*k)) for x in u]
    assert sum(fractions)==1
    minimum=dot(u,fractions)+k*dot(fractions,fractions)
    boundary=512*(u[2]+k-minimum);boundary=I(boundary.numerator)/boundary.denominator
    p=I(P.numerator)/P.denominator;eta=I(1)/5
    initial=512*((1-p)*(-eta*p).exp()+p*(eta*(1-p)).exp()).ln()
    comb=I(math.comb(8192,1974)).ln()
    regional=comb+1974*p.ln()+6218*(1-p).ln()
    return ((I(1)/10)*371901-comb+256*regional-boundary-initial)/4095


def positive_logdet(matrix):
    a=[row.copy() for row in matrix];result=I()
    for j in range(len(a)):
        pivot=a[j][j];assert pivot.lo>0,('nonpositive pivot',j,pivot.lo)
        result+=pivot.ln()
        for i in range(j+1,len(a)):
            for k in range(j+1,len(a)):a[i][k]-=a[i][j]*a[j][k]/pivot
    return result


def interval_solve(matrix,vector):
    a=[row.copy()+[x] for row,x in zip(matrix,vector)];n=len(vector)
    for j in range(n):
        pivot=a[j][j];assert pivot.lo>0
        a[j]=[x/pivot for x in a[j]]
        for i in range(n):
            if i!=j:
                factor=a[i][j];a[i]=[x-factor*y for x,y in zip(a[i],a[j])]
    return [row[-1] for row in a]


def spectral_boundary(matrix,lower):
    eigen,_=eigensystem([[x.midfloat() for x in row] for row in matrix])
    uncertainty=max(sum(float(x.hi-x.lo)/2 for x in row) for row in matrix)
    for padding in (1e-6,1e-4,.001,.01,.05,.1):
        proposal=min(eigen)-2*uncertainty-padding if lower else max(eigen)+2*uncertainty+padding
        rounded=math.floor(proposal*10**6) if lower else math.ceil(proposal*10**6)
        value=I(str(rounded/10**6));assert value.lo>0,'spectral boundary nonpositive'
        difference=[[matrix[i][j]-value*int(i==j) for j in range(4)] for i in range(4)]
        if not lower:difference=[[-x for x in row] for row in difference]
        try:positive_logdet(difference)
        except AssertionError:continue
        return value
    raise AssertionError('could not certify spectral boundary')


def gaussian_shell_mass(precision,mean,radius,lower):
    mass=I(1)
    if lower:
        least=spectral_boundary(precision,True)
        distance=I(max(D(0),(I(lower)-length(mean)).lo))
        x=least*distance.sq()/2;mass=(1+x)*(-x).exp()
    if radius is None:return I(mass.hi)
    radius=I(radius)
    # The triangle inequality and the two spectral bounds compare both
    # radial events with the exact chi-square-four survival function.
    largest=spectral_boundary(precision,False)
    y=largest*(radius+length(mean).upper()).sq()/2
    mass-=(1+y)*(-y).exp();assert mass.hi>0
    return I(min(D(1),mass.hi))


def interval_data(h,center,log_weights=LOG_WEIGHT_NUMERATORS):
    n=sum(h);theta=I(1)/10;eta=I(1)/5;kappa=I(2)/25
    u=[I(x)/10**9 for x in log_weights];beta=(2*kappa/n).sqrt()
    basis=[[I() for _ in range(4)] for _ in range(5)]
    for j in range(4):
        factor=1/I((j+1)*(j+2)).sqrt()
        for i in range(j+1):basis[i][j]=factor
        basis[j+1][j]=-(j+1)*factor
    projected={key:[beta*sum(basis[c][i]*out[c] for c in range(5)) for i in range(4)] for key,out in OUTPUT.items()}
    base=[u[c]+2*kappa*center[c]/n for c in range(5)]
    entries={};core=I();gradient=[I() for _ in range(4)]
    for a,b in TYPES:
        weights=[I(LAW[j].numerator)/LAW[j].denominator*(theta*LOSS[a,b,j]+eta*j+dot(base,OUTPUT[a,b,j])).exp() for j in range(3)]
        mass=sum(weights);weights=[x/mass for x in weights];f=mass.ln()
        g=[sum(weights[j]*projected[a,b,j][i] for j in range(3)) for i in range(4)]
        centered=[([projected[a,b,j][i]-g[i] for i in range(4)],weights[j]) for j in range(3)]
        covariance=[[sum(weight*v[i]*v[k] for v,weight in centered) for k in range(4)] for i in range(4)]
        diameter=I(max(length([x-y for x,y in zip(projected[a,b,j],projected[a,b,k])]).hi for j in range(3) for k in range(3)))
        numerator,denominator=probability_colors(h,(a,b));count=I((n//2)*numerator)/denominator
        entries[a,b]=(f,g,covariance,diameter,count,centered);core+=count*f
        for i in range(4):gradient[i]+=count*g[i]
    core-=eta*n*I(P.numerator)/P.denominator+dot(u,h)+kappa*(dot(h,h)+dot(center,center))/n
    for i in range(4):gradient[i]-=beta*sum(basis[c][i]*center[c] for c in range(5))
    return entries,core,gradient,beta,kappa


def checked_switches(floating,entries,subset,matrix):
    results=[]
    for _,_,_,key in switches(floating,subset):
        v=[sum(coef*entries[pair][0] for pair,coef in key)]
        v += [sum(coef*entries[pair][1][i] for pair,coef in key) for i in range(4)]
        positive_logdet([[matrix[i][j]-v[i]*v[j] for j in range(5)] for i in range(5)])
        hessian=[[sum(coef*entries[pair][2][i][j] for pair,coef in key) for j in range(4)] for i in range(4)]
        frobenius=length([x for row in hessian for x in row]).upper()
        results.append((key,frobenius,v[0].absupper()))
    return results


def shell(h,entries,core,gradient,matrices,options,radius,lower,slack,priority,quotas=None,box_delta=None,detail=None,box_polynomials=None,
          matrix_mode='interval',control_degree=2,control_values=None,control_vertices=None,control_pairs=None):
    n=sum(h);radius=I(radius);slack=I(slack)
    curvature=[[sum(value[4]*weight*(length(v)*radius).exp()*v[i]*v[j]
                     for value in entries.values() for v,weight in value[5]) for j in range(4)] for i in range(4)]
    curvature_polys={}
    if box_delta is not None:
        for i in range(4):
            for j in range(4):
                scores={key:sum(weight*(length(v)*radius).exp()*v[i]*v[j] for v,weight in value[5]) for key,value in entries.items()}
                base=I()
                for key,score in scores.items():
                    numerator,denominator=probability_colors(h,key)
                    base+=I((n//2)*numerator)/denominator*score
                polynomial=polynomial_coefficients(base,scores,h,axes=box_polynomials[4])
                curvature_polys[i,j]=polynomial
                curvature[i][j]=evaluate_polynomial(polynomial,box_delta)
    budget=[[I() for _ in range(5)] for _ in range(5)];remaining=n//2-1
    phase_matrices=[]
    for index,a in enumerate(priority):
        if quotas is None and remaining<=0:break
        if quotas is not None and index>=len(quotas):break
        subset=tuple(sorted(priority[index:]));steps=min(h[a],remaining);remaining-=steps
        if quotas is not None:steps=quotas[index]
        bounds=[]
        for key,frobenius,_ in options[subset]:
            residual=sum(abs(coef)*((entries[pair][3]*radius).exp()-1)*sum(entries[pair][2][i][i] for i in range(4)) for pair,coef in key)
            bounds.append((frobenius+residual).hi)
        remainder=I(max(bounds,default=D(0)))/2
        for i in range(5):
            for j in range(5):budget[i][j]+=steps*(1+slack)*matrices[subset][i][j]/8
        isotropic=steps*(1+1/slack)*remainder.sq()*radius.sq()/8
        for i in range(1,5):budget[i][i]+=isotropic
        phase=[[(1+slack)*matrices[subset][i][j]/8 for j in range(5)] for i in range(5)]
        for i in range(1,5):phase[i][i]+=(1+1/slack)*remainder.sq()*radius.sq()/8
        phase_matrices.append(phase)
    if box_delta is not None:
        quota_polynomials=box_polynomials[2]
        assert len(quota_polynomials)==len(phase_matrices)
        if box_polynomials[3]:
            # A nested quadratic majorant makes the greedy exposure budget
            # matrix-concave in prefix populations. Its reference affine
            # branch is then an upper bound beyond the branch boundary too.
            for k in range(len(phase_matrices)-2,-1,-1):
                current=phase_matrices[k];following=phase_matrices[k+1]
                difference=[[following[i][j]-current[i][j] for j in range(5)] for i in range(5)]
                shift=I(max(D(0),max((difference[i][i]+sum(I(difference[i][j].absupper()) for j in range(5) if j!=i)).hi for i in range(5))))+I('1e-20')
                for i in range(5):current[i][i]+=shift
                positive_logdet([[current[i][j]-following[i][j] for j in range(5)] for i in range(5)])
        budget_polys={}
        for i in range(5):
            for j in range(5):
                polynomial={}
                for quota,matrix in zip(quota_polynomials,phase_matrices):
                    polynomial=poly_add(polynomial,quota,matrix[i][j])
                budget_polys[i,j]=polynomial
                budget[i][j]=evaluate_polynomial(polynomial,box_delta)
    precision=[[I(int(i==j))-curvature[i][j]-2*budget[i+1][j+1] for j in range(4)] for i in range(4)]
    if box_polynomials is not None and matrix_mode in ('common','bernstein','polytope','triangulated'):
        precision_polys={}
        for i in range(4):
            for j in range(4):
                polynomial=poly_add({(0,0,0,0):I(int(i==j))},curvature_polys[i,j],I(-1))
                precision_polys[i,j]=poly_add(polynomial,budget_polys[i+1,j+1],I(-2))
        if matrix_mode in ('bernstein','polytope','triangulated'):
            core_poly,gradient_poly,_,_,_=box_polynomials
            constant_poly=poly_add(core_poly,budget_polys[0,0])
            linear_poly=[poly_add(gradient_poly[i],budget_polys[0,i+1],I(2)) for i in range(4)]
            return bernstein_gaussian_shell(constant_poly,linear_poly,precision_polys,box_delta,radius,lower,detail,
                                            control_degree,control_values,control_vertices,control_pairs)
        precision=common_precision_bound(precision_polys,box_delta)
    determinant=positive_logdet(precision)
    linear=[gradient[i]+2*budget[0][i+1] for i in range(4)];mean=interval_solve(precision,linear)
    value=core+budget[0][0]-determinant/2+dot(linear,mean)/2
    if box_polynomials is not None:
        core_poly,gradient_poly,_,_,_=box_polynomials
        linear_poly=[poly_add(gradient_poly[i],budget_polys[0,i+1],I(2)) for i in range(4)]
        inverse_columns=[interval_solve(precision,[I(int(i==j)) for i in range(4)]) for j in range(4)]
        energy=poly_add(core_poly,budget_polys[0,0])
        for i in range(4):
            for j in range(4):energy=poly_add(energy,poly_multiply(linear_poly[i],linear_poly[j]),inverse_columns[j][i]/2)
        value=bernstein_polynomial(energy,box_delta)-determinant/2
        if matrix_mode=='common':
            mean=[]
            for i in range(4):
                polynomial={}
                for j in range(4):polynomial=poly_add(polynomial,linear_poly[j],inverse_columns[j][i])
                mean.append(bernstein_polynomial(polynomial,box_delta))
    if detail is not None:
        detail.update(core_upper=str(core.hi),budget_constant_upper=str(budget[0][0].hi),
                      determinant_cost_upper=str((-determinant/2).hi),
                      mean_cost_upper=str((dot(linear,mean)/2).hi),gradient_norm_upper=str(length(gradient).hi))
    value+=gaussian_shell_mass(precision,mean,radius,lower).ln()
    return value


def poly_add(a,b,scale=I(1)):
    result=a.copy()
    for key,value in b.items():result[key]=result.get(key,I())+scale*value
    return result


def poly_multiply(a,b):
    result={}
    for x,v in a.items():
        for y,w in b.items():
            key=tuple(i+j for i,j in zip(x,y));result[key]=result.get(key,I())+v*w
    return result


def evaluate_polynomial(coefficients,delta):
    result=I()
    for powers,coefficient in coefficients.items():
        term=coefficient
        for power,x in zip(powers,delta):
            if power==1:term*=x
            elif power==2:term*=x.sq()
            elif power==3:term*=I.rounded(x.lo**3,x.hi**3)
            elif power==4:term*=x.sq().sq()
            else:assert power==0
        result+=term
    return result


def unit_polynomial(coefficients,box):
    unit={}
    for powers,coefficient in coefficients.items():
        for exponents in product(*(range(power+1) if box[i].hi>box[i].lo else (0,) for i,power in enumerate(powers))):
            term=coefficient
            for power,degree,x in zip(powers,exponents,box):
                low=I(x.lo);width=I(x.hi)-I(x.lo)
                term*=math.comb(power,degree)
                for _ in range(power-degree):term*=low
                for _ in range(degree):term*=width
            unit[exponents]=unit.get(exponents,I())+term
    return unit


def bernstein_coefficients_reference(coefficients,box,degrees=None):
    """Tensor coefficients; an explicit degree supports common matrix weights."""
    unit=unit_polynomial(coefficients,box)
    natural=[max(key[i] for key in unit) for i in range(4)]
    degrees=natural if degrees is None else degrees
    assert all(n<=d for n,d in zip(natural,degrees))
    keys=tuple(product(*(range(degree+1) for degree in degrees)))
    values={key:unit.get(key,I()) for key in keys}
    for axis,degree in enumerate(degrees):
        output={}
        for key in keys:
            value=I()
            for k in range(key[axis]+1):
                source=list(key);source[axis]=k
                value+=values[tuple(source)]*math.comb(key[axis],k)/math.comb(degree,k)
            output[key]=value
        values=output
    return values


def bernstein_coefficients(coefficients,box,degrees=None):
    """Same tensor transform, with shared directed binomial ratios.

    The ratios depend only on the degree, not on the coefficient or axis.
    Rebuild them per call so changing Decimal precision cannot reuse stale
    intervals. Exact zeros and factors of one need no rounded arithmetic.
    """
    unit=unit_polynomial(coefficients,box)
    natural=[max(key[i] for key in unit) for i in range(4)]
    degrees=natural if degrees is None else degrees
    assert all(n<=d for n,d in zip(natural,degrees))
    keys=tuple(product(*(range(degree+1) for degree in degrees)))
    values={key:unit.get(key,I()) for key in keys}
    ratios={degree:tuple(tuple(None if math.comb(row,k)==math.comb(degree,k)
                              else I(math.comb(row,k))/math.comb(degree,k)
                              for k in range(row+1)) for row in range(degree+1))
            for degree in set(degrees)}
    for axis,degree in enumerate(degrees):
        output={}
        for key in keys:
            value=None
            for k,ratio in enumerate(ratios[degree][key[axis]]):
                source=list(key);source[axis]=k;term=values[tuple(source)]
                if term.lo==0 and term.hi==0:continue
                if ratio is not None:term=term*ratio
                value=term if value is None else value+term
            output[key]=I() if value is None else value
        values=output
    return values


def bernstein_polynomial(coefficients,box):
    """Enclose a polynomial on a box by its tensor Bernstein coefficients."""
    values=bernstein_coefficients(coefficients,box)
    return I(min(x.lo for x in values.values()),max(x.hi for x in values.values()))


def common_precision_bound(polynomials,box):
    """A constant P_* <= P(h), checked on every matrix Bernstein coefficient."""
    coefficients={key:bernstein_coefficients(poly,box,(2,2,2,2)) for key,poly in polynomials.items()}
    keys=tuple(coefficients[0,0]);vertices=[]
    for key in keys:vertices.append([[coefficients[i,j][key] for j in range(4)] for i in range(4)])
    reference=[]
    for i in range(4):
        row=[]
        for j in range(4):
            value=sum(vertex[i][j].midfloat() for vertex in vertices)/len(vertices)
            row.append(I(str(value)))
        reference.append(row)
    for i in range(4):
        for j in range(i):reference[i][j]=reference[j][i]
    shift=max(0.,max(-min(eigensystem([[vertex[i][j].midfloat()-reference[i][j].midfloat() for j in range(4)] for i in range(4)])[0])
                     for vertex in vertices))
    for padding in (1e-9,1e-7,1e-5):
        candidate=[[reference[i][j]-I(str(shift+padding))*int(i==j) for j in range(4)] for i in range(4)]
        try:
            positive_logdet(candidate)
            for vertex in vertices:positive_logdet([[vertex[i][j]-candidate[i][j] for j in range(4)] for i in range(4)])
        except AssertionError:continue
        return candidate
    raise AssertionError('common precision bound failed')


def weak_compositions(total,parts):
    if parts==1:
        yield (total,);return
    for first in range(total+1):
        for suffix in weak_compositions(total-first,parts-1):yield (first,)+suffix


def polar_coefficients_reference(polynomial,vertices,degree=2):
    """Quadratic coefficients on barycentric coordinates of any vertex hull."""
    result={}
    for i,v in enumerate(vertices):
        for j in range(i,len(vertices)):
            w=vertices[j];value=I()
            for powers,coefficient in polynomial.items():
                factors=[a for a,k in enumerate(powers) for _ in range(k)]
                assert len(factors)<=2
                if not factors:weight=I(1)
                elif len(factors)==1:weight=(I(v[factors[0]])+I(w[factors[0]]))/2
                else:
                    a,b=factors;weight=(I(v[a])*w[b]+I(w[a])*v[b])/2
                value+=coefficient*weight
            result[i,j]=value
    return elevate_polar(result,len(vertices),degree)


def polar_coefficients(polynomial,vertices,degree=2,control_pairs=None):
    """Use exact integer vertex products before directed coefficient arithmetic."""
    if any(not isinstance(x,int) for v in vertices for x in v):
        result=polar_coefficients_reference(polynomial,vertices,2)
        if control_pairs is not None:result={key:result[key] for key in control_pairs}
        return elevate_polar(result,len(vertices),degree,control_pairs)
    terms=[]
    for powers,coefficient in polynomial.items():
        factors=tuple(a for a,k in enumerate(powers) for _ in range(k))
        assert len(factors)<=2
        terms.append((factors,coefficient))
    result={}
    for i,v in enumerate(vertices):
        for j in range(i,len(vertices)):
            if control_pairs is not None and (i,j) not in control_pairs:continue
            w=vertices[j];value=None
            for factors,coefficient in terms:
                if not factors:numerator=2
                elif len(factors)==1:numerator=v[factors[0]]+w[factors[0]]
                else:
                    a,b=factors;numerator=v[a]*w[b]+w[a]*v[b]
                if not numerator:continue
                if numerator%2:term=coefficient*numerator/2
                else:
                    scale=numerator//2
                    term=coefficient if scale==1 else -coefficient if scale==-1 else coefficient*scale
                value=term if value is None else value+term
            result[i,j]=I() if value is None else value
    return elevate_polar(result,len(vertices),degree,control_pairs)


@lru_cache(maxsize=16)
def restricted_elevation_plan(size,degree,pairs):
    """Retain multi-indices whose active vertices form an allowed-edge clique.

    Every simplex's vertex set is such a clique. Extra cliques are harmless
    upper-bound controls. This retains a common basis across all shells.
    """
    assert 2<degree<=4
    assert math.comb(size+degree-1,degree)<=100000,'elevation enumeration guard'
    allowed=set(pairs);plan=[]
    for indices in combinations_with_replacement(range(size),degree):
        active=tuple(sorted(set(indices)))
        if any((i,j) not in allowed for i in active for j in active if i<=j):continue
        alpha=tuple(indices.count(i) for i in range(size));terms=[]
        for i in active:
            if alpha[i]>1:terms.append(((i,i),alpha[i]*(alpha[i]-1)))
            for j in active:
                if j>i:terms.append(((i,j),2*alpha[i]*alpha[j]))
        assert sum(weight for _,weight in terms)==degree*(degree-1)
        plan.append((alpha,tuple(terms)))
        assert len(plan)<=2000,'restricted polar control count guard'
    return tuple(plan)


def elevate_polar(result,size,degree,control_pairs=None):
    if degree==2:return result
    if control_pairs is not None:
        plan=restricted_elevation_plan(size,degree,tuple(sorted(control_pairs)))
        return {alpha:sum((result[pair]*weight for pair,weight in terms),I())/(degree*(degree-1))
                for alpha,terms in plan}
    assert 2<degree<=4 and math.comb(size+degree-1,degree)<=2000,'polar control count guard'
    elevated={}
    for alpha in weak_compositions(degree,size):
        active=[(i,a) for i,a in enumerate(alpha) if a]
        value=sum(result[i,i]*a*(a-1) for i,a in active)
        for k,(i,a) in enumerate(active):
            for j,b in active[k+1:]:value+=2*a*b*result[i,j]
        elevated[alpha]=value/(degree*(degree-1))
    return elevated


def bernstein_gaussian_shell(constant,linear,precision,box,radius,lower,detail=None,degree=2,control_values=None,vertices=None,control_pairs=None):
    """Convexity bounds the shell integral by its Bernstein control integrals.

    For each h the exponent is a convex combination, with the same Bernstein
    weights for its constant, linear and quadratic coefficients. Jensen's
    inequality applies before integration over this fixed Gaussian shell.
    """
    polynomials=[constant]+linear+[precision[i,j] for i in range(4) for j in range(4)]
    assert degree>=2
    degrees=[degree if x.lo!=x.hi else 0 for x in box]
    assert control_pairs is None or vertices is not None
    coefficients=[bernstein_coefficients(poly,box,degrees) if vertices is None else polar_coefficients(poly,vertices,degree,control_pairs)
                  for poly in polynomials]
    worst=None;worst_key=None;largest_normalizer=None
    for key in coefficients[0]:
        assert time.monotonic()<VERIFY_DEADLINE,'region verification deadline'
        c=coefficients[0][key];v=[coefficients[1+i][key] for i in range(4)]
        matrix=[[coefficients[5+4*i+j][key] for j in range(4)] for i in range(4)]
        determinant=positive_logdet(matrix);mean=interval_solve(matrix,v)
        normalizer=c-determinant/2+dot(v,mean)/2
        value=normalizer+gaussian_shell_mass(matrix,mean,radius,lower).ln()
        if control_values is not None:control_values[key]=I(value.hi)
        if worst is None or value.hi>worst:worst=value.hi;worst_key=key
        largest_normalizer=normalizer.hi if largest_normalizer is None else max(largest_normalizer,normalizer.hi)
    if detail is not None:detail.update(bernstein_controls=len(coefficients[0]),worst_control=worst_key,
                                      largest_log_normalizer=str(largest_normalizer))
    return I(worst)


def simplex_polynomial(coefficients,box):
    """Box coordinates x=lo+width*u with u>=0 and sum(u)<=1."""
    unit=unit_polynomial(coefficients,box);degree=max(sum(key) for key in unit)
    lower=None;upper=None
    def falling(n,k):return math.prod(range(n-k+1,n+1))
    for alpha in product(range(degree+1),repeat=4):
        if sum(alpha)>degree:continue
        value=I()
        for powers,coefficient in unit.items():
            if any(k>a for k,a in zip(powers,alpha)):continue
            numerator=math.prod(falling(a,k) for a,k in zip(alpha,powers))
            value+=coefficient*numerator/falling(degree,sum(powers))
        lower=value.lo if lower is None else min(lower,value.lo)
        upper=value.hi if upper is None else max(upper,value.hi)
    return I(lower,upper)


def polynomial_coefficients(base_value,scores,h,linear=None,diagonal=I(),axes=(0,1,2,4)):
    """Exact quadratic expansion with one coordinate absorbing the sum."""
    n=sum(h);absorber=next(a for a in range(5) if a not in axes)
    basis=[[int(a==axis)-int(a==absorber) for a in range(5)] for axis in axes]
    matrix=[[scores[min(a,b),max(a,b)] for b in range(5)] for a in range(5)]
    linear=linear or [I() for _ in range(5)]
    grad=[sum(matrix[a][b]*h[b] for b in range(5))/(n-1)-matrix[a][a]/(2*(n-1))
          +linear[a]+2*diagonal*h[a] for a in range(5)]
    result={(0,0,0,0):base_value}
    for i,vector in enumerate(basis):
        key=tuple(int(j==i) for j in range(4));result[key]=dot(vector,grad)
    for i in range(4):
        for j in range(4):
            coefficient=sum(basis[i][a]*basis[j][b]*(matrix[a][b]/(2*(n-1))+diagonal*int(a==b)) for a in range(5) for b in range(5))
            key=tuple(int(k==i)+int(k==j) for k in range(4))
            result[key]=result.get(key,I())+coefficient
    return result


def polynomial_box(base_value,scores,h,delta,linear=None,diagonal=I(),axes=(0,1,2,4)):
    return evaluate_polynomial(polynomial_coefficients(base_value,scores,h,linear,diagonal,axes),delta)


def phase_quotas(h,delta,priority,axes=(0,1,2,4),final_phase=None,vertices=None):
    """Affine exposure quotas on a box with a fixed final exposure phase."""
    absorber=next(a for a in range(5) if a not in axes);zero=(0,0,0,0)
    state=[]
    for a in range(5):
        polynomial={zero:I(h[a])}
        for i,axis in enumerate(axes):
            polynomial[tuple(int(j==i) for j in range(4))]=I(int(a==axis)-int(a==absorber))
        state.append(polynomial)
    quotas=[];prefix={zero:I()};total=0
    if final_phase is not None:
        assert 0<=final_phase<5 and vertices
        assert all(sum(v[a] for a in priority[:final_phase])<=255<=
                   sum(v[a] for a in priority[:final_phase+1]) for v in vertices)
    for index,a in enumerate(priority):
        next_prefix=poly_add(prefix,state[a]);total+=h[a]
        if (total>=255 if final_phase is None else index==final_phase):
            before=evaluate_polynomial(prefix,delta);after=evaluate_polynomial(next_prefix,delta)
            # Histograms are integral, so floor/ceil removes outward noise
            # at exact phase boundaries without weakening the condition.
            stable=final_phase is not None or (math.floor(before.hi)<=255 and math.ceil(after.lo)>=255)
            quotas.append(poly_add({zero:I(255)},prefix,I(-1)))
            return quotas,stable
        quotas.append(state[a]);prefix=next_prefix
    raise AssertionError('missing final exposure phase')


def phase_and_shell_tests():
    checks=0;weights=(25,16,9,4,1)
    for prefix in product(range(9),repeat=4):
        if sum(prefix)>8:continue
        populations=prefix+(8-sum(prefix),);remaining=3;actual=0
        for count,weight in zip(populations,weights):
            amount=min(remaining,count);actual+=amount*weight;remaining-=amount
        for final in range(5):
            extension=3*weights[final]
            for i in range(final):extension+=sum(populations[:i+1])*(weights[i]-weights[i+1])
            assert extension>=actual;checks+=1
        h=tuple(x*64 for x in populations);priority=tuple(sorted(range(5),key=lambda a:(h[a],a)))
        polynomials,stable=phase_quotas(h,[I()]*4,priority);assert stable
        remaining=255
        for a,polynomial in zip(priority,polynomials):
            expected=min(remaining,h[a]);remaining-=expected
            value=evaluate_polynomial(polynomial,[I()]*4)
            assert value.lo<=expected<=value.hi;checks+=1
        assert remaining==0
    identity=[[I(int(i==j)) for j in range(4)] for i in range(4)]
    for lower,radius in ((0,2),(1,2),(2,4),(4,6)):
        upper=gaussian_shell_mass(identity,[I()]*4,radius,lower)
        x=I(lower*lower)/2;y=I(radius*radius)/2
        exact=(1+x)*(-x).exp()-(1+y)*(-y).exp()
        assert upper.lo>=exact.hi;checks+=1
    return checks


def polynomial_tests():
    h=(185,188,83,56,0);scores={key:F((key[0]+2*key[1])**2-7,11) for key in TYPES}
    linear=[F(a-2,13) for a in range(5)];diagonal=F(-2,25*512)
    def exact(v):
        return sum(F(v[a]*(v[a]-1),2*511)*value if a==b else F(v[a]*v[b],511)*value
                   for (a,b),value in scores.items())+dot(linear,v)+diagonal*dot(v,v)
    def enclose(x):return I(x.numerator)/x.denominator
    coeff=polynomial_coefficients(enclose(exact(h)),{key:enclose(x) for key,x in scores.items()},h,
                                  [enclose(x) for x in linear],enclose(diagonal))
    box=[I(-1,1)]*3+[I(0,1)];cover=bernstein_polynomial(coeff,box)
    square_cover=bernstein_polynomial(poly_multiply(coeff,coeff),box);checks=0
    for delta in product((-1,0,1),(-1,0,1),(-1,0,1),(0,1)):
        v=[h[i]+delta[i] for i in range(3)]+[h[3]-sum(delta),delta[3]]
        truth=exact(v);point=evaluate_polynomial(coeff,[I(x) for x in delta])
        assert F(point.lo)<=truth<=F(point.hi)
        assert F(cover.lo)<=truth<=F(cover.hi)
        squared=evaluate_polynomial(poly_multiply(coeff,coeff),[I(x) for x in delta])
        assert F(squared.lo)<=truth*truth<=F(squared.hi)
        assert F(square_cover.lo)<=truth*truth<=F(square_cover.hi);checks+=4
    return checks


def tighten_region(bounds,n=512):
    """Preserve exactly the integer vectors in a box with coordinate sum n."""
    assert len(bounds)==5
    bounds=tuple((max(0,lo),min(n,hi)) for lo,hi in bounds)
    while True:
        if any(lo>hi for lo,hi in bounds) or sum(lo for lo,_ in bounds)>n or sum(hi for _,hi in bounds)<n:return None
        lows=sum(lo for lo,_ in bounds);highs=sum(hi for _,hi in bounds)
        updated=tuple((max(lo,n-highs+hi),min(hi,n-lows+lo)) for lo,hi in bounds)
        if updated==bounds:return bounds
        bounds=updated


def region_count(bounds,n=512):
    bounds=tighten_region(bounds,n)
    if bounds is None:return 0
    residual=n-sum(lo for lo,_ in bounds);widths=[hi-lo+1 for lo,hi in bounds];total=0
    for mask in range(32):
        leftover=residual-sum(widths[i] for i in range(5) if mask>>i&1)
        if leftover>=0:total+=(-1 if mask.bit_count()%2 else 1)*math.comb(leftover+4,4)
    assert total>=0
    return total


def region_vertices(bounds,n=512):
    """Vertices of a coordinate box intersected with sum h=n."""
    vertices=set()
    for free in range(5):
        fixed=[a for a in range(5) if a!=free]
        for endpoints in product((0,1),repeat=4):
            vertex=[0]*5
            for a,endpoint in zip(fixed,endpoints):vertex[a]=bounds[a][endpoint]
            vertex[free]=n-sum(vertex)
            if bounds[free][0]<=vertex[free]<=bounds[free][1]:vertices.add(tuple(vertex))
    assert vertices
    return tuple(sorted(vertices))


def exact_rank(rows):
    matrix=[[F(x) for x in row] for row in rows];rank=0
    if not matrix:return 0
    for col in range(len(matrix[0])):
        pivot=next((i for i in range(rank,len(matrix)) if matrix[i][col]),None)
        if pivot is None:continue
        matrix[rank],matrix[pivot]=matrix[pivot],matrix[rank]
        scale=matrix[rank][col];matrix[rank]=[x/scale for x in matrix[rank]]
        for i in range(rank+1,len(matrix)):
            scale=matrix[i][col]
            if scale:matrix[i]=[x-scale*y for x,y in zip(matrix[i],matrix[rank])]
        rank+=1
        if rank==len(matrix):break
    return rank


def pulling_simplices(vertices,facet_normals=None):
    """Exact pulling triangulation with a complete list of facet normals.

    Default facets are coordinate bounds intersected with sum h=n.
    Additional cuts must supply all their normals. Cone the least vertex
    over all facets not containing it, using the same order on every face.
    Returns vertex indices, not a sampled or approximate hull.
    """
    assert vertices and all(len(v)==5 and sum(v)==sum(vertices[0]) for v in vertices)
    normals=(tuple(tuple(int(i==j) for i in range(5)) for j in range(5))
             if facet_normals is None else tuple(facet_normals))
    cache={}
    def dimension(ids):
        base=vertices[ids[0]]
        return exact_rank([[x-y for x,y in zip(vertices[i],base)] for i in ids[1:]])
    def visit(ids):
        if ids in cache:return cache[ids]
        d=dimension(ids)
        if len(ids)==d+1:return (ids,)
        pivot=ids[0];facets=set()
        for normal in normals:
            values={i:sum(a*b for a,b in zip(normal,vertices[i])) for i in ids}
            for endpoint in (min(values.values()),max(values.values())):
                face=tuple(i for i in ids if values[i]==endpoint)
                if pivot not in face and dimension(face)==d-1:facets.add(face)
        assert facets,'missing facets in pulling triangulation'
        result=tuple(sorted({tuple(sorted((pivot,)+s)) for face in sorted(facets) for s in visit(face)}))
        assert all(len(s)==d+1 and dimension(s)==d for s in result)
        cache[ids]=result
        return result
    return visit(tuple(range(len(vertices))))


def triangulation_pairs(simplices):
    return frozenset((i,j) for simplex in simplices for i in simplex for j in simplex if i<=j)


def bisect_simplex_edge(vertices,simplices,edge):
    """Cover the same hull by splitting every simplex incident to this edge."""
    i,j=sorted(edge);assert i!=j and (i,j) in triangulation_pairs(simplices)
    midpoint=tuple(F(a+b,2) for a,b in zip(vertices[i],vertices[j]))
    midpoint=tuple(x.numerator if x.denominator==1 else x for x in midpoint)
    assert midpoint not in vertices,'edge already contains a mesh vertex'
    new_id=len(vertices);result=[]
    for simplex in simplices:
        if i in simplex and j in simplex:
            result.extend(tuple(sorted(tuple(k for k in simplex if k!=drop)+(new_id,))) for drop in (i,j))
        else:result.append(simplex)
    return vertices+(midpoint,),tuple(sorted(set(result)))


def population_triangulation(bounds,n=512,star=False):
    vertices=region_vertices(bounds,n)
    if not star:return vertices,pulling_simplices(vertices)
    center=region_reference(bounds,n)
    # An integral center need not be relative-interior in a thin region.
    # In that case retain the complete pulling triangulation.
    if center in vertices or not all(lo==hi or lo<center[a]<hi for a,(lo,hi) in enumerate(bounds)):
        return vertices,pulling_simplices(vertices)
    d=exact_rank([[x-y for x,y in zip(v,vertices[0])] for v in vertices[1:]])
    facets=set()
    for a,(lo,hi) in enumerate(bounds):
        for endpoint in (lo,hi):
            ids=tuple(i for i,v in enumerate(vertices) if v[a]==endpoint)
            if ids and exact_rank([[x-y for x,y in zip(vertices[i],vertices[ids[0]])] for i in ids[1:]])==d-1:
                facets.add(ids)
    assert facets
    center_id=len(vertices);simplices=[]
    for ids in sorted(facets):
        for simplex in pulling_simplices(tuple(vertices[i] for i in ids)):
            simplices.append(tuple(sorted(tuple(ids[i] for i in simplex)+(center_id,))))
    return vertices+(center,),tuple(sorted(set(simplices)))


def region_reference(bounds,n=512):
    """Round the mean of at most80 feasible vertices to an integer point."""
    vertices=region_vertices(bounds,n)
    totals=[sum(v[a] for v in vertices) for a in range(5)];denominator=len(vertices)
    result=[value//denominator for value in totals]
    order=sorted(range(5),key=lambda a:(-(totals[a]%denominator),a))
    for a in order:
        if sum(result)==n:break
        if result[a]<bounds[a][1]:result[a]+=1
    assert sum(result)==n and all(lo<=value<=hi for value,(lo,hi) in zip(result,bounds))
    return tuple(result)


def polar_square_controls(polynomials,vertices):
    """Exact degree-four Bernstein coefficients of a sum of quadratic squares.

    Multiply degree-two barycentric representations. Keeping signed cross
    terms avoids the Jensen loss from squaring controls independently.
    """
    size=len(vertices)
    assert math.comb(size+3,4)<=2000,'polar control count guard'
    controls=[polar_coefficients(poly,vertices) for poly in polynomials]
    pairs=list(controls[0]);result={}
    for index,ab in enumerate(pairs):
        for cd in pairs[index:]:
            alpha=tuple((ab+cd).count(i) for i in range(size))
            multinomial=24//math.prod(math.factorial(x) for x in alpha)
            numerator=(1 if ab[0]==ab[1] else 2)*(1 if cd[0]==cd[1] else 2)*(1 if ab==cd else 2)
            value=sum(poly[ab]*poly[cd] for poly in controls)*numerator/multinomial
            result[alpha]=result.get(alpha,I())+value
    return result


def coupled_radial_bound(core_poly,gradient_poly,budget_constant,a,linear,delta,range_bound,lower=0,
                         vertices=None,degree=2):
    """Integrate the core's linear Gaussian term before bounding its h-dependence.

    Only the matching correction's radial linear term uses Young's inequality.
    The degree-four energy retains cancellation between core and gradient.
    """
    norm_square={}
    for polynomial in gradient_poly:norm_square=poly_add(norm_square,poly_multiply(polynomial,polynomial))
    square_upper=max(D(0),range_bound(norm_square,delta).hi)
    controls=None
    if vertices is not None:
        # The same nonnegative polar weights represent both C(h) and g(h).
        # Convexity of ||g||^2 preserves their common population dependence.
        constants=polar_coefficients(poly_add(core_poly,budget_constant),vertices,degree)
        if degree==4:
            squares=polar_square_controls(gradient_poly,vertices)
            assert squares.keys()==constants.keys()
            controls={key:(c,squares[key]) for key,c in constants.items()}
        else:
            gradients=[polar_coefficients(poly,vertices,degree) for poly in gradient_poly]
            controls={key:(c,sum(g[key].sq() for g in gradients)) for key,c in constants.items()}
        square_upper=min(square_upper,max(D(0),max(value[1].hi for value in controls.values())))
    # A floating scalar proposal is fixed and checked below. It is not evidence.
    aa=float(a.lo);cc=float(linear.hi)**2;gg=float(square_upper)
    lo,hi=0.,aa
    for _ in range(50):
        middle=(lo+hi)/2;gamma=aa-middle
        derivative=-cc/(4*middle*middle)+gg/(4*gamma*gamma)+2/gamma
        if derivative>0:hi=middle
        else:lo=middle
    epsilon=I(str(max(1e-15,min(.999*aa,(lo+hi)/2))))
    assert epsilon.lo>0 and epsilon.hi<a.lo
    gamma=a-epsilon
    energy=poly_add(poly_add(core_poly,budget_constant),norm_square,1/(4*gamma))
    energy_bound=range_bound(energy,delta)
    if controls is not None:
        polar_upper=max((c+square/(4*gamma)).hi for c,square in controls.values())
        energy_bound=I(min(energy_bound.hi,polar_upper))
    bound=energy_bound+linear.sq()/(4*epsilon)-(4*gamma.sq()).ln()
    if lower:
        mean_norm=I(square_upper).sqrt()/(2*gamma)
        distance=I(max(D(0),(I(lower)-mean_norm).lo))
        x=I(gamma.lo)*distance.sq()
        bound+=((1+x)*(-x).exp()).ln()
    return bound


def radial_budget(ranges,quota_polynomials,delta,range_bound,polynomial_terms=None):
    """Monotone scalar coefficients justify the affine phase extension."""
    for index in range(len(ranges)-2,-1,-1):
        for j in range(2):ranges[index][j]=I(max(ranges[index][j].hi,ranges[index+1][j].hi))
    polynomials=[{}, {}, {}]
    for quota,(offset,slope) in zip(quota_polynomials,ranges):
        for index,coefficient in enumerate((offset.sq()/8,offset*slope/4,slope.sq()/8)):
            polynomials[index]=poly_add(polynomials[index],quota,coefficient)
    coefficients=[range_bound(poly,delta) for poly in polynomials]
    full_offset,full_slope=ranges[0]
    caps=(255*full_offset.sq()/8,255*full_offset*full_slope/4,255*full_slope.sq()/8)
    constant,linear,quadratic=[I(max(D(0),min(value.hi,cap.hi))) for value,cap in zip(coefficients,caps)]
    if polynomial_terms is not None:polynomial_terms.update(constant=polynomials[0],linear=polynomials[1],quadratic=polynomials[2])
    return constant,linear,quadratic,polynomials[0]


def global_radial_budget(floating,entries,beta,priority,quota_polynomials,delta,range_bound,polynomial_terms=None,
                         ball_switch_bounds=None):
    lipschitz=4*I(2).sqrt()*beta
    ranges=[]
    for index in range(len(quota_polynomials)):
        subset=tuple(sorted(priority[index:]));largest=D(0)
        for _,_,_,key in switches(floating,subset):
            difference=sum(coef*entries[pair][0] for pair,coef in key)
            largest=max(largest,difference.absupper())
        slope=lipschitz if len(subset)>1 else I()
        if ball_switch_bounds is not None:
            from random_pairing_switch_ball import subset_squared_bound
            slope=I(min(slope.hi,(beta*subset_squared_bound(ball_switch_bounds,subset).sqrt()).hi))
        ranges.append([I(largest),slope])
    return radial_budget(ranges,quota_polynomials,delta,range_bound,polynomial_terms)


def global_covariance_polynomials(h,entries,axes):
    """Outcome-geometry covariance bounds valid for every Gaussian tilt.

    For categorical probabilities p, diag(p)-p p^T <= I/2 by Gershgorin.
    Centering the distinct possible outcome vectors gives the stated bound,
    independent of their tilted probabilities. Duplicate outcomes are removed.
    """
    matrices={}
    for pair,value in entries.items():
        unique={}
        for j in range(3):unique.setdefault(OUTPUT[pair[0],pair[1],j],value[5][j][0])
        vectors=list(unique.values());center=[sum(v[i] for v in vectors)/len(vectors) for i in range(4)]
        centered=[[v[i]-center[i] for i in range(4)] for v in vectors]
        matrices[pair]=[[sum(v[i]*v[j] for v in centered)/2 for j in range(4)] for i in range(4)]
    result={}
    for i in range(4):
        for j in range(4):
            scores={key:matrix[i][j] for key,matrix in matrices.items()}
            base=sum(entries[key][4]*score for key,score in scores.items())
            result[i,j]=polynomial_coefficients(base,scores,h,axes=axes)
    return result


def coarse_region_bound(h,floating,entries,core,gradient,beta,kappa,priority,quota_polynomials,delta,range_bound,
                        coupled_polynomials=None,detail=None,vertices=None,degree=2):
    """All-Gaussian radial bound retaining rare-color exposure counts."""
    constant,linear,quadratic,constant_poly=global_radial_budget(floating,entries,beta,priority,quota_polynomials,delta,range_bound)
    a=(I('0.5')-kappa-quadratic);b=length(gradient).upper()+linear
    assert a.lo>0
    # Young's inequality b*r <= epsilon*r^2+b^2/(4*epsilon).
    b_upper=float(b.hi);a_lower=float(a.lo)
    proposal=2*b_upper*a_lower/(math.sqrt(b_upper*b_upper+32*a_lower)+b_upper) if b_upper else 0
    epsilon=I(str(max(1e-15,min(.999*a_lower,proposal))))
    assert epsilon.lo>0 and epsilon.hi<a.lo
    result=core+constant+b.sq()/(4*epsilon)-(4*(a-epsilon).sq()).ln()
    if detail is not None:
        detail.update(core_upper=str(core.hi),gradient_norm_upper=str(length(gradient).hi),
                      exposure_constant_upper=str(constant.hi),exposure_linear_upper=str(linear.hi),
                      exposure_quadratic_upper=str(quadratic.hi),radial_log_upper=str(result.hi))
    if coupled_polynomials is not None:
        core_poly,gradient_poly=coupled_polynomials
        coupled=coupled_radial_bound(core_poly,gradient_poly,constant_poly,a,linear,delta,range_bound,
                                    vertices=vertices,degree=degree)
        if detail is not None:detail['coupled_log_upper']=str(coupled.hi)
        result=I(min(result.hi,coupled.hi))
    return result


def local_radial_region_bound(h,floating,entries,core,gradient,beta,kappa,priority,quota_polynomials,
                              delta,range_bound,core_poly,gradient_poly,axes,radius=16,lower=0,vertices=None,degree=2,
                              ball_switch_bounds=None):
    """A ball-local scalar curvature bound plus an independent global tail."""
    radius=I(radius);curvature=[[I() for _ in range(4)] for _ in range(4)]
    for i in range(4):
        for j in range(4):
            scores={key:sum(weight*(length(v)*radius).exp()*v[i]*v[j] for v,weight in value[5])
                    for key,value in entries.items()}
            base=sum(entries[key][4]*value for key,value in scores.items())
            curvature[i][j]=range_bound(polynomial_coefficients(base,scores,h,axes=axes),delta)
    gershgorin=max((curvature[i][i]+sum(I(curvature[i][j].absupper()) for j in range(4) if j!=i)).hi for i in range(4))
    curvature_bound=min(2*kappa.hi,gershgorin)
    ranges=[];global_ranges=[];global_lip=4*I(2).sqrt()*beta
    for index in range(len(quota_polynomials)):
        subset=tuple(sorted(priority[index:]));offset=D(0);slope=D(0)
        for _,_,_,key in switches(floating,subset):
            difference=sum(coef*entries[pair][0] for pair,coef in key)
            grad=[sum(coef*entries[pair][1][i] for pair,coef in key) for i in range(4)]
            hessian=[[sum(coef*entries[pair][2][i][j] for pair,coef in key) for j in range(4)] for i in range(4)]
            residual=sum(abs(coef)*((entries[pair][3]*radius).exp()-1)*sum(entries[pair][2][i][i] for i in range(4))
                         for pair,coef in key)
            hessian_norm=length([x for row in hessian for x in row])+residual
            offset=max(offset,difference.absupper())
            slope=max(slope,(length(grad)+radius*hessian_norm/2).hi)
        if ball_switch_bounds is not None:
            from random_pairing_switch_ball import subset_squared_bound
            slope=min(slope,(beta*subset_squared_bound(ball_switch_bounds,subset).sqrt()).hi)
        ranges.append([I(offset),I(min(slope,global_lip.hi))])
        global_ranges.append([I(offset),global_lip if len(subset)>1 else I()])
    constant,linear,quadratic,constant_poly=radial_budget(ranges,quota_polynomials,delta,range_bound)
    a=I('0.5')-I(curvature_bound)/2-quadratic
    assert a.lo>0,'local radial coefficient nonpositive'
    ball=coupled_radial_bound(core_poly,gradient_poly,constant_poly,a,linear,delta,range_bound,lower,vertices,degree)
    c0,c1,c2,global_constant_poly=radial_budget(global_ranges,quota_polynomials,delta,range_bound)
    tail_a=I('0.5')-kappa-c2;tail_b=length(gradient).upper()+c1
    assert tail_a.lo>0
    epsilon=I(str(max(1e-15,min(float(tail_a.lo)/2,float(tail_b.hi)/(2*float(radius.lo))))))
    assert epsilon.lo>0 and epsilon.hi<tail_a.lo
    gamma=tail_a-epsilon
    tail=core+c0+tail_b.sq()/(4*epsilon)+(1+gamma*radius.sq()).ln()-(4*gamma.sq()).ln()-gamma*radius.sq()
    coupled_tail=coupled_radial_bound(core_poly,gradient_poly,global_constant_poly,tail_a,c1,delta,range_bound,radius.lo,
                                     vertices,degree)
    tail=I(min(tail.hi,coupled_tail.hi))
    return (ball.exp()+tail.exp()).ln(),dict(ball_log_upper=str(ball.hi),tail_log_upper=str(tail.hi),
        core_curvature_upper=str(curvature_bound),local_c0=str(constant.hi),local_c1=str(linear.hi),local_c2=str(quadratic.hi))


def verify_region(bounds,coarse_only=False,coupled_coarse=False,local_radial=False,matrix_mode='interval',control_degree=2,
                  exposure_priority=None,log_weights=LOG_WEIGHT_NUMERATORS,radial_controls=False,matrix_gate=None,simplex_splits=(),
                  bounded_switch=False,exposure_phase=None,force_nested_phase=False,center_override=None,integrated_degree=0):
    """Bound every integer histogram in a five-coordinate clipped box."""
    star=matrix_mode=='star'
    if star:matrix_mode='triangulated'
    assert not integrated_degree or (bounded_switch and matrix_mode=='triangulated' and control_degree<=integrated_degree<=12)
    bounds=tighten_region(bounds);assert bounds is not None
    phase_vertices=None;phase_simplices=None;histograms=region_count(bounds)
    assert not force_nested_phase or exposure_phase is not None
    if exposure_phase is None:h=region_reference(bounds)
    else:
        assert exposure_priority is not None and matrix_mode=='triangulated'
        assert not (coarse_only or coupled_coarse or local_radial or radial_controls or matrix_gate)
        from random_pairing_flow_phases import phase_geometry,phase_count
        phase_vertices,phase_simplices,h=phase_geometry(bounds,tuple(exposure_priority),exposure_phase,star=star)
        histograms=phase_count(bounds,tuple(exposure_priority),exposure_phase)
        assert phase_vertices and histograms>0,'empty owned exposure phase'
    u=[x/10**9 for x in log_weights]
    floating,_,_=local_data(h,u,.1,.2,.08)
    _,proposal,_=variational(h,u,.1,.2,.08);center,_=moments(h,tables(proposal))
    center=[I(str(round(x,12))) for x in center[:4]];center.append(I(512)-sum(center))
    if center_override is not None:
        # The Gaussian identity permits any fixed center. This option
        # isolates center loss from population-control loss in diagnostics.
        assert len(center_override)==4
        center=[I(str(x)) for x in center_override];center.append(I(512)-sum(center))
        floating,_,_=local_data(h,u,.1,.2,.08,center=[float(x.lo) for x in center])
    entries,core,gradient,beta,kappa=interval_data(h,center,log_weights)
    absorber=max(range(5),key=lambda a:(bounds[a][1]-bounds[a][0],h[a],a==3,-a))
    axes=tuple(a for a in range(5) if a!=absorber)
    delta=[I(bounds[a][0]-h[a],bounds[a][1]-h[a]) for a in axes]
    variables=[I(lo,hi) for lo,hi in bounds]
    rational_u=[I(x)/10**9 for x in log_weights]
    core_poly=polynomial_coefficients(core,{key:value[0] for key,value in entries.items()},h,
                                     [-x for x in rational_u],-kappa/512,axes)
    gradient_poly=[polynomial_coefficients(gradient[i],{key:value[1][i] for key,value in entries.items()},h,axes=axes) for i in range(4)]
    residual=512-sum(lo for lo,_ in bounds)
    pure_simplex=all(hi-lo==residual for lo,hi in bounds)
    range_bound=simplex_polynomial if pure_simplex else bernstein_polynomial
    core=range_bound(core_poly,delta)
    gradient=[range_bound(poly,delta) for poly in gradient_poly]
    priority=tuple(sorted(range(5),key=lambda a:(h[a],a))) if exposure_priority is None else tuple(exposure_priority)
    assert sorted(priority)==list(range(5))
    matrices={};options={}
    quota_polynomials,stable_phase=phase_quotas(h,delta,priority,axes,exposure_phase,phase_vertices)
    if force_nested_phase:stable_phase=False
    quotas=[evaluate_polynomial(poly,delta) for poly in quota_polynomials]
    coarse_detail={}
    radial_vertices=[tuple(v[a]-h[a] for a in axes) for v in region_vertices(bounds)] if radial_controls else None
    coarse=coarse_region_bound(h,floating,entries,core,gradient,beta,kappa,priority,quota_polynomials,delta,range_bound,
                              (core_poly,gradient_poly) if coupled_coarse else None,coarse_detail,radial_vertices,control_degree)
    allowance=allowance_interval(log_weights)
    if local_radial and (local_radial!='auto' or coarse.hi>=allowance.lo):
        def ball_switches(radius):
            if not bounded_switch:return None
            from random_pairing_switch_ball import switch_squared_bounds
            return switch_squared_bounds(radius,log_weights,center)
        if local_radial in ('shells','shells32'):
            total=I();local_detail={}
            annuli=((0,4),(4,8),(8,16))
            if local_radial=='shells32':annuli+=((16,32),)
            for lower,radius in annuli:
                _,part=local_radial_region_bound(h,floating,entries,core,gradient,beta,kappa,priority,quota_polynomials,
                                               delta,range_bound,core_poly,gradient_poly,axes,radius,lower,radial_vertices,control_degree,
                                               ball_switches(radius))
                total+=I(part['ball_log_upper']).exp();local_detail[str((lower,radius))]=part
            total+=I(part['tail_log_upper']).exp();local=total.ln()
        else:
            try:
                local,local_detail=local_radial_region_bound(h,floating,entries,core,gradient,beta,kappa,priority,quota_polynomials,
                                                            delta,range_bound,core_poly,gradient_poly,axes,
                                                            vertices=radial_vertices,degree=control_degree,
                                                            ball_switch_bounds=ball_switches(16))
            except AssertionError as error:
                if local_radial!='auto' or str(error)!='local radial coefficient nonpositive':raise
                local=None;local_detail={'local_radial_inconclusive':str(error)}
        coarse_detail.update(local_detail)
        coarse_detail['bounded_switch']=bounded_switch
        if local is not None:
            coarse_detail['local_radial_log_upper']=str(local.hi)
            coarse=I(min(coarse.hi,local.hi))
    gate_rejects=matrix_gate is not None and (max(hi-lo for lo,hi in bounds)>matrix_gate[0] or coarse.hi>D(matrix_gate[1]))
    if exposure_phase is None and center_override is None and (coarse_only or coarse.hi<allowance.lo or gate_rejects):
        return dict(bounds=bounds,reference=h,absorber=absorber,integer_histograms=region_count(bounds),
                    log_weight_numerators=tuple(log_weights),log_weight_denominator=10**9,
                    radial_controls=radial_controls,control_degree=control_degree,
                    matrix_gate_rejected=gate_rejects,
                    directed_log_upper=str(coarse.hi),passes=coarse.hi<allowance.lo,
                    proof_method='radial_combination' if local_radial else 'coarse_radial',pure_simplex=pure_simplex,
                    coarse_detail=coarse_detail,whole_code_certificate=False)
    for index in range(len(quotas)):
        subset=tuple(sorted(priority[index:]));candidate=affine_ellipsoid(floating,subset)
        matrix=[[I(str(candidate[i][j]))*I('1.000001')+I('0.000000000001')*int(i==j) for j in range(5)] for i in range(5)]
        matrices[subset]=matrix;options[subset]=checked_switches(floating,entries,subset,matrix)
    adjusted={}
    for (a,b),value in entries.items():
        count=variables[a]*(variables[a]-1)/(2*511) if a==b else variables[a]*variables[b]/511
        # Counts are nonnegative on integer histograms, including h4=0 or1.
        count=I(max(D(0),count.lo),count.hi)
        adjusted[a,b]=value[:4]+(count,)+value[5:]
    total=I();detail={};control_totals={}
    control_vertices=None;control_pairs=None
    if matrix_mode in ('polytope','triangulated'):
        if matrix_mode=='triangulated':
            if exposure_phase is None:vertices,simplices=population_triangulation(bounds,star=star)
            else:vertices,simplices=phase_vertices,phase_simplices
            for edge in simplex_splits:vertices,simplices=bisect_simplex_edge(vertices,simplices,edge)
        else:vertices=region_vertices(bounds)
        control_vertices=[tuple(v[a]-h[a] for a in axes) for v in vertices]
        if matrix_mode=='triangulated':
            control_pairs=triangulation_pairs(simplices)
            retained=(len(control_pairs) if control_degree==2 else
                      len(restricted_elevation_plan(len(vertices),control_degree,tuple(sorted(control_pairs)))))
            detail.update(simplices=len(simplices),retained_controls=retained,
                          unrestricted_controls=math.comb(len(vertices)+control_degree-1,control_degree),
                          star_center=star,simplex_splits=simplex_splits)
    for lower,outer,slack in ((0,4,'0.2'),(4,5,'0.2'),(5,6,'0.4'),(6,8,'0.4'),(8,12,'0.4'),(12,16,'0.4')):
        controls={} if matrix_mode in ('bernstein','polytope','triangulated') else None
        total+=shell(h,adjusted,core,gradient,matrices,options,outer,lower,slack,priority,quotas,delta,
                     detail if lower==0 else None,(core_poly,gradient_poly,quota_polynomials,not stable_phase,axes),
                     matrix_mode,control_degree,controls,control_vertices,control_pairs).exp()
        if controls is not None:
            if control_totals:assert control_totals.keys()==controls.keys()
            for key,value in controls.items():control_totals[key]=control_totals.get(key,I())+value.exp()
    if control_totals:
        # The histogram's Bernstein weights are the same in every annulus.
        # Sum each control's contributions before maximizing over controls.
        detail['separate_shell_max_log_upper']=str(total.ln().hi)
        total=I(max(value.hi for value in control_totals.values()))
        detail['joint_shell_max_log_upper']=str(total.ln().hi)
        if control_vertices is not None:
            worst=max(control_totals,key=lambda key:control_totals[key].hi)
            alpha=([sum(i==j for j in worst) for i in range(len(vertices))] if control_degree==2 else worst)
            assert sum(alpha)==control_degree
            detail['joint_worst_control']=worst
            detail['joint_control_barycenter']=tuple(str(F(sum(a*v[i] for a,v in zip(alpha,vertices)),control_degree)) for i in range(5))
            pure=[value.hi for key,value in control_totals.items()
                  if (key[0]==key[1] if control_degree==2 else max(key)==control_degree)]
            detail['pure_vertex_shell_log_upper']=str(I(max(pure)).ln().hi)
    full=tuple(sorted(priority));deltamax=I(max(x[2] for x in options[full]));lip=4*I(2).sqrt()*beta
    a=I('0.5')-kappa-255*lip.sq()/8
    b=length(gradient)+255*deltamax*lip/4;d=255*deltamax.sq()/8
    tail_polynomials={}
    c0,c1,c2,constant_poly=global_radial_budget(floating,entries,beta,priority,quota_polynomials,delta,range_bound,tail_polynomials)
    tail=coupled_radial_bound(core_poly,gradient_poly,constant_poly,I('0.5')-kappa-c2,c1,delta,range_bound,16)
    effective=a-b/16
    if effective.lo>0:
        old_tail=core+d+(1+effective*256).ln()-(4*effective.sq()).ln()-effective*256
        tail=I(min(tail.hi,old_tail.hi))
    detail['tail_log_upper']=str(tail.hi)
    total+=tail.exp()
    if control_totals:
        tail_a=I('0.5')-kappa-c2;epsilon=I(str(tail_a.lo/4))
        assert epsilon.lo>0 and epsilon.hi<tail_a.lo
        gamma=tail_a-epsilon
        tail_constant=poly_add(poly_add(core_poly,constant_poly),
                               poly_multiply(tail_polynomials['linear'],tail_polynomials['linear']),1/(4*epsilon))
        tail_precision={(i,j):{(0,0,0,0):2*gamma*int(i==j)} for i in range(4) for j in range(4)}
        tail_controls={}
        control_tail=bernstein_gaussian_shell(tail_constant,gradient_poly,tail_precision,delta,None,16,
                                              degree=control_degree,control_values=tail_controls,vertices=control_vertices,control_pairs=control_pairs)
        assert tail_controls.keys()==control_totals.keys()
        joint=I(max((control_totals[key]+value.exp()).hi for key,value in tail_controls.items()))
        detail['control_tail_log_upper']=str(control_tail.hi)
        detail['joint_total_log_upper']=str(joint.ln().hi)
        total=I(min(total.hi,joint.hi))
        geometry=global_covariance_polynomials(h,entries,axes)
        geometry_precision={(i,j):poly_add({(0,0,0,0):(I(1)-2*c2-2*epsilon)*int(i==j)},geometry[i,j],I(-1))
                            for i in range(4) for j in range(4)}
        geometry_controls={}
        try:
            geometry_tail=bernstein_gaussian_shell(tail_constant,gradient_poly,geometry_precision,delta,None,16,
                                                   degree=control_degree,control_values=geometry_controls,vertices=control_vertices,control_pairs=control_pairs)
        except AssertionError as error:
            if not str(error).startswith("('nonpositive pivot'"):raise
            detail['geometry_tail_inconclusive']=str(error)
        else:
            assert geometry_controls.keys()==control_totals.keys()
            geometry_joint=I(max((control_totals[key]+value.exp()).hi for key,value in geometry_controls.items()))
            detail['geometry_tail_log_upper']=str(geometry_tail.hi)
            detail['geometry_total_log_upper']=str(geometry_joint.ln().hi)
            total=I(min(total.hi,geometry_joint.hi))
        if bounded_switch:
            from random_pairing_switch_ball import switch_squared_bounds
            # Same population controls for the central shells, bounded
            # annulus 16..32, and unrestricted tail beyond32.
            local_polys={}
            _,lc1,lc2,local_constant=global_radial_budget(
                floating,entries,beta,priority,quota_polynomials,delta,range_bound,local_polys,
                switch_squared_bounds(32,log_weights,center))
            local_a=I('0.5')-kappa-lc2;local_epsilon=I(str(local_a.lo/4))
            assert 0<local_epsilon.lo and local_epsilon.hi<local_a.lo
            local_c=poly_add(poly_add(core_poly,local_constant),
                             poly_multiply(local_polys['linear'],local_polys['linear']),1/(4*local_epsilon))
            local_p={(i,j):poly_add({(0,0,0,0):(I(1)-2*lc2-2*local_epsilon)*int(i==j)},geometry[i,j],I(-1))
                     for i in range(4) for j in range(4)}
            finite_controls={};far_controls={}
            try:
                finite=bernstein_gaussian_shell(local_c,gradient_poly,local_p,delta,32,16,
                    degree=control_degree,control_values=finite_controls,vertices=control_vertices,control_pairs=control_pairs)
                far=bernstein_gaussian_shell(tail_constant,gradient_poly,geometry_precision,delta,None,32,
                    degree=control_degree,control_values=far_controls,vertices=control_vertices,control_pairs=control_pairs)
            except AssertionError as error:
                if not str(error).startswith("('nonpositive pivot'"):raise
                detail['bounded_control_tail_inconclusive']=str(error)
            else:
                assert finite_controls.keys()==far_controls.keys()==control_totals.keys()
                combined=I(max((control_totals[key]+finite_controls[key].exp()+far_controls[key].exp()).hi
                               for key in control_totals))
                detail.update(bounded_annulus_log_upper=str(finite.hi),beyond32_log_upper=str(far.hi),
                              bounded_control_total_log_upper=str(combined.ln().hi))
                total=I(min(total.hi,combined.hi))
                if integrated_degree:
                    from random_pairing_flow_integrated_controls import integrated_upper
                    integrated={key:control_totals[key]+finite_controls[key].exp()+far_controls[key].exp() for key in control_totals}
                    refined,metadata=integrated_upper(integrated,len(vertices),simplices,control_degree,integrated_degree)
                    detail.update(metadata,integrated_log_upper=str(refined.ln().hi))
                    total=I(min(total.hi,refined.hi))
    result=total.ln();allowance=allowance_interval(log_weights)
    return dict(bounds=bounds,reference=h,absorber=absorber,integer_histograms=histograms,exposure_phase=exposure_phase,
                gaussian_center=tuple(str(x.lo) for x in center[:4]),
                log_weight_numerators=tuple(log_weights),log_weight_denominator=10**9,
                radial_controls=radial_controls,control_degree=control_degree,
                directed_log_upper=str(result.hi),passes=result.hi<allowance.lo,
                proof_method='local_quadratic',matrix_mode=matrix_mode,exposure_priority=priority,phase_extension=not stable_phase,
                core_width=str(core.hi-core.lo),central_terms=detail,whole_code_certificate=False)


def verify_box(radius,fifth,fifth_min=0,first=(185,188,83)):
    assert len(first)==3 and 0<=fifth_min<=fifth
    bounds=tuple((x-radius,x+radius) for x in first)+((0,512),(fifth_min,fifth))
    result=verify_region(bounds)
    result.update(first=first,radius=radius,fifth_min=fifth_min,fifth_max=fifth)
    return result


def verify_histogram(h,log_weights=LOG_WEIGHT_NUMERATORS):
    assert len(h)==5 and sum(h)==512 and min(h)>=0
    u=[x/10**9 for x in log_weights]
    floating,_,_=local_data(h,u,.1,.2,.08)
    _,proposal,_=variational(h,u,.1,.2,.08);center,_=moments(h,tables(proposal))
    center=[I(str(round(x,12))) for x in center[:4]];center.append(I(512)-sum(center))
    entries,core,gradient,beta,kappa=interval_data(h,center,log_weights)
    priority=tuple(sorted((a for a in range(5) if h[a]),key=lambda a:(h[a],a)))
    matrices={};options={};remaining=255
    for index,a in enumerate(priority):
        if remaining<=0:break
        remaining-=min(h[a],remaining)
        subset=tuple(sorted(priority[index:]));candidate=affine_ellipsoid(floating,subset)
        matrix=[[I(str(candidate[i][j]))*I('1.000001')+I('0.000000000001')*int(i==j) for j in range(5)] for i in range(5)]
        matrices[subset]=matrix;options[subset]=checked_switches(floating,entries,subset,matrix)
    total=I();records=[]
    for lower,radius,slack in ((0,4,'0.2'),(4,5,'0.2'),(5,6,'0.4'),(6,8,'0.4'),(8,12,'0.4'),(12,16,'0.4')):
        value=shell(h,entries,core,gradient,matrices,options,radius,lower,slack,priority)
        total+=value.exp();records.append((lower,radius,str(value.hi)))
    # Global tail bound, independent of the local Taylor envelopes.
    full=tuple(sorted(priority));delta=I(max((x[2] for x in options[full]),default=D(0)))
    lip=4*I(2).sqrt()*beta;steps=255
    a=I('0.5')-kappa-steps*lip.sq()/8
    b=length(gradient)+steps*delta*lip/4;d=steps*delta.sq()/8
    effective=a-b/16;assert effective.lo>0
    tail=core+d+(1+effective*256).ln()-(4*effective.sq()).ln()-effective*256
    total+=tail.exp();result=total.ln();allowance=allowance_interval(log_weights)
    hypothetical_exponent=4095*(result-allowance)/I(2).ln()
    return dict(histogram=h,shell_log_upper_bounds=records,tail_log_upper=str(tail.hi),
               log_weight_numerators=tuple(log_weights),log_weight_denominator=10**9,
               directed_log_moment_upper=str(result.hi),allowance_lower=str(allowance.lo),
               exponent_if_this_bound_were_uniform=str(hypothetical_exponent.hi),
               passes=result.hi<allowance.lo,uniform_histogram_certificate=False)


def main():
    global VERIFY_DEADLINE
    assert shutil.disk_usage(Path.cwd()).free>=20*(1<<30),'disk guard'
    if '--precision80' in sys.argv:getcontext().prec=80
    started=time.monotonic();VERIFY_DEADLINE=started+25
    checks=arithmetic_tests()+polynomial_tests()+phase_and_shell_tests();base=(185,188,83,56,0)
    result=verify_histogram(base);print(result,flush=True);assert result['passes']
    if '--neighbors' in sys.argv:
        count=0;failures=[];worst=result
        for a in range(4):
            for b in range(5):
                if a==b:continue
                assert time.monotonic()-started<24
                h=list(base);h[a]-=1;h[b]+=1
                current=verify_histogram(tuple(h));count+=1
                if D(current['directed_log_moment_upper'])>D(worst['directed_log_moment_upper']):worst=current
                if not current['passes']:failures.append((h,current['directed_log_moment_upper']))
        print(dict(neighbors_checked=count,failures=failures,worst_histogram=worst['histogram'],
                   worst_directed_log_upper=worst['directed_log_moment_upper']),flush=True)
        assert not failures
    if '--boxes' in sys.argv:
        for radius,fifth,fifth_min in ((1,0,0),(2,0,0),(2,2,2),(4,8,4),(8,16,8)):
            assert time.monotonic()-started<24
            result=verify_box(radius,fifth,fifth_min);print(result,flush=True)
            if (radius,fifth)!=(2,0):assert result['passes']
    if '--boundary-cover' in sys.argv:
        bounds=[];failures=[]
        for offset in product((-1,1),repeat=3):
            assert time.monotonic()-started<24
            first=tuple(x+y for x,y in zip(base[:3],offset))
            result=verify_box(1,0,0,first);bounds.append(D(result['directed_log_upper']))
            if not result['passes']:failures.append(first)
            print(dict(first=first,log_upper=result['directed_log_upper'],passes=result['passes']),flush=True)
        print(dict(boundary_histograms=125,boxes=8,failed_boxes=failures,
                   uniform_log_upper=str(max(bounds))),flush=True)
        assert not failures
    if '--phase-cover' in sys.argv:
        for radius,fifth,first in ((2,48,(48,48,48)),(2,32,(64,320,64)),(4,100,(100,100,100)),(4,64,(64,64,64))):
            assert time.monotonic()-started<24
            result=verify_box(radius,fifth,fifth,first);print(result,flush=True);assert result['passes']
    assert time.monotonic()-started<25
    print(dict(seconds=time.monotonic()-started,precision=getcontext().prec,arithmetic_checks=checks))


if __name__=='__main__':main()
