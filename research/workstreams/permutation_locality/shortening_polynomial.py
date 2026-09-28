"""Exact positive-polynomial witnesses for binary shortening dimensions.

Floating roots propose a polynomial only. Its coefficients, sign on all
allowed distances, and final bound are verified over the rationals.
"""
from fractions import Fraction as Q
from math import comb, log2
import numpy as np
from scipy.linalg import eigvalsh_tridiagonal


def values(n, x, degree):
    k = [Q(1)]
    if degree:
        k.append(n-2*x)
    for j in range(1,degree):
        k.append(((n-2*x)*k[-1]-(n-j+1)*k[-2])/(j+1))
    return k


def witness(n, d, degree, root):
    """f(x)=(d-x)p(x)^2; p and (d-x)p have nonnegative K coefficients."""
    assert 0<=degree<n and 0<d<=n
    kernel = values(n,root,degree)
    p = [v/comb(n,j) for j,v in enumerate(kernel)]
    if any(v<0 for v in p):
        return None
    g = []
    for j in range(degree+2):
        value = (d-Q(n,2))*p[j] if j<=degree else Q(0)
        if j:
            value += Q(j,2)*p[j-1]
        if j+1<=degree:
            value += Q(n-j,2)*p[j+1]
        g.append(value)
    if any(v<0 for v in g):
        return None
    # Krawtchouk products have nonnegative integer linearization
    # coefficients. Thus p*g has nonnegative K coefficients.
    constant = sum(p[j]*g[j]*comb(n,j) for j in range(degree+1))
    if constant<=0:
        return None
    at_zero = sum(p[j]*comb(n,j) for j in range(degree+1))
    upper = d*at_zero**2/constant
    return upper,p,g


def verify(n,d,result):
    """Independently expand f by orthogonality, at every binary distance."""
    upper,p,g = result
    rows = [values(n,Q(x),n) for x in range(n+1)]
    evaluations = []
    for x,row in enumerate(rows):
        px = sum(c*k for c,k in zip(p,row))
        gx = sum(c*k for c,k in zip(g,row))
        assert gx == (d-x)*px
        evaluations.append((d-x)*px*px)
    assert all(v<=0 for v in evaluations[d:])
    coefficients = [sum(comb(n,x)*evaluations[x]*rows[x][j] for x in range(n+1)) /
                    (2**n*comb(n,j)) for j in range(n+1)]
    assert all(v>=0 for v in coefficients) and coefficients[0]>0
    assert upper == evaluations[0]/coefficients[0]


def bound(n,d=38,maximum_degree=40):
    if n<d:
        return Q(1),None
    best = None
    for degree in range(min(n-1,maximum_degree)+1):
        off = -.5*np.sqrt(np.arange(1,degree+1)*(n-np.arange(degree)))
        root = float(eigvalsh_tridiagonal(np.full(degree+1,n/2),off,
                     select='i',select_range=(0,0))[0]) if degree else n/2
        if root>=d:
            continue
        candidate = witness(n,d,degree,Q(root).limit_denominator(10**10))
        if candidate is not None and (best is None or candidate[0]<best[0]):
            best = candidate
    if best is None:
        return Q(1<<n),None
    return best[0],best


def self_test():
    for n,d,size in ((4,2,8),(8,4,16),(7,3,16),(12,4,128),(8,8,2)):
        upper,result = bound(n,d,maximum_degree=n-1)
        assert upper>=size
        assert result is not None
        verify(n,d,result)
    print('Positive-polynomial witnesses: exact full expansions passed on small codes',flush=True)


def improve_dimensions(dimensions,d=38,last=192):
    """Retain prior valid bounds; accept only exact positive witnesses."""
    result = dimensions[:]
    accepted = 0
    for n in range(1,len(result)):
        result[n] = min(result[n],result[n-1]+1)
        if d<=n<=last:
            upper,witness_data = bound(n,d)
            if witness_data is not None:
                accepted += 1
                while 1<<result[n]>upper:
                    result[n] -= 1
    print('Positive shortening polynomials:',accepted,'exact positive witnesses;',
          sum(a>b for a,b in zip(dimensions,result)),'dimension caps improved',flush=True)
    return result


if __name__=='__main__':
    self_test()
    for n in (80,96,104,112,128,144,160,176,192):
        upper,result = bound(n)
        if result is not None:
            verify(n,38,result)
        dimension = 0
        while 1<<(dimension+1)<=upper:
            dimension += 1
        print(n,'log2 bound',log2(upper.numerator)-log2(upper.denominator),
              'dimension cap',dimension,'degree',len(result[1])-1 if result else None,flush=True)
