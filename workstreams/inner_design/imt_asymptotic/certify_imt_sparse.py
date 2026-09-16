"""Exact rational polynomial certificate for the IMT sparse Collatz inequality.

This is a local transfer certificate, not the complete asymptotic theorem.
The polynomial variable x=10000*alpha ranges over [0,1].
"""
import argparse
from fractions import Fraction as F
import hashlib
import json
import math
from pathlib import Path
import sys

from flint import fmpq, fmpq_poly, ctx
import screen
import sparse
import certify_imt_dense as dense


def rat(x):
    x = F(x)
    return fmpq(x.numerator,x.denominator)


def polynomial(x):
    return fmpq_poly([rat(x)])


def sign_certificate(poly,strict=False):
    """Certify p(x)<=0 on [0,1], or p(x)<0 on (0,1] when strict.

    Remove an exact zero at x=0 before bounding. First try a positive-tail
    power bound, then the exact Bernstein convex-hull bound.
    """
    coefficients = list(poly.coeffs())
    if not coefficients:
        assert not strict
        return dict(method='zero',zero_order=0,upper='0')
    order = 0
    while coefficients and coefficients[0]==0:
        coefficients.pop(0);order+=1
    assert coefficients
    upper = coefficients[0]+sum((max(c,fmpq(0)) for c in coefficients[1:]),fmpq(0))
    if upper<0 or (not strict and upper<=0):
        return dict(method='positive_power_tail',zero_order=order,upper=str(upper))
    degree = len(coefficients)-1
    bernstein = []
    for k in range(degree+1):
        bernstein.append(sum((coefficients[i]*fmpq(math.comb(k,i),math.comb(degree,i)) for i in range(k+1)),fmpq(0)))
    upper = max(bernstein)
    assert upper<0 if strict else upper<=0,('polynomial sign unresolved',order,degree,float(upper))
    return dict(method='bernstein',zero_order=order,upper=str(upper))


def best_upper(polys,checks):
    """Choose a polynomial majorizing every candidate on the entire interval."""
    choice = max(polys,key=lambda p:p(fmpq(1,2)))
    for p in polys:
        checks.append(sign_certificate(p-choice))
    return choice


def smallest_bound(polys):
    # Each input is already a valid upper bound. No claim that the selected
    # one is minimal throughout the interval is needed.
    return min(polys,key=lambda p:p(fmpq(1,2)))


def build(model,gamma=96):
    e = model.engine
    n,m = e.n,e.m
    v,h,_ = sparse.witness(model)
    x = fmpq_poly([0,1])
    alpha = x/10000
    z = 1-rat(F(8,5))*alpha
    beta = rat(F(4,5))*alpha
    zp = [z**i for i in range(129)]
    bp = [beta**i for i in range(129)]
    cp = [(1-beta)**i for i in range(129)]
    vector = [polynomial(1)]+[rat(v)*(1+rat(hi)*alpha) for hi in h]
    for w in vector:
        sign_certificate(-w,strict=True)
    live_average = sum((vector[i+2]*int(count) for i,count in enumerate(model.counts)),fmpq_poly())/m
    assert live_average==polynomial(v)
    rows = [fmpq_poly() for _ in range(n)]
    dominance = []
    for j in range(129):
        total = math.comb(128,j)
        nonzero = total-e.kernel[j]
        live = polynomial(F(nonzero,total))
        moments = []
        for w in e.levels:
            moment = sum((math.comb(w,k)*math.comb(128-w,j-k)*zp[w+j-2*k]
                          for k in range(max(0,j-128+w),min(w,j)+1)),fmpq_poly())/total
            moments.append(moment)
        arbitrary = best_upper(moments,dominance)
        cap = int(e.caps[j]['cap'])
        bounds = [arbitrary,live,rat(F(cap,total))*zp[min(abs(w-j) for w in e.levels)]]
        if j in e.low:
            patterns = [sum((count*zp[w] for w,count in pattern),fmpq_poly())/total
                        for pattern in e.low[j]['patterns']]
            bounds.append(best_upper(patterns,dominance) if patterns else fmpq_poly())
        cancel = smallest_bound(bounds)
        fresh_cancel = smallest_bound([arbitrary,live])
        action = [rat(F(e.kernel[j],total))*zp[j]+live*zp[j]*vector[1],
                  cancel/2+fresh_cancel/(2*m)+arbitrary*(vector[1]+live_average)/2]
        for i,w in enumerate(e.levels):
            moment = moments[i]
            mass = min(nonzero,e.spectrum[w]*cap)
            bounds = [moment,rat(F(mass,e.spectrum[w]*total))*zp[abs(w-j)]]
            if j in e.low:
                bounds.append(sum((count*zp[out] for out,count in e.low[j]['by_weight'].get(w,{}).items()),fmpq_poly())/(e.spectrum[w]*total))
            cancel = smallest_bound(bounds)
            fresh_cancel = smallest_bound([moment,live])
            lazy_target = vector[1] if nonzero else vector[i+2]
            action.append(cancel/2+fresh_cancel/(2*m)+moment*(live_average+lazy_target)/2)
        probability = total*bp[j]*cp[128-j]
        for i in range(n):
            rows[i] += probability*action[i]
        if j%32==0:
            print('polynomial input weight',j,flush=True)
    residuals = [r-(1-gamma*alpha)*w for r,w in zip(rows,vector)]
    return vector,rows,residuals,dominance


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output',type=Path,required=True)
    p.add_argument('--verify',type=Path)
    args = p.parse_args()
    if args.output.exists():
        raise FileExistsError(args.output)
    saved = json.loads(args.verify.read_text()) if args.verify else None
    if saved:
        for name,digest in saved['source_sha256'].items():
            assert screen.sha(screen.ROOT/name)==digest,name
    model = screen.Model()
    vector,rows,residuals,dominance = build(model)
    checks = []
    for i,poly in enumerate(residuals):
        assert poly[0]==0
        proof = sign_certificate(poly,strict=True)
        record = dict(row=i,degree=poly.degree(),coefficients_sha256=hashlib.sha256(str(poly).encode()).hexdigest(),**proof)
        checks.append(record)
        print('row',i,'degree',poly.degree(),proof['method'],'upper',float(F(proof['upper'])),flush=True)
    ctx.prec=256
    num=dense.number
    # A deliberately conservative uniform bound on the natural exponent
    # coefficient, including the existing outer likelihood and support costs.
    coefficient = num(2).log()/2+num(F(5,8)).log()+num(F(3,5)+F(11,100)*F(8,5))/num(1-F(8,5)*F(1,10000))-num(F(96,128))+num(F(1281,100000))+num(2).log()/num(F(39,4))
    coefficient_upper = dense.upper(coefficient)
    assert coefficient_upper<0
    if saved:
        assert checks==saved['row_checks']
        assert saved['inner']==model.engine.identity()['inner']
    paths=[Path(__file__),Path(screen.__file__),Path(sparse.__file__),Path(dense.__file__),screen.HERE.parent/'NO_CONSTANT_MAP.json']
    paths += [Path(mod.__file__).resolve() for mod in list(sys.modules.values()) if getattr(mod,'__file__',None)
              and Path(mod.__file__).resolve().is_relative_to(screen.ROOT) and Path(mod.__file__).suffix=='.py']
    result=dict(status='EXACT_IMT_SPARSE_COLLATZ_POLYNOMIAL_CERTIFICATE',
                alpha_interval=['0','1/10000'],strict_for_positive_alpha=True,
                contraction='1-96*alpha',beta='(4/5)*alpha',z='1-(8/5)*alpha',
                inner=model.engine.identity()['inner'],vector_coefficients=[[str(c) for c in p.coeffs()] for p in vector],
                row_checks=checks,dominance_comparisons=len(dominance),
                conditional_structured_exponent_coefficient_upper=str(coefficient_upper),
                source_sha256={q.relative_to(screen.ROOT).as_posix():screen.sha(q) for q in set(paths)},
                full_asymptotic_theorem=False,
                scope='Uniform iid-input transfer inequality; the uniform sparse route reduction and fixed-occupancy limits remain separate obligations.')
    with args.output.open('x',encoding='utf-8') as stream:
        json.dump(result,stream,indent=2);stream.write('\n')
    print('PASS: uniform sparse Collatz; conditional exponent coefficient',float(coefficient_upper))


if __name__=='__main__':
    main()
