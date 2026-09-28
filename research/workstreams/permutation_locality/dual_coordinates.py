"""Exact change of BCH sandwich coordinates for dual-shell LP witnesses.

p is A(P^perp); a is (A(Q^perp)-A(P^perp))/255, the average of
nonzero coset spectra. Then A(C^perp)=p+7a. No claim of dual-coset
transitivity is needed. Every accepted bound uses exact residual repair.
"""
import argparse
from fractions import Fraction as Q
from math import comb, log2
import numpy as np
from scipy.optimize import linprog
from flint import fmpz_mat
from dual_sandwich import raw_model, repaired_upper, price_upper, exact_solver
from dual_moments import dual_shell_caps
from shortened_bound import krawtchouk


def folded_transform(n):
    assert n%4==0
    weights=list(range(0,n//2+1,2))
    columns={w:krawtchouk(n,w) for w in weights}
    return weights,[[columns[w][j]*(2 if w<n//2 else 1) for w in weights] for j in weights]


def transform_checks():
    checks=0
    for n in (4,8,12,16,256):
        weights,matrix=folded_transform(n)
        size=len(weights);square=fmpz_mat(matrix)**2
        assert square==fmpz_mat([[int(i==j)*(1<<n) for j in range(size)] for i in range(size)])
        checks+=size*size
    # Block transform identities hold independently of a spectrum's values.
    # [[1,255],[1,-1]] squares to 256*I, with 131+133=256+8.
    mixing=fmpz_mat([[1,255],[1,-1]])
    assert mixing**2==fmpz_mat([[256,0],[0,256]])
    # C=q+31h and D=2^-128 K C=p+7a under the inverse map.
    assert 1+31==32 and 255-31==7*32
    # Exhaustive actual-code check: repetition Q inside first-order RM P.
    # Its nonzero Q-cosets all have the same (middle-weight) spectrum.
    from joint_support import span
    for n in (8,16):
        ones=(1<<n)-1
        linear=[sum(((j>>bit)&1)<<j for j in range(n)) for bit in range(n.bit_length()-1)]
        qb=[ones];pb=qb+linear;cb=pb[:3]
        words=lambda basis:span(basis)
        dual=lambda basis:[x for x in range(1<<n) if all((x&b).bit_count()%2==0 for b in basis)]
        spectrum=lambda xs:[sum(x.bit_count()==w for x in xs) for w in range(n+1)]
        q,p,c=map(lambda b:spectrum(words(b)),(qb,pb,cb))
        dq,dp,dc=map(lambda b:spectrum(dual(b)),(qb,pb,cb))
        quotient=1<<(len(pb)-len(qb));selected=1<<(len(cb)-len(qb))
        h=[Q(y-x,quotient-1) for x,y in zip(q,p)]
        a=[Q(y-x,quotient-1) for x,y in zip(dp,dq)]
        assert c==[x+(selected-1)*y for x,y in zip(q,h)]
        assert dc==[x+(quotient//selected-1)*y for x,y in zip(dp,a)]
        assert all(x>=0 for x in a)
        ws,mat=folded_transform(n)
        for original,dual_spectrum,dimension in ((q,dq,len(qb)),(p,dp,len(pb)),(c,dc,len(cb))):
            result=[sum(mat[i][j]*original[w] for j,w in enumerate(ws)) for i in range(len(ws))]
            assert result==[dual_spectrum[w]*(1<<dimension) for w in ws]
            checks+=len(ws)
    print('Exact folded/inverse-transform and actual-code checks:',checks,'entries passed',flush=True)


def model():
    oldvars,raw,oldfixed,caps,table=raw_model()
    from bch_quotient import generator_polynomial,binary_poly_divmod
    assert binary_poly_divmod((1<<255)-1,generator_polynomial(39))[1]==0
    assert binary_poly_divmod((1<<255)-1,generator_polynomial(37))[1]==0
    weights,matrix=folded_transform(256)
    allvars=[f'{p}_{w}' for p in ('p','a') for w in weights]
    # q=2^-133 K(p+255a), h=2^-133 K(p-a).
    inverse=[]
    for old in oldvars:
        prefix,weight=old.split('_');row=matrix[weights.index(int(weight))]
        inverse.append(row+[x*(255 if prefix=='q' else -1) for x in row])
    # Old nonnegativity must survive the coordinate change.
    raw=list(raw)+[dict(name=f'retain_nonnegative_{v}',sense='ge',rhs=0,coeffs={v:1}) for v in oldvars]
    raw += [dict(name=f'retain_fixed_{v}',sense='eq',rhs=x,coeffs={v:1}) for v,x in oldfixed.items()]
    transformed=fmpz_mat([[int(r['coeffs'].get(v,0)) for v in oldvars] for r in raw])*fmpz_mat(inverse)
    newraw=[]
    for i,row in enumerate(raw):
        newraw.append((row['name'],row['sense'],Q(int(row['rhs'])),
                       [Q(int(transformed[i,j]),1<<133) for j in range(len(allvars))]))
    dcaps=dual_shell_caps()
    for w in weights:
        coeff=[Q(int(v==f'p_{w}')+7*int(v==f'a_{w}')) for v in allvars]
        newraw.append((f'dual_cap_{w}','le',Q(dcaps[w]),coeff))
    # Q^perp has minimum 30, P^perp is contained in it. Q and P
    # contain ones and are even, so both dual spectra are even/symmetric.
    fixed={f'{p}_{w}':int(p=='p' and w==0) for p in ('p','a') for w in weights if w<30}
    variables=[v for v in allvars if v not in fixed]
    scales=[max(1,comb(256,int(v.split('_')[1]))//(1<<130)) for v in variables]
    upper=[Q(dcaps[int(v.split('_')[1])],s*(7 if v.startswith('a_') else 1)) for v,s in zip(variables,scales)]
    positions=[allvars.index(v) for v in variables]
    inequalities=[];equalities=[];seen=set()
    for name,sense,rhs,row in newraw:
        rhs-=sum((row[allvars.index(v)]*x for v,x in fixed.items()),Q(0))
        row=[row[j]*s for j,s in zip(positions,scales)]
        if sense=='ge':row=[-x for x in row];rhs=-rhs;sense='le'
        if not any(row):
            assert rhs==0 if sense=='eq' else rhs>=0,(name,rhs)
            continue
        scale=max([abs(rhs)]+[abs(x) for x in row])
        row=[x/scale for x in row];rhs/=scale
        # Normalize equality signs before removing exact duplicate rows.
        if sense=='eq' and next(x for x in row if x)<0:row=[-x for x in row];rhs=-rhs
        key=(sense,tuple(row),rhs)
        if key in seen:continue
        seen.add(key)
        (equalities if sense=='eq' else inequalities).append((row,rhs,name))
    print('Dual coordinates:',len(variables),'variables,',len(inequalities),'inequalities,',len(equalities),'equalities',flush=True)
    return variables,scales,upper,inequalities,equalities


def shell(data,w,old,direct=False,limit=100,soplex=False,tight=False):
    variables,scales,upper,inequalities,equalities=data
    unit=1<<(old.bit_length()-1)
    objective=[Q(s*(int(v==f'p_{w}')+7*int(v==f'a_{w}')),unit) for v,s in zip(variables,scales)]
    if soplex or tight:
        exact_solver(objective,Q(0),upper,inequalities,equalities,variables,tight)
        return old
    if direct:
        rows=inequalities+equalities;ni=len(inequalities);ne=len(equalities);nv=len(variables)
        matrix=np.array([[float(x) for x in row] for row,_,_ in rows])
        fit=linprog(np.array([float(rhs) for _,rhs,_ in rows]+list(map(float,upper))),
                    A_ub=np.hstack((-matrix.transpose(),-np.eye(nv))),
                    b_ub=-np.array(list(map(float,objective))),
                    bounds=[(0,limit)]*ni+[(-limit,limit)]*ne+[(0,None)]*nv,
                    method='highs',options={'time_limit':30})
        if not fit.success:
            print('No direct coordinate witness',w,fit.message,flush=True);return old
        prices=[max(Q(0),Q(float(x))) for x in fit.x[:ni]]
        eqprices=[Q(float(x)) for x in fit.x[ni:ni+ne]]
        bound,correction=price_upper(objective,Q(0),upper,inequalities,equalities,prices,eqprices)
        numerical=fit.fun
    else:
        bound,correction,numerical=primal_witness(objective,upper,inequalities,equalities,w)
        if bound is None:return old
    assert bound>=0
    value=bound*unit
    checked=min(old,value.numerator//value.denominator)
    print('COORDINATE SHELL',w,'log2 old/checked',log2(old),log2(checked) if checked else '-inf',
          'numerical scaled',numerical,'repair scaled',float(correction),'checked scaled',float(bound),flush=True)
    return checked


def primal_witness(objective,upper,inequalities,equalities,w):
    result=linprog(-np.array(list(map(float,objective))),
                   A_ub=np.array([[float(x) for x in row] for row,_,_ in inequalities]),
                   b_ub=np.array([float(rhs) for _,rhs,_ in inequalities]),
                   A_eq=np.array([[float(x) for x in row] for row,_,_ in equalities]),
                   b_eq=np.array([float(rhs) for _,rhs,_ in equalities]),
                   bounds=[(0,float(x)) for x in upper],method='highs',
                   options={'time_limit':30,'dual_feasibility_tolerance':1e-9,'primal_feasibility_tolerance':1e-9})
    if not result.success:
        print('No coordinate witness',w,result.message,flush=True);return None,None,None
    bound,correction=repaired_upper(objective,Q(0),upper,inequalities,equalities,result)
    return bound,correction,-result.fun


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--weights',nargs='+',type=int,default=[30,38,56,64])
    parser.add_argument('--self-test-only',action='store_true')
    parser.add_argument('--direct-dual',action='store_true')
    parser.add_argument('--multiplier-limit',type=float,default=100)
    parser.add_argument('--soplex',action='store_true')
    parser.add_argument('--soplex-tight',action='store_true')
    args=parser.parse_args();transform_checks()
    if not args.self_test_only:
        data=model();old=dual_shell_caps()
        for w in args.weights:
            assert 30<=w<=128 and w%2==0
            shell(data,w,old[w],args.direct_dual,args.multiplier_limit,args.soplex,args.soplex_tight)
