"""GF16 shell bounds for the shared-four-row packet ensemble.

The four binary rows form the scalar extension of C to GF16. Symbol weight
is their union weight, and its dual is the scalar extension of C-perp.
The bounds here concern genuine shells, unlike differences of CDF caps.
Floating LPs propose dual multipliers; rational residual checks certify them.
"""
import argparse
import json
from fractions import Fraction as Q
from math import comb, isfinite, log2
from pathlib import Path

import numpy as np
from scipy.optimize import linprog
from scipy.sparse import csr_matrix

import shared_support
from constraint_counts import normalize, verify_dual


def krawtchouk(n, u, degree, q=16):
    if (any(type(v) is not int for v in (n, u, degree, q))
            or q < 2 or q & (q-1) or not 0 <= u <= n or not 0 <= degree <= n):
        raise ValueError('binary-extension alphabet and valid integer weights required')
    values = [1]
    if degree:
        values.append((q-1)*n-q*u)
    for j in range(1, degree):
        numerator = ((q-1)*n-q*u-(q-2)*j)*values[j]-(q-1)*(n-j+1)*values[j-1]
        value, remainder = divmod(numerator, j+1)
        assert remainder == 0
        values.append(value)
    return values


def shell_caps(n, k, dual_distance, q=16):
    """Christoffel bound from degree<dual_distance uniform coordinate moments."""
    if (any(type(v) is not int for v in (n, k, dual_distance, q))
            or q < 2 or q & (q-1) or not 0 <= k <= n or not 1 <= dual_distance <= n+1):
        raise ValueError('valid scalar-extension code parameters required')
    degree = (dual_distance-1)//2
    result = []
    for u in range(n+1):
        values = krawtchouk(n, u, degree, q)
        kernel = sum((Q(v*v, (q-1)**j*comb(n,j)) for j,v in enumerate(values)), Q(0))
        cap = Q(q**k)/kernel
        result.append(cap.numerator//cap.denominator)
    return result


class Constraints:
    def __init__(self, cdfs, shells, k, q=16, last_transform=None):
        if (len(cdfs) != 2 or len(shells) != 2 or not cdfs[0]
                or any(len(row) != len(cdfs[0]) for row in (*cdfs,*shells))):
            raise ValueError('two matching CDF and shell-cap arrays required')
        n = len(cdfs[0])-1
        if type(k) is not int or not 0 <= k <= n or type(q) is not int or q<2 or q&(q-1):
            raise ValueError('valid dimension and binary-extension alphabet required')
        if last_transform is None:last_transform=n
        if type(last_transform) is not int or not 0<=last_transform<=n:
            raise ValueError('valid maximum transform degree required')
        for side, dimension in enumerate((k,n-k)):
            cdf, caps = cdfs[side], shells[side]
            if (cdf[0] != 1 or caps[0] != 1 or cdf[-1] != q**dimension
                    or any(type(v) is not int or v < 0 for row in (cdf,caps) for v in row)
                    or any(a > b for a,b in zip(cdf,cdf[1:]))):
                raise ValueError('integer CDF including zero with exact final mass required')
        self.n,self.k,self.q = n,k,q
        self.cdfs = cdfs
        self.caps = [[1]+[min(a-1,b) for a,b in zip(cdf[1:],row[1:])]
                     for cdf,row in zip(cdfs,shells)]
        self.keys = [(side,u) for side in range(2) for u in range(1,n+1) if self.caps[side][u]]
        self.indices = {key:i for i,key in enumerate(self.keys)}
        self.rows = []
        self.kw = [krawtchouk(n,u,n,q) for u in range(n+1)]

        def term(row,side,u,c):
            if (side,u) in self.indices:
                i = self.indices[side,u]
                row[i] = row.get(i,0)+c*self.caps[side][u]

        def equality(row,rhs):
            self.rows.append(normalize(row,rhs))
            self.rows.append(normalize({i:-v for i,v in row.items()},-rhs))

        for j in range(last_transform+1):
            row = {}
            for u in range(1,n+1):
                term(row,0,u,self.kw[u][j])
            term(row,1,j,-q**k)
            equality(row,q**k*(j==0)-self.kw[0][j])
        for side,dimension in enumerate((k,n-k)):
            row = {}
            for u in range(1,n+1):
                term(row,side,u,1)
                self.rows.append(normalize(row,cdfs[side][u]-1))
            equality(row,q**dimension-1)
        values,ri,ci = [],[],[]
        for j,(row,_) in enumerate(self.rows):
            for i,value in row.items():
                values.append(float(value));ri.append(j);ci.append(i)
        self.matrix = csr_matrix((values,(ri,ci)),shape=(len(self.rows),len(self.keys)))
        self.rhs = np.array([float(b) for _,b in self.rows])

    def check_values(self,enumerators):
        values = [Q(enumerators[side][u],self.caps[side][u]) for side,u in self.keys]
        assert all(0<=v<=1 for v in values)
        for side in range(2):
            assert all(not enumerators[side][u] for u in range(1,self.n+1) if (side,u) not in self.indices)
        for row,rhs in self.rows:
            assert sum(v*values[i] for i,v in row.items()) <= rhs

    def verify(self,support,entries):
        if type(support) is not int or not 1<=support<=self.n:
            raise ValueError('valid nonzero support required')
        if (0,support) not in self.indices:
            return 0,Q(0)
        multipliers = [Q(0)]*len(self.rows)
        seen = set()
        for i,value in entries:
            if type(i) is not int or not 0<=i<len(self.rows) or i in seen or Q(value)<0:
                raise ValueError('distinct valid indices and nonnegative multipliers required')
            seen.add(i);multipliers[i] = Q(value)
        upper,correction = verify_dual(self.rows,[],{self.indices[0,support]:Q(1)},multipliers,[],len(self.keys))
        if upper<0:
            raise ArithmeticError('dual contradicts coding premises')
        cap = min(Q(1),upper)*self.caps[0][support]
        return cap.numerator//cap.denominator,correction

    def transform_cap(self,support):
        """Inverse transform with a nonnegative decreasing CDF majorant."""
        if type(support) is not int or not 1<=support<=self.n:
            raise ValueError('valid nonzero support required')
        weights=[self.kw[j][support] for j in range(self.n+1)]
        high=0;majorant=[0]*(self.n+1)
        for j in range(self.n,0,-1):
            high=max(high,weights[j]);majorant[j]=high
        cdf=self.cdfs[1]
        by_cdf=weights[0]+sum((cdf[j]-cdf[j-1])*majorant[j] for j in range(1,self.n+1))
        by_shell=weights[0]+sum(self.caps[1][j]*max(0,weights[j]) for j in range(1,self.n+1))
        return min(self.caps[0][support],min(by_cdf,by_shell)//self.q**(self.n-self.k))

    def solve(self,support,seconds=5.):
        if not isfinite(seconds) or not 0<seconds<=120:
            raise ValueError('positive bounded LP time required')
        prior,_ = self.verify(support,[])
        result = dict(support=support,prior=prior,cap=prior,attempts=[])
        if not prior:return result
        objective = np.zeros(len(self.keys));objective[self.indices[0,support]]=-1
        for cutoff in (0,1e-9):
            matrix = self.matrix.copy();rhs = self.rhs.copy()
            if cutoff:
                removed = matrix.copy()
                removed.data = np.where(np.abs(removed.data)<cutoff,np.minimum(removed.data,0),0)
                rhs -= np.asarray(removed.sum(axis=1)).ravel()
                matrix.data[np.abs(matrix.data)<cutoff]=0
                matrix.eliminate_zeros();rhs+=1e-8
            fit = linprog(objective,A_ub=matrix,b_ub=rhs,bounds=(0,1),method='highs',options={'time_limit':seconds})
            attempt = dict(cutoff=cutoff,status=int(fit.status),message=fit.message)
            result['attempts'].append(attempt)
            if not fit.success:continue
            entries = [[i,str(Q(round(max(0.,-float(v))*(1<<48)),1<<48))]
                       for i,v in enumerate(fit.ineqlin.marginals) if v<0]
            cap,correction = self.verify(support,entries)
            attempt.update(proposal=float(-fit.fun),cap=cap,correction=str(correction))
            if cap<result['cap']:result.update(cap=cap,dual=entries)
        return result


def actual(coupled=True):
    from count_refinements import refine_pair
    cdf,ranks = shared_support.shared_counts(True,coupled)
    primal,dual,_,_ = refine_pair(ranks,shared_support.authenticated_caps())
    cdfs = [[1+sum(row[u] for row in rs) for u in range(257)] for rs in (primal,dual)]
    shells = [shell_caps(256,128,d) for d in (30,38)]
    for row in shells:row[0]=1
    return Constraints(cdfs,shells,128)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--supports',type=int,nargs='+',default=[112,120,128,144,160,176,192,224,240])
    parser.add_argument('--seconds',type=float,default=5.)
    parser.add_argument('--cached-model',type=Path,help='Exploratory re-use of an earlier authenticated CDF receipt')
    parser.add_argument('--last-transform',type=int)
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args()
    if args.output.exists():parser.error('refusing to overwrite an existing receipt')
    if args.cached_model:
        record=json.loads(args.cached_model.read_text())
        model=Constraints(record['cdfs'],record['shell_caps'],128,last_transform=args.last_transform)
    else:
        model=actual()
        if args.last_transform is not None:
            model=Constraints(model.cdfs,model.caps,128,last_transform=args.last_transform)
    print('GF16 MACWILLIAMS MODEL',len(model.keys),'variables',len(model.rows),'inequalities',flush=True)
    rows=[]
    for u in args.supports:
        row=model.solve(u,args.seconds);rows.append(row)
        row['transform_cap']=model.transform_cap(u)
        original=model.cdfs[0][u]-1
        print('SHELL',u,'CDF bits',log2(original) if original else None,'prior bits',log2(row['prior']) if row['prior'] else None,
              'verified bits',log2(row['cap']) if row['cap'] else None,'statuses',[v['status'] for v in row['attempts']],flush=True)
        args.output.parent.mkdir(parents=True,exist_ok=True)
        args.output.write_text(json.dumps(dict(schema='shared-relaxed-gf16-shell-lp-1',q=16,
            cdfs=model.cdfs,shell_caps=model.caps,rows=rows,last_transform=args.last_transform,
            cached_model=str(args.cached_model) if args.cached_model else None,
            proof_status='exact shell caps only; not a distance certificate'),indent=2)+'\n')


if __name__=='__main__':main()
