"""Exact-checked quaternary MacWilliams constraints for BCH row pairs.

A pair x,y in a binary linear code C identifies x+omega*y in its scalar
extension to GF4. Its symbol weight is the pair union weight. The dual is
the scalar extension of C-perp. Floating LPs propose nonnegative dual
multipliers; exact rational box-residual verification supplies every bound.
"""
import argparse
import json
from fractions import Fraction as Q
from math import isfinite, log2
from pathlib import Path

import numpy as np
from scipy.optimize import linprog
from scipy.sparse import csr_matrix

import pairwise_support
import pairwise_moments
from constraint_counts import normalize, verify_dual


class Constraints:
    def __init__(self, cdfs, shells, k):
        if (len(cdfs) != 2 or len(shells) != 2 or not cdfs[0]
                or any(len(row) != len(cdfs[0]) for row in (*cdfs,*shells))):
            raise ValueError('two matching CDF and shell-cap arrays required')
        n = len(cdfs[0])-1
        if type(k) is not int or not 0 <= k <= n:
            raise ValueError('valid binary dimension required')
        for side, dimension in enumerate((k,n-k)):
            cdf, caps = cdfs[side], shells[side]
            if (cdf[0] != 1 or caps[0] != 1 or cdf[-1] != 1 << (2*dimension)
                    or any(type(v) is not int or v < 0 for row in (cdf,caps) for v in row)
                    or any(a > b for a,b in zip(cdf,cdf[1:]))):
                raise ValueError('valid pair CDF including zero and exact total required')
        self.n, self.k = n,k
        self.cdfs = cdfs
        self.caps = [[1]+[min(a-1,b) for a,b in zip(cdf[1:],row[1:])]
                     for cdf,row in zip(cdfs,shells)]
        self.keys = [(side,u) for side in range(2) for u in range(1,n+1) if self.caps[side][u]]
        self.indices = {key:i for i,key in enumerate(self.keys)}
        self.rows = []
        self.kw = [pairwise_moments.krawtchouk4(n,u,n) for u in range(n+1)]
        def term(row,side,u,c):
            if (side,u) in self.indices:
                i = self.indices[side,u]
                row[i] = row.get(i,0)+c*self.caps[side][u]
        def equality(row,rhs):
            self.rows.append(normalize(row,rhs))
            self.rows.append(normalize({i:-v for i,v in row.items()},-rhs))
        for j in range(n+1):
            row = {}
            for u in range(1,n+1):
                term(row,0,u,self.kw[u][j])
            term(row,1,j,-(1 << (2*k)))
            equality(row,(1 << (2*k))*(j == 0)-self.kw[0][j])
        for side,dimension in enumerate((k,n-k)):
            row = {}
            for u in range(1,n+1):
                term(row,side,u,1)
                self.rows.append(normalize(row,cdfs[side][u]-1))
            equality(row,(1 << (2*dimension))-1)
        values, ri, ci = [],[],[]
        for j,(row,_) in enumerate(self.rows):
            for i,value in row.items():
                values.append(float(value));ri.append(j);ci.append(i)
        self.matrix = csr_matrix((values,(ri,ci)),shape=(len(self.rows),len(self.keys)))
        self.rhs = np.array([float(b) for _,b in self.rows])

    def check_values(self, enumerators):
        """Exact small-code audit of every row, independent of the LP."""
        values = [Q(enumerators[side][u],self.caps[side][u]) for side,u in self.keys]
        assert all(0 <= v <= 1 for v in values)
        for side in range(2):
            assert all(not enumerators[side][u] for u in range(1,self.n+1) if (side,u) not in self.indices)
        for row,rhs in self.rows:
            assert sum(v*values[i] for i,v in row.items()) <= rhs

    def transform_cap(self, support):
        """Inverse transform with a decreasing majorant; no floating LP.

        B_0=1 is exact. For j>0 replace K_support(j) by its nonnegative
        suffix maximum. Summation by parts then bounds the functional
        using the dual CDF, without treating CDF differences as shell caps.
        A separate positive-coefficient sum uses actual shell caps.
        """
        if type(support) is not int or not 1 <= support <= self.n:
            raise ValueError('nonzero support within block required')
        weights=[self.kw[j][support] for j in range(self.n+1)]
        high=0; majorant=[0]*(self.n+1)
        for j in range(self.n,0,-1):
            high=max(high,weights[j]);majorant[j]=high
        cdf=self.cdfs[1]
        by_cdf=weights[0]+sum((cdf[j]-cdf[j-1])*majorant[j]
                             for j in range(1,self.n+1))
        by_shell=weights[0]+sum(self.caps[1][j]*max(0,weights[j])
                               for j in range(1,self.n+1))
        return min(self.caps[0][support],min(by_cdf,by_shell)//(1 << (2*(self.n-self.k))))

    def verify(self, support, entries):
        if type(support) is not int or not 1 <= support <= self.n:
            raise ValueError('nonzero support within block required')
        multipliers = [Q(0)]*len(self.rows)
        seen = set()
        for i,value in entries:
            if type(i) is not int or not 0 <= i < len(self.rows) or i in seen or Q(value) < 0:
                raise ValueError('distinct valid indices and nonnegative dual multipliers required')
            seen.add(i);multipliers[i] = Q(value)
        if (0,support) not in self.indices:
            return 0,Q(0)
        upper, correction = verify_dual(self.rows,[],{self.indices[0,support]:Q(1)},
                                       multipliers,[],len(self.keys))
        if upper < 0:
            raise ArithmeticError('dual contradicts the coding premises')
        value = min(Q(1),upper)*self.caps[0][support]
        return value.numerator//value.denominator,correction

    def solve(self, support, seconds=3.):
        if not isfinite(seconds) or not 0 < seconds <= 60:
            raise ValueError('bounded positive LP time limit required')
        prior,_ = self.verify(support,[])
        result = dict(support=support,prior=prior,cap=prior,attempts=[])
        if not prior:
            return result
        objective = np.zeros(len(self.keys));objective[self.indices[0,support]]=-1
        for cutoff in (0,1e-9):
            matrix = self.matrix.copy();rhs = self.rhs.copy()
            if cutoff:
                removed = matrix.copy()
                removed.data = np.where(np.abs(removed.data)<cutoff,np.minimum(removed.data,0),0)
                rhs -= np.asarray(removed.sum(axis=1)).ravel()
                matrix.data[np.abs(matrix.data)<cutoff]=0
                matrix.eliminate_zeros();rhs += 1e-8
            fit = linprog(objective,A_ub=matrix,b_ub=rhs,bounds=(0,1),method='highs',
                          options={'time_limit':seconds})
            attempt = dict(coefficient_cutoff=cutoff,status=int(fit.status),message=fit.message)
            result['attempts'].append(attempt)
            if not fit.success:
                continue
            entries = [[i,str(Q(round(max(0.,-float(v))*(1<<48)),1<<48))]
                       for i,v in enumerate(fit.ineqlin.marginals) if v < 0]
            cap,correction = self.verify(support,entries)
            attempt.update(proposal=float(-fit.fun),verified_cap=cap,correction=str(correction))
            if cap < result['cap']:
                result.update(cap=cap,dual=entries)
        return result


def actual():
    from shared_support import shared_counts,authenticated_caps
    from count_refinements import refine_pair
    # Authenticate primal/dual distances, dimensions and production maps.
    _, ranks = shared_counts(True)
    primal,dual,_,_ = refine_pair(ranks,authenticated_caps())
    cdfs = [pairwise_support.pair_cdf(rows) for rows in (primal,dual)]
    shells = [pairwise_moments.shell_caps(256,128,d) for d in (30,38)]
    for row in shells:
        row[0]=1
    return Constraints(cdfs,shells,128)


def replay_shell_witnesses(model,record,caps):
    """Recheck saved multipliers; never import a saved numerical bound."""
    if (record.get('schema')!='pairwise-quaternary-shell-lp-1'
            or not isinstance(record.get('rows'),list) or len(caps)!=model.n+1):
        raise ValueError('matching quaternary shell witness record required')
    result=list(caps);checked=0
    for row in record['rows']:
        u=row.get('support')
        if type(u) is not int or not 1<=u<=model.n:
            raise ValueError('valid target shell required')
        if 'dual' not in row:continue
        upper,_=model.verify(u,row['dual'])
        result[u]=min(result[u],upper);checked+=1
    if not checked:raise ValueError('no exact shell multipliers to replay')
    return result,checked


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--supports',type=int,nargs='+',default=[112,128,144,160,176,192])
    parser.add_argument('--seconds',type=float,default=3.)
    parser.add_argument('--output',type=Path)
    args = parser.parse_args()
    model = actual()
    print('QUATERNARY model',len(model.keys),'variables',len(model.rows),'exact inequalities',flush=True)
    results = []
    for u in args.supports:
        row = model.solve(u,args.seconds);results.append(row)
        row['transform_cap']=model.transform_cap(u)
        gain = log2(row['prior'])-log2(row['cap']) if row['cap'] else None
        print('QUATERNARY support',u,'prior log2',log2(row['prior']) if row['prior'] else None,
              'gain bits',gain,'statuses',[v['status'] for v in row['attempts']],flush=True)
        print('  exact CDF transform gain',log2(row['prior'])-log2(row['transform_cap'])
              if row['transform_cap'] else None,flush=True)
        if args.output:
            args.output.parent.mkdir(parents=True,exist_ok=True)
            args.output.write_text(json.dumps(dict(schema='pairwise-quaternary-shell-lp-1',
                seconds=args.seconds,rows=results,note='Exact checked shell caps only; not a distance certificate.'),indent=2)+'\n')


if __name__ == '__main__':
    main()
