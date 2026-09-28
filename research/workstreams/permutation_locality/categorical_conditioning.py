"""Exact multitype conditioning identity for independently shuffled groups.

Finite toy transfer matrices verify the normalization. No SPIN local
operators or outer count bounds are supplied here; no certificate follows.
"""
from itertools import product
from math import factorial,prod
from fractions import Fraction as Q
from flint import fmpq_mat


def identity():return fmpq_mat([[1,0],[0,1]])


def self_test():
    checks=0
    for groups,length,categories in ((1,4,3),(2,3,2),(2,4,3)):
        probabilities=[tuple(Q(c+1+i,sum(range(1+i,categories+1+i))) for c in range(categories))
                       for i in range(groups)]
        labels=list(product(range(categories),repeat=groups))
        # Deliberately noncommuting, nonnegative rational local operators.
        operators={xs:fmpq_mat([[1+sum(xs),1+xs[0]],[xs[-1],2]]) for xs in labels}
        average=fmpq_mat(2,2)
        for xs,t in operators.items():
            p=prod(probabilities[i][x] for i,x in enumerate(xs))
            average+=t*(p.numerator)/p.denominator
        by_profile={};ways={}
        for sequence in product(labels,repeat=length):
            profile=tuple(tuple(sum(xs[i]==c for xs in sequence) for c in range(categories))
                          for i in range(groups))
            t=identity()
            for xs in sequence:t=t*operators[xs]
            by_profile[profile]=by_profile.get(profile,fmpq_mat(2,2))+t
            ways[profile]=ways.get(profile,0)+1
        recovered=fmpq_mat(2,2);iid=average**length
        probability_sum=Q(0)
        for profile,total in by_profile.items():
            count=prod(factorial(length)//prod(factorial(v) for v in counts) for counts in profile)
            assert count==ways[profile]
            one_sequence=prod(probabilities[i][c]**n for i,row in enumerate(profile) for c,n in enumerate(row))
            probability=count*one_sequence
            conditional=total/count
            contribution=conditional*probability.numerator/probability.denominator
            assert all(contribution[i,j]<=iid[i,j] for i in range(2) for j in range(2))
            recovered+=contribution;probability_sum+=probability;checks+=1
        assert recovered==iid and probability_sum==1
    print('Categorical conditioning:',checks,'exact profile/matrix identities and upper bounds passed',flush=True)


if __name__=='__main__':self_test()
