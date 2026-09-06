"""Fixed-row-weight bound from positive exponential region modes.

Given R_j <= sum_i P_i (a+i*d)^j entrywise, Maclaurin's elementary
symmetric mean inequality bounds each permuted row by its enumerator
at a+d*b/n, where b is the sum of the mode labels across n regions.
Matrix polynomial coefficients retain the noncommutative region order.
"""
from fractions import Fraction as F
import activation_bridge as q1


def polynomial_power(modes,positions,number=F):
    current=[tuple(number(int(i==j)) for i in range(3) for j in range(3))]
    for _ in range(positions):
        updated=[[number(0)]*9 for _ in range(len(current)+len(modes)-1)]
        for degree,left in enumerate(current):
            for label,right in enumerate(modes):
                product=q1.positive_mul(left,right)
                updated[degree+label]=[x+y for x,y in zip(updated[degree+label],product)]
        current=[tuple(row) for row in updated]
    return current


def numerator(modes,a,delta,positions,occupancy,caps):
    coefficients=polynomial_power(modes,positions)
    return sum((sum(matrix[:3])*sum((count*(a+delta*F(b,positions))**w for w,count in caps.items()),F(0))**occupancy
                for b,matrix in enumerate(coefficients)),F(0))
