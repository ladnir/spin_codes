"""Rank-potential bound retaining the joint rank of required-live forms.

Fix U<=F_q^m of dimension d and a numerical matrix B of rank r. Sample
the m diagonal entries independently from F_q^*. This bounds
E[tau^(-rank(B D |_U)); ker(B D |_U) contains a full-support vector].
All equations in B must be charged jointly before these diagonals are
exposed. It is not a bound after conditioning on some of those equations.
"""
from fractions import Fraction as F


def live_kernel_bound(dimension, coordinates, equation_rank, tau, q=256,
                      sharpen=True):
    d, m, r = dimension, coordinates, equation_rank
    if not (1 <= d <= m < q and 0 <= r <= m and 1 <= tau <= q-1):
        raise ValueError('require 1<=d<=m<q, 0<=r<=m, 1<=tau<=q-1')
    if isinstance(tau, int):
        tau = F(tau)
    maximum_kernel = min(d, m-r)
    if maximum_kernel == 0:
        return tau*0
    expected_live = F(q**d-1, (q-1)**r)
    if not sharpen:
        return min(tau**0, expected_live*tau**(1-d)/(q-m))
    # All coordinate kernels contain zero. If a full-support vector exists,
    # each coordinate kernel is a proper hyperplane of the j-dimensional
    # kernel. Counting their shared zero improves the elementary union bound.
    factor = max(tau**(j-d)/((q**j-1)-m*(q**(j-1)-1))
                 for j in range(1, maximum_kernel+1))
    return min(tau**0, expected_live*factor)
