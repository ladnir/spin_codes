"""Retain even row weights after GF randomization, before regional routing.

Conditional on the packet activity support, independent GF multipliers
erase the lane values. The remaining row-parity factor has 16 characters.
"""
from fractions import Fraction as Q
from math import prod
import scalar_cover as sc


def parity_factor(probabilities,minimum_support):
    """Uniform upper on Pr[all row parities even | a support of size >=d]."""
    probabilities=list(map(Q,probabilities))
    if (len(probabilities)!=4 or any(not 0<=p<=1 for p in probabilities)
            or type(minimum_support) is not int or minimum_support<1):
        raise ValueError('four probabilities and positive integer support floor required')
    empty=prod(1-p for p in probabilities);activity=1-empty
    if not activity:return Q(1)
    characters=[(prod(1-2*p for i,p in enumerate(probabilities) if mask>>i&1)-empty)/activity
                for mask in range(16)]
    if any(abs(v)>1 for v in characters):raise ArithmeticError('invalid conditional character')
    return sum((abs(v)**minimum_support for v in characters),Q(0))/16


def components(rows,minimum_support):
    """Positive iid activity comparisons for nonzero groups with support >=d."""
    if len(rows)>10:raise ValueError('single-digit row component names required')
    result=[]
    for name,coefficient,law,active in sc.group_components(rows):
        factor=parity_factor([rows[int(i)][1] for i in name],minimum_support) if active else Q(1)
        result.append((name,coefficient*factor,law,active))
    return result


def actual_components(central_bits=129,theta=Q(21,50),central_scale=Q(1001,1000)):
    caps=sc.authenticated_caps()
    if any(caps[w] for w in range(1,len(caps),2)):
        raise ValueError('authenticated caps must exclude every odd row weight')
    support=min(w for w in range(1,len(caps)) if caps[w])
    rows=sc.envelope(caps,Q(central_scale)*(1<<central_bits),Q(theta))
    return components(rows,support)
