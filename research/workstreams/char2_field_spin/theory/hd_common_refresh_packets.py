"""Exact common-refresh operators for direct sums of RS[4h,h] blocks.

The local construction is specified in HD_COMMON_REFRESH.md. Each entry
is a rational majorant before explicit upward dyadic rounding. This module
does not claim whole-code coverage or include the diagonal density loss.
"""
from fractions import Fraction as F
from functools import lru_cache
from math import comb, prod

from hd_bucket_feedback import profiles


@lru_cache(None)
def principal_defects(h=4):
    from hd_common_refresh_principal import nullities
    values = nullities(h)
    return tuple(max(int(values[mask]) for mask in range(1 << (4*h))
                     if mask.bit_count() == c) for c in range(4*h+1))


def local_return_cap(l, z, q=65536, h=4):
    z = F(z)
    if h not in (1,2,4):raise ValueError('Unsupported local rank')
    length=4*h
    nu = principal_defects(h)
    rho = min(h, l)
    result = F(0)
    for c in range(l+1):
        for e in range(min(h, length+1-l)):
            rank = max(0, c-min(nu[c], max(0, h-l+c)))
            count = F(q**min(rho, h-e)-1, q**(l-max(0, l-c-h)))
            probability = min(1-F(1, q**rho), F(1, q**rank), count)
            result += comb(l, c)*comb(length-l, e)*z**(length-c-e)*(1-z)**(c+e)*probability
    return result


def local_moments(z, q=65536, h=4):
    z = F(z)
    length=4*h
    quiet = tuple(sum(F(comb(l, c), q**(c+min(h, l-c)))*z**(l-c)*(1-z)**c
                      for c in range(l+1)) for l in range(length+1))
    flat = tuple(sum(F(comb(n, c), q**min(h, c))*z**(n-c)*(1-z)**c
                     for c in range(n+1)) for n in range(length+1))
    ret = tuple(min((z+(1-z)/q)**l-quiet[l],
                    (1-F(1, q**min(h, l)))*z**max(0, length-h+1-l),
                    local_return_cap(l, z, q, h)) for l in range(length+1))
    return quiet, flat, ret


def matrix(occupancies, z, local=None, q=65536, h=4):
    if h not in (1,2,4) or len(occupancies) != 16//h or any(not 0 <= l <= 4*h for l in occupancies):
        raise ValueError('Invalid common-refresh occupancies')
    z = F(z)
    quiet, flat, ret = local_moments(z, q, h) if local is None else local
    d = q**16-1
    k = prod(quiet[l] for l in occupancies)
    t = (z+(1-z)/q)**sum(occupancies)
    v = t*(q**16*prod(flat[4*h-l] for l in occupancies)-1)/d
    j = min(prod(quiet[l]+ret[l] for l in occupancies)-k, d*v)
    return k, t-k, j/d, v


def grouped_packets(z, bits=12, maximum=64, h=4):
    from hd_sumrank_paired_certificate import round_up
    if not 0 <= maximum <= 64:
        raise ValueError('Invalid packet degree cutoff')
    local = local_moments(z,h=h)
    groups = {}
    for occupancies, count in profiles(16, 64, h):
        j = sum(occupancies)
        if j > maximum:
            continue
        entries = tuple(round_up(x, bits) for x in matrix(occupancies, z, local,h=h))
        key = (j, entries)
        groups[key] = groups.get(key, 0)+count
    result = tuple((j, count, entries) for (j, entries), count in groups.items())
    for j in range(maximum+1):
        assert sum(count for k, count, _ in result if k == j) == comb(64, j)
    return result
