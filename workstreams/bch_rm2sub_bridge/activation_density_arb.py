"""Outward four-class transfer with exact character bounds on activation.

Classes Z,D,U,L and the inequalities are those of the landscape's
occupation_refresh_v1 and syndrome_density_v1. No binary64 moment is used.
"""
import math

from flint import arb, arb_poly

import certificate_search_core as core
import syndrome_density_v1 as reference


def mul(a, b):
    return tuple(sum((a[4*i+k]*b[4*k+j] for k in range(4)), arb(0))
                 for i in range(4) for j in range(4))


def power(a, count):
    if count < 0:
        raise ValueError('Negative matrix exponent')
    result = tuple(arb(int(i == j)) for i in range(4) for j in range(4))
    while count:
        if count & 1:
            result = mul(result, a)
        count >>= 1
        if count:
            a = mul(a, a)
    return result


def epoch(t, s, spectrum, kernel, z, maximum=None):
    """Point upper bounds on every entry; exact integers select the density rule."""
    if maximum is None:
        maximum = t
    core.require(0 < z < 1 and s >= 2 and 0 <= maximum <= t, 'Invalid epoch parameters')
    field = 1 << s
    M = field-1
    core.require(sum(spectrum.values()) == M and sum(kernel.values()) == 1 << (t-s), 'Wrong spectrum mass')
    chars = {w: reference.krawtchouk_row(t, w) for w in spectrum}
    d = min(spectrum)
    powers = [z**i for i in range(t+1)]
    output = []
    for j in range(maximum+1):
        choose = math.comb(t, j)
        K = kernel.get(j, 0)
        signed = choose+sum(n*chars[w][j] for w, n in spectrum.items())
        core.require(signed == field*K, 'Kernel/character mismatch')
        row = [arb(0) for _ in range(16)]
        if j == 0:
            arbitrary = powers[d].upper()
            uniform = (sum((n*powers[w] for w, n in spectrum.items()), arb(0))/M).upper()
            row[0] = arb(1)
            row[6] = arbitrary
            row[10] = uniform
            row[14] = min(arbitrary, (arb(M)/(M-1)*uniform).upper())
        else:
            moments = {}
            for w in spectrum:
                moments[w] = (sum((math.comb(w, v)*math.comb(t-w, j-v)*powers[w+j-2*v]
                    for v in range(max(0, j-t+w), min(w, j)+1)), arb(0))/choose).upper()
            arbitrary = max(moments.values())
            uniform = (sum((spectrum[w]*moments[w] for w in spectrum), arb(0))/M).upper()
            live = (arbitrary, uniform, min(arbitrary, (arb(M)/(M-1)*uniform).upper()))
            beta = arb(K)/choose
            nonkernel = arb(choose-K)/choose
            row[0] = (beta*powers[j]).upper()
            row[1] = (nonkernel*powers[j]).upper()
            numerator = min(choose+sum(n*abs(chars[w][j]) for w, n in spectrum.items()), field*(choose-K))
            mass_numerator = max(field*(choose-K), (field-2)*numerator)
            if choose > K and mass_numerator <= 4*field*(choose-K):
                row[1] = arb(0)
                row[3] = (arb(mass_numerator)/(field*choose)*powers[j]).upper()
            distance = powers[max(0, d-j)]
            for state, moment in enumerate(live, 1):
                kernel_upper = min(moment, (beta*distance).upper())
                nonkernel_upper = min(moment, (nonkernel*distance).upper())
                row[4*state] = (nonkernel_upper/M).upper()
                row[4*state+3] = (arb(M-1)/M*moment+kernel_upper/M).upper()
        output.append(tuple(row))
    return output


def regions(t, s, spectrum, kernel, z, maximum, length):
    core.require(0 <= maximum <= length and length >= t and length % t == 0, 'Invalid region geometry')
    coefficients = epoch(t, s, spectrum, kernel, z, min(maximum, t))
    a = tuple(arb_poly([row[k]*math.comb(t, j) for j, row in enumerate(coefficients)]) for k in range(16))
    result = tuple(arb_poly([int(i == j)]) for i in range(4) for j in range(4))
    def product(left, right):
        return tuple(sum((left[4*i+k]*right[4*k+j] for k in range(4)), arb_poly()).truncate(maximum+1)
                     for i in range(4) for j in range(4))
    count = length//t
    while count:
        if count & 1:
            result = product(result, a)
        count >>= 1
        if count:
            a = product(a, a)
    return [tuple(max(arb(0), (p[j]/math.comb(length, j)).upper()) for p in result)
            for j in range(maximum+1)]
