"""Exact interval anti-concentration for the original adaptive quartet law.

For each of the15 old types0<=w<=u<=4, normalize the five label-count
coefficients of the hole and selected factors at a rational x>0. Let
 a=min_(all30laws) sum_(j=0)^3 p_j*p_(j+1).
Then a>0 and a<=1/4. A product of an even number m of these factors has
every coefficient probability at most
 eta=[z^0](a*z^-1+1-2*a+a*z)^(m/2).
Indeed each characteristic function has squared modulus at most
1-2*a+2*a*cos(t); Fourier inversion gives the displayed central coefficient.
The same bound holds for positive mixtures, including each physical
population transition after normalization at x. Any integer interval
of length ell therefore has probability at most min(1,ell*eta).

This factor applies to a current-update label-count interval only. It
must be applied before minimizing quota evaluations. It does not apply
automatically to arbitrary positive polynomial test matrices, and does
not replace sixteen-update quota conditioning. The eta and local laws
here are exact Fractions; callers using FLOAT arithmetic remain FLOAT
diagnostics. There is no numerical work on import. Parent owns execution.
"""

import argparse
from fractions import Fraction
from functools import lru_cache
from math import comb
from operator import index
import os
import shutil
import threading


def _qa_parameters(x, quartets):
    x, quartets = Fraction(x), index(quartets)
    if x <= 0 or quartets < 2 or quartets > 128 or quartets % 2:
        raise ValueError('Use rational x>0 and an even integer quartet count2..128')
    return x, quartets


@lru_cache(maxsize=32)
def _qa_local_laws(x):
    """Thirty exact tilted probability vectors, with old type/branch labels."""
    laws = []
    for u in range(5):
        denominator = (1 << (4+u))-1
        for w in range(u+1):
            hole, selected = [], []
            for j in range(5):
                t, total = Fraction(0), Fraction(0)
                for r in range(max(0, j-(4-w)), min(w, j)+1):
                    multiplicity = comb(w, r)*comb(4-w, j-r)*(1 << (w-r))
                    loss = max(0, 4+u-2*j+r)
                    g = Fraction((1 << loss)-1, denominator)
                    t += multiplicity*g
                    total += multiplicity
                hole.append(t*x**j)
                selected.append((total-t)*x**j)
            for name, coefficients in (('hole', hole), ('selected', selected)):
                mass = sum(coefficients)
                if mass <= 0 or any(value < 0 for value in coefficients):
                    raise ArithmeticError('Invalid original-network local polynomial')
                probabilities = tuple(value/mass for value in coefficients)
                adjacency = sum(probabilities[j]*probabilities[j+1] for j in range(4))
                if adjacency <= 0 or 4*adjacency > 1:
                    raise ArithmeticError('Missing adjacent support or invalid probability law')
                laws.append((u, w, name, probabilities, adjacency))
    if len(laws) != 30:
        raise AssertionError('Expected fifteen old types and two branches')
    return tuple(laws)


@lru_cache(maxsize=128)
def _qa_exact_data(x, quartets):
    laws = _qa_local_laws(x)
    adjacency = min(law[4] for law in laws)
    a, denominator = adjacency.numerator, adjacency.denominator
    stay = denominator-2*a
    half = quartets//2
    numerator = sum(comb(half, 2*k)*comb(2*k, k)*a**(2*k)*stay**(half-2*k)
                    for k in range(quartets//4+1))
    eta = Fraction(numerator, denominator**half)
    if not 0 < eta <= 1:
        raise ArithmeticError('Invalid exact anti-concentration constant')
    minimizing_types = tuple((u, w, name) for u, w, name, _, value in laws if value == adjacency)
    return adjacency, eta, minimizing_types


def quota_anticoncentration_eta(x, quartets=128):
    """Return the exact Fraction eta for one update with m=quartets."""
    x, quartets = _qa_parameters(x, quartets)
    return _qa_exact_data(x, quartets)[1]


def quota_anticoncentration_data(x, quartets=128):
    """Return fresh metadata; all numerical constants remain exact Fractions."""
    x, quartets = _qa_parameters(x, quartets)
    adjacency, eta, types = _qa_exact_data(x, quartets)
    return dict(x=x, quartets=quartets, width=4*quartets, a=adjacency, eta=eta,
                minimizing_types=types, local_laws=30, exact_fraction=True,
                physical_entry_mixtures=True, regional_quota_replaced=False)


def _qa_product(left, right):
    result = [Fraction(0)]*(len(left)+len(right)-1)
    for i, a in enumerate(left):
        for j, b in enumerate(right):
            result[i+j] += a*b
    return result


def _qa_factor_product(factors):
    result = [Fraction(1)]
    for factor in factors:
        result = _qa_product(result, factor)
    if sum(result) != 1 or any(value < 0 for value in result):
        raise AssertionError('Exact factor product is not a probability law')
    return result


def _qa_intervals(probabilities, eta, all_lengths):
    prefix = [Fraction(0)]
    for value in probabilities:
        prefix.append(prefix[-1]+value)
    lengths = range(1, len(probabilities)+1) if all_lengths else (1, 2, len(probabilities))
    checks = 0
    for length in lengths:
        bound = min(Fraction(1), length*eta)
        for lower in range(len(probabilities)-length+1):
            assert prefix[lower+length]-prefix[lower] <= bound
            checks += 1
    return checks


def quota_anticoncentration_checks():
    counts = dict(local_laws=0, exact_central_coefficients=0,
                  homogeneous_products=0, mixed_pairs=0, mixed_four_products=0,
                  positive_mixtures=0, interval_bounds=0, production_constants=0,
                  argument_guards=0)
    for x in (Fraction(1), Fraction(2), Fraction(4)):
        laws = _qa_local_laws(x)
        probabilities = [law[3] for law in laws]
        for u, w, _, p, adjacency in laws:
            assert sum(p) == 1 and adjacency == sum(p[j]*p[j+1] for j in range(4))
            even, odd = sum(p[::2]), sum(p[1::2])
            assert 0 < adjacency <= even*odd <= Fraction(1, 4)
            counts['local_laws'] += 1
        for quartets in (2, 4):
            data = quota_anticoncentration_data(x, quartets)
            a, eta = data['a'], data['eta']
            walk = _qa_factor_product([(a, 1-2*a, a)]*(quartets//2))
            assert walk[quartets//2] == eta
            counts['exact_central_coefficients'] += 1
            for p in probabilities:
                product = _qa_factor_product([p]*quartets)
                assert max(product) <= eta
                counts['homogeneous_products'] += 1
                counts['interval_bounds'] += _qa_intervals(product, eta, True)
        eta2, eta4 = quota_anticoncentration_eta(x, 2), quota_anticoncentration_eta(x, 4)
        for i, first in enumerate(probabilities):
            product = _qa_factor_product((first, probabilities[(i+1) % 30]))
            assert max(product) <= eta2
            counts['mixed_pairs'] += 1
            counts['interval_bounds'] += _qa_intervals(product, eta2, False)
            four = _qa_factor_product([probabilities[j % 30] for j in (i, i+1, i+7, i+17)])
            assert max(four) <= eta4
            counts['mixed_four_products'] += 1
            counts['interval_bounds'] += _qa_intervals(four, eta4, True)
            second_four = _qa_factor_product([probabilities[(i+9) % 30]]*4)
            mixture = [(a+2*b)/3 for a, b in zip(four, second_four)]
            assert sum(mixture) == 1 and max(mixture) <= eta4
            counts['positive_mixtures'] += 1
            counts['interval_bounds'] += _qa_intervals(mixture, eta4, False)
        production = quota_anticoncentration_data(x, 128)
        assert 0 < production['eta'] <= eta4 <= eta2 <= 1
        counts['production_constants'] += 1
    for x, quartets in ((0, 2), (-1, 2), (1, 0), (1, 1), (1, 3), (1, 130), (1, 2.0)):
        try:
            quota_anticoncentration_eta(x, quartets)
        except (ValueError, TypeError):
            counts['argument_guards'] += 1
        else:
            raise AssertionError('Invalid anti-concentration parameters were accepted')
    print(dict(**counts, all_arithmetic_exact=True, bounded_width8_and16_products=True,
               full_operator_or_distance_certificate=False), flush=True)


def quota_anticoncentration_main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--checks', required=True, action='store_true')
    parser.parse_args()
    import resource
    if shutil.disk_usage('/tmp').free < 20*(1 << 30):
        raise RuntimeError('20GiB /tmp free-space guard')
    resource.setrlimit(resource.RLIMIT_AS, (2*(1 << 30), 2*(1 << 30)))
    timer = threading.Timer(175, lambda: os._exit(124))
    timer.daemon = True
    timer.start()
    try:
        quota_anticoncentration_checks()
    finally:
        timer.cancel()


if __name__ == '__main__':
    quota_anticoncentration_main()
