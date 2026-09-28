"""Binary64 scalar Holder envelopes for independent-row category families.

No inner construction or distance bound is evaluated. The printed costs
measure only the positive product envelope, not certified message counts.
"""
import argparse
from collections import Counter
from itertools import combinations_with_replacement
from math import comb, factorial, log, prod
from pathlib import Path
import sys

import numpy as np
from scipy.optimize import minimize_scalar
from scipy.special import logsumexp

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from bch_joint_support import authenticated_caps
from two_group_screen import log_binomial_mass

NAMES = ('0', 'L', 'M', 'H', '1')
INTERVALS = ((0,0), (38,62), (64,192), (194,218), (256,256))
SCORES = (-2,-1,0,1,2)


def category_gamma(caps, lo, hi, p):
    if p == 0:
        assert lo == hi == 0
        return log(caps[0])
    if p == 1:
        assert lo == hi == 256
        return log(caps[256])
    return max(log(caps[w])-log_binomial_mass(256,w,p)
               for w in range(lo,hi+1) if caps[w])


def types(caps, p):
    probabilities = (0,p,.5,1-p,1)
    gammas = [category_gamma(caps,lo,hi,x) for (lo,hi),x in zip(INTERVALS,probabilities)]
    result = []
    for labels in combinations_with_replacement(range(5),4):
        if labels == (0,0,0,0):
            continue
        histogram = Counter(labels)
        multiplicity = factorial(4)//prod(factorial(n) for n in histogram.values())
        logh = log(multiplicity)+sum(gammas[label] for label in labels)
        law = np.array([1.])
        for label in labels:
            x = probabilities[label]
            law = np.convolve(law,[1-x,x])
        assert abs(law.sum()-1) < 1e-12 and np.min(law) >= 0
        result.append((labels,logh,law))
    assert len(result) == 69
    return gammas,result


def families(rows):
    exception = [row for row in rows if row[0] != (2,2,2,2)]
    side = lambda row: sum(SCORES[label] for label in row[0])
    result = {'middle':[row for row in rows if row[0] == (2,2,2,2)],
              'all_nonzero':rows, 'all_exceptions':exception,
              'exceptions_low_or_tied':[row for row in exception if side(row) <= 0],
              'exceptions_high':[row for row in exception if side(row) > 0],
              'exceptions_low_strict':[row for row in exception if side(row) < 0],
              'exceptions_balanced':[row for row in exception if side(row) == 0]}
    assert len(result['exceptions_low_or_tied'])+len(result['exceptions_high']) == 68
    for count in range(1,5):
        selected = [row for row in exception if sum(label != 2 for label in row[0]) == count]
        result[f'exactly_{count}_exceptions'] = selected
        result[f'exactly_{count}_low_or_tied'] = [row for row in selected if side(row) <= 0]
        result[f'exactly_{count}_high'] = [row for row in selected if side(row) > 0]
        result[f'exactly_{count}_low_strict'] = [row for row in selected if side(row) < 0]
        balanced = [row for row in selected if side(row) == 0]
        if balanced:
            result[f'exactly_{count}_balanced'] = balanced
    assert sum(len(result[f'exactly_{count}_exceptions']) for count in range(1,5)) == 68
    return result


def envelope(rows):
    logh = np.array([row[1] for row in rows])
    law = np.array([row[2] for row in rows])
    loglaw = np.full(law.shape,-np.inf)
    np.log(law,out=loglaw,where=law > 0)
    lognu = logsumexp(logh[:,None]+256*loglaw,axis=0)/256
    logmass = logsumexp(lognu)
    normalized = np.exp(lognu-logmass)
    return 256*logmass/log(2), logsumexp(logh)/log(2), normalized, np.exp(lognu)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--p', type=float, default=.25)
    parser.add_argument('--optimize', action='store_true')
    args = parser.parse_args()
    if not 0 < args.p < .5:
        parser.error('0 < p < .5 required')
    caps = authenticated_caps()
    covered = [w for lo,hi in INTERVALS for w in range(lo,hi+1)]
    assert len(set(covered)) == len(covered)
    assert all(w in covered for w,a in enumerate(caps) if a)
    gammas,rows = types(caps,args.p)
    family = families(rows)
    middle = envelope(family['middle'])[0]
    print('BINARY64 category witness p/log2gammas:',args.p,[g/log(2) for g in gammas],flush=True)
    print('BINARY64 middle product cost:',middle,flush=True)
    for name,selected in family.items():
        cost,count,law,nu = envelope(selected)
        print('FAMILY',name,'types',len(selected),'cost',cost,'log2mass',count,
              'Holder inflation',cost-count,'gap_below_middle',middle-cost,
              'normalized_packet',list(law),'nu',list(nu),flush=True)
    if args.optimize:
        targets = ['all_exceptions','exceptions_low_or_tied','exceptions_high']
        targets += [f'exactly_{n}_{side}' for n in range(1,5) for side in ('low_or_tied','high')]
        targets += [name for name in family if name.startswith('exactly_') and
                    (name.endswith('low_strict') or name.endswith('balanced'))]
        for name in targets:
            def objective(p):
                return envelope(families(types(caps,p)[1])[name])[0]
            fit = minimize_scalar(objective,bounds=(.12,.46),method='bounded',options={'xatol':1e-7})
            optimum = min((fit.fun,fit.x),(objective(args.p),args.p))
            _,selected = types(caps,optimum[1])
            cost,count,law,nu = envelope(families(selected)[name])
            print('OPTIMIZED FAMILY',name,'p',optimum[1],'cost',cost,'log2mass',count,
                  'Holder inflation',cost-count,'gap_below_middle',middle-cost,
                  'normalized_packet',list(law),'nu',list(nu),flush=True)
    print('Scalar binary64 envelopes only: no inner moment, certificate, or occupancy cover.',flush=True)


if __name__ == '__main__':
    main()
