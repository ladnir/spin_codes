"""Support covers using both cumulative and individual shell bounds.

For nonnegative shell counts x[u], callers supply C[u] >= sum_{v<=u} x[v]
and optionally S[u] >= x[u].  IntervalFolds bounds the weighted sum on one
support interval.  Floating arithmetic proposes witnesses; outward replay
uses only exact counts and Arb enclosures.

The optional prefix-rank bound combines both constraints.  For a set A of
indices in [lo, hi], sum_{u in A} x[u] is at most each of

    sum_{u in A} S[u],
    C[k] + sum_{u in A, u>k} S[u]  (lo <= k <= hi).

Sort the nonnegative weights in descending order and express their weighted
sum as a positive combination of these nested-set sums.  Applying the
displayed bounds gives a valid upper bound without an LP solver.  This
argument does not treat differences of CDF bounds as shell bounds.

The separate cover implementation below preserves the shared-route driver's
box geometry, multiplicities, retained parents, and outward replay.  It does
not modify that driver or change any of its module globals.
"""
from collections import Counter
from fractions import Fraction
from heapq import heappop, heappush
from math import comb, factorial, log
from pathlib import Path
import sys

import numpy as np
from flint import arb, arb_mat
from scipy.optimize import minimize, minimize_scalar
from scipy.special import gammaln, logsumexp

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from group_rank_one_verify import up
from occupancy_adaptive import volume, split
from occupancy_cdf_cover import RetainedCover
from occupancy_memory_verify import DENOMINATOR
from occupancy_screen import matrix_for_probabilities
from two_group_screen import log_power_moment


def _log_exact(value):
    value = Fraction(value)
    return log(value.numerator) - log(value.denominator) if value else -np.inf


def _arb_exact(value):
    value = Fraction(value)
    return arb(value.numerator) / value.denominator


def prefix_rank_bound(cdf, shells, weights):
    """Bound sum x[i]*weights[i] from prefix caps and shell caps.

    All three sequences refer to the same interval.  cdf retains its global
    prefix cap at the interval's left endpoint.  Inputs may be exact rational
    numbers, or exact nonnegative Arb points for weights.  With Arb weights,
    the result is an enclosure whose upper endpoint is a valid bound.
    """
    assert len(cdf) == len(shells) == len(weights) and len(cdf) > 0
    assert all(x >= 0 for x in cdf) and all(x >= 0 for x in shells)
    assert all(x >= 0 for x in weights)
    order = sorted(range(len(weights)), key=lambda i: weights[i], reverse=True)
    # Cut -1 has no prefix allowance; it uses only individual shell caps.
    costs = [0, *cdf]
    total = 0
    for j, index in enumerate(order):
        for cut in range(index + 1):
            costs[cut] += shells[index]
        next_weight = weights[order[j + 1]] if j + 1 < len(order) else 0
        rank = min(costs)
        total += (weights[index] - next_weight) * (_arb_exact(rank) if isinstance(weights[index], arb) else rank)
    return total


class IntervalFolds:
    """Reusable interval objectives; n is the outer constituent length."""
    def __init__(self, cdf, shells=None, *, n=256, prefix_rank=False):
        self.n = n
        self.cdf = tuple(Fraction(x) for x in cdf)
        self.shells = None if shells is None else tuple(Fraction(x) for x in shells)
        assert len(self.cdf) == n + 1
        assert all(x >= 0 for x in self.cdf)
        assert all(a <= b for a, b in zip(self.cdf, self.cdf[1:]))
        if self.shells is not None:
            assert len(self.shells) == n + 1 and all(x >= 0 for x in self.shells)
        self.prefix_rank = prefix_rank and self.shells is not None
        self._functions = {}

    def function(self, lo, hi):
        """Return a binary64 proposal objective in log space."""
        assert 0 <= lo <= hi <= self.n
        key = lo, hi
        if key in self._functions:
            return self._functions[key]
        supports = np.arange(lo, hi + 1)
        choose = gammaln(self.n + 1) - gammaln(supports + 1) - gammaln(self.n - supports + 1)
        counts = [self.cdf[lo], *(self.cdf[u] - self.cdf[u - 1] for u in range(lo + 1, hi + 1))]
        cdf_logs = np.array([_log_exact(x) for x in counts])
        shell_logs = None if self.shells is None else np.array([_log_exact(x) for x in self.shells[lo:hi + 1]])
        if self.prefix_rank:
            cap_logs = np.array([_log_exact(x) for x in self.cdf[lo:hi + 1]])
            scale = max(float(np.max(cap_logs)), float(np.max(shell_logs)))
            if np.isfinite(scale):
                scaled_cdf = np.exp(cap_logs - scale)
                scaled_shells = np.exp(shell_logs - scale)
            else:
                scaled_cdf = scaled_shells = np.zeros(len(supports))

        def evaluate(p):
            assert 0 < p < 1
            weights = -choose - supports * np.log(p) - (self.n - supports) * np.log1p(-p)
            majorant = np.maximum.accumulate(weights[::-1])[::-1]
            result = float(logsumexp(cdf_logs + majorant))
            if shell_logs is not None:
                result = min(result, float(logsumexp(shell_logs + weights)))
            if self.prefix_rank and np.isfinite(scale):
                order = np.argsort(-weights)
                normalized = np.exp(weights[order] - weights[order[0]])
                differences = normalized - np.r_[normalized[1:], 0.]
                costs = np.r_[0., scaled_cdf].copy()
                ranked_sum = 0.
                for j, index in enumerate(order):
                    costs[:index + 1] += scaled_shells[index]
                    ranked_sum += differences[j] * np.min(costs)
                rank_log = -np.inf if ranked_sum == 0 else log(ranked_sum) + weights[order[0]] + scale
                result = min(result, float(rank_log))
            return result

        self._functions[key] = evaluate
        return evaluate

    def log(self, lo, hi, p):
        return self.function(lo, hi)(p)

    def outward(self, lo, hi, p):
        """Return an exact Arb point above the weighted interval sum."""
        assert 0 <= lo <= hi <= self.n and 0 < p < 1
        weights = []
        for u in range(lo, hi + 1):
            mass = arb(comb(self.n, u)) * p**u * (1 - p)**(self.n - u)
            lower = arb(mass.lower())
            assert lower > 0
            weights.append(up(1 / lower))
        majorant = weights.copy()
        for j in range(len(majorant) - 2, -1, -1):
            majorant[j] = max(majorant[j], majorant[j + 1])
        counts = [self.cdf[lo], *(self.cdf[u] - self.cdf[u - 1] for u in range(lo + 1, hi + 1))]
        result = up(sum((_arb_exact(c) * w for c, w in zip(counts, majorant)), arb(0)))
        if self.shells is not None:
            shells = self.shells[lo:hi + 1]
            shell_bound = up(sum((_arb_exact(c) * w for c, w in zip(shells, weights)), arb(0)))
            result = min(result, shell_bound)
            if self.prefix_rank:
                # Counts stay exact when constructing all nested-set caps.
                rank_bound = up(prefix_rank_bound(self.cdf[lo:hi + 1], shells, weights))
                result = min(result, rank_bound)
        return result


def cover(args, operators, count_sets, terminal, shell_sets=None, *, prefix_rank=False):
    """Complete support cover for args.groups; optional exact shell caps.

    The arguments match occupancy_cdf_cover.cover.  count_sets and shell_sets
    use the same penalty-string keys.  prefix_rank enables the combined
    nested-prefix bound; otherwise each fold takes min(CDF, shell).
    """
    size = len(terminal)
    assert all(len(exact) > args.groups and len(region) > args.groups
               for exact, region in operators.values())
    folds = {penalty: IntervalFolds(counts, None if shell_sets is None else shell_sets[penalty],
                                    prefix_rank=prefix_rank)
             for penalty, counts in count_sets.items()}
    cache = {}

    def proposal(choice, interval):
        key = choice, interval
        if key not in cache:
            lo, hi = interval
            region = operators[choice][1]
            fold = folds[choice[1]].function(lo, hi)
            def objective(z):
                p = 1 / (1 + np.exp(-z))
                return log_power_moment(matrix_for_probabilities(region, [p] * args.groups), 256, terminal) + args.groups * fold(p)
            fit = minimize_scalar(objective, bounds=(-8., 16.), method='bounded', options={'xatol': 1e-7})
            p = 1 / (1 + np.exp(-fit.x))
            numerator = max(1, min(DENOMINATOR - 1, round(p * DENOMINATOR)))
            cache[key] = numerator, fold(numerator / DENOMINATOR)
        return cache[key]

    def witness(box):
        best = None
        candidates = []
        for choice, (_, region) in operators.items():
            tilt, penalty = choice
            items = [proposal(choice, interval) for interval in box]
            nums = [n for n, _ in items]
            value = log_power_moment(matrix_for_probabilities(region, [n / DENOMINATOR for n in nums]), 256, terminal)
            value += float(tilt) * 209715 + sum(w for _, w in items)
            if best is None or value < best[0]:
                best = value, choice, nums
            candidates.append((value, choice, nums))
        if args.joint_witness:
            intervals = sorted(set(box))
            indices = [intervals.index(interval) for interval in box]
            for _, choice, nums in sorted(candidates)[:args.joint_top]:
                region = operators[choice][1]
                functions = [folds[choice[1]].function(lo, hi) for lo, hi in intervals]
                start = np.array([nums[box.index(interval)] / DENOMINATOR for interval in intervals])
                def objective(logits):
                    ps = 1 / (1 + np.exp(-logits))
                    value = log_power_moment(matrix_for_probabilities(region, [ps[i] for i in indices]), 256, terminal)
                    terms = [f(p) for f, p in zip(functions, ps)]
                    return value + float(choice[0]) * 209715 + sum(terms[i] for i in indices)
                fit = minimize(objective, np.log(start / (1 - start)), method='L-BFGS-B', bounds=[(-8., 16.)] * len(start),
                               options={'maxiter': 60, 'ftol': 1e-11})
                ps = 1 / (1 + np.exp(-fit.x))
                numbers = [max(1, min(DENOMINATOR - 1, round(p * DENOMINATOR))) for p in ps]
                quantized = np.array(numbers) / DENOMINATOR
                value = objective(np.log(quantized / (1 - quantized)))
                if value < best[0]:
                    best = value, choice, [numbers[i] for i in indices]
        return best

    if args.probe_supports or args.probe_vector:
        vectors = [(u,) * args.groups for u in args.probe_supports] + [tuple(v) for v in args.probe_vector]
        for vector in vectors:
            assert len(vector) == args.groups and all(38 <= u <= 256 for u in vector)
            value, choice, nums = witness(tuple((u, u) for u in sorted(vector)))
            multiplicity = factorial(args.groups)
            for n in Counter(vector).values():
                multiplicity //= factorial(n)
            score = (value + log(comb(2048, args.groups) * multiplicity)) / log(2)
            print('POINT support vector', vector, 'log2 upper', score, 'choice', choice, 'numerators', nums, flush=True)
        print('Selected-point binary64 diagnostic; no complete coverage or outward certificate.')
        return
    heap = []
    serial = 0
    tree = RetainedCover() if args.retain_parents else None
    locations = comb(2048, args.groups)
    def push(box, mult, parent=None):
        nonlocal serial
        score, choice, ps = witness(box)
        serial += 1
        item = (-score - log(locations * mult), serial, box, mult, choice, ps)
        heappush(heap, item)
        if tree is not None:
            tree.add(item, parent)
    push(((38, 256),) * args.groups, 1)
    threshold = -(args.target_bits + 2) * log(2)
    score = -heap[0][0]
    steps = 0
    while heap and steps < args.max_splits and score > threshold:
        item = heappop(heap)
        _, _, box, mult, _, _ = item
        if all(lo == hi for lo, hi in box):
            if tree is not None:
                continue
            heappush(heap, item)
            print('Dominating singleton', box, flush=True)
            break
        for child, child_mult in split(box, mult):
            push(child, child_mult, item[1] if tree is not None else None)
        if tree is not None:
            tree.update(item[1])
        steps += 1
        if steps % 25 == 0:
            score = tree.nodes[1]['best'] if tree is not None else float(logsumexp([-x[0] for x in heap]))
            if steps % 250 == 0:
                print('Shell/CDF-cover splits', steps, 'leaves', len(heap), 'log2 union', score / log(2), flush=True)
    if tree is not None:
        heap = tree.selected()
    score = float(logsumexp([-x[0] for x in heap]))
    assert sum(volume(x[2], x[3]) for x in heap) == 219**args.groups
    print('Full-domain shell/CDF cover:', steps, 'splits;', len(heap), 'leaves; log2 union', score / log(2), flush=True)
    print('Largest leaves:', [(-x[0] / log(2), x[2], x[4]) for x in sorted(heap)[:3]], flush=True)
    if score > threshold or args.screen_only:
        print('No outward certificate from this run.', flush=True)
        return
    total = arb(0)
    for index, (_, _, box, mult, choice, nums) in enumerate(heap, 1):
        tilt, penalty = choice
        ps = [arb(n) / DENOMINATOR for n in nums]
        masses = [arb(1)]
        for p in ps:
            updated = [arb(0)] * (len(masses) + 1)
            for j, m in enumerate(masses):
                updated[j] += m * (1 - p)
                updated[j + 1] += m * p
            masses = updated
        matrix = sum((m * r for m, r in zip(masses, operators[choice][0])), arb_mat(size, size))**256
        term = sum((matrix[0, j] for j in range(size) if terminal[j]), arb(0)) * (arb(tilt) * 209715).exp() * locations * mult
        for (lo, hi), p in zip(box, ps):
            term = up(term * folds[penalty].outward(lo, hi, p))
        total = up(total + term)
        if index % 1000 == 0:
            print('Outward shell/CDF-cover leaves', index, '/', len(heap), flush=True)
    assert 0 < total < arb(2)**-args.target_bits
    print('VERIFIED all supports at occupancy', args.groups, 'precision', args.precision, 'upper', total,
          'margin', -total.log() / arb(2).log(), flush=True)
    print('Other occupancies remain; not a full-code certificate.')
    return total
