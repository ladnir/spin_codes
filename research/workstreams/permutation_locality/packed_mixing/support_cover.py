"""Sparse support-box covers for rational expected CDFs, with explicit domain.

This is a scoped adaptation of occupancy_cdf_cover. Legacy BCH-only covers
retain their original domain. Every accepted result here is freshly evaluated
with outward arithmetic, including supports below the BCH distance.
"""
from collections import Counter
from fractions import Fraction as Q
from heapq import heappop, heappush
from math import comb, log

import numpy as np
from scipy.optimize import minimize, minimize_scalar
from scipy.special import logsumexp
from flint import arb, arb_mat, ctx

import occupancy_cdf_cover as old
from gf16_packets.cdf_gradient import JointObjective
from occupancy_rank import aq, up


def domain(counts, groups, support_min):
    if (len(counts) != 257 or counts[0] != 0 or counts[-1] <= 0
            or any(c < 0 for c in counts) or any(a > b for a, b in zip(counts, counts[1:]))
            or type(groups) is not int or not 1 <= groups <= 2048
            or type(support_min) is not int or not 1 <= support_min <= 256
            or any(counts[:support_min]) or not counts[support_min]):
        raise ValueError('full rational CDF and exact first-positive support floor required')
    return ((support_min, 256),)*groups


def fold_arb(counts, lo, hi, p):
    """Directed interval CDF fold, retaining fractional expected counts."""
    if not 0 < p < 1 or not 0 <= lo <= hi <= 256:
        raise ValueError('interior probability and valid support interval required')
    weights = []
    for u in range(lo, hi+1):
        mass = arb(comb(256, u))*p**u*(1-p)**(256-u)
        lower = arb(mass.lower())
        if not lower > 0:
            raise ArithmeticError('strictly positive binomial denominator required')
        weights.append(up(1/lower))
    for i in range(len(weights)-2, -1, -1):
        weights[i] = max(weights[i], weights[i+1])
    return up(sum((aq(Q(c))*w for c, w in zip(old.increments(counts, lo, hi), weights)), arb(0)))


def checked_partition(tree, selected, root_box):
    """Authenticate the split tree, not merely its aggregate volume."""
    chosen = {item[1] for item in selected}
    if len(chosen) != len(selected) or tree.nodes[1]['item'][2:4] != (root_box, 1):
        raise ValueError('unique leaves and exact root domain required')
    reached = set()
    visited = set()
    def visit(identifier):
        if identifier in visited or identifier not in tree.nodes:
            raise ValueError('split tree contains a cycle, duplicate edge, or absent node')
        visited.add(identifier)
        node = tree.nodes[identifier]
        item = node['item']
        if identifier in chosen:
            reached.add(identifier)
            return
        expected = list(old.split(item[2], item[3]))
        actual = [tree.nodes[c]['item'][2:4] for c in node['children']]
        if not expected or actual != expected:
            raise ValueError('selected leaves do not exhaust a checked binary split')
        for c in node['children']:
            if tree.nodes[c]['parent'] != identifier:
                raise ValueError('invalid split ancestry')
            visit(c)
    visit(1)
    if reached != chosen or sum(old.volume(x[2], x[3]) for x in selected) != old.volume(root_box, 1):
        raise ValueError('support cover is not exhaustive')
    return True


def cover(args, operators, counts, terminal, *, cutoff, support_min):
    root_box = domain(counts, args.groups, support_min)
    if type(cutoff) is not int or not 0 <= cutoff < (1 << 21):
        raise ValueError('integer output cutoff in the code length required')
    if type(args.target_bits) is not int or args.target_bits < 1:
        raise ValueError('positive integer target bits required')
    size = len(terminal)
    if any(len(exact) <= args.groups or len(region) <= args.groups for exact, region in operators.values()):
        raise ValueError('operators must include the requested occupancy')
    cache = {}
    def proposal(choice, interval):
        key = choice, interval
        if key not in cache:
            region = operators[choice][1]
            fold = old.fold_function(counts, *interval)
            def objective(z):
                p = 1/(1+np.exp(-z))
                return old.log_power_moment(old.matrix_for_probabilities(region, [p]*args.groups), 256, terminal)+args.groups*fold(p)
            fit = minimize_scalar(objective, bounds=(-8., 16.), method='bounded', options={'xatol': 1e-7})
            p = 1/(1+np.exp(-fit.x))
            numerator = max(1, min(old.DENOMINATOR-1, round(p*old.DENOMINATOR)))
            cache[key] = numerator, fold(numerator/old.DENOMINATOR)
        return cache[key]
    def witness(box):
        candidates = []
        for choice, (_, region) in operators.items():
            if choice[1] != '1':
                raise ValueError('only unweighted expected counts supported')
            items = [proposal(choice, interval) for interval in box]
            nums = [n for n, _ in items]
            value = old.log_power_moment(old.matrix_for_probabilities(region, [n/old.DENOMINATOR for n in nums]), 256, terminal)
            value += float(choice[0])*cutoff+sum(w for _, w in items)
            candidates.append((value, choice, nums))
        best = min(candidates)
        if args.joint_witness:
            intervals = sorted(set(box))
            indices = [intervals.index(interval) for interval in box]
            for _, choice, nums in sorted(candidates)[:args.joint_top]:
                region = operators[choice][1]
                start = np.array([nums[box.index(interval)]/old.DENOMINATOR for interval in intervals])
                objective = JointObjective(region, box, counts, terminal, float(choice[0])*cutoff)
                fit = minimize(objective, np.log(start/(1-start)), jac=True, method='L-BFGS-B',
                    bounds=[(-8., 16.)]*len(start), options={'maxiter': 60, 'ftol': 1e-11})
                ps = 1/(1+np.exp(-fit.x))
                numbers = [max(1, min(old.DENOMINATOR-1, round(p*old.DENOMINATOR))) for p in ps]
                quantized = np.array(numbers)/old.DENOMINATOR
                value = objective(np.log(quantized/(1-quantized)))[0]
                if value < best[0]:
                    best = value, choice, [numbers[i] for i in indices]
        return best
    heap = []
    tree = old.RetainedCover()
    serial = 0
    locations = comb(2048, args.groups)
    def push(box, mult, parent=None):
        nonlocal serial
        score, choice, nums = witness(box)
        serial += 1
        item = (-score-log(locations*mult), serial, box, mult, choice, nums)
        heappush(heap, item)
        tree.add(item, parent)
    push(root_box, 1)
    threshold = -(args.target_bits+2)*log(2)
    score = -heap[0][0]
    steps = 0
    while heap and steps < args.max_splits and score > threshold:
        item = heappop(heap)
        if all(lo == hi for lo, hi in item[2]):
            continue
        for child, mult in old.split(item[2], item[3]):
            push(child, mult, item[1])
        tree.update(item[1])
        score = tree.nodes[1]['best']
        steps += 1
        if steps % 25 == 0:
            print('PACKED cover q', args.groups, 'splits', steps, 'log2 upper proposal', score/log(2), flush=True)
    selected = tree.selected()
    checked_partition(tree, selected, root_box)
    score = float(logsumexp([-x[0] for x in selected]))
    print('PACKED complete support domain q', args.groups, 'floor', support_min,
          'splits', steps, 'leaves', len(selected), 'log2 proposal', score/log(2), flush=True)
    details = dict(method='full support-box CDF cover', support_min=support_min,
        splits=steps, domain_volume=(257-support_min)**args.groups,
        proposal_log2=score/log(2), leaves=[dict(box=item[2], multiplicity=item[3],
            tilt=item[4][0], numerators=item[5], proposal_log2=-item[0]/log(2)) for item in selected])
    # Retain the full split geometry so an external replay can authenticate
    # coverage without rerunning floating optimization.
    details['tree'] = [dict(id=i, parent=n['parent'], children=n['children'],
        box=n['item'][2], multiplicity=n['item'][3]) for i, n in tree.nodes.items()]
    details['selected_ids'] = [x[1] for x in selected]
    if score > threshold:
        return None, details
    return replay(details, operators, counts, terminal, groups=args.groups,
        cutoff=cutoff, precision=args.precision, target_bits=args.target_bits), details


def replay(details, operators, counts, terminal, *, groups, cutoff, precision, target_bits):
    """Verify saved rational witnesses and split geometry without optimizing."""
    root_box = domain(counts, groups, details['support_min'])
    if (type(cutoff) is not int or not 0 <= cutoff < (1 << 21)
            or type(precision) is not int or precision < 128
            or type(target_bits) is not int or target_bits < 1):
        raise ValueError('valid cutoff, precision and target required')
    size = len(terminal)
    if (size < 1 or any(len(exact) <= groups
            or any(matrix.nrows() != size or matrix.ncols() != size for matrix in exact[:groups+1])
            for exact, _ in operators.values())):
        raise ValueError('replay needs every occupancy coefficient with consistent square dimension')
    tree = old.RetainedCover()
    if details['domain_volume'] != old.volume(root_box, 1):
        raise ValueError('recorded volume differs from the complete support domain')
    for node in details['tree']:
        identifier = node['id']
        if type(identifier) is not int or identifier < 1 or identifier in tree.nodes:
            raise ValueError('unique positive integer tree identifiers required')
        box = tuple(tuple(pair) for pair in node['box'])
        mult = node['multiplicity']
        if (len(box) != groups or box != tuple(sorted(box)) or type(mult) is not int or mult < 1
                or any(len(pair) != 2 or any(type(x) is not int for x in pair)
                       or not details['support_min'] <= pair[0] <= pair[1] <= 256 for pair in box)):
            raise ValueError('valid sorted support box and multiplicity required')
        tree.nodes[identifier] = dict(item=(0., identifier, box, mult),
            parent=node['parent'], children=node['children'])
    ids = details['selected_ids']
    if len(ids) != len(details['leaves']) or any(type(i) is not int or i not in tree.nodes for i in ids):
        raise ValueError('selected leaves must match tree nodes')
    selected = []
    for identifier, leaf in zip(ids, details['leaves']):
        node = tree.nodes[identifier]
        box = tuple(tuple(pair) for pair in leaf['box'])
        mult = leaf['multiplicity']
        nums = leaf['numerators']
        choice = leaf['tilt'], '1'
        if (node['item'][2:4] != (box, mult) or choice not in operators
                or Q(choice[0]) <= 0 or len(nums) != groups
                or any(type(n) is not int or not 0 < n < old.DENOMINATOR for n in nums)):
            raise ValueError('matching geometry and valid exact probabilities required')
        selected.append((0., identifier, box, mult, choice, nums))
    checked_partition(tree, selected, root_box)
    ctx.prec = precision
    locations = comb(2048, groups)
    total = arb(0)
    folds = {}
    for _, _, box, mult, choice, nums in selected:
        masses = old.bernoulli_masses(nums)
        matrix = sum((m*r for m, r in zip(masses, operators[choice][0])), arb_mat(size, size))**256
        term = sum((matrix[0, j] for j in range(size) if terminal[j]), arb(0))
        term *= (aq(Q(choice[0]))*cutoff).exp()*locations*mult
        for (lo, hi), numerator in zip(box, nums):
            key = lo, hi, numerator
            if key not in folds:
                folds[key] = fold_arb(counts, lo, hi, arb(numerator)/old.DENOMINATOR)
            term = up(term*folds[key])
        total = up(total+term)
    if not 0 < total < arb(2)**(-target_bits):
        raise ArithmeticError('outward replay failed the requested margin')
    return total
