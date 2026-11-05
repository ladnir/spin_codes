r"""Exact bounded checks for QUARTET_SHORT_MATCHING_HALL_COUNT.md.

J is the non-H input set. Under2|J|=nT, the short-matching deficit equals
max_(C subset J)(2|C|-|N^-(C)|). Canonical saturation fixes C1=J1 and,
at later updates, fills every non-H input in every selected quartet.
Selected quartets must meet J. Their future unions give the reverse
old-label state A subset B and the positive15-composition local kernel.

These checks independently enumerate input subsets, predecessor arcs,
canonical selector bits, integer matching cuts, and rational y2/y4
weights. They verify both flip normalizations and the asymmetric random
cover identity under exact balance. Unbalanced controls retain the exact
offset, rather than silently applying the balanced identity.

Width4,T2/3 uses every balanced label pattern and every input subset.
Six width8 instances use fresh varying partitions and at most six genuine
quartet selector bits. A tiny three-state polynomial action also checks
the exact total-quota coefficient and the forced first-update endpoint.
No large operator, probability certificate, or numerical job on import.
"""

import argparse
from fractions import Fraction
from itertools import combinations, product
from math import comb
import os
import shutil
import threading
import time


def _dh_dependency(name):
    if name in globals():
        return globals()[name]
    if name == 'short_matching_solve':
        import quartet_short_matching as module
    else:
        import quartet_fixed_row_flow as module
    return getattr(module, name)


def _dh_geometry(J, partitions, width):
    full = (1 << width)-1
    H = tuple(full ^ mask for mask in J)
    _dh_dependency('_fr_validate_instance')(H, partitions, width)
    quartet_masks, by_coordinate = [], []
    for partition in partitions:
        masks = tuple(sum(1 << i for i in quartet) for quartet in partition)
        lookup = [0]*width
        for quartet, mask in zip(partition, masks):
            for coordinate in quartet:
                lookup[coordinate] = mask
        quartet_masks.append(masks)
        by_coordinate.append(tuple(lookup))
    return H, tuple(quartet_masks), tuple(by_coordinate)


def _dh_predecessors(C, by_coordinate):
    """Expand physical input predecessor arcs, with no formal endpoint input."""
    predecessors = list(C)
    for t in range(1, len(C)):
        mask, earlier = C[t], 0
        while mask:
            bit = mask & -mask
            earlier |= by_coordinate[t][bit.bit_length()-1]
            mask -= bit
        predecessors[t-1] |= earlier
        if t >= 2:
            predecessors[t-2] |= earlier
    return tuple(predecessors)


def _dh_score(C, by_coordinate):
    predecessors = _dh_predecessors(C, by_coordinate)
    return 2*sum(mask.bit_count() for mask in C)-sum(mask.bit_count() for mask in predecessors)


def _dh_saturate(C, J, quartet_masks):
    result, selectors = [J[0]], [0]
    for t in range(1, len(C)):
        selected = 0
        for quartet in quartet_masks[t]:
            if C[t] & quartet:
                selected |= quartet
        selectors.append(selected)
        result.append(selected & J[t])
    return tuple(result), tuple(selectors)


def _dh_input_subsets(J):
    """Gray enumeration changes one input bit, without rebuilding every subset."""
    positions = []
    for t, mask in enumerate(J):
        while mask:
            bit = mask & -mask
            positions.append((t, bit))
            mask -= bit
    current = [0]*len(J)
    yield tuple(current)
    for counter in range(1, 1 << len(positions)):
        position = (counter & -counter).bit_length()-1
        t, bit = positions[position]
        current[t] ^= bit
        yield tuple(current)


def _dh_canonical_witnesses(J, quartet_masks, by_coordinate):
    all_quartets = [(t, quartet, (quartet & J[t]).bit_count())
                    for t in range(1, len(J)) for quartet in quartet_masks[t]]
    candidates = [entry for entry in all_quartets if entry[2]]
    witnesses = []
    for bits in range(1 << len(candidates)):
        U, C = [0]*len(J), [J[0]]+[0]*(len(J)-1)
        for i, (t, quartet, _) in enumerate(candidates):
            if bits >> i & 1:
                U[t] |= quartet
                C[t] |= quartet & J[t]
        predecessors = _dh_predecessors(C, by_coordinate)
        formula = []
        for t in range(len(J)):
            future = (U[t+1] if t+1 < len(J) else 0) | (U[t+2] if t+2 < len(J) else 0)
            formula.append(C[t] | future)
        assert tuple(formula) == predecessors
        assert C[0] == J[0]
        assert all(not (U[t] & quartet) or (quartet & J[t])
                   for t in range(1, len(J)) for quartet in quartet_masks[t])
        total_inputs = sum(mask.bit_count() for mask in C)
        predecessor_count = sum(mask.bit_count() for mask in predecessors)
        witnesses.append((tuple(C), 2*total_inputs-predecessor_count, predecessor_count))
    # Enumerating all formal bits independently must reject every all-H selection.
    accepted = set()
    for bits in range(1 << len(all_quartets)):
        C, valid = [J[0]]+[0]*(len(J)-1), True
        for i, (t, quartet, j) in enumerate(all_quartets):
            if bits >> i & 1:
                if j == 0:
                    valid = False
                    break
                C[t] |= quartet & J[t]
        if valid:
            accepted.add(tuple(C))
    assert accepted == {entry[0] for entry in witnesses}
    assert len(accepted) == 1 << len(candidates)
    return tuple(candidates), tuple(witnesses)


def _dh_weights(y, j, mode):
    if j == 0:
        return Fraction(1), Fraction(0)
    if mode == 'raw':
        return Fraction(1), Fraction(1)
    a, b = 8-j, 2*j
    if mode == 'symmetric':
        weight = 1/(1+y**(-max(a, b)))
        return weight, weight
    if mode != 'asymmetric':
        raise ValueError('Unknown normalization')
    A, B = y**(-a), y**(-b)
    return (1-A)/(1-A*B), (1-B)/(1-A*B)


def _dh_normalizations(J, candidates, witnesses, y, counters):
    scores = [entry[1] for entry in witnesses]
    delta = max(scores)
    raw = sum((y**score for score in scores), Fraction(0))
    denominator = Fraction(1)
    for _, _, j in candidates:
        denominator *= 1+y**(-max(8-j, 2*j))
    totals = {'raw': raw}
    for mode in ('symmetric', 'asymmetric'):
        weights = [_dh_weights(y, entry[2], mode) for entry in candidates]
        weighted = Fraction(0)
        probability_total = Fraction(0)
        cover_moment = Fraction(0)
        for bits, (_, score, predecessors) in enumerate(witnesses):
            weight = probability = Fraction(1)
            for i, (_, _, j) in enumerate(candidates):
                chosen = (bits >> i) & 1
                weight *= weights[i][chosen]
                if mode == 'asymmetric':
                    probability *= weights[i][1] if chosen else y**(-2*j)*weights[i][0]
            weighted += weight*y**score
            if mode == 'asymmetric':
                probability_total += probability
                cover_moment += probability*y**(len(J)*counters['width']-predecessors)
        assert y**delta <= weighted
        totals[mode] = weighted
        if mode == 'symmetric':
            assert weighted == raw/denominator
        else:
            assert probability_total == 1
            imbalance = len(J)*counters['width']-2*sum(mask.bit_count() for mask in J)
            assert cover_moment == weighted*y**imbalance
            if imbalance == 0:
                assert cover_moment == weighted
                counters['balanced_cover_identities'] += 1
            else:
                assert cover_moment != weighted
                counters['unbalanced_cover_controls'] += 1
        counters['normalized_majorants'] += 1
    return totals


def _dh_instance(J, partitions, width, counters, quota_totals=None):
    H, quartet_masks, by_coordinate = _dh_geometry(J, partitions, width)
    candidates, witnesses = _dh_canonical_witnesses(J, quartet_masks, by_coordinate)
    maximum = -width*len(J)
    for C in _dh_input_subsets(J):
        score = _dh_score(C, by_coordinate)
        maximum = max(maximum, score)
        saturated, U = _dh_saturate(C, J, quartet_masks)
        added = sum((new & ~old).bit_count() for old, new in zip(C, saturated))
        assert all(not (old & ~new) for old, new in zip(C, saturated))
        assert _dh_score(saturated, by_coordinate) >= score+added
        again, again_U = _dh_saturate(saturated, J, quartet_masks)
        assert again == saturated and again_U == U
        assert saturated[0] == J[0]
        counters['dual_input_subsets'] += 1
    canonical_maximum = max(entry[1] for entry in witnesses)
    assert maximum == canonical_maximum
    matching = _dh_dependency('short_matching_solve')(H, partitions, width)
    offset = 2*sum(mask.bit_count() for mask in J)-width*len(J)
    assert maximum == matching['integer_deficit']+offset
    counters['dual_matching_maxima'] += 1
    counters['canonical_witnesses'] += len(witnesses)
    # Check additions and all bounded-loss paths around arbitrary reference witnesses.
    scores = [entry[1] for entry in witnesses]
    for bits, score in enumerate(scores):
        for i, (_, _, j) in enumerate(candidates):
            if not (bits >> i) & 1:
                difference = scores[bits | (1 << i)]-score
                assert j-8 <= difference <= 2*j
                counters['single_flip_bounds'] += 1
        for target, target_score in enumerate(scores):
            loss = 0
            changed = bits ^ target
            for i, (_, _, j) in enumerate(candidates):
                if changed >> i & 1:
                    loss += 2*j if bits >> i & 1 else 8-j
            assert target_score >= score-loss
            counters['multiple_flip_bounds'] += 1
    counters['width'] = width
    for y in (Fraction(2), Fraction(4)):
        totals = _dh_normalizations(J, candidates, witnesses, y, counters)
        if quota_totals is not None:
            for mode, value in totals.items():
                key = (len(J), y, mode)
                quota_totals[key] = quota_totals.get(key, Fraction(0))+value


def _dh_local_formula(y, mode):
    """Positive coefficient form of the homogeneous old/new-label polynomial."""
    result = {}
    for a in range(5):
        for e in range(5-a):
            r = a+e
            multiplicity = comb(4, a)*comb(4-a, e)
            for j in range(5):
                absent, present = _dh_weights(y, j, mode)
                result[(a, e, 0, a, j)] = multiplicity*comb(4, j)*y**(-r)*absent
                if j:
                    overlap = sum((Fraction(comb(r, ell)*comb(4-r, j-ell))*y**ell
                                   for ell in range(max(0, j-(4-r)), min(r, j)+1)), Fraction(0))
                    result[(a, e, 4, 0, j)] = multiplicity*y**(j-r)*overlap*present
    return result


def _dh_local_physical(y, mode):
    """Enumerate actual coordinate labels and J masks; compute union sizes."""
    result = {}
    for labels in product(range(3), repeat=4):
        A = sum(1 << i for i, label in enumerate(labels) if label == 0)
        E = sum(1 << i for i, label in enumerate(labels) if label == 1)
        B = A | E
        for J in range(16):
            j = J.bit_count()
            for selected in range(2 if j else 1):
                U, C = (15, J) if selected else (0, 0)
                new_A, new_E = U, A & ~U
                score = 2*C.bit_count()-(C | B).bit_count()
                key = (A.bit_count(), E.bit_count(), new_A.bit_count(), new_E.bit_count(), j)
                weight = _dh_weights(y, j, mode)[selected]*y**score
                result[key] = result.get(key, Fraction(0))+weight
    return result


def _dh_terminal_coefficients(y, b):
    return [sum((Fraction(comb(b, ell)*comb(4-b, j-ell))*y**(j-b+ell)
                 for ell in range(max(0, j-(4-b)), min(b, j)+1)), Fraction(0))
            for j in range(5)]


def _dh_poly_add_product(target, left, right):
    for i, a in enumerate(left):
        if a:
            for j, b in enumerate(right):
                if b:
                    target[i+j] += a*b


def _dh_toy_quota_coefficient(updates, y, mode):
    """Width4 exact reverse-time action; no dense production operator."""
    states = ((0, 0), (0, 4), (4, 4))
    lookup = {state: i for i, state in enumerate(states)}
    local = _dh_local_formula(y, mode)
    vector = [[Fraction(1)], [Fraction(0)], [Fraction(0)]]
    for _ in range(updates-1):
        degree = len(vector[0])+4
        after = [[Fraction(0)]*degree for _ in states]
        for source, (a, b) in enumerate(states):
            multiplicity = comb(4, a)*comb(4-a, b-a)
            absent = [local[(a, b-a, 0, a, j)]/multiplicity for j in range(5)]
            selected = [Fraction(0)]+[local[(a, b-a, 4, 0, j)]/multiplicity for j in range(1, 5)]
            _dh_poly_add_product(after[lookup[(0, a)]], vector[source], absent)
            _dh_poly_add_product(after[lookup[(4, 4)]], vector[source], selected)
        vector = after
    output = [Fraction(0)]*(4*updates+1)
    for source, (_, b) in enumerate(states):
        _dh_poly_add_product(output, vector[source], _dh_terminal_coefficients(y, b))
    return output[2*updates]


def dual_hall_checks():
    counters = dict(dual_input_subsets=0, dual_matching_maxima=0, canonical_witnesses=0,
                    single_flip_bounds=0, multiple_flip_bounds=0, normalized_majorants=0,
                    balanced_cover_identities=0, unbalanced_cover_controls=0, width=4)
    quota_totals = {}
    quartet = ((0, 1, 2, 3),)
    for updates in (2, 3):
        for selected in combinations(range(4*updates), 2*updates):
            packed = sum(1 << slot for slot in selected)
            J = tuple((packed >> (4*t)) & 15 for t in range(updates))
            _dh_instance(J, (quartet,)*updates, 4, counters, quota_totals)
    for updates, seeds in ((3, (0, 1, 7, 19)), (4, (0, 7))):
        for seed in seeds:
            region = 3 if updates == 3 else 2
            H, partitions = _dh_dependency('fixed_row_flow_sample')(8, updates, region, 4*region, seed)
            J = tuple(255 ^ mask for mask in H)
            _dh_instance(J, partitions, 8, counters)
    for J in ((0, 0), (15, 15), (1, 0), (15, 7)):
        _dh_instance(J, (quartet, quartet), 4, counters)
    coefficients = terminal = weight_identities = quota_coefficients = 0
    for y in (Fraction(2), Fraction(4)):
        for j in range(1, 5):
            a, b = 8-j, 2*j
            w0, w1 = _dh_weights(y, j, 'asymmetric')
            assert w0 > 0 and w1 > 0
            assert w0+y**(-a)*w1 == 1
            assert y**(-b)*w0+w1 == 1
            weight_identities += 2
        for mode in ('raw', 'symmetric', 'asymmetric'):
            formula, physical = _dh_local_formula(y, mode), _dh_local_physical(y, mode)
            assert formula == physical
            assert all(coefficient > 0 for coefficient in formula.values())
            coefficients += len(formula)
            for updates in (2, 3):
                coefficient = _dh_toy_quota_coefficient(updates, y, mode)
                assert coefficient == quota_totals[(updates, y, mode)]
                quota_coefficients += 1
        for a in range(5):
            for e in range(5-a):
                b = a+e
                B = (1 << b)-1
                physical = [Fraction(0)]*5
                for J in range(16):
                    physical[J.bit_count()] += y**(2*J.bit_count()-(J | B).bit_count())
                assert physical == _dh_terminal_coefficients(y, b)
                terminal += 5
        # At T1, C1=J1 is the only canonical witness; there is no U1 bit.
        assert _dh_toy_quota_coefficient(1, y, 'asymmetric') == comb(4, 2)*y**2
        terminal += 1
    counters.pop('width')
    print(dict(**counters, old_label_compositions=15,
               positive_local_coefficient_checks=coefficients,
               exact_terminal_coefficient_checks=terminal,
               asymmetric_weight_identities=weight_identities,
               exact_quota_coefficient_checks=quota_coefficients,
               rational_markers=(2, 4), large_operator_built=False,
               full_code_certificate=False), flush=True)


def dual_hall_main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--checks', action='store_true', required=True)
    parser.parse_args()
    import resource
    if shutil.disk_usage('/tmp').free < 20*(1 << 30):
        raise RuntimeError('20GiB /tmp free-space guard')
    resource.setrlimit(resource.RLIMIT_AS, (2*(1 << 30), 2*(1 << 30)))
    timer = threading.Timer(175, lambda: os._exit(124))
    timer.daemon = True
    timer.start()
    try:
        dual_hall_checks()
    finally:
        timer.cancel()


if __name__ == '__main__':
    dual_hall_main()
