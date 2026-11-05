"""Exact small-quotient tests for the relaxed coupled rank transfer.

No K20 certificate, simulation, or benchmark. Enumerates bounded GF(16)
linear-form fixtures and tiny pair permutations. Writes no files; use -B.
"""

from collections import Counter, defaultdict
from fractions import Fraction as F
from itertools import permutations, product

from coupled_rank_transfer import (symbolic_transitions, route_profile_normalizer,
                                   local_row, paired_route_row, physical_pair_weights,
                                   pair_antipodal_cycle)
from test_local_algebra import MUL, COEFFICIENTS, require


INV = {a: next(b for b in range(1, 16) if MUL[a][b] == 1) for a in range(1, 16)}


def scale(a, form):
    return MUL[a][form & 15] | (MUL[a][form >> 4] << 4)


def reduce_form(form, basis):
    for pivot, row in enumerate(basis):
        coefficient = (form >> (4 * pivot)) & 15
        if row and coefficient:
            form ^= scale(coefficient, row)
    return form


def zero_span(forms):
    basis = [0, 0]
    for form in forms:
        form = reduce_form(form, basis)
        if form:
            pivot = 0 if (form & 15) else 1
            basis[pivot] = scale(INV[(form >> (4 * pivot)) & 15], form)
    return tuple(basis)


def support(forms):
    return sum(1 << j for j, form in enumerate(forms) if form)


def branch_potential(inputs, state, outputs, following, out_mask, next_mask, history, tau=2):
    zeros = [form for j, form in enumerate(outputs) if not (out_mask >> j) & 1]
    zeros += [form for j, form in enumerate(following) if not (next_mask >> j) & 1]
    basis = zero_span(zeros)
    required = [form for form in inputs + state if form] + list(history)
    required += [form for j, form in enumerate(outputs) if (out_mask >> j) & 1]
    required += [form for j, form in enumerate(following) if (next_mask >> j) & 1]
    if any(reduce_form(form, basis) == 0 for form in required):
        return F(0)
    return F(1, tau**sum(bool(row) for row in basis))


def check_quotients(tau=2):
    # Forms on F^2: zero, two independent bases, their sum, and theta times e0.
    # The last form exercises exceptional-coefficient branches.
    forms = (0, 1, 16, 17, 2)
    fixtures = rank_branches = dependent_branches = 0
    y = F(1, 3)
    for values in product(forms, repeat=4):
        inputs, state = values[:2], values[2:]
        input_mask, state_mask = support(inputs), support(state)
        symbolic = symbolic_transitions(2, input_mask, state_mask)
        upper = {
            outgoing: sum((multiplicity * y**weight / min(tau, 2)**charge
                           for weight, charge, multiplicity in terms), F(0))
            for outgoing, terms in symbolic.items()
        }
        # A pre-existing nonzero requirement tests that old branches are retained,
        # rather than conditioning and renormalizing the coefficient distribution.
        for history in ((), (33,)):
            actual = defaultdict(F)
            for coefficients in product(COEFFICIENTS, repeat=2):
                outputs = tuple(x ^ scale(a, q)
                                for x, a, q in zip(inputs, coefficients, state))
                following = (state[1] ^ outputs[0], state[0] ^ outputs[1])
                for out_mask, next_mask in product(range(4), repeat=2):
                    potential = branch_potential(inputs, state, outputs, following,
                                                 out_mask, next_mask, history, tau)
                    actual[next_mask] += potential * y**out_mask.bit_count() / 4
                    if potential and potential < 1:
                        rank_branches += 1
                    if potential == 1 and state_mask and out_mask != 3:
                        dependent_branches += 1
            for next_mask in range(4):
                require(actual[next_mask] <= upper.get(next_mask, F(0)),
                        ("quotient domination", inputs, state, history, next_mask,
                         actual[next_mask], upper.get(next_mask, F(0))))
            fixtures += 1
    require(fixtures == 1250, "fixture coverage")
    require(rank_branches > 0 and dependent_branches > 0, "missing quotient case")
    return fixtures, rank_branches, dependent_branches


def check_pair_profiles():
    profiles = sequences = 0
    # Enumerate labeled pair permutations and all independent swaps.
    for count in range(1, 5):
        for n0 in range(count + 1):
            for n1 in range(count - n0 + 1):
                n2 = count - n0 - n1
                labels = [0] * n0 + [1] * n1 + [3] * n2
                outcomes = Counter()
                for order in permutations(range(count)):
                    for swaps in product((0, 1), repeat=count):
                        word = tuple(2 if labels[index] == 1 and flip else labels[index]
                                     for index, flip in zip(order, swaps))
                        outcomes[word] += 1
                normalizer = route_profile_normalizer(count, n0, n1, n2)
                require(len(outcomes) == normalizer * 2**n1, "pair profile support size")
                total = sum(outcomes.values())
                for multiplicity in outcomes.values():
                    require(F(multiplicity, total) == F(1, normalizer * 2**n1),
                            "pair profile probability")
                profiles += 1
                sequences += total
    return profiles, sequences


def check_generating_rows():
    cases = 0
    y, tau = F(2, 5), F(7, 3)
    x = (F(3, 7), F(2, 9), F(5, 11))
    for width in (2, 4, 8):
        weights = physical_pair_weights(width, *x)
        require(sum(weights) == sum(x)**(width // 2), "pair generating total")
        if width == 8:
            continue  # Do not materialize a costly width-eight generating matrix.
        for cycle in (tuple(range(width)), pair_antipodal_cycle(width)):
            for state in (0, 1, (1 << width) - 1):
                expected = [F(0)] * (1 << width)
                for inputs, weight in enumerate(weights):
                    row = local_row(width, inputs, state, y, tau, cycle)
                    for outgoing, value in enumerate(row):
                        expected[outgoing] += weight * value
                actual = paired_route_row(width, state, y, tau, *x, cycle)
                require(actual == tuple(expected), "paired generating row")
                cases += 1
    return cases


def check_shortening_bridge():
    # Binary [3,2] parity code, extended to GF(16). The full support exists
    # only after extension: (1,2,3), whereas every nonzero binary word has weight2.
    columns = (1, 16, 17)
    binary_words = ((1, 0, 1), (0, 1, 1), (1, 1, 0))
    require(support((1, 2, 3)) == 7 and all(support(c) != 7 for c in binary_words),
            "exact-support extension negative control")
    u = F(2, 7)
    lhs = F(0)
    active = []
    for mask in range(8):
        zeros = [form for j, form in enumerate(columns) if not (mask >> j) & 1]
        dimension = 2 - sum(bool(row) for row in zero_span(zeros))
        within = sum((support(c) & ~mask) == 0 for c in binary_words)
        require(within == 2**dimension - 1, "shortening cardinality")
        lhs += u**mask.bit_count() * within
        if dimension:
            require(2**(dimension - 1) <= within, "single-group potential bridge")
            active.append((dimension, within))
    rhs = sum((u**support(c).bit_count() * (1 + u)**(3 - support(c).bit_count())
               for c in binary_words), F(0))
    require(lhs == rhs, "support-containment transform")
    for (d1, count1), (d2, count2) in product(active, repeat=2):
        require(2**(d1 + d2 - 1) <= 2 * count1 * count2, "two-group bridge factor")
    return len(active)**2


def main():
    # Explicit linear-span sanity checks independent of the transfer code.
    require(reduce_form(17, zero_span((1, 16))) == 0, "rank-two span")
    require(reduce_form(2, zero_span((1,))) == 0, "dependent span")
    require(reduce_form(16, zero_span((1,))) != 0, "independent form")
    fixtures, rank_branches, dependent_branches = check_quotients()
    large_tau_fixtures, _, _ = check_quotients(tau=8)
    profiles, sequences = check_pair_profiles()
    rows = check_generating_rows()
    bridge_cases = check_shortening_bridge()
    print(f"quotient fixtures={fixtures}; retained rank-gaining branches={rank_branches}; "
          f"retained dependent/zero-rank branches={dependent_branches}")
    print(f"tiny pair profiles={profiles}; labeled permutation/swap samples={sequences}")
    print(f"tau=8 > q=2 fixed-mask fixtures={large_tau_fixtures}; coefficient charge capped at q")
    print(f"paired generating rows={rows}; two-group shortening fixtures={bridge_cases}")
    print("PASS: local common-quotient domination and paired normalization; not a certificate")


if __name__ == "__main__":
    main()
