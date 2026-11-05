"""Exact width-two/GF16 checks of the finite-state-rank Bellman bound.

Uses the degree-four theta with coefficient set {theta,theta+1}; its two
roots give q=2 in this toy test. No K20 certificate or numerical files.
"""

from fractions import Fraction as F
from itertools import product

from empty_span_rank import empty_step, terminal_ones
from test_local_algebra import COEFFICIENTS, require
from test_rank_transfer import zero_span, reduce_form, scale, support


def actual(state, history, length, terminal, z, tau):
    if not length:
        rank = sum(bool(row) for row in zero_span(state))
        return terminal[rank][support(state)]
    total = F(0)
    for coefficients in product(COEFFICIENTS, repeat=2):
        following = (state[1] ^ scale(coefficients[0], state[0]),
                     state[0] ^ scale(coefficients[1], state[1]))
        for mask in range(1, 4):
            zeros = [form for j, form in enumerate(following) if not (mask >> j) & 1]
            basis = zero_span(zeros)
            required = tuple(history) + tuple(form for form in state if form)
            required += tuple(form for j, form in enumerate(following) if (mask >> j) & 1)
            if any(reduce_form(form, basis) == 0 for form in required):
                continue
            child = tuple(reduce_form(form, basis) for form in following)
            requirements = tuple(set(reduce_form(form, basis) for form in required))
            gained = sum(bool(row) for row in basis)
            total += F(1, 4 * tau**gained) * actual(child, requirements, length - 1, terminal, z, tau)
    return z**support(state).bit_count() * total


def main():
    fixtures = 0
    z = F(3, 5)
    terminals = [terminal_ones(2, F(1)), terminal_ones(2, F(1))]
    for rank in (1, 2):
        for mask in range(1, 4):
            if rank <= mask.bit_count():
                terminals[1][rank][mask] = F(1 + 3 * rank + 2 * mask, 11)
    forms = (0, 1, 2, 16, 17, 33)
    for terminal, tau in product(terminals, (1, 2, 8)):
        for length in (1, 2):
            bound = terminal
            for _ in range(length):
                bound = empty_step(bound, 2, z, F(tau), coefficient_count=2)
            for state in product(forms, repeat=2):
                if not any(state):
                    continue
                rank = sum(bool(row) for row in zero_span(state))
                mask = support(state)
                for history in ((), (19,)):
                    value = actual(state, history, length, terminal, z, tau)
                    require(value <= bound[rank][mask],
                            ("rank Bellman", tau, terminal, length, state, history, value, bound[rank][mask]))
                    fixtures += 1
    print(f"PASS: {fixtures} exact empty-span fixtures; tau=1,2,8; q=2; lengths1/2; rank1/2.")


if __name__ == "__main__":
    main()
