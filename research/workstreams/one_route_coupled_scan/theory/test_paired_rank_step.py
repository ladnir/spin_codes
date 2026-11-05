"""Bounded exact GF16 checks of the paired active-step rank bound.

Uses an independent Gaussian oracle on up to six-variable nibble forms.
The coefficient set is {2,3} in GF(2)[t]/(t^4+t+1), so q=2 here.
No distance certificate, performance claim, file output, or random sweep.
Run serially with Python -B. There are 144 terminal/price comparisons.
"""

from collections import Counter
from fractions import Fraction as F
from itertools import product

from paired_rank_step import active_step


def require(condition, detail):
    if not condition:
        raise AssertionError(detail)


def field_product(left, right):
    answer = 0
    for _ in range(4):
        if right & 1:
            answer ^= left
        right >>= 1
        left <<= 1
        if left & 16:
            left ^= 0x13
    return answer


MUL = tuple(tuple(field_product(a, b) for b in range(16)) for a in range(16))
INV = (0,) + tuple(next(b for b in range(1, 16) if MUL[a][b] == 1)
                   for a in range(1, 16))
COEFFICIENTS = (2, 3)


def packed(*coefficients):
    return sum(value << (4 * j) for j, value in enumerate(coefficients))


def scale(coefficient, word, variables):
    return sum(MUL[coefficient][(word >> (4 * j)) & 15] << (4 * j)
               for j in range(variables))


def reduce_word(word, basis, variables):
    for pivot, row in basis:
        coefficient = (word >> (4 * pivot)) & 15
        if coefficient:
            word ^= scale(coefficient, row, variables)
    return word


def span(words, variables):
    """Canonical RREF, using highest nonzero coordinates as pivots."""
    rows = {}
    for word in words:
        require(0 <= word < (1 << (4 * variables)), "form exceeds ambient space")
        ordered = tuple(sorted(rows.items(), reverse=True))
        word = reduce_word(word, ordered, variables)
        if not word:
            continue
        pivot = (word.bit_length() - 1) // 4
        word = scale(INV[(word >> (4 * pivot)) & 15], word, variables)
        for old_pivot, old_row in tuple(rows.items()):
            coefficient = (old_row >> (4 * pivot)) & 15
            if coefficient:
                rows[old_pivot] = old_row ^ scale(coefficient, word, variables)
        rows[pivot] = word
    return tuple(sorted(rows.items(), reverse=True))


def reduced(words, basis, variables):
    return tuple(reduce_word(word, basis, variables) for word in words)


def support(words):
    return sum(1 << j for j, word in enumerate(words) if word)


def rotated(words, cycle):
    result = [0] * len(words)
    for index, source in enumerate(cycle):
        result[cycle[(index + 1) % len(cycle)]] = words[source]
    return tuple(result)


def image_data(state, inputs, zero_basis, variables):
    q = reduced(state, zero_basis, variables)
    x = reduced(inputs, zero_basis, variables)
    old = span(q, variables)
    augmented = span(q + x, variables)
    outside = support(reduced(x, old, variables))
    return len(old), len(augmented), outside


def check_pivots(state, inputs, zero_equations, variables, coverage):
    """Check each independent equation's old/augmented/outside update."""
    imposed = []
    basis = ()
    for coordinate, equation in zero_equations:
        residual = reduce_word(equation, basis, variables)
        if not residual:
            continue
        a, b, outside = image_data(state, inputs, basis, variables)
        old = span(reduced(state, basis, variables), variables)
        is_outside = bool(reduce_word(residual, old, variables))
        imposed.append(equation)
        basis = span(imposed, variables)
        aa, bb, pp = image_data(state, inputs, basis, variables)
        require(bb == b - 1, ("augmented rank update", a, b, aa, bb))
        if is_outside:
            bit = 1 << coordinate
            require(inputs[coordinate] and outside & bit, "outside pivot must use its active input")
            expected = 0 if b - a == 1 else outside & ~bit
            require(aa == a and pp == expected,
                    ("outside-image update", a, b, outside, aa, bb, pp, expected))
            coverage["outside pivots"] += 1
        else:
            require(aa == a - 1 and pp == outside,
                    ("old-span update", a, b, outside, aa, bb, pp))
            coverage["old-span pivots"] += 1


def exact_profile(fixture, coverage):
    name, variables, state, inputs, history, initial_zeros, fresh = fixture
    width = len(state)
    cycle = (0, 1) if width == 2 else (0, 2, 1, 3)
    initial_basis = span(initial_zeros, variables)
    state = reduced(state, initial_basis, variables)
    inputs = reduced(inputs, initial_basis, variables)
    history = reduced(history, initial_basis, variables)
    require(all(history), (name, "initial history is inconsistent"))
    input_mask, state_mask = support(inputs), support(state)
    rank = len(span(state, variables))
    augmented = len(span(state + inputs, variables))
    k = input_mask.bit_count()
    require(1 <= k <= 2, (name, "input-mask size"))
    if fresh:
        require(augmented == rank + k, (name, "freshness precondition"))
    profile = Counter()
    old_requirements = tuple(word for word in state + inputs if word) + history
    for coefficients in product(COEFFICIENTS, repeat=width):
        outputs = tuple(x ^ scale(a, q, variables)
                        for x, a, q in zip(inputs, coefficients, state))
        prior = rotated(state, cycle)
        following = tuple(v ^ y for v, y in zip(prior, outputs))
        require(span(following + inputs, variables) == span(state + inputs, variables),
                (name, "unquotiented augmented-span identity", coefficients))
        for output_mask, next_mask in product(range(1 << width), repeat=2):
            # Identically zero coordinates cannot meet a nonzero requirement.
            if output_mask & ~support(outputs) or next_mask & ~support(following):
                continue
            zeros = [(j, word) for j, word in enumerate(outputs) if not output_mask & (1 << j)]
            zeros += [(j, word) for j, word in enumerate(following) if not next_mask & (1 << j)]
            basis = span((word for _, word in zeros), variables)
            required = old_requirements
            required += tuple(word for j, word in enumerate(outputs) if output_mask & (1 << j))
            required += tuple(word for j, word in enumerate(following) if next_mask & (1 << j))
            if any(not reduce_word(word, basis, variables) for word in required):
                continue
            gain = len(basis)
            q_next = reduced(following, basis, variables)
            next_rank = len(span(q_next, variables))
            a, b, outside = image_data(state, inputs, basis, variables)
            require(b == augmented - gain, (name, "rank conservation"))
            require(span(q_next + reduced(inputs, basis, variables), variables)
                    == span(reduced(state + inputs, basis, variables), variables),
                    (name, "quotiented augmented-span identity"))
            require(max(b - k, b - a, 0) <= next_rank <= min(b, next_mask.bit_count()),
                    (name, "outgoing rank interval", a, b, k, next_rank, next_mask))
            require((next_rank == 0) == (next_mask == 0), (name, "quiet rank/support"))
            require((outside == 0) == (b == a), (name, "outside image rank"))
            # Match the operator's coordinate order, but not its classification.
            check_pivots(state, inputs, sorted(zeros, key=lambda pair: pair[0]), variables, coverage)
            profile[gain, next_rank, next_mask, output_mask.bit_count()] += 1
            coverage["retained branches"] += 1
            coverage["positive rank gain" if gain else "zero rank gain"] += 1
            if not next_mask:
                coverage["resets"] += 1
    require(profile, (name, "empty exact profile"))
    return profile, rank, state_mask, input_mask, cycle


def fixtures():
    e = tuple(1 << (4 * j) for j in range(6))
    return (
        ("quiet singleton", 1, (0, 0), (e[0], 0), (), (), True),
        ("inside span and history", 2, (e[0], e[1]), (e[0] ^ e[1], 0),
         (packed(1, 2),), (), False),
        ("dependent coefficient roots", 1, (e[0], packed(2)), (0, packed(3)), (), (), False),
        ("two fresh inputs", 4, (e[0], e[1]), (e[2] ^ e[3], packed(0, 0, 2, 1)), (), (), True),
        ("shared outside image", 4, (e[0], e[1]), (e[0] ^ e[2], e[1] ^ e[3]),
         (packed(1, 2),), (e[2] ^ e[3],), False),
        ("rank-one state, two input images", 2, (e[0], e[0]), (e[1], e[0] ^ e[1]),
         (packed(1, 2),), (), False),
        ("quiet physical pair", 2, (0, 0, 0, 0), (e[0], e[1], 0, 0), (), (), True),
        ("six-variable fresh pair", 6, (e[0], e[1], e[2], e[3]),
         (packed(0, 0, 0, 0, 1, 2), packed(0, 0, 0, 0, 3, 1), 0, 0), (), (), True),
        ("four coordinates, old rank two", 2,
         (e[0], e[1], e[0] ^ e[1], packed(2, 3)),
         (packed(1, 2), packed(3, 1), 0, 0), (packed(1, 3),), (), False),
        ("paired reset geometry", 2, (e[0], 0, e[1], 0),
         (packed(2), e[1], 0, 0), (), (), False),
        ("cross pair, quotient-correlated inputs", 4,
         (e[0], e[1], packed(2, 1), packed(1, 3)),
         (e[0] ^ e[2], 0, e[1] ^ e[3], 0), (packed(1, 2, 3),), (e[2] ^ e[3],), False),
        ("fresh singleton, old rank three", 4, (e[0], e[1], e[2], e[0] ^ e[1]),
         (0, 0, e[3], 0), (packed(1, 2),), (), True),
    )


def terminals(width):
    tables = [[[F(0) for _ in range(1 << width)] for _ in range(width + 1)]
              for _ in range(4)]
    tables[0][0][0] = F(1)
    tables[2][0][0] = F(3, 7)
    for mask in range(1, 1 << width):
        for rank in range(1, mask.bit_count() + 1):
            tables[1][rank][mask] = F(1)
            tables[2][rank][mask] = F(2 + 5 * rank + 3 * mask, 11 + width)
    tables[3][1][1] = F(1)
    tables[3][1][(1 << width) - 1] = F(7, 5)
    tables[3][2][(1 << width) - 1] = F(5, 4)
    return tables


def main():
    require(MUL[2][8] == 3 and MUL[3][3] == 5, "GF16 polynomial convention")
    for variables in (2, 6):
        e0, e1 = 1, 16
        require(span((e0, e1), variables) == span((e0 ^ e1, scale(3, e1, variables)), variables),
                "Gaussian basis independence")
        require(len(span((e0, scale(2, e0, variables)), variables)) == 1, "dependent forms")
    coverage = Counter()
    comparisons = 0
    z = F(3, 5)
    for fixture in fixtures():
        profile, rank, state_mask, input_mask, cycle = exact_profile(fixture, coverage)
        name, _, state, _, _, _, fresh = fixture
        width = len(state)
        for tau, terminal in product((1, 2, 8), terminals(width)):
            bound = active_step(terminal, width, input_mask, z, F(tau),
                                coefficient_count=2, cycle=cycle, fresh=fresh)
            actual = sum((F(count, len(COEFFICIENTS)**width * tau**gain)
                          * z**weight * terminal[out_rank][out_mask]
                          for (gain, out_rank, out_mask, weight), count in profile.items()), F(0))
            require(actual <= bound[rank][state_mask],
                    (name, "paired active domination", tau, terminal, actual, bound[rank][state_mask]))
            comparisons += 1
    require(comparisons == 144, ("fixture count", comparisons))
    for key in ("outside pivots", "old-span pivots", "positive rank gain", "zero rank gain", "resets"):
        require(coverage[key] > 0, ("missing coverage", key))
    print(f"PASS: {comparisons} exact paired-active comparisons across {len(fixtures())} form fixtures; "
          "widths2/4, tau1/2/8, q2, up to6 variables, fresh and ordinary inputs.")
    print("Retained exact branches and span checks:", dict(sorted(coverage.items())))
    print("No whole-code certificate or performance claim.")


if __name__ == "__main__":
    main()
