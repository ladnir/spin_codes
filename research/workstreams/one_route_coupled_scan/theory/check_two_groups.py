"""Bounded two-adjacent-group gate using the actual paired regional route.

No independent-with-replacement approximation. Two-state active envelopes
are uniform in all common quotients. Final arithmetic is positive/upward.
This first gate restricts both groups to the SAME full-column support.
"""

from decimal import Decimal as D, ROUND_CEILING, localcontext
from fractions import Fraction as F
from math import comb

from check_empty_span_potential import check
from check_sparse_full_gate import multiply, power, upward
from coupled_rank_transfer import local_matrix, pair_antipodal_cycle, rotate_state_mask


def fresh_terms(mask, state, z, tau):
    """Fresh source pivots are independent of all old forms and constraints."""
    rotated = rotate_state_mask(4, state, pair_antipodal_cycle(4))
    terms = []
    for output in range(16):
        if (mask ^ state) & ~output or output & ~(mask | state):
            continue
        for following in range(16):
            if (output ^ rotated) & ~following or following & ~(output | rotated):
                continue
            active_charges = (mask & ~(output & following)).bit_count()
            old_charges = (state & ~mask & ~(output & following)).bit_count()
            weight = z**output.bit_count() / (tau**active_charges * min(tau, F(256))**old_charges)
            terms.append((following, weight, active_charges))
    return terms


def active(vector, z, tau, masks, *, fresh=False):
    """Uniform envelope for each mask, then average its independent choice."""
    valid = [(r, h) for h in range(1, 16) for r in range(1, h.bit_count() + 1)]
    result = [[F(0), F(0)], [F(0), F(0)]]
    for mask in masks:
        weight = mask.bit_count()
        if fresh:
            terms = {h: fresh_terms(mask, h, z, tau) for h in range(1, 16)}
            a = z**weight * vector[weight][mask]
            b = max(sum((term for following, term, _ in terms[h] if not following), F(0))
                    / vector[r][h] for r, h in valid)
            d = max(sum((term * max(vector[s][following] for s in range(1, min(r + weight - gained, following.bit_count()) + 1))
                         for following, term, gained in terms[h] if following), F(0)) / vector[r][h]
                    for r, h in valid)
            result[0][1] += a / len(masks)
            result[1][0] += b / len(masks)
            result[1][1] += d / len(masks)
            continue
        matrix = local_matrix(4, mask, z, min(tau, F(256)), pair_antipodal_cycle(4))
        a = z**weight * max(vector[r][mask] for r in range(1, weight + 1))
        b = max(matrix[h][0] / vector[r][h] for r, h in valid)
        d = max(sum((matrix[h][following]
                     * max(vector[s][following] for s in range(1, min(r + weight, following.bit_count()) + 1))
                     for following in range(1, 16)), F(0)) / vector[r][h] for r, h in valid)
        result[0][1] += a / len(masks)
        result[1][0] += b / len(masks)
        result[1][1] += d / len(masks)
    return result


def poly_multiply(left, right, degree):
    zero = D(0)
    result = [[[zero, zero], [zero, zero]] for _ in range(degree + 1)]
    for i, a in enumerate(left):
        for j, b in enumerate(right[:degree + 1 - i]):
            product = multiply(a, b)
            for r in range(2):
                for c in range(2):
                    result[i + j][r][c] += product[r][c]
    return result


def insertion_moment(empty, event, steps, degree):
    """[u^degree](empty+u event)^steps, preserving chronological order."""
    result = [[[D(1), D(0)], [D(0), D(1)]]]
    base = [empty, event]
    while steps:
        if steps & 1:
            result = poly_multiply(result, base, degree)
        steps >>= 1
        if steps:
            base = poly_multiply(base, base, degree)
    return result[degree]


def regional_matrices(contraction, singleton, paired, cross, steps=2048):
    """Caller supplies Decimal matrices under upward rounding."""
    empty = [[D(1), D(0)], [D(0), upward(contraction)]]
    pair_region = insertion_moment(empty, paired, steps, 1)
    collision = insertion_moment(empty, cross, steps, 1)
    distinct = insertion_moment(empty, singleton, steps, 2)
    for i in range(2):
        for j in range(2):
            pair_region[i][j] /= D(steps)
            collision[i][j] /= D(steps * (2 * steps - 1))
            # Probability of distinct batches, divided by C(steps,2).
            distinct[i][j] *= D(2 * steps - 2)
            distinct[i][j] /= D((2 * steps - 1) * comb(steps, 2))
    split_region = [[collision[i][j] + distinct[i][j] for j in range(2)] for i in range(2)]
    return pair_region, split_region


def evaluate(precision, vector, z, tau, contraction, *, fresh_prefix=False):
    singleton = active(vector, z, tau, (1, 2, 4, 8))
    paired = active(vector, z, tau, (3, 12))
    cross = active(vector, z, tau, (5, 6, 9, 10))
    if fresh_prefix:
        fresh_singleton = active(vector, z, tau, (1, 2, 4, 8), fresh=True)
        fresh_paired = active(vector, z, tau, (3, 12), fresh=True)
        fresh_cross = active(vector, z, tau, (5, 6, 9, 10), fresh=True)
    minimum = min(vector[r][h] for h in range(1, 16) for r in range(1, h.bit_count() + 1))
    with localcontext() as context:
        context.prec = precision
        context.rounding = ROUND_CEILING
        convert = lambda matrix: [[upward(x) for x in row] for row in matrix]
        pair_region, split_region = regional_matrices(contraction, convert(singleton), convert(paired), convert(cross))
        if fresh_prefix:
            fresh_pair, fresh_split = regional_matrices(contraction, convert(fresh_singleton), convert(fresh_paired), convert(fresh_cross))
            return evaluate_fresh_columns(pair_region, split_region, fresh_pair, fresh_split,
                                          minimum, z, tau, contraction)
        empty_column = [[D(1), D(0)], [D(0), upward(contraction**8192)]]
        moments = []
        for a, b in ((pair_region, split_region), (split_region, pair_region)):
            occupied = multiply(a, b)
            occupied = multiply(occupied, occupied)
            rows = [[D(0), D(0)] for _ in range(65)]
            rows[0][0] = D(1)
            for column in range(64):
                following = [[D(0), D(0)] for _ in rows]
                for h in range(column + 2):
                    for j in range(2):
                        following[h][j] = sum((rows[h][k] * empty_column[k][j] for k in range(2)), D(0))
                        if h:
                            following[h][j] += sum((rows[h - 1][k] * occupied[k][j] for k in range(2)), D(0))
                rows = following
            moments.append([row[0] + row[1] * upward(1 / minimum) for row in rows])
        common = D(4096) * power(upward(1 / z), (7 * (1 << 21)) // 20)
        # Two shortened groups, each dimension4(h-32); 4096 adjacent pairs per parity.
        return {h: common * upward(tau**(8 * (h - 32) - 1)) * (moments[0][h] + moments[1][h])
                for h in range(33, 65)}


def evaluate_fresh_columns(pair_region, split_region, fresh_pair, fresh_split, minimum, z, tau, contraction,
                           late_prefix=None):
    """First h-32 OCCUPIED columns are fresh; zero columns do not advance it."""
    empty = [[D(1), D(0)], [D(0), upward(contraction**8192)]]
    parities = []
    for a, b, fa, fb in ((pair_region, split_region, fresh_pair, fresh_split),
                         (split_region, pair_region, fresh_split, fresh_pair)):
        occupied = multiply(a, b)
        occupied = multiply(occupied, occupied)
        if late_prefix is not None:
            occupied = multiply(late_prefix, occupied)
        new_occupied = multiply(fa, fb)
        new_occupied = multiply(new_occupied, new_occupied)
        parities.append((occupied, new_occupied))
    common = D(4096) * power(upward(1 / z), (7 * (1 << 21)) // 20)
    bounds = {}
    for target in range(33, 65):
        total = D(0)
        for occupied, new_occupied in parities:
            rows = [[D(0), D(0)] for _ in range(target + 1)]
            rows[0][0] = D(1)
            for column in range(64):
                following = [[D(0), D(0)] for _ in rows]
                for h in range(min(target, column + 1) + 1):
                    for j in range(2):
                        following[h][j] = sum((rows[h][k] * empty[k][j] for k in range(2)), D(0))
                        if h:
                            matrix = new_occupied if h - 1 < target - 32 else occupied
                            following[h][j] += sum((rows[h - 1][k] * matrix[k][j] for k in range(2)), D(0))
                rows = following
            total += rows[target][0] + rows[target][1] * upward(1 / minimum)
        bounds[target] = common * upward(tau**(8 * (target - 32) - 1)) * total
    return bounds


def main():
    vector, z, tau, contraction = check(F(999, 1000), F(4981, 5000), tau=F(65536))
    for precision in (80, 120):
        bounds = evaluate(precision, vector, z, tau, contraction, fresh_prefix=True)
        with localcontext() as context:
            context.prec = precision
            context.rounding = ROUND_CEILING
            total = sum(bounds.values(), D(0))
            print(f"UPWARD precision={precision}: two adjacent groups, same full columns <= {total}")
            print(f"h33={bounds[33]:.12E}; h64={bounds[64]:.12E}; worst h={max(bounds,key=bounds.get)}")
            if total * D(1 << 40) >= 1:
                print("OPEN: this envelope does not close40 bits; not a counterexample.")
    print("Scope excludes different column supports and partial columns in either group.")


if __name__ == "__main__":
    main()
