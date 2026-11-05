"""Keep quiet/rank1..4 at regional boundaries, with full tables inside.

Uses the online route proof and the actual two-adjacent/full-column class.
The source dimension, paired collisions, invertibility conditioning, and
free terminal remain unchanged. No files, jobs, or caches at import.
Float mode is diagnostic; directed mode uses positive Decimal arithmetic.
"""

from decimal import Decimal as D, ROUND_CEILING, localcontext
from fractions import Fraction as F

from coupled_rank_transfer import pair_antipodal_cycle
from empty_span_rank import empty_step
from online_route import VALID, region, zeros


def basis(vector, rank):
    result = zeros(vector[0][0])
    if rank == 0:
        result[0][0] = vector[0][0] * 0 + 1
    else:
        for r, h in VALID[1:]:
            if r == rank:
                result[r][h] = vector[r][h]
    return result


def collapse(vector, table):
    return [table[0][0]] + [max(table[r][h] / vector[r][h]
                              for rr, h in VALID if rr == r) for r in range(1, 5)]


def pivot_step(table, tau):
    """One first-nonzero source partition, preserving the residual rank."""
    out = zeros(tau)
    one = tau * 0 + 1
    weights = [one / tau**i for i in range(4)]
    out[0][0] = sum(weights) * table[0][0]
    for r, h in VALID[1:]:
        out[r][h] = sum(weight * max(table[s][h] for s in range(max(1, r - i), r + 1))
                        for i, weight in enumerate(weights))
    return out


def empty_interval(table, z, tau, steps):
    if not 0 <= steps <= 8192:
        raise ValueError("bounded empty interval")
    for _ in range(steps):
        table = empty_step(table, 4, z, tau, cycle=pair_antipodal_cycle(4))
    return table


def envelopes(vector, z, tau, steps=2048, progress=None, *, rank_floor=False):
    """Nonnegative 5x5 matrices, indexed by exact input/output residual rank."""
    zero = z * 0
    matrices = {name: [[zero for _ in range(5)] for _ in range(5)]
                for name in ("empty", "pivot", "late_pair", "late_split", "fresh_pair", "fresh_split")}
    for b in range(5):
        terminal = basis(vector, b)
        columns = {
            "empty": collapse(vector, empty_interval(terminal, z, tau, 4 * steps)),
            "pivot": collapse(vector, pivot_step(terminal, tau)),
        }
        for kind, mode in (("late", dict(outer_randomized=True)), ("fresh", dict(fresh=True))):
            pair, split = region(terminal, z, tau, steps=steps, rank_floor=rank_floor, **mode)
            columns[kind + "_pair"] = collapse(vector, pair)
            columns[kind + "_split"] = collapse(vector, split)
        for name, column in columns.items():
            for a, value in enumerate(column):
                matrices[name][a][b] = value
        if progress:
            progress(b)
    return matrices


def multiply(a, b):
    zero = a[0][0] * 0
    return [[sum((a[i][k] * b[k][j] for k in range(5)), zero) for j in range(5)] for i in range(5)]


def row_multiply(a, b):
    zero = a[0] * 0
    return [sum((a[k] * b[k][j] for k in range(5)), zero) for j in range(5)]


def column_matrices(matrices):
    result = []
    p2 = multiply(matrices["pivot"], matrices["pivot"])
    for first, second in (("pair", "split"), ("split", "pair")):
        late = multiply(matrices["late_" + first], matrices["late_" + second])
        late = multiply(p2, multiply(late, late))
        fresh = multiply(matrices["fresh_" + first], matrices["fresh_" + second])
        fresh = multiply(fresh, fresh)
        result.append((late, fresh))
    return result


def moment(empty, late, fresh, terminal, occupied, columns=64, fresh_count=None):
    """Ordered column choices, with fresh_count initial OCCUPIED columns."""
    if fresh_count is None:
        fresh_count = occupied - 32
    zero = terminal[0] * 0
    rows = [[zero] * 5 for _ in range(occupied + 1)]
    rows[0][0] = zero + 1
    for col in range(columns):
        following = [[zero] * 5 for _ in rows]
        for h in range(min(occupied, col + 1) + 1):
            following[h] = row_multiply(rows[h], empty)
            if h:
                added = row_multiply(rows[h - 1], fresh if h - 1 < fresh_count else late)
                following[h] = [a + b for a, b in zip(following[h], added)]
        rows = following
    return sum((a * b for a, b in zip(rows[occupied], terminal)), zero)


def evaluate(matrices, vector, z, tau):
    """Caller sets upward Decimal context; no float matrix becomes certified."""
    from check_sparse_full_gate import power, upward
    terminal = [D(1)] + [upward(1 / min(vector[r][h] for rr, h in VALID if rr == r))
                         for r in range(1, 5)]
    parities = column_matrices(matrices)
    common = D(4096) * power(upward(1 / z), (7 * (1 << 21)) // 20) * upward(F(256, 255)**64)
    return {h: common * upward(tau**(8 * (h - 32) - 1))
            * sum((moment(matrices["empty"], late, fresh, terminal, h) for late, fresh in parities), D(0))
            for h in range(33, 65)}


def main():
    import argparse
    import math
    from check_empty_span_potential import check
    from check_sparse_full_gate import upward

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--precision", type=int, default=0, choices=(0, 80, 120))
    parser.add_argument("--rank-floor", action="store_true")
    args = parser.parse_args()
    vector, z, tau, _ = check(F(499, 500), F(993, 1000), tau=F(1024))
    with localcontext() as context:
        context.prec = args.precision or 80
        context.rounding = ROUND_CEILING
        convert = upward if args.precision else float
        work_vector = [[convert(x) for x in row] for row in vector]
        matrices = envelopes(work_vector, convert(z), convert(tau),
                             progress=lambda b: print(f"completed boundary rank {b}/4", flush=True),
                             rank_floor=args.rank_floor)
        if not args.precision:
            matrices = {name: [[D(str(x)) for x in row] for row in matrix] for name, matrix in matrices.items()}
        bounds = evaluate(matrices, vector, z, tau)
        print("UPWARD" if args.precision else "FLOAT DIAGNOSTIC", "precision", args.precision,
              "rank_floor", args.rank_floor, "log2 union", math.log2(sum(bounds.values(), D(0))))
        print("selected log2 bounds", {h: math.log2(bounds[h]) for h in (33, 48, 64)})
        print("h64 late/fresh column diagonal log2", [[math.log2(m[r][r]) if m[r][r] else None
              for r in range(5)] for m in column_matrices(matrices)[0]])
        print("Only shared-full-column adjacent pairs; no whole-code certificate.")


if __name__ == "__main__":
    main()
