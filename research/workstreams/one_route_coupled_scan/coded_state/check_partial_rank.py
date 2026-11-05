"""Exhaustive toy verification of two-tap partial-support rank lemmas.

This checks local linear algebra, not an all-message distance certificate.
Run serially with check_maps.py available in the same directory.
"""
import itertools
import time
from check_maps import Field, Maps


def subsets(values):
    for r in range(len(values) + 1):
        yield from itertools.combinations(values, r)


def rank_or_zero(f, matrix):
    return f.rank(matrix) if matrix else 0


def main():
    start = time.monotonic()
    f, s = Field(4, 0x13), 4
    maps = Maps(f, s)
    a2 = [maps.a2([int(i == j) for i in range(s)]) for j in range(s)]
    cases = strict = max_gap = 0
    for active in subsets(tuple(range(s))):
        r = len(active)
        if not r:
            continue
        matrix = [[maps.a_matrix[i][j] for j in active] + [a2[j][i] for j in active]
                  for i in range(2 * s)]
        for zeros in subsets(tuple(range(2 * s))):
            z_left = sum(i < s for i in zeros)
            z_right = len(zeros) - z_left
            left_rank, right_rank = min(r, z_left), min(r, z_right)
            output_rows = [matrix[i] for i in zeros]
            assert rank_or_zero(f, output_rows) == left_rank + right_rank
            for state_zeros in subsets(active):
                # B[A,A2]=[I,I]. A state-coordinate equation has this row.
                rows = output_rows + [[int(j == i) for j in active] * 2
                                      for i in state_zeros]
                actual = rank_or_zero(f, rows)
                lower = (left_rank + right_rank +
                         min(len(state_zeros), r - min(left_rank, right_rank)))
                assert actual >= lower, (active, zeros, state_zeros, actual, lower)
                cases += 1
                strict += actual > lower
                max_gap = max(max_gap, actual - lower)
    # Full-state feedback law in a smaller model, enumerated independently.
    tiny = Field(2, 0x7)
    tm = Maps(tiny, 1)
    q = 2
    seen = {}
    for d, e in itertools.product(range(4), repeat=2):
        if d ^ e == 1:
            continue
        left = tm.a([tiny.table[d][q]])
        right = tm.a2([tiny.table[e][q]])
        y = tuple(a ^ b for a, b in zip(left, right))
        seen[y] = seen.get(y, 0) + 1
    expected = {y for y in itertools.product(range(4), repeat=2) if tm.bmap(list(y))[0] != q}
    assert set(seen) == expected and set(seen.values()) == {1}
    separated = {}
    for h in range(1, 4):
        for e in range(4):
            # d=1+h+e: the state multiplier h is uniform nonzero,
            # while the kernel multiplier e is independently uniform.
            aq = tm.a([tiny.table[1 ^ h][q]])
            eq = tiny.table[e][q]
            kq = [a ^ b for a, b in zip(tm.a([eq]), tm.a2([eq]))]
            y = tuple(a ^ b for a, b in zip(aq, kq))
            separated[y] = separated.get(y, 0) + 1
            assert q ^ tm.bmap(list(y))[0] == tiny.table[h][q]
    assert separated == seen
    print(f"Partial-support ranks: {cases} cases PASS; bound strict in {strict}; maximum gap {max_gap}")
    print(f"Full-state conditional-uniform law: {len(seen)} distinct permitted outputs PASS")
    print("Independent nonzero state multiplier / uniform kernel multiplier law PASS")
    print(f"Elapsed {time.monotonic()-start:.3f}s; no whole-code certificate")


if __name__ == '__main__':
    main()
