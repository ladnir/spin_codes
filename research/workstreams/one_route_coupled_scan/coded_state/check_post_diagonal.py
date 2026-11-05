"""Exact checks for Y=X+D*A*Q, Q'=Q+B*Y. Not a full certificate."""
import itertools
import random
import time
from check_maps import Field, Maps, weight


def main():
    start = time.monotonic()
    f, s = Field(4, 0x13), 4
    maps = Maps(f, s)
    masks = set()
    least = 2 * s
    for state in itertools.product(range(f.size), repeat=s):
        if not any(state):
            continue
        expanded = maps.a(list(state))
        support = tuple(i for i, x in enumerate(expanded) if x)
        least = min(least, len(support))
        masks.add(support)
    assert least == s + 1
    for support in masks:
        assert f.rank([[row[j] for j in support] for row in maps.b_matrix]) == s
    print(f"All {f.size**s-1} nonzero toy states: min expanded weight={least}; "
          f"all {len(masks)} resulting compression restrictions have full rank")

    # Enumerate the exact refresh distribution in the smallest nontrivial model.
    tiny = Field(2, 0x7)
    tm = Maps(tiny, 1)
    fixtures = 0
    for q in range(1, 4):
        aq = tm.a([q])
        for x in itertools.product(range(4), repeat=2):
            counts = [0] * 4
            for d in itertools.product(range(4), repeat=2):
                y = [xx ^ tiny.table[dd][a] for xx, dd, a in zip(x, d, aq)]
                next_q = q ^ tm.bmap(y)[0]
                counts[next_q] += 1
                recovered = [yy ^ tiny.table[dd][a] for yy, dd, a in zip(y, d, aq)]
                assert recovered == list(x)
            assert counts == [4] * 4
            fixtures += 1
    print(f"Exact uniform next-state / inverse check: {fixtures} GF4 input-state fixtures PASS")

    f, s = Field(8, 0x11b), 64
    maps = Maps(f, s)
    # Construct an expanded word with the minimum weight, using 63 explicit roots in Y.
    roots = f.points[:s-1]
    def evaluate(x):
        result = 1
        for root in roots:
            result = f.table[result][x ^ root]
        return result
    state = [evaluate(x) for x in f.points[2*s:3*s]]
    aq = maps.a(state)
    assert aq == [evaluate(x) for x in f.points[:2*s]]
    assert weight(aq) == 65
    support = [i for i, x in enumerate(aq) if x]
    assert f.rank([[row[j] for j in support] for row in maps.b_matrix]) == 64
    rng = random.Random(20261021)
    state = [0] * s
    for _ in range(32):
        aq = maps.a(state)
        x = [rng.randrange(256) for _ in range(2*s)]
        d = [rng.randrange(256) for _ in range(2*s)]
        y = [xx ^ f.table[dd][a] for xx, dd, a in zip(x, d, aq)]
        assert x == [yy ^ f.table[dd][a] for yy, dd, a in zip(y, d, aq)]
        state = [q ^ by for q, by in zip(state, maps.bmap(y))]
    print("GF256 minimum-weight expanded state, full refresh rank, and 32-batch inverse PASS")
    # Diagnose a structural zero pattern in the original projection A*B.
    # Row scaling B preserves its MDS property but can remove this alignment.
    map_rng = random.Random(2026102101)
    coefficients = [map_rng.randrange(1, 256) for _ in range(s)]
    plain_zeros = scaled_zeros = 0
    for col in range(2*s):
        column = [row[col] for row in maps.b_matrix]
        plain = maps.a(column)
        scaled = maps.a([f.table[r][x] for r, x in zip(coefficients, column)])
        plain_zeros += sum(x == 0 for x in plain)
        scaled_zeros += sum(x == 0 for x in scaled)
    assert plain_zeros == 2*s*(s-1)
    print(f"Structural diagnostic A*B zeros={plain_zeros}, A*diag(r)*B zeros={scaled_zeros} "
          "(fixed test coefficients; NOT a rank certificate)")
    k, n = 1 << 20, 1 << 21
    w, holes = 49*n//100, 4587
    margin = -(n + 128 * (w-k+holes))
    assert margin == 128
    print(f"Exact bounded-hole class: H<={holes}, delta=49/100, failure <=2^-{margin}")
    batches = n // 128
    sparse_margin = 128 * ((65*(batches-1)+1)//2 - 31*batches) - n
    assert sparse_margin == 1044480
    print(f"Exact full-scan live / <=31 output symbols per batch class: <=2^-{sparse_margin}")
    print(f"Elapsed {time.monotonic()-start:.3f}s; no whole-code certificate")


if __name__ == '__main__':
    main()
