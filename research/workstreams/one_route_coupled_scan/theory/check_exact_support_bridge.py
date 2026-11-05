"""Small-field exhaustive support checks and exact bridge constant; no files.

Enumerates systematic GF4 [n,k] codes for (n,k)=(1,1),(2,1),(3,2),(4,2),
including every parity matrix. Checks their scalar extensions to GF16.
This tests the algebraic bridge, not the RS/Toeplitz probability estimate.
"""

from collections import Counter
from fractions import Fraction as F
from itertools import product

from test_local_algebra import MUL, require
from test_rank_transfer import zero_span, reduce_form


E = (0, 1, 6, 7)  # The unique GF4 subfield of this GF16 representation.


def words(columns, alphabet, k):
    counts = Counter()
    for values in product(alphabet, repeat=k):
        mask = 0
        for index, column in enumerate(columns):
            value = MUL[values[0]][column & 15]
            if k == 2:
                value ^= MUL[values[1]][column >> 4]
            if value:
                mask |= 1 << index
        counts[mask] += 1
    return counts


def check_small_codes():
    require(all(MUL[x][MUL[x][MUL[x][x]]] == x for x in E), "GF4 subfield")
    codes = supports = 0
    for n, k in ((1, 1), (2, 1), (3, 2), (4, 2)):
        identity = (1,) if k == 1 else (1, 16)
        for parity in product(E, repeat=k * (n - k)):
            columns = identity + tuple(
                parity[k * j] | ((parity[k * j + 1] << 4) if k == 2 else 0)
                for j in range(n - k)
            )
            base = words(columns, E, k)
            extension = words(columns, range(16), k)
            require(set(base) == set(extension), ("support preservation", columns))
            for support, count in base.items():
                if not support:
                    continue
                # Inactive coordinate forms generate the shortening constraints.
                basis = zero_span([c for j, c in enumerate(columns) if not (support >> j) & 1])
                dimension = k - sum(bool(row) for row in basis)
                size = support.bit_count()
                require(dimension >= 1, ("nonzero shortened space", columns, support))
                require(all(reduce_form(c, basis) for j, c in enumerate(columns) if (support >> j) & 1),
                        ("realizable coordinate forms", columns, support))
                if dimension == 1:
                    require(count == 3, ("one-dimensional exact support", columns, support))
                else:
                    require(count >= (5 - size) * 3 * 4**(dimension - 2), "hyperplane lower bound")
                require(count >= 3**(dimension - 1) * (3 - size + dimension), "grid lower bound")
                for tau in (1, 2, 4):
                    require(F(tau**(dimension - 1)) <= F(4, 5 - size) * F(count, 3),
                            ("hyperplane bridge", tau, columns, support))
                supports += 1
            codes += 1
    # The length<=q condition is necessary: a binary parity code supplies a
    # negative control, with support111 present only after field extension.
    parity = (1, 16, 17)
    require(7 not in words(parity, (0, 1), 2), "binary parity base control")
    require(7 in words(parity, E, 2), "binary parity extension control")
    return codes, supports


def check_actual_constant():
    bound = F(255, 127) * F(256, 255)**127
    require(bound < F(331, 100), "uniform exact-support factor")
    # Under the same point cap and uniform factor, containment is tighter.
    require(F(256, 255)**256 < bound, "comparison with containment envelope")
    return bound


if __name__ == "__main__":
    codes, supports = check_small_codes()
    check_actual_constant()
    print(f"PASS: {codes} GF4 codes, {supports} nonzero exact supports and GF16 extensions.")
    print("PASS: exact C_* < 3.31; out-of-scope binary counterexample retained.")
