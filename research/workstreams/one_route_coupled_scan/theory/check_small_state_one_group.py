"""Integer-only checks for SMALL_STATE_PROOF_GATE.md; no files or dependencies.

This checks finite arithmetic and tiny interval geometry, not the general
field/rank lemmas or a whole-code certificate. The theorems concern only
messages confined to one or two retained historical outer groups.
"""

from math import comb


def intervals(mask, length):
    """Return the maximal nonempty intervals of a small binary mask."""
    result = []
    start = None
    for index in range(length + 1):
        bit = index < length and ((mask >> index) & 1)
        if bit and start is None:
            start = index
        elif not bit and start is not None:
            result.append((start, index))
            start = None
    return result


def check_tiny_geometry():
    fixtures = 0
    for column_length in range(1, 5):
        for length in range(1, 11):
            for mask in range(1, 1 << length):
                pieces = intervals(mask, length)
                selected = mask.bit_count()
                touched = len({i // column_length for i in range(length)
                               if (mask >> i) & 1})
                assert column_length * touched <= (
                    selected + 2 * column_length * len(pieces)
                )
                # An interval representation may split adjacent pieces;
                # increasing c only weakens the required inequality.
                for start, end in pieces:
                    count = ((end - 1) // column_length
                             - start // column_length + 1)
                    assert column_length * count <= (
                        end - start + 2 * column_length
                    )
                fixtures += 1
    return fixtures


def main():
    q = 1 << 32
    groups = 1 << 13
    columns = 64
    batch = 64
    state = 16
    total = 1 << 21
    steps = total // batch
    column_batches = steps // columns
    threshold = (7 * total) // 20
    assert steps == 1 << 15
    assert column_batches == 512
    assert threshold == 734003

    # Bernoulli: (1-1/q)^m >= 1-m/q, hence the reciprocal is
    # <= q/(q-m) < 2. No enormous rational powers are constructed.
    assert 2 * state * steps < q
    assert q < 2 * (q - state * steps)
    results = []
    expected = {
        1: (14846, 10, 20, -412, -742042, 411),
        2: (14318, 18, 24, -407, -501389, 406),
    }
    class_count = 0
    for active_groups in (1, 2):
        rows = []
        for occupied in range(33, 65):
            dimension = 4 * active_groups * (occupied - 32)
            cutoff = dimension // state + 2
            minimum_live = (column_batches * (occupied - 2 * cutoff)
                            - cutoff)
            assert dimension <= 128 * active_groups
            assert 0 < minimum_live < steps
            assert batch * minimum_live > threshold
            rows.append((occupied, dimension, cutoff, minimum_live))
        class_count += len(rows)
        minimum_live = min(row[3] for row in rows)
        maximum_returns = max(row[2] for row in rows)
        minimum_charge = min(state * row[2] - row[1] for row in rows)
        assert sum(comb(columns, row[0]) for row in rows) <= 1 << columns
        assert 2 * state * maximum_returns < q
        group_bits = active_groups * (groups.bit_length() - 1)
        assert comb(groups, active_groups) <= 1 << group_bits
        step_bits = steps.bit_length() - 1
        high_exponent = (group_bits + columns + step_bits * maximum_returns
                         - 32 * minimum_charge + 1)

        # C(64a,64a-e) <= 2^(64a), while the rank bound costs
        # q^(-ceil((64a-e)/4)) <= 2^(-8*(64a-e)).
        live_exponent = 8 * threshold - 7 * batch * minimum_live
        assert -7 * batch < 0  # Decreasing as a increases.
        low_exponent = (live_exponent + 32 * 128 * active_groups + steps
                        + columns + group_bits + 1)
        assert low_exponent < high_exponent
        margin = -high_exponent - 1
        observed = (minimum_live, maximum_returns, minimum_charge,
                    high_exponent, low_exponent, margin)
        assert observed == expected[active_groups]
        results.append((active_groups, high_exponent, low_exponent, margin))
    combined_margin = min(row[3] for row in results) - 1
    assert combined_margin == 405

    # Record the explicit boundary of this route-independent argument.
    three_group_dimension = 3 * 4 * (64 - 32)
    three_group_cutoff = three_group_dimension // state + 2
    three_group_live = (column_batches * (64 - 2 * three_group_cutoff)
                        - three_group_cutoff)
    assert (three_group_dimension, three_group_cutoff, three_group_live) == (
        384, 26, 6118
    )
    assert batch * three_group_live < threshold

    geometry = check_tiny_geometry()
    print(f"PASS: {class_count} column classes; {geometry} tiny interval fixtures")
    for count, high, low, margin in results:
        print(f"{count} group(s): high < 2^({high}), low < 2^({low}), "
              f"union < 2^(-{margin})")
    print(f"Combined one-/two-group union < 2^(-{combined_margin})")
    print("Scope: at most two historical outer groups; not a whole-code certificate.")


if __name__ == "__main__":
    main()
