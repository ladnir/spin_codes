"""Independent small-timeline and full integer-coverage checks."""
from sparse_certificate import cell, N, TERM_BITS


def timeline_check(regions, regional_batches):
    n = regions * regional_batches
    cases = 0
    for support in range(1, 1 << n):
        h = len({i // regional_batches for i in range(n) if support >> i & 1})
        # Merge trajectories only when all quantities used by the lemma agree.
        histories = {(False, 0, 0, 0)}  # live state, resets, live batches, activations
        for i in range(n):
            active_input = bool(support >> i & 1)
            next_histories = set()
            for live, resets, length, activations in histories:
                if live:
                    next_histories.add((True, resets, length+1, activations))
                    next_histories.add((False, resets+1, length+1, activations))
                else:
                    next_histories.add((active_input, resets, length,
                                        activations+int(active_input)))
            histories = next_histories
        for live, resets, length, activations in histories:
            assert activations <= resets+1
            lower = (regional_batches*(h-1) if resets == 0 else
                     regional_batches*(h-2*resets-1)-(resets+1))
            assert length >= lower, (support, live, resets, length, activations, h, lower)
            cases += 1
    return cases


def main():
    cases = timeline_check(4, 2) + timeline_check(3, 3)
    print(f"Exhaustive small timeline summaries: {cases} PASS (free final state)")
    count = 0
    threshold = 491*N//1000
    least = N
    for p in range(1, 32):
        for h in range(129, 257):
            c = cell(p, h)
            # Recompute thresholds without the certificate's ceildiv helper.
            count_bound = 13*p+257+128*(p*(h-128)-1)
            assert count_bound == c['log_count']
            assert count_bound-8178*(c['r']+1) <= -TERM_BITS
            assert count_bound-107*c['m'] <= -TERM_BITS
            assert c['output_lower'] > threshold
            least = min(least, c['output_lower'])
            count += 1
    assert count == 3968 and least == 1039505
    assert 3*count < 1 << 14
    assert TERM_BITS-14 == 50
    failure = cell(32, 256)
    assert failure['output_lower'] == 1018159 < threshold
    print(f"All {count} p<=31 cells PASS; minimum weight lower bound={least}; union <2^-50")
    print("p=32,h=256 envelope failure reproduced; not a code counterexample")


if __name__ == '__main__':
    main()
