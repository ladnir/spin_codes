"""Integer coverage gate for sparse outer supports in post-diagonal SPIN.

See SPARSE_PROOF.md. No floating point, generated output files, or benchmark.
Failed cells mean this envelope is inconclusive, not that the code is bad.
"""
import argparse

K = 1 << 20
N = 2 * K
OUTER_BLOCKS = N // 256
BATCH = 128
STATE = 64
REGIONAL_BATCHES = OUTER_BLOCKS // BATCH
LOG_Q = 128
LOG_BATCHES = 14
LOG_COORDS = 21
TERM_BITS = 64


def ceildiv(a, b):
    return (a + b - 1) // b


def cell(p, h):
    assert 1 <= p <= 64 and 129 <= h <= 256
    dimension = p * (h - 128)
    # Choose source blocks, a containing h-region set, and a projective word.
    log_count = 13*p + 257 + LOG_Q*(dimension-1)
    # At most r resets outside an event of union probability <=2^-TERM_BITS.
    r = ceildiv(log_count + TERM_BITS, STATE*LOG_Q - LOG_BATCHES) - 1
    # Each of the two total-defect events has the same union bound.
    m = ceildiv(log_count + TERM_BITS, LOG_Q - LOG_COORDS)
    assert log_count - (STATE*LOG_Q - LOG_BATCHES)*(r+1) <= -TERM_BITS
    assert log_count - (LOG_Q - LOG_COORDS)*m <= -TERM_BITS
    if r == 0:
        live = REGIONAL_BATCHES*(h-1)
    else:
        live = max(0, REGIONAL_BATCHES*(h-2*r-1) - (r+1))
    activation_holes = (BATCH-(BATCH-STATE+1))*(r+1)
    output_lower = BATCH*live - activation_holes - 2*(m-1)
    return dict(p=p, h=h, dimension=dimension, r=r, m=m, live=live,
                output_lower=output_lower, log_count=log_count)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--numerator', type=int, default=49)
    parser.add_argument('--denominator', type=int, default=100)
    args = parser.parse_args()
    assert 0 < args.numerator < args.denominator
    threshold = N*args.numerator // args.denominator
    prefix = 0
    worst_prefix = None
    first_failure = None
    for p in range(1, 65):
        rows = [cell(p, h) for h in range(129, 257)]
        worst = min(rows, key=lambda c: c['output_lower'])
        if worst['output_lower'] > threshold and first_failure is None:
            prefix = p
            if worst_prefix is None or worst['output_lower'] < worst_prefix['output_lower']:
                worst_prefix = worst
        elif first_failure is None:
            first_failure = worst
    assert prefix > 0
    terms = 3 * prefix * 128
    union_bits = TERM_BITS - (terms-1).bit_length()
    assert union_bits >= 40
    print(f"Target >{args.numerator}/{args.denominator}, W={threshold}")
    print(f"Complete coverage p=1..{prefix}, h=129..256: {prefix*128} cells, "
          f"failure <=2^-{union_bits}")
    print(f"Worst certified cell: {worst_prefix}; weight slack={worst_prefix['output_lower']-threshold}")
    print(f"First unsupported p, worst cell: {first_failure}")
    print("All larger source supports remain outside this certificate.")


if __name__ == '__main__':
    main()
