"""Bounded exact outer-CDF comparison across canonical block widths.

Regenerate one authenticated refined dimension table, then reuse it for every
width. No inner computation, benchmarks, or whole-code certificates.
"""
from canonical_counts import authenticated_bch_dimensions, canonical_rank_caps, total_caps, transport_cdf, display_log
from local_models import full_block, mean


def compare(dimensions, widths=(1, 2, 4, 8)):
    result = {}
    for width in widths:
        caps = total_caps(canonical_rank_caps(dimensions, width))
        result[width] = dict(canonical_cdf=caps,
                            output_cdf=transport_cdf(caps, full_block(width)),
                            local_mean=mean(full_block(width)))
    return result


def main():
    dimensions = authenticated_bch_dimensions(104, True)
    rows = compare(dimensions)
    print('Same dimension bounds only: width1 is NOT the stronger existing BCH-spectrum certificate.')
    print('width,minimum_possible_support,mean_per_nonempty_block')
    for width, row in rows.items():
        first = next(u for u, bound in enumerate(row['output_cdf']) if bound)
        print(width, first, float(row['local_mean']), sep=',')
    print('support,' + ','.join(f'width{w}_log2_expected_CDF' for w in rows))
    for u in (5, 10, 19, 32, 38, 40, 48, 64, 80, 96, 104, 112, 115, 120, 128, 160, 192, 224, 240, 256):
        print(u, *(display_log(row['output_cdf'][u]) for row in rows.values()), sep=',')
    print('The support moments/CDF bounds are exact; logarithms are display only. No distance margin.')


if __name__ == '__main__':
    main()
