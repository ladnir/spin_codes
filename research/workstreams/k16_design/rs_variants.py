"""Screen two exact RS outer ensembles with identical K16 packet geometry.

The GL32/GF256 [8,4] and GL16/GF16 [16,8] variants both encode 128 message
bits to 64 four-bit packets per group. This driver changes only their exact
expected support counts; the retained physical inner and routing diagnostic
are unchanged. An occupancy-one result is not a whole-code certificate.

The GF16 parity check is an algebraic circuit check, not a runtime prototype.
At evaluation points 0..15, its parity block is XOR convolution on 0..7.
In the characteristic-two group algebra, writing y_i = g_i + 1 gives y_i^2=0.
The subset-zeta change of basis therefore turns the fixed convolution into
the disjoint-subset product used below. Both changes of basis use XORs only.
"""

import argparse
from pathlib import Path

from packet_q1 import evaluate_q1
from rs_outer import expected_group_support_counts
from rs16_maps import check_parity


def run(variant, tilts, precision=192, output=None):
    if variant == 'rs16-gf16':
        n, k, packets_per_symbol, map_bits = 16, 8, 4, 16
        parity_check = check_parity()
    elif variant == 'rs8-gf256':
        n, k, packets_per_symbol, map_bits = 8, 4, 8, 32
        parity_check = None
    else:
        raise ValueError('unknown RS ensemble')
    counts = expected_group_support_counts(n=n, k=k, packet_bits=4,
                                           packets_per_symbol=packets_per_symbol)
    return evaluate_q1(counts, group_count=512, regions=64,
                       epochs_per_region=16, group_dimension=128,
                       count_kind='shells', tilts=tilts, precision=precision,
                       output=output,
                       metadata=dict(outer=f'four parallel GF{1 << packets_per_symbol} RS[{n},{k}] rows',
                           label_mixing=f'independent uniform GL{map_bits} per symbol and per group',
                           exact_expected_shell_counts=True,
                           variant=variant, parity_check=parity_check))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--variant', choices=('rs16-gf16', 'rs8-gf256'),
                        default='rs16-gf16')
    parser.add_argument('--tilts', nargs='+',
                        default=['.00512', '.01024', '.0256', '.0512'])
    parser.add_argument('--precision', type=int, default=192)
    parser.add_argument('--output', type=Path)
    parser.add_argument('--parity-only', action='store_true')
    args = parser.parse_args()
    if args.parity_only:
        print(check_parity())
    else:
        run(args.variant, args.tilts, args.precision, args.output)


if __name__ == '__main__':
    main()
