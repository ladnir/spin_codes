"""Reverse the route's materialization: contiguous inner stores, outer gathers.

Every packet is still read once at its original outer position. The caller
provides inverse-route block offsets, relative to the full contiguous scratch
array. This is an exact-map experiment, not a new permutation distribution.
"""
import argparse
from pathlib import Path
import subprocess
import sys

import outer_layout_codegen
from packed_coeff_codegen import definition


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('header', type=Path)
    parser.add_argument('--prefetch', type=int, choices=(0, 1, 2, 4, 8), default=0,
                        help='prefetch this many outer groups ahead, within each tile')
    args = parser.parse_args()
    result = subprocess.run([sys.executable, str(Path(outer_layout_codegen.__file__).resolve()),
                             str(args.header), '--mode', '2', '--tile-mode', '5'],
                            capture_output=True, text=True, check=True)
    source = result.stdout
    old = definition(source, 'SPIN_NOINLINE void bchPackedCoeffCompact(')
    new = old.replace('bchPackedCoeffCompact(', 'bchGatherCoeffCompact(').replace(
        'const std::uint64_t* __restrict coeff)',
        'const std::uint64_t* __restrict coeff,const std::uint32_t* __restrict packetOffsets)')
    assert new != old and 'packetOffsets)' in new
    for j in range(8):
        original = f'a+4*(8*group+{j})'
        assert new.count(original) == 1
        new = new.replace(original, f'a+packetOffsets[8*group+{j}]')
    if args.prefetch:
        prefix = 'for(unsigned group=0;group<32;++group) {'
        assert new.count(prefix) == 1
        hint = f'if(group+{args.prefetch}<32) {{\n' + '\n'.join(
            f'__builtin_prefetch(a+packetOffsets[8*(group+{args.prefetch})+{j}],0,3);'
            for j in range(8)) + '\n}'
        new = new.replace(prefix, prefix + '\n' + hint)
    # Only keep the new public export. Helpers have internal linkage so the
    # object can coexist with the unchanged control in a single executable.
    for storage in ('Original', 'Aligned'):
        source = source.replace(definition(source, f'SPIN_NOINLINE void bchPackedCoeff{storage}('), '')
    source = source.replace(old, new)
    for name in ('outerLayoutMode', 'packedCoeffTileMode'):
        source = source.replace(f'unsigned {name}()', f'static unsigned gather_{name}()')
    print(source)
    print(result.stderr, file=sys.stderr, end='')
    print('Gather route: eight indexed full-packet loads replace eight direct loads; '
          'all packed arithmetic unchanged.', file=sys.stderr)


if __name__ == '__main__':
    main()
