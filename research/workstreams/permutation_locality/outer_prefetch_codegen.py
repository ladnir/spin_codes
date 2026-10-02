"""Exact retained outer kernel with within-tile prefetch experiments.

All addresses are inside the declared input, output, or coefficient tile.
Only the compact entry point changes; arithmetic and external layouts do not.
"""
import argparse
from pathlib import Path
import subprocess
import sys

import packed_coeff_codegen as retained


def generate(header, mode):
    result = subprocess.run([sys.executable, str(Path(retained.__file__)),
        str(header), '--tile-mode', '1'], check=True, capture_output=True, text=True)
    source = result.stdout
    marker = 'SPIN_NOINLINE void bchPackedCoeffCompact('
    original = retained.definition(source, marker)
    modified = original
    loop = 'for(unsigned group=0;group<32;++group) {'
    if original.count(loop) != 1:
        raise ValueError('unique compact preparation loop required')
    prefetch = []
    if mode in (1, 2, 3):
        write = 0 if mode == 3 else 1
        prefetch += [f'__builtin_prefetch(out+{128*row}+4*group,{write},3);'
                     for row in range(4)]
    if mode == 2:
        prefetch += ['if(group+4<32) {']
        prefetch += [f'__builtin_prefetch(a+4*(8*(group+4)+{col}),0,3);' for col in range(8)]
        prefetch += ['__builtin_prefetch(coeff+16*(group+4),0,3);',
                     '__builtin_prefetch(coeff+16*(group+4)+8,0,3);', '}']
    if prefetch:
        modified = modified.replace(loop, loop+'\n'+'\n'.join(prefetch))
    for group in range(32):
        for row in range(4):
            assert 0 <= 128*row+4*group < 512
        if group+4 < 32:
            assert 4*(8*(group+4)+7) < 1024
            assert 16*(group+4)+8 < 512
    source = source.replace(original, modified)
    source = source.replace('unsigned packedCoeffTileMode() {return 1;}',
                            f'unsigned packedCoeffTileMode() {{return {300+mode};}}')
    return source


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('header', type=Path)
    parser.add_argument('--mode', type=int, choices=range(4), required=True)
    args = parser.parse_args()
    print(generate(args.header, args.mode))
    print(f'Prefetch mode {args.mode}: exact retained arithmetic; all prefetch addresses in-tile', file=sys.stderr)


if __name__ == '__main__':
    main()
