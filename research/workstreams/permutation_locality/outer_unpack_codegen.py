"""Generate exact tile-5/layout-2 BCH with one fused VBMI output network.

The eight final bit transposes are unchanged.  The remaining byte transpose,
coordinate assembly, and row-major output transpose become three stages of
eight two-source byte permutes.  No setup or mathematical map changes.
"""
import argparse
from contextlib import redirect_stdout
from io import StringIO
from pathlib import Path
import subprocess
import sys

import outer_layout_codegen as layout
import packed_bch_tune_codegen as tune
from packed_coeff_codegen import definition


def unpack_bytes(a, b, width, high=False):
    result = []
    for lane in range(4):
        start = 16 * lane + (8 if high else 0)
        for byte in range(start, start + 8, width):
            result.extend(a[byte:byte + width])
            result.extend(b[byte:byte + width])
    return result


def shuffle_lanes(a, b, immediate):
    lanes_a = [a[16*i:16*i+16] for i in range(4)]
    lanes_b = [b[16*i:16*i+16] for i in range(4)]
    return sum(tune.shuffle128(lanes_a, lanes_b, immediate), [])


def retained_byte_outputs():
    """Literal symbolic copy of the retained post-GFNI byte operations."""
    x = [list(range(64*i, 64*i+64)) for i in range(8)]
    t = [unpack_bytes(x[i], x[i+1], 1, high)
         for i in (0, 2, 4, 6) for high in (False, True)]
    u = [unpack_bytes(t[i], t[j], 2, high)
         for i, j in ((0, 2), (1, 3), (4, 6), (5, 7))
         for high in (False, True)]
    x = [unpack_bytes(u[i], u[i+4], 4, high)
         for i in range(4) for high in (False, True)]
    y = [unpack_bytes(x[i], x[i+4], 1, high)
         for i in range(4) for high in (False, True)]
    outputs = []
    for first in (0, 4):
        a = shuffle_lanes(y[first], y[first+1], 0x44)
        b = shuffle_lanes(y[first], y[first+1], 0xee)
        c = shuffle_lanes(y[first+2], y[first+3], 0x44)
        d = shuffle_lanes(y[first+2], y[first+3], 0xee)
        outputs.extend((shuffle_lanes(a, c, 0x88),
                        shuffle_lanes(a, c, 0xdd),
                        shuffle_lanes(b, d, 0x88),
                        shuffle_lanes(b, d, 0xdd)))
    if sorted(sum(outputs, [])) != list(range(512)):
        raise ArithmeticError("retained byte operations are not a permutation")
    return outputs


def network():
    targets = retained_byte_outputs()
    for target in targets:
        if [sum(byte // 64 == plane for byte in target) for plane in range(8)] != [8]*8:
            raise ArithmeticError("unbalanced output cannot use the compact network")
    nodes = {f'x{i}': list(range(64*i, 64*i+64)) for i in range(8)}
    gates = []

    def gate(name, left, right, wanted):
        source = nodes[left] + nodes[right]
        if len(wanted) != 64 or len(set(source)) != 128:
            raise ArithmeticError("invalid network node")
        lookup = {value: i for i, value in enumerate(source)}
        index = tuple(lookup[value] for value in wanted)
        got = [source[i] for i in index]
        if got != wanted:
            raise ArithmeticError("network gate changed a byte")
        nodes[name] = got
        gates.append((name, left, right, index))

    # Four plane pairs, each pair split between the low/high four targets.
    for pair in range(4):
        for half in range(2):
            wanted = [byte for target in targets[4*half:4*half+4]
                      for byte in target if byte // 128 == pair]
            gate(f'a{2*pair+half}', f'x{2*pair}', f'x{2*pair+1}', wanted)
    # Two four-plane sets, each split into four adjacent target pairs.
    for half in range(2):
        for pair in range(4):
            wanted = [byte for target in targets[2*pair:2*pair+2]
                      for byte in target if byte // 256 == half]
            gate(f'b{4*half+pair}', f'a{4*half+pair//2}',
                 f'a{4*half+2+pair//2}', wanted)
    for target in range(8):
        gate(f'c{target}', f'b{target//2}', f'b{4+target//2}', targets[target])
    if len(gates) != 24 or any(nodes[f'c{i}'] != targets[i] for i in range(8)):
        raise ArithmeticError("fused unpack changed output bytes")
    # Every bit of every input byte follows that byte without modification.
    # This checks all4096 post-GFNI input-bit bases symbolically at once.
    for output, wanted in enumerate(targets):
        for byte in range(64):
            for bit in range(8):
                if 8*nodes[f'c{output}'][byte]+bit != 8*wanted[byte]+bit:
                    raise ArithmeticError("fused output changed a payload bit")
    return gates


def emit_fused_output(parallel, gates):
    # Indices are immutable constants, shared by both output groups.  Express
    # them as qwords to avoid compiler-specific 64-argument byte constructors.
    indices = list(dict.fromkeys(gate[3] for gate in gates))
    for i, index in enumerate(indices):
        words = [sum(index[8*j+k] << (8*k) for k in range(8)) for j in range(8)]
        print(f'const auto fusedIndex{i}=_mm512_setr_epi64('
              + ','.join(f'0x{word:016x}ULL' for word in words) + ');')
    for output in range(parallel):
        print('{')
        print('const auto basis=_mm512_set1_epi64(0x8040201008040201ULL);')
        for plane in range(8):
            print(f'const auto x{plane}=_mm512_gf2p8affine_epi64_epi8('
                  f'basis,y{8*output+plane},0);')
        for name, left, right, index in gates:
            print(f'const auto {name}=_mm512_permutex2var_epi8('
                  f'{left},fusedIndex{indices.index(index)},{right});')
            if name.startswith('c'):
                target = int(name[1:])
                row, first = target % 4, 4*(target//4)
                print(f'_mm512_storeu_si512(out+{128*row}+8*(output+{output})+{first},{name});')
        print('}')


def generate(header):
    gates = network()
    reference = subprocess.run([sys.executable, str(Path(layout.__file__).resolve()),
                                str(header), '--mode', '2', '--tile-mode', '5'],
                               check=True, capture_output=True, text=True)
    if reference.stderr:
        print(reference.stderr, file=sys.stderr, end='')
    source = reference.stdout
    old_tile = definition(source, 'template<unsigned output> static SPIN_NOINLINE void packedCoeffTile(')
    old_output = StringIO()
    with redirect_stdout(old_output):
        tune.emit_output(2, packed_layout=True)
    if old_tile.count(old_output.getvalue().rstrip()) != 1:
        raise ArithmeticError("retained tile output no longer matches authenticated emitter")
    new_output = StringIO()
    with redirect_stdout(new_output):
        emit_fused_output(2, gates)
    new_tile = old_tile.replace(old_output.getvalue().rstrip(), new_output.getvalue().rstrip())
    return ('// Exact post-GFNI output layout: 24 VBMI permutes replace48 shuffles.\n'
            '#ifndef __AVX512VBMI__\n'
            '#if !defined(SPIN_PACKET_ENABLE_VBMI) || !SPIN_PACKET_ENABLE_VBMI\n'
            '#error Fused BCH output requires AVX512VBMI\n#endif\n#endif\n'
            + source.replace(old_tile, new_tile))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('header', type=Path, nargs='?')
    parser.add_argument('--output', type=Path)
    parser.add_argument('--verify-only', action='store_true')
    args = parser.parse_args()
    gates = network()
    if args.verify_only:
        tune.verify()
        print('All4096 fused byte-network bit bases and retained pack/unpack bases verified; '
              f'24 VBMI gates, {len(set(gate[3] for gate in gates))} unique index vectors.')
        return
    if args.header is None:
        parser.error('header is required unless --verify-only is supplied')
    source = generate(args.header)
    if args.output:
        args.output.write_text(source, encoding='utf-8')
    else:
        print(source)
    print('Fused BCH output: eight unchanged GFNI +24 VBMI permutes per8-coordinate group; '
          'all original/aligned/compact exports retained.', file=sys.stderr)


if __name__ == '__main__':
    main()
