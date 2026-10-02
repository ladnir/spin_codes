"""Generate the independent RS16/GF16+GL16 outer transpose candidate.

Each four-row GF16 symbol is16 binary coordinates. Packing concatenates two
row nibbles per byte, so a GL16 is four8-by8 blocks rather than sixteen4-by4
blocks. Four ZMM vectors carry all128 payload bits of one symbol. The exact
19-product parity circuit is lifted with duplicated-nibble binary adjoints.
No compilation, remote work, or timing is performed by this generator.
"""
from contextlib import redirect_stdout
from io import StringIO
from pathlib import Path
import argparse
import hashlib
import random
import sys

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
import rs16_maps as maps


def affine(value, matrix):
    return sum((((matrix >> (8 * (7 - bit))) & value).bit_count() & 1) << bit
               for bit in range(8))


def adjoint_matrix(coefficient):
    return sum((maps.multiply(coefficient, 1 << (bit % 4)) << (4 * (bit // 4)))
               << (8 * (7 - bit)) for bit in range(8))


def pack_targets():
    # Input vectors are four columns, each containing four128-bit row records.
    # Before bit transpose, each qword holds8 reversed coordinate bytes.
    return [[64 * (coordinate % 4) + 16 * (2 * group + coordinate // 4)
             + 8 * half + payload_byte
             for payload_byte in range(8) for coordinate in range(7, -1, -1)]
            for group in range(2) for half in range(2)]


def pack_network():
    targets = pack_targets()
    nodes = {f"raw{i}": list(range(64 * i, 64 * (i + 1))) for i in range(4)}
    gates = []

    def gate(name, left, right, wanted):
        source = nodes[left] + nodes[right]
        lookup = {value: i for i, value in enumerate(source)}
        index = tuple(lookup[value] for value in wanted)
        assert len(index) == 64
        nodes[name] = [source[i] for i in index]
        assert nodes[name] == wanted
        gates.append((name, left, right, index))

    for pair in range(2):
        for group in range(2):
            wanted = [byte for target in targets[2 * group:2 * group + 2]
                      for byte in target if byte // 128 == pair]
            gate(f"a{2 * pair + group}", f"raw{2 * pair}", f"raw{2 * pair + 1}", wanted)
    for output, wanted in enumerate(targets):
        gate(f"b{output}", f"a{output // 2}", f"a{2 + output // 2}", wanted)
    assert sorted(sum(targets, [])) == list(range(256))
    return gates


def unpack_indices(row_in_pair):
    # After inverse bit transpose, byte8*b+c is coordinate c's payload byte b.
    return tuple(64 * (payload_byte // 8) + 8 * (payload_byte % 8)
                 + 4 * row_in_pair + column
                 for column in range(4) for payload_byte in range(16))


def transpose_bits(vector, descending):
    return [sum(((vector[base + 7 - row] >> (7 - byte if descending else byte)) & 1) << row
                for row in range(8))
            for base in range(0, 64, 8) for byte in range(8)]


def verify():
    gates = pack_network()
    # Complete2048-bit input-symbol basis: coordinate order, payload-bit order,
    # and inverse row-major stores, including both64-bit halves of every record.
    for bit in range(2048):
        nodes = {f"raw{i}": [0] * 64 for i in range(4)}
        nodes[f"raw{bit // 512}"][(bit % 512) // 8] = 1 << (bit % 8)
        for name, left, right, index in gates:
            source = nodes[left] + nodes[right]
            nodes[name] = [source[i] for i in index]
        packed = [transpose_bits(nodes[f"b{i}"], True) for i in range(4)]
        row = (bit % 512) // 128
        column = bit // 512
        payload = bit % 128
        expected = [[0] * 64 for _ in range(4)]
        expected[2 * (row // 2) + payload // 64][8 * ((payload // 8) % 8) + 7 - payload % 8] = 1 << (4 * (row % 2) + column)
        assert packed == expected
        before_store = [transpose_bits(value, False) for value in packed]
        restored = []
        for output_row in range(4):
            source = before_store[2 * (output_row // 2)] + before_store[2 * (output_row // 2) + 1]
            restored.append([source[i] for i in unpack_indices(output_row % 2)])
        wanted = [[0] * 64 for _ in range(4)]
        wanted[row][16 * column + payload // 8] = 1 << (payload % 8)
        assert restored == wanted
    # Every coefficient entry of the two-by-two8-bit GL16 blocking.
    for output in range(16):
        for source in range(16):
            coefficients = [0] * 4
            coefficients[2 * (output // 8) + source // 8] = 1 << (8 * (7 - output % 8) + source % 8)
            value = 1 << source
            lo = affine(value & 255, coefficients[0]) ^ affine(value >> 8, coefficients[1])
            hi = affine(value & 255, coefficients[2]) ^ affine(value >> 8, coefficients[3])
            assert lo | (hi << 8) == 1 << output
    assert maps.check_parity()["factored_nontrivial_products"] == 19
    for coefficient in maps.SQUARE_ZERO_COEFFICIENTS:
        for byte in range(256):
            expected = maps.adjoint_multiply(coefficient, byte & 15) | (maps.adjoint_multiply(coefficient, byte >> 4) << 4)
            assert affine(byte, adjoint_matrix(coefficient)) == expected
    # Every complete-group coordinate basis under arbitrary GL16 maps, followed
    # by the factored binary-adjoint parity, compared with the full RS matrix.
    rng = random.Random(0x52533136)
    matrices = [[rng.getrandbits(16) for _ in range(16)] for _ in range(16)]
    for basis in range(256):
        mixed = [[0] * 4 for _ in range(16)]
        symbol, bit = divmod(basis, 16)
        for output, mask in enumerate(matrices[symbol]):
            mixed[symbol][output // 4] |= ((mask >> bit) & 1) << (output % 4)
        for row in range(4):
            got = maps.transpose_row([mixed[s][row] for s in range(16)])
            expected = [0] * 8
            for output in range(8):
                for source in range(16):
                    expected[output] ^= maps.adjoint_multiply(maps.GENERATOR[source][output], mixed[source][row])
            assert got == tuple(expected)
    assert sorted(32 * row + 4 * symbol + column for row in range(4)
                  for symbol in range(8) for column in range(4)) == list(range(128))


def emit_index(name, index):
    words = [sum(index[8 * j + k] << (8 * k) for k in range(8)) for j in range(8)]
    print(f"const auto {name}=_mm512_setr_epi64("
          + ",".join(f"0x{word:016x}ULL" for word in words) + ");")


def emit_mix():
    print("static SPIN_FORCEINLINE void mixSymbol(const Block* input,__m512i* packed,const std::uint64_t* coeff){")
    gates = pack_network()
    indices = list(dict.fromkeys(gate[3] for gate in gates))
    for i, index in enumerate(indices):
        emit_index(f"index{i}", index)
    for i in range(4):
        print(f"const auto raw{i}=_mm512_loadu_si512(input+{4 * i});")
    for name, left, right, index in gates:
        print(f"const auto {name}=_mm512_permutex2var_epi8({left},index{indices.index(index)},{right});")
    print("const auto basis=_mm512_set1_epi64(0x0102040810204080ULL);")
    for i in range(4):
        print(f"const auto v{i}=_mm512_gf2p8affine_epi64_epi8(basis,b{i},0);")
        print(f"const auto m{i}=_mm512_set1_epi64(coeff[{i}]);")
    for group in range(2):
        for half in range(2):
            print(f"_mm512_store_si512(packed+{2 * group + half},_mm512_xor_si512("
                  f"_mm512_gf2p8affine_epi64_epi8(v{half},m{2 * group},0),"
                  f"_mm512_gf2p8affine_epi64_epi8(v{2 + half},m{2 * group + 1},0)));")
    print("}")


def emit_factor():
    print("template<unsigned plane> static SPIN_FORCEINLINE void parityPlane(__m512i* packed){")
    print("static_assert(plane<4);")
    for coefficient in sorted(set(maps.SQUARE_ZERO_COEFFICIENTS) - {0, 1}):
        print(f"const auto m{coefficient}=_mm512_set1_epi64(0x{adjoint_matrix(coefficient):016x}ULL);")
    for symbol in range(8):
        print(f"auto x{symbol}=_mm512_load_si512(packed+{32 + 4 * symbol}+plane);")
    for bit in (1, 2, 4):
        for mask in range(8):
            if not mask & bit:
                print(f"x{mask}=_mm512_xor_si512(x{mask},x{mask | bit});")
    for output in range(8):
        print(f"auto y{output}=x{output};")
        subset = output
        while subset:
            coefficient = maps.SQUARE_ZERO_COEFFICIENTS[subset]
            if coefficient == 1:
                term = f"x{output ^ subset}"
            elif coefficient:
                term = f"_mm512_gf2p8affine_epi64_epi8(x{output ^ subset},m{coefficient},0)"
            else:
                subset = (subset - 1) & output
                continue
            print(f"y{output}=_mm512_xor_si512(y{output},{term});")
            subset = (subset - 1) & output
    for bit in (1, 2, 4):
        for mask in range(8):
            if not mask & bit:
                print(f"y{mask}=_mm512_xor_si512(y{mask},y{mask | bit});")
    for symbol in range(8):
        print(f"_mm512_store_si512(packed+{4 * symbol}+plane,"
              f"_mm512_xor_si512(_mm512_load_si512(packed+{4 * symbol}+plane),y{symbol}));")
    print("}")


def emit_output():
    print("template<unsigned symbol> static SPIN_FORCEINLINE void storeSymbol(const __m512i* packed,Block* out){")
    print("static_assert(symbol<8);")
    emit_index("index0", unpack_indices(0))
    emit_index("index1", unpack_indices(1))
    print("const auto basis=_mm512_set1_epi64(0x8040201008040201ULL);")
    for i in range(4):
        print(f"const auto x{i}=_mm512_gf2p8affine_epi64_epi8(basis,_mm512_load_si512(packed+4*symbol+{i}),0);")
    for row in range(4):
        print(f"_mm512_storeu_si512(out+{32 * row}+4*symbol,"
              f"_mm512_permutex2var_epi8(x{2 * (row // 2)},index{row % 2},x{2 * (row // 2) + 1}));")
    print("}")


def generate():
    output = StringIO()
    with redirect_stdout(output):
        print("// Generated by rs16_outer_codegen.py; regenerate rather than hand-edit.")
        print(f"// Source rs16_maps.py: {hashlib.sha256(Path(maps.__file__).read_bytes()).hexdigest()}")
        print("#pragma once\n#include \"RsPrototype.h\"\nnamespace spin::research::rs::fast16 {")
        emit_mix()
        emit_factor()
        emit_output()
        print("static SPIN_NOINLINE void outerGroup(const Block* __restrict input,Block* __restrict output,const std::uint64_t* __restrict coeff){")
        print("alignas(64) __m512i packed[64];")
        print("for(unsigned symbol=0;symbol<16;++symbol)mixSymbol(input+16*symbol,packed+4*symbol,coeff+4*symbol);")
        for plane in range(4):
            print(f"parityPlane<{plane}>(packed);")
        for symbol in range(8):
            print(f"storeSymbol<{symbol}>(packed,output);")
        print("}\n} // namespace spin::research::rs::fast16")
    return output.getvalue()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    verify()
    result = generate()
    target = HERE / "Rs16OuterFast.h"
    if args.check:
        if not target.exists() or target.read_text() != result:
            raise SystemExit("Rs16OuterFast.h differs from the verified generator")
    else:
        target.write_text(result, encoding="utf-8")
    print("PASS:2048 pack/unpack bits,256 GL16 entries,GF16 nibble adjoints,256 complete outer bases; generated RS16 header current")


if __name__ == "__main__":
    main()
