"""Generate the fixed-width RS outer transpose, with exact composition checks.

The selected inner is included unchanged by RsFast.cpp. This generator reuses
only authenticated pack/GL32/unpack instruction layouts, not the BCH map.
Run without arguments to regenerate RsOuterFast.h, or --check to compare it.
Neither mode compiles or times a kernel.
"""
from contextlib import redirect_stdout
from io import StringIO
from pathlib import Path
import argparse
import hashlib
import random
import sys

HERE = Path(__file__).resolve().parent
LOCALITY = HERE.parents[1] / "permutation_locality"
sys.path.insert(0, str(LOCALITY))
sys.path.insert(0, str(HERE.parent))
import outer_layout_codegen as layout
import outer_unpack_codegen as unpack
import packed_bch_tune_codegen as tune
from packed_coeff_codegen import definition
import rs_maps


def adjoint_matrix(coefficient):
    # Row b of M_c^T is column b of M_c, the byte c*(1<<b).
    return sum(rs_maps.multiply(coefficient, 1 << bit) << (8 * (7 - bit))
               for bit in range(8))


def factor(values):
    a, b, c, d = values
    total = a ^ b ^ c ^ d
    u = layout.affine(total, adjoint_matrix(6))
    v = layout.affine(total, adjoint_matrix(8))
    common = (layout.affine(total, adjoint_matrix(26))
              ^ layout.affine(b ^ d, adjoint_matrix(6))
              ^ layout.affine(c ^ d, adjoint_matrix(8)))
    return (a ^ common, b ^ common ^ u, c ^ common ^ v, d ^ common ^ u ^ v)


def verify():
    tune.verify()              # Every retained 8-coordinate/128-payload pack bit.
    layout.verify()            # Every arbitrary GL32 matrix entry and lane route.
    unpack.network()           # Every bit of the fused byte-output permutation.
    generator = rs_maps.systematic_generator()
    assert generator[:4] == tuple(tuple(int(i == j) for j in range(4)) for i in range(4))
    assert generator[4:] == ((27, 28, 18, 20), (28, 27, 20, 18),
                             (18, 20, 27, 28), (20, 18, 28, 27))
    # All32 binary basis vectors distinguish field-symbol transpose from
    # the actual binary transpose in the polynomial basis.
    for index in range(32):
        value = tuple(1 << (index % 8) if j == index // 8 else 0 for j in range(4))
        expected = [0] * 4
        for output in range(4):
            for source in range(4):
                expected[output] ^= layout.affine(value[source], adjoint_matrix(generator[4 + source][output]))
        assert factor(value) == tuple(expected)
    # Compose arbitrary GL32s with RS on every one-payload-bit basis of a
    # complete group. GL32 acts across row lanes; RS acts within each row.
    rng = random.Random(0x52533834)
    matrices = [[rng.getrandbits(32) for _ in range(32)] for _ in range(8)]
    for basis in range(256):
        source = [[0] * 4 for _ in range(8)]
        source[basis // 32][(basis % 32) // 8] = 1 << (basis % 8)
        mixed = [[0] * 4 for _ in range(8)]
        for symbol in range(8):
            word = sum(source[symbol][row] << (8 * row) for row in range(4))
            for bit, mask in enumerate(matrices[symbol]):
                mixed[symbol][bit // 8] |= ((mask & word).bit_count() & 1) << (bit % 8)
        for row in range(4):
            parity = factor(tuple(mixed[4 + j][row] for j in range(4)))
            got = tuple(mixed[j][row] ^ parity[j] for j in range(4))
            expected = [0] * 4
            for output in range(4):
                for symbol in range(8):
                    expected[output] ^= layout.affine(mixed[symbol][row], adjoint_matrix(generator[symbol][output]))
            assert got == tuple(expected)
    # New row stride32 exactly partitions the128 output records.
    stores = [32 * row + 8 * symbol + first + offset
              for symbol in range(4) for row in range(4)
              for first in (0, 4) for offset in range(4)]
    assert sorted(stores) == list(range(128))


def emit_indices(gates):
    indices = list(dict.fromkeys(gate[3] for gate in gates))
    for i, index in enumerate(indices):
        words = [sum(index[8 * j + k] << (8 * k) for k in range(8)) for j in range(8)]
        print(f"const auto index{i}=_mm512_setr_epi64("
              + ",".join(f"0x{word:016x}ULL" for word in words) + ");")
    return indices


def emit_mix():
    print("static SPIN_FORCEINLINE void mixSymbol(const Block* a,__m512i* packed,const std::uint64_t* coeff){")
    for plane in range(8):
        print(f"auto v{plane}=_mm512_loadu_si512(a+4*{plane});")
    print("packedTunePack5(" + ",".join(f"v{p}" for p in range(8)) + ");")
    print("const auto duplicate=_mm512_setr_epi64(0,0,1,1,2,2,3,3);")
    print("const auto duplicateOdd=_mm512_setr_epi64(3,3,0,0,1,1,2,2);")
    for diagonal in range(4):
        print(f"const auto c{diagonal}=_mm256_load_si256(reinterpret_cast<const __m256i*>(coeff+{4 * diagonal}));")
        index = "duplicateOdd" if diagonal % 2 else "duplicate"
        print(f"const auto m{diagonal}=_mm512_permutexvar_epi64({index},_mm512_castsi256_si512(c{diagonal}));")
    for plane in range(8):
        print("{")
        print(f"const auto half=_mm512_shuffle_i32x4(v{plane},v{plane},0x4e);")
        print(f"const auto odd=_mm512_xor_si512(_mm512_gf2p8affine_epi64_epi8(v{plane},m1,0),"
              "_mm512_gf2p8affine_epi64_epi8(half,m3,0));")
        print(f"const auto mixed=_mm512_ternarylogic_epi64(_mm512_gf2p8affine_epi64_epi8(v{plane},m0,0),"
              "_mm512_gf2p8affine_epi64_epi8(half,m2,0),_mm512_shuffle_i32x4(odd,odd,0x39),0x96);")
        print(f"_mm512_store_si512(packed+{plane},mixed);")
        print("}")
    print("}")


def emit_factor():
    print("template<unsigned plane> static SPIN_FORCEINLINE void parityPlane(__m512i* packed){")
    print("static_assert(plane<8);")
    for name, coefficient in (("m6", 6), ("m8", 8), ("m26", 26)):
        print(f"const auto {name}=_mm512_set1_epi64(0x{adjoint_matrix(coefficient):016x}ULL);")
    for symbol, name in enumerate(("a", "b", "c", "d")):
        print(f"const auto {name}=_mm512_load_si512(packed+{32 + 8 * symbol}+plane);")
    print("const auto total=_mm512_ternarylogic_epi64(a,b,_mm512_xor_si512(c,d),0x96);")
    print("const auto u=_mm512_gf2p8affine_epi64_epi8(total,m6,0);")
    print("const auto v=_mm512_gf2p8affine_epi64_epi8(total,m8,0);")
    print("const auto common=_mm512_ternarylogic_epi64("
          "_mm512_gf2p8affine_epi64_epi8(total,m26,0),"
          "_mm512_gf2p8affine_epi64_epi8(_mm512_xor_si512(b,d),m6,0),"
          "_mm512_gf2p8affine_epi64_epi8(_mm512_xor_si512(c,d),m8,0),0x96);")
    print("const auto commonU=_mm512_xor_si512(common,u);")
    for symbol, (name, term) in enumerate((("a", "common"), ("b", "commonU"),
                                         ("c", "_mm512_xor_si512(common,v)"),
                                         ("d", "_mm512_xor_si512(commonU,v)"))):
        print(f"_mm512_store_si512(packed+{8 * symbol}+plane,"
              f"_mm512_ternarylogic_epi64(_mm512_load_si512(packed+{8 * symbol}+plane),{name},{term},0x96));")
    print("}")


def emit_output():
    gates = unpack.network()
    print("template<unsigned symbol> static SPIN_FORCEINLINE void storeSymbol(const __m512i* packed,Block* out){")
    print("static_assert(symbol<4);")
    indices = emit_indices(gates)
    print("const auto basis=_mm512_set1_epi64(0x8040201008040201ULL);")
    for plane in range(8):
        print(f"const auto x{plane}=_mm512_gf2p8affine_epi64_epi8(basis,_mm512_load_si512(packed+8*symbol+{plane}),0);")
    for name, left, right, index in gates:
        print(f"const auto {name}=_mm512_permutex2var_epi8({left},index{indices.index(index)},{right});")
        if name.startswith("c"):
            target = int(name[1:])
            print(f"_mm512_storeu_si512(out+{32 * (target % 4)}+8*symbol+{4 * (target // 4)},{name});")
    print("}")


def generate():
    output = StringIO()
    with redirect_stdout(output):
        print("// Generated by rs_outer_codegen.py; regenerate rather than hand-edit.")
        for path in (Path(tune.__file__), Path(layout.__file__), Path(unpack.__file__), Path(rs_maps.__file__)):
            print(f"// Source {path.name}: {hashlib.sha256(path.read_bytes()).hexdigest()}")
        print("#pragma once\n#include \"RsPrototype.h\"\nnamespace spin::research::rs::fast {")
        helpers = StringIO()
        with redirect_stdout(helpers):
            tune.emit_byte_transpose_helpers()
        for marker in ("static SPIN_FORCEINLINE void packedTuneByteTranspose5(",
                       "static SPIN_FORCEINLINE void packedTunePack5("):
            print(definition(helpers.getvalue(), marker))
        emit_mix()
        emit_factor()
        emit_output()
        print("static SPIN_NOINLINE void outerGroup(const Block* __restrict a,Block* __restrict out,const std::uint64_t* __restrict coeff){")
        print("alignas(64) __m512i packed[64];")
        print("for(unsigned symbol=0;symbol<8;++symbol)mixSymbol(a+32*symbol,packed+8*symbol,coeff+16*symbol);")
        for plane in range(8):
            print(f"parityPlane<{plane}>(packed);")
        for symbol in range(4):
            print(f"storeSymbol<{symbol}>(packed,out);")
        print("}\n} // namespace spin::research::rs::fast")
    return output.getvalue()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    verify()
    result = generate()
    target = HERE / "RsOuterFast.h"
    if args.check:
        if not target.exists() or target.read_text() != result:
            raise SystemExit("RsOuterFast.h differs from the verified generator")
    else:
        target.write_text(result, encoding="utf-8")
    print("PASS: pack/GL32/unpack basis checks,32 RS-adjoint basis checks,256 complete outer coordinate bases; generated header current")


if __name__ == "__main__":
    main()
