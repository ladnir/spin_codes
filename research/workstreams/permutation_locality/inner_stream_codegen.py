"""Exact four-packet streaming expansion, routing, and feedback moments.

The state uses InnerPacketMaps.h Maps<true>'s basis.  Original certified
sources are read and authenticated, never rewritten.  Emitters receive
complete four-element packets in descending physical index order.
"""
import argparse
from pathlib import Path
import sys

import inner_packet_codegen as maps


HIGH = maps.PACKET_MONOMIALS


def vx(*values):
    return tuple(__import__('functools').reduce(int.__xor__, lane, 0) for lane in zip(*values))


def prepare():
    rows, columns, _, _, _, digest, _ = maps.retained.prepare()
    monomial_rows = [maps.anf(row) for row in rows]
    old_coefficients = {m: sum(((row >> j) & 1) << i for i, row in enumerate(monomial_rows))
                        for j, m in enumerate(maps.MONOMIALS)}
    basis = [1 << j for j in range(7)] + [old_coefficients[m] for m in (3, 5, 9, 17, 33, 6, 10, 34, 12)]
    inverse = maps.inverse(basis)
    coefficients = {m: maps.compose(form, inverse) for m, form in old_coefficients.items()}
    expected = [maps.compose(column, inverse) for column in columns]
    maps.validate_packets(coefficients, expected)
    return coefficients, expected, digest


def verify(coefficients, expected):
    d = maps.packet_coefficient_forms(coefficients)
    for pruned in (False, True):
        cache = (vx(d[0], d[8]), vx(d[1], d[9]), vx(d[2], d[10]))
        packets, groups = {}, {}
        for group in (3, 2, 1, 0):
            if group == 3:
                if pruned:
                    l0, l1, l2 = vx(cache[0], d[4], d[12]), vx(cache[1], d[5]), vx(cache[2], d[6])
                else:
                    l0, l1, l2 = vx(d[0], d[4], d[8], d[12]), vx(d[1], d[5], d[9]), vx(d[2], d[6], d[10])
            elif group == 2:
                l0, l1, l2 = cache if pruned else (vx(d[0], d[8]), vx(d[1], d[9]), vx(d[2], d[10]))
            elif group == 1:
                l0, l1, l2 = vx(d[0], d[4]), vx(d[1], d[5]), vx(d[2], d[6])
            else:
                l0, l1, l2 = d[0], d[1], d[2]
            expansion = (l0, vx(l0, l1), vx(l0, l2), vx(l0, l1, l2, d[3]))
            ys = []
            for lane in range(4):
                h = 4 * group + lane
                if expansion[lane] != tuple(expected[4 * h + k] for k in range(4)):
                    raise ArithmeticError("streaming expansion changed coordinates")
                # 64 independent input words plus 16 independent state words.
                y = tuple((1 << (4 * h + k)) ^ (form << 64) for k, form in enumerate(expansion[lane]))
                packets[h] = y
                ys.append(y)
            groups[group] = (vx(*ys), vx(ys[1], ys[3]), vx(ys[2], ys[3]), ys[3])
        incremental = {}
        for group in (3, 2, 1, 0):
            for h in HIGH:
                if group & (h >> 2) == (h >> 2):
                    term = groups[group][h & 3]
                    incremental[h] = vx(incremental[h], term) if h in incremental else term
        a, b, c, dgroup = (groups[g] for g in range(4))
        grouped = {3: vx(a[3], b[3], c[3], dgroup[3]), 12: dgroup[0]}
        for low in range(3):
            grouped[low] = vx(a[low], b[low], c[low], dgroup[low])
            grouped[4 + low] = vx(b[low], dgroup[low])
            grouped[8 + low] = vx(c[low], dgroup[low])
        for h in HIGH:
            wanted = vx(*(packets[p] for p in range(16) if p & h == h))
            if incremental[h] != wanted or grouped[h] != wanted:
                raise ArithmeticError("streaming high moments changed the linear map")
        for mask in maps.MONOMIALS:
            h, low = mask >> 2, mask & 3
            got = 0
            for lane in range(4):
                if lane & low == low:
                    got ^= incremental[h][lane]
            want = 0
            for p in range(64):
                if p & mask == mask:
                    want ^= packets[p >> 2][p & 3]
            if got != want:
                raise ArithmeticError("streaming low-lane fold changed the linear map")


def emit_coefficients(coefficients):
    code, cost = maps.emit_cse([coefficients[m] for m in maps.MONOMIALS],
        [f"state[{j}]" for j in range(16)], "cs", [f"c{m}" for m in maps.MONOMIALS])
    code.insert(0, f"// State to 22 ANF coefficients: {cost} XMM XORs.")
    for h in HIGH:
        m = 4 * h
        if h.bit_count() == 2:
            code.append(f"const auto d{h}=_mm512_broadcast_i32x4(c{m});")
        else:
            code.extend((f"const auto l{h}1=vx(c{m},c{m+1});", f"const auto l{h}2=vx(c{m},c{m+2});"))
            last = f"vx(l{h}1,c{m+2})"
            if h == 0:
                last = f"vx({last},c3)"
            code.extend((f"const auto l{h}3={last};", f"const auto d{h}=join(c{m},l{h}1,l{h}2,l{h}3);"))
    code.extend((
        "// Shared group-2 coefficients reduce the block evaluator from29 to26 XORs.",
        "const auto cache0=Pruned?_mm512_xor_si512(d0,d8):_mm512_setzero_si512();",
        "const auto cache1=Pruned?_mm512_xor_si512(d1,d9):_mm512_setzero_si512();",
        "const auto cache2=Pruned?_mm512_xor_si512(d2,d10):_mm512_setzero_si512();",
    ))
    return code


def xor_expression(*terms):
    result = terms[0]
    for term in terms[1:]:
        result = f"_mm512_xor_si512({result},{term})"
    return result


def emit_group(group):
    code = [f"__m512i g{group}0,g{group}1,g{group}2,g{group}3;", "{"]
    if group == 3:
        code.extend((
            "const auto b0=Pruned?" + xor_expression("cache0", "d4", "d12") + ":" + xor_expression("d0", "d4", "d8", "d12") + ";",
            "const auto b1=Pruned?" + xor_expression("cache1", "d5") + ":" + xor_expression("d1", "d5", "d9") + ";",
            "const auto b2=Pruned?" + xor_expression("cache2", "d6") + ":" + xor_expression("d2", "d6", "d10") + ";",
        ))
    elif group == 2:
        code.extend(f"const auto b{j}=Pruned?cache{j}:_mm512_xor_si512(d{j},d{8+j});" for j in range(3))
    elif group == 1:
        code.extend(f"const auto b{j}=_mm512_xor_si512(d{j},d{4+j});" for j in range(3))
    else:
        code.extend(f"const auto b{j}=d{j};" for j in range(3))
    code.extend((
        "const auto e1=_mm512_xor_si512(b0,b1);",
        "const auto e2=_mm512_xor_si512(b0,b2);",
        "const auto e3=_mm512_xor_si512(_mm512_xor_si512(e1,b2),d3);",
    ))
    for low in (3, 2, 1, 0):
        h = 4 * group + low
        expansion = f"e{low}" if low else "b0"
        code.append(f"const auto y{low}=_mm512_xor_si512({expansion},_mm512_loadu_si512(raw+{4*h}));")
        code.append(f"emit(packetBase+{h},y{low});")
    code.extend((
        f"g{group}1=_mm512_xor_si512(y1,y3);",
        f"g{group}2=_mm512_xor_si512(y2,y3);",
        f"g{group}0=_mm512_ternarylogic_epi64(y0,y2,g{group}1,0x96);",
        f"g{group}3=y3;",
        "}",
    ))
    return code


def emit_incremental():
    result, nodes, number = [], {}, 0
    for group in (3, 2, 1, 0):
        result.extend(emit_group(group))
        for h in HIGH:
            if group & (h >> 2) != (h >> 2):
                continue
            term = f"g{group}{h&3}"
            if h in nodes:
                name = f"s{number}"
                number += 1
                result.append(f"const auto {name}=_mm512_xor_si512({nodes[h]},{term});")
                nodes[h] = name
            else:
                nodes[h] = term
    if number != 18:
        raise ArithmeticError("unexpected incremental high-moment cost")
    result.extend(f"high[{h}]={nodes[h]};" for h in HIGH)
    return result


def emit_grouped():
    result = emit_group(3) + emit_group(2)
    # Fold C/D before generating B/A. Only seven moments survive this point,
    # not all eight C/D values. The four output packets have already died.
    result.extend(f"const auto cd{j}=_mm512_xor_si512(g2{j},g3{j});" for j in range(4))
    result.extend(emit_group(1))
    result.extend(emit_group(0))
    for j in range(3):
        result.append(f"high[{j}]=_mm512_ternarylogic_epi64(g0{j},g1{j},cd{j},0x96);")
        result.append(f"high[{4+j}]=_mm512_xor_si512(g1{j},g3{j});")
        result.append(f"high[{8+j}]=cd{j};")
    result.extend(("high[3]=_mm512_ternarylogic_epi64(g03,g13,cd3,0x96);", "high[12]=g30;"))
    return result


def generate():
    coefficients, expected, digest = prepare()
    verify(coefficients, expected)
    source = [f"""// Generated by inner_stream_codegen.py. Exact packet-streaming t64/s16 maps.
// Authenticated selected map SHA256: {digest}
#pragma once
#include "InnerPacketMaps.h"
#include "InnerPacketFeedback.h"
#include <cstddef>
namespace spin::research::packet_inner {{
// Pruned shares quadratic subexpressions:26 vs29 source ZMM XORs.
// Incremental immediately accumulates eleven high moments:30 feedback
// vector operations. Grouped uses a pruned C/D tree:23 feedback operations.
// Both schedules retain at most four completed output packets at a time.
// state is in Maps<true>'s basis. raw contains64 input elements. high[h]
// receives sum(packet[j] : (j&h)==h) for popcount(h)<=2; other slots untouched.
// State, raw input, and high output must not overlap. emit must not mutate
// them; the routing emitter writes a separate scratch allocation.
template<bool Pruned,bool Incremental=true,class Emit>
static SPIN_FORCEINLINE void streamStepHigh(const __m128i* state,
const spin::detail::kernel::block* raw,std::size_t packetBase,Emit& emit,__m512i* high){{"""]
    source.extend(emit_coefficients(coefficients))
    source.append("if constexpr(Incremental){")
    source.extend(emit_incremental())
    source.append("}else{")
    source.extend(emit_grouped())
    source.append("}\n}")
    source.append("""template<bool Pruned,bool Incremental=true,class Emit>
static SPIN_FORCEINLINE void streamStep(const __m128i* state,
const spin::detail::kernel::block* raw,std::size_t packetBase,Emit& emit,__m128i* moments){
alignas(64) __m512i high[16];
streamStepHigh<Pruned,Incremental>(state,raw,packetBase,emit,high);
using feedback_detail::storeLowMoments;""")
    for h in HIGH:
        source.append(f"storeLowMoments<{4*h},{2-h.bit_count()}>(high[{h}],moments);")
    source.append("}\n}\n")
    print("Verified original and shared packet evaluation, incremental and grouped high moments, "
          "and all22 low moments on the joint64-input/16-state symbolic basis; no permutation changed.", file=sys.stderr)
    return "\n".join(source)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    header = generate()
    if args.output:
        args.output.write_text(header, encoding="utf-8", newline="\n")
    else:
        print(header, end="")
