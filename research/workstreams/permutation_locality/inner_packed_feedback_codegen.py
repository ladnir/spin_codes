"""Emit exact P A^T feedback directly into the retained Packed<16> layout.

The pruned high-index zeta and low-lane moments are retained. The final
16-word syndrome is replaced by column-wise accumulation of the byte
matrices consumed by four GFNI instructions. This is a bounded experiment:
its additional byte shuffles and constants may outweigh removed packing.
Generated source is printed to stdout; this script does not change files.
"""

import argparse
from collections import defaultdict
import sys

import inner_packet_codegen as packet


def check(condition, message):
    if not condition:
        raise ArithmeticError(message)


def xor(values):
    result = 0
    for value in values:
        result ^= value
    return result


def unpack_pair(a, b, width, high):
    start = 8 if high else 0
    return [value for offset in range(start, start + 8, width)
            for source in (a, b) for value in source[offset:offset + width]]


def byte_transpose(words):
    t = [unpack_pair(words[i], words[i + 1], 1, high)
         for i in range(0, 8, 2) for high in (False, True)]
    u = [unpack_pair(t[a], t[b], 2, high)
         for a, b in ((0, 2), (1, 3), (4, 6), (5, 7))
         for high in (False, True)]
    return [unpack_pair(u[a], u[b], 4, high)
            for a, b in ((0, 4), (1, 5), (2, 6), (3, 7))
            for high in (False, True)]


def prepare():
    rows, _, _, _, _, digest, source_digest = packet.retained.prepare()
    monomial_rows = [packet.anf(row) for row in rows]
    columns = {m: sum(((row >> j) & 1) << i for i, row in enumerate(monomial_rows))
               for j, m in enumerate(packet.MONOMIALS)}
    basis = [1 << j for j in range(7)] + [columns[m] for m in (3, 5, 9, 17, 33, 6, 10, 34, 12)]
    packet.inverse(basis)
    feedback = [packet.compose(row, monomial_rows) for row in basis]
    physical = [packet.compose(row, rows) for row in basis]
    groups = []
    for group in range(2):
        grouped = defaultdict(int)
        for j in range(len(packet.MONOMIALS)):
            column = sum(((feedback[8 * group + r] >> j) & 1) << r for r in range(8))
            if column:
                grouped[column] |= 1 << j
        groups.append(sorted(grouped.items()))

    # Independently authenticate every high-node/lane and low moment as a
    # linear form in the 64 input words, including pruned boundary masks.
    high = {h: [sum(1 << (4 * p + lane) for p in range(16) if p & h == h)
                for lane in range(4)] for h in packet.PACKET_MONOMIALS}
    moment_forms = []
    for m in packet.MONOMIALS:
        form = xor(high[m >> 2][lane] for lane in range(4) if lane & (m & 3) == (m & 3))
        check(form == sum(1 << p for p in range(64) if p & m == m), "incorrect high/low zeta")
        moment_forms.append(form)
    for group, terms in enumerate(groups):
        for row in range(8):
            reconstructed = xor(form for column, form in terms if column >> row & 1)
            check(reconstructed == feedback[8 * group + row], "column grouping changed feedback")
            check(packet.compose(reconstructed, moment_forms) == physical[8 * group + row],
                  "feedback changed the physical map")

    # Symbolic bytes carry all eight bit expressions at once. Compare the
    # candidate's entire 256-byte GFNI input against the exact unpack-based
    # byteTranspose network in gfni_r4::pack, not merely an assumed layout.
    def payload_byte(form, byte):
        return tuple(sum(1 << (128 * p + 8 * byte + bit) for p in range(64) if form >> p & 1)
                     for bit in range(8))

    zero = (0,) * 8
    symbolic_bits = 0
    for group, terms in enumerate(groups):
        state = [[payload_byte(physical[8 * group + row], byte) for byte in range(16)]
                 for row in range(8)]
        reference = sum(byte_transpose(list(reversed(state))), [])
        for half in range(2):
            candidate = [zero] * 64
            for column, form in terms:
                raw_form = packet.compose(form, moment_forms)
                for pos in range(64):
                    if column >> (7 - pos % 8) & 1:
                        source = payload_byte(raw_form, 8 * half + pos // 8)
                        candidate[pos] = tuple(a ^ b for a, b in zip(candidate[pos], source))
            check(candidate == reference[64 * half:64 * half + 64],
                  "composed byte matrix differs from retained pack")
            symbolic_bits += 512
    check(symbolic_bits == 2048, "incomplete packed-state bit coverage")
    return groups, basis, digest, source_digest


def scalar_expression(form):
    names = [f"z{m}" for j, m in enumerate(packet.MONOMIALS) if form >> j & 1]
    check(bool(names), "empty direct-packing term")
    value = names[0]
    position = 1
    while position + 1 < len(names):
        value = f"_mm_ternarylogic_epi64({value},{names[position]},{names[position + 1]},0x96)"
        position += 2
    if position < len(names):
        value = f"_mm_xor_si128({value},{names[position]})"
    return value


def generate(groups, basis, digest, source_digest):
    lines = [f"""// Generated by inner_packed_feedback_codegen.py; exact transformed feedback.
// Selected map SHA256: {digest}
// Selected JSON SHA256: {source_digest}
// Direct byte-matrix accumulation: an experimental schedule, not a speedup claim.
#pragma once
#include "InnerPacketMaps.h"
#include "InnerPacketFeedback.h"
#include <array>
#include <cstdint>
#include <cstring>
namespace spin::research::packet_inner {{
namespace direct_feedback_detail {{
template<unsigned Column,unsigned Half>
static SPIN_FORCEINLINE __m512i matrixBytes(__m512i repeated) {{
    static_assert(Column>0 && Column<256 && Half<2);
    // Each qword is a GFNI matrix: reversed output-coordinate bytes.
    // Shuffle controls include zeroing, so no separate coefficient mask is needed.
    alignas(64) static constexpr auto controls=[] {{
        std::array<std::uint8_t,64> result{{}};
        for(unsigned pos=0;pos<64;++pos)
            result[pos]=(Column>>(7-pos%8))&1 ? 8*Half+pos/8 : 0x80;
        return result;
    }}();
    return _mm512_shuffle_epi8(repeated,_mm512_load_si512(controls.data()));
}}
// The same 23-operation high-index zeta as packetMoments, retaining only
// the eleven needed masks. All physical 128-bit lanes remain in order.
static SPIN_FORCEINLINE void highMoments(const __m512i* packets,__m512i* high) {{
    const auto a=feedback_detail::firstStage<0>(packets);
    const auto b=feedback_detail::firstStage<4>(packets);
    const auto c=feedback_detail::firstStage<8>(packets);
    const auto d=feedback_detail::firstStage<12>(packets);
    high[4]=_mm512_xor_si512(b.z0,d.z0);
    high[8]=_mm512_xor_si512(c.z0,d.z0);
    high[0]=_mm512_ternarylogic_epi64(a.z0,c.z0,high[4],0x96);
    high[12]=d.z0;
    high[5]=_mm512_xor_si512(b.z1,d.z1);
    high[9]=_mm512_xor_si512(c.z1,d.z1);
    high[1]=_mm512_ternarylogic_epi64(a.z1,c.z1,high[5],0x96);
    high[6]=_mm512_xor_si512(b.z2,d.z2);
    high[10]=_mm512_xor_si512(c.z2,d.z2);
    high[2]=_mm512_ternarylogic_epi64(a.z2,c.z2,high[6],0x96);
    high[3]=_mm512_ternarylogic_epi64(a.z3,b.z3,_mm512_xor_si512(c.z3,d.z3),0x96);
}}
}} // namespace direct_feedback_detail
// high[h]=XOR packets[j] over (j&h)==h, with four physical 128-bit words
// per packet. Only popcount(h)<=2 is read. The output is the exact existing
// Packed<16> representation of Maps<true>::finish; no syndrome-word array
// or generic pack call occurs here. Requires AVX2, AVX-512F/BW/VL and GFNI.
static SPIN_FORCEINLINE void directFeedbackHigh(
    const __m512i* high,spin::research::gfni_r4::Packed<16>& out) {{
    using direct_feedback_detail::matrixBytes;"""]
    for j, row in enumerate(basis):
        lines.append(f"static_assert(Maps<true>::basisRows[{j}]=={hex(row)});")
    for h in packet.PACKET_MONOMIALS:
        m = 4 * h
        lines += [f"const auto upper{h}=_mm512_extracti64x4_epi64(high[{h}],1);",
                  f"const auto folded{h}=_mm256_xor_si256(_mm512_castsi512_si256(high[{h}]),upper{h});",
                  f"const auto odd{h}=_mm256_extracti128_si256(folded{h},1);",
                  f"const auto z{m}=_mm_xor_si128(_mm256_castsi256_si128(folded{h}),odd{h});"]
        if h.bit_count() <= 1:
            lines += [f"const auto last{h}=_mm256_extracti128_si256(upper{h},1);",
                      f"const auto z{m+1}=odd{h};",
                      f"const auto z{m+2}=_mm_xor_si128(_mm256_castsi256_si128(upper{h}),last{h});"]
        if h == 0:
            lines.append("const auto z3=last0;")
    for group, terms in enumerate(groups):
        lines.append(f"{{ // Packed coordinate group {group}: {len(terms)} distinct coefficient columns.")
        # Each moment combination is broadcast once and consumed by both
        # payload halves. Streaming the pair bounds accumulator liveness.
        for term, (column, form) in enumerate(terms):
            lines.append(f"const auto r{term}=_mm512_broadcast_i32x4({scalar_expression(form)});")
            for half in range(2):
                lines.append(f"const auto t{half}_{term}=matrixBytes<{column},{half}>(r{term});")
            if term == 0:
                lines.extend(f"auto a{half}=t{half}_0;" for half in range(2))
            elif term % 2 == 0:
                lines.extend(f"a{half}=_mm512_ternarylogic_epi64(a{half},t{half}_{term-1},t{half}_{term},0x96);"
                             for half in range(2))
        if len(terms) % 2 == 0:
            lines.extend(f"a{half}=_mm512_xor_si512(a{half},t{half}_{len(terms)-1});" for half in range(2))
        lines.append("const auto basis=_mm512_set1_epi64(0x0102040810204080ULL);")
        lines.extend(f"out.v[{2*group+half}]=_mm512_gf2p8affine_epi64_epi8(basis,a{half},0);" for half in range(2))
        lines.append("}")
    lines.append("""}
// Ordinary-function scratch only; do not inline this array into a coroutine.
static SPIN_FORCEINLINE void directFeedback(const __m512i* packets,
    spin::research::gfni_r4::Packed<16>& out) {
    alignas(64) __m512i high[16];
    direct_feedback_detail::highMoments(packets,high);
    directFeedbackHigh(high,out);
}
// Setup-only native check against the retained word finisher and pack.
// All 8192 individual input-bit bases are checked, including both payload
// halves and every packed state coordinate. Checks stay active under NDEBUG.
[[nodiscard]] inline bool directFeedbackSelfCheck() {
    alignas(64) std::uint64_t words[128]{};
    alignas(64) __m512i packets[16];
    alignas(64) __m128i moments[64],syndrome[16];
    spin::research::gfni_r4::Packed<16> expected,actual;
    for(unsigned p=0;p<64;++p) for(unsigned bit=0;bit<128;++bit) {
        words[2*p+(bit>>6)]=std::uint64_t(1)<<(bit&63);
        for(unsigned j=0;j<16;++j)packets[j]=_mm512_load_si512(words+8*j);
        packetMoments(packets,moments);
        Maps<true>::finish(moments,syndrome);
        spin::research::gfni_r4::pack<16>(syndrome,expected);
        directFeedback(packets,actual);
        if(std::memcmp(&expected,&actual,sizeof(actual))!=0)return false;
        words[2*p+(bit>>6)]=0;
    }
    return true;
}
} // namespace spin::research::packet_inner
""")
    return "\n".join(lines)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--verify-only", action="store_true")
    args = parser.parse_args()
    groups, basis, digest, source_digest = prepare()
    broadcasts = sum(len(group) for group in groups)
    accum = sum(2 * (len(group) // 2) for group in groups)
    print(f"Verified all 2048 packed output bits as linear forms in all 8192 input bits; "
          f"retained transpose and exact P A^T agree. Columns {[len(g) for g in groups]}; "
          f"after low folding: {broadcasts} broadcasts, {2*broadcasts} byte shuffles, "
          f"{accum} accumulator logical ops, 4 final GFNI. No speedup is assumed.", file=sys.stderr)
    if not args.verify_only:
        print(generate(groups, basis, digest, source_digest), end="")


if __name__ == "__main__":
    main()
