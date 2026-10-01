"""Generate exact t64/s16 packet-native expansion and conjugated feedback.

No setup distribution or encoded coordinate is changed.  The optional basis
is an invertible reparameterization of the sixteen persistent state words.
All generated XOR circuits are authenticated symbolically before emission.
Run this script to stdout; the generated header is InnerPacketMaps.h.
"""
import argparse
from collections import Counter
from hashlib import sha256
from itertools import combinations
from pathlib import Path
import sys

import gfni_t64_codegen as retained


MONOMIALS = tuple(m for m in range(64) if m.bit_count() <= 2)
PACKET_MONOMIALS = tuple(m for m in range(16) if m.bit_count() <= 2)
REFERENCE_SHA256 = "8165d10c7cc752b5640749bb4aa31e3d295e552e22ff4426f468668bec73c49e"


def compose(form, rows):
    result = 0
    for j, row in enumerate(rows):
        if (form >> j) & 1:
            result ^= row
    if form >> len(rows):
        raise ArithmeticError("linear form exceeds its input dimension")
    return result


def inverse(rows):
    n = len(rows)
    work = [row | (1 << (n + j)) for j, row in enumerate(rows)]
    for j in range(n):
        pivot = next((k for k in range(j, n) if (work[k] >> j) & 1), None)
        if pivot is None:
            raise ArithmeticError("singular proposed state basis")
        work[j], work[pivot] = work[pivot], work[j]
        for k in range(n):
            if k != j and ((work[k] >> j) & 1):
                work[k] ^= work[j]
    result = [row >> n for row in work]
    if any(compose(row, result) != 1 << j for j, row in enumerate(rows)):
        raise ArithmeticError("state basis inverse failed")
    if any(compose(row, rows) != 1 << j for j, row in enumerate(result)):
        raise ArithmeticError("state inverse basis failed")
    return result


def anf(row):
    values = [(row >> p) & 1 for p in range(64)]
    for bit in (1, 2, 4, 8, 16, 32):
        for p in range(64):
            if p & bit:
                values[p] ^= values[p ^ bit]
    if any(values[m] for m in range(64) if m.bit_count() > 2):
        raise ArithmeticError("selected map is not quadratic")
    return sum(values[m] << j for j, m in enumerate(MONOMIALS))


def cse(forms, dimension):
    """Deterministic common-pair extraction, with full linear-form checks."""
    terms = [set(j for j in range(dimension) if (form >> j) & 1) for form in forms]
    if any(form >> dimension for form in forms):
        raise ArithmeticError("CSE form exceeds dimension")
    nodes = [1 << j for j in range(dimension)]
    gates = []
    while True:
        counts = Counter(pair for term in terms for pair in combinations(sorted(term), 2))
        if not counts:
            break
        pair, count = max(counts.items(), key=lambda item: (item[1], tuple(-j for j in item[0])))
        if count < 2:
            break
        replacement = len(nodes)
        gates.append(pair)
        nodes.append(nodes[pair[0]] ^ nodes[pair[1]])
        for term in terms:
            if pair[0] in term and pair[1] in term:
                term.difference_update(pair)
                term.add(replacement)
    for wanted, term in zip(forms, terms):
        got = 0
        for node in term:
            got ^= nodes[node]
        if got != wanted:
            raise ArithmeticError("CSE changed a linear form")
    count = len(gates) + sum(max(0, len(term) - 1) for term in terms)
    return gates, [sorted(term) for term in terms], count


def expression(terms, names):
    if not terms:
        return "_mm_setzero_si128()"
    value = names[terms[0]]
    for term in terms[1:]:
        value = f"vx({value},{names[term]})"
    return value


def emit_cse(forms, inputs, prefix, outputs):
    gates, terms, count = cse(forms, len(inputs))
    names = list(inputs)
    lines = []
    for j, (a, b) in enumerate(gates):
        name = f"{prefix}{j}"
        lines.append(f"const auto {name}=vx({names[a]},{names[b]});")
        names.append(name)
    for name, term in zip(outputs, terms):
        lines.append(f"const auto {name}={expression(term, names)};")
    return lines, count


def packet_coefficient_forms(coefficients):
    result = {}
    for h in PACKET_MONOMIALS:
        low = [coefficients.get(4 * h + b, 0) for b in range(4)]
        result[h] = (low[0], low[0] ^ low[1], low[0] ^ low[2],
                     low[0] ^ low[1] ^ low[2] ^ low[3])
    return result


def pruned_evaluator():
    """A 26-XOR vector-valued quadratic evaluator on four index bits."""
    gates = []

    def recurse(n, degree, base):
        if n == 0 or degree == 0:
            if base not in PACKET_MONOMIALS:
                raise ArithmeticError("unexpected evaluator coefficient")
            return [f"d{base}"] * (1 << n)
        lo = recurse(n - 1, degree, base)
        hi = recurse(n - 1, degree - 1, base | (1 << (n - 1)))
        result = list(lo)
        for a, b in zip(lo, hi):
            name = f"e{len(gates)}"
            gates.append((name, a, b))
            result.append(name)
        return result

    outputs = recurse(4, 2, 0)
    if len(gates) != 26 or len(outputs) != 16:
        raise ArithmeticError("unexpected pruned evaluation cost")
    return gates, outputs


def validate_packets(coefficients, expected):
    packets = packet_coefficient_forms(coefficients)
    # Authenticate the ordinary four-stage zeta evaluator.
    full = [packets.get(h, (0, 0, 0, 0)) for h in range(16)]
    for bit in (1, 2, 4, 8):
        for h in range(16):
            if h & bit:
                full[h] = tuple(a ^ b for a, b in zip(full[h], full[h ^ bit]))
    # Authenticate the separately generated pruned straight-line evaluator.
    nodes = {f"d{h}": value for h, value in packets.items()}
    gates, outputs = pruned_evaluator()
    for name, a, b in gates:
        nodes[name] = tuple(x ^ y for x, y in zip(nodes[a], nodes[b]))
    for h in range(16):
        for lane in range(4):
            wanted = expected[4 * h + lane]
            if full[h][lane] != wanted or nodes[outputs[h]][lane] != wanted:
                raise ArithmeticError("packet evaluator changed an encoded coordinate")


def emit_map(basis, rows, columns, monomial_rows, original_body):
    inv = inverse(basis)
    coefficients = {
        m: compose(sum(((row >> j) & 1) << i for i, row in enumerate(monomial_rows)), inv)
        for j, m in enumerate(MONOMIALS)
    }
    expected = [compose(column, inv) for column in columns]
    validate_packets(coefficients, expected)
    feedback = [compose(p, monomial_rows) for p in basis]
    # Reconstruct each feedback row directly over all 64 physical inputs.
    monomial_sums = [sum(1 << p for p in range(64) if p & m == m) for m in MONOMIALS]
    for p, form in zip(basis, feedback):
        if compose(form, monomial_sums) != compose(p, rows):
            raise ArithmeticError("folded feedback changed the physical input map")
    transformed = basis != [1 << j for j in range(16)]
    result = [f"template<> struct Maps<{str(transformed).lower()}> {{"]
    result.append("static constexpr std::array<std::uint16_t,16> basisRows{" + ",".join(map(hex, basis)) + "};")
    result.append("static constexpr std::array<std::uint16_t,16> inverseRows{" + ",".join(map(hex, inv)) + "};")
    result.append("template<bool Pruned=false> static SPIN_FORCEINLINE void emission(const __m128i* state,__m512i* packets){")
    generated, coeff_cost = emit_cse([coefficients[m] for m in MONOMIALS],
        [f"state[{j}]" for j in range(16)], "cse", [f"c{m}" for m in MONOMIALS])
    result.append(f"// All 22 scalar ANF coefficients: {coeff_cost} XMM XORs.")
    result.extend(generated)
    # Five nonconstant-lane coefficients need 16 XMM XORs and 15 inserts.
    # The remaining six are a single 128-bit broadcast each.
    for h in PACKET_MONOMIALS:
        if h.bit_count() == 2:
            result.append(f"const auto d{h}=_mm512_broadcast_i32x4(c{4*h});")
        else:
            m = 4 * h
            result.append(f"const auto l{h}1=vx(c{m},c{m+1});")
            result.append(f"const auto l{h}2=vx(c{m},c{m+2});")
            last = f"vx(l{h}1,c{m+2})"
            if h == 0:
                last = f"vx({last},c3)"
            result.append(f"const auto l{h}3={last};")
            result.append(f"const auto d{h}=join(c{m},l{h}1,l{h}2,l{h}3);")
    result.append("if constexpr(Pruned){")
    gates, outputs = pruned_evaluator()
    for name, a, b in gates:
        result.append(f"const auto {name}=_mm512_xor_si512({a},{b});")
    for h, name in enumerate(outputs):
        result.append(f"packets[{h}]={name};")
    result.append("}else{")
    for h in range(16):
        value = f"d{h}" if h in PACKET_MONOMIALS else "_mm512_setzero_si512()"
        result.append(f"packets[{h}]={value};")
    for bit in (1, 2, 4, 8):
        for h in range(16):
            if h & bit:
                result.append(f"packets[{h}]=_mm512_xor_si512(packets[{h}],packets[{h^bit}]);")
    result.append("}\n}")
    result.append("static SPIN_FORCEINLINE void finish(const __m128i* z,__m128i* out){")
    if not transformed:
        # Retain the hand-optimized 36-XOR original finisher exactly.
        result.append(original_body)
        feedback_cost = original_body.count("vx(")
    else:
        generated, feedback_cost = emit_cse(feedback, [f"z[{m}]" for m in MONOMIALS],
            "fse", [f"f{j}" for j in range(16)])
        result.append(f"// P A^T on the same 22 zeta moments: {feedback_cost} XMM XORs.")
        result.extend(generated)
        result.extend(f"out[{j}]=f{j};" for j in range(16))
    result.append("}")
    result.append("static SPIN_FORCEINLINE void finishPacked(const __m128i* z,spin::research::gfni_r4::Packed<16>& packed){")
    if not transformed:
        result.append("alignas(32) __m128i values[16];")
        result.append("finish(z,values);spin::research::gfni_r4::pack<16>(values,packed);")
        packed_feedback_cost = feedback_cost
    else:
        packed_feedback_cost = 0
        # Independent CSE circuits deliberately bound live scalar outputs to
        # eight before converting to packed bytes.  Each half is checked by
        # emit_cse and the parent 16x64 feedback map was checked above.
        for half in range(2):
            generated, cost = emit_cse(feedback[8 * half:8 * half + 8],
                [f"z[{m}]" for m in MONOMIALS], "hse", [f"h{j}" for j in range(8)])
            packed_feedback_cost += cost
            result.append(f"{{ // Feedback coordinates {8*half}..{8*half+7}: {cost} XMM XORs.")
            result.append("alignas(32) __m128i half[8];")
            result.extend(generated)
            result.extend(f"half[{j}]=h{j};" for j in range(8))
            result.append(f"packEight(half,packed.v+{2*half});")
            result.append("}")
    result.append("}\n};")
    return "\n".join(result), coeff_cost, feedback_cost, packed_feedback_cost


def generate():
    rows, columns, _, _, body, digest, source_digest = retained.prepare()
    monomial_rows = [anf(row) for row in rows]
    coefficients = {m: sum(((row >> j) & 1) << i for i, row in enumerate(monomial_rows))
                    for j, m in enumerate(MONOMIALS)}
    basis = [1 << j for j in range(7)] + [coefficients[m] for m in (3, 5, 9, 17, 33, 6, 10, 34, 12)]
    inverse(basis)
    if any((a & b).bit_count() & 1 for a in rows for b in rows):
        raise ArithmeticError("selected expansion no longer has A^T A=0")
    original, c0, f0, p0 = emit_map([1 << j for j in range(16)], rows, columns, monomial_rows, body)
    transformed, c1, f1, p1 = emit_map(basis, rows, columns, monomial_rows, body)
    header = f"""// Generated by inner_packet_codegen.py; exact t64/s16 packet maps.
// Selected map SHA256: {digest}
// Selected JSON SHA256: {source_digest}
// Fixed state basis affects only representation. For the reverse encoder,
// use P M^T P^-1 updates and P A^T feedback with the transformed map.
#pragma once
#include "FusedR4Gfni.h"
#include <array>
#include <cstdint>
namespace spin::research::packet_inner {{
static SPIN_FORCEINLINE __m128i vx(__m128i a,__m128i b){{return _mm_xor_si128(a,b);}}
static SPIN_FORCEINLINE __m512i join(__m128i a,__m128i b,__m128i c,__m128i d){{
auto x=_mm512_castsi128_si512(a);
x=_mm512_inserti32x4(x,b,1);x=_mm512_inserti32x4(x,c,2);
return _mm512_inserti32x4(x,d,3);
}}
// The retained Packed<S> API accepts only S=16/19. This is precisely one
// eight-coordinate group of its pack function, without instantiating S=8.
static SPIN_FORCEINLINE void packEight(const __m128i* state,__m512i* out){{
auto x0=state[0],x1=state[1],x2=state[2],x3=state[3];
auto x4=state[4],x5=state[5],x6=state[6],x7=state[7];
spin::research::gfni_r4::byteTranspose(x7,x6,x5,x4,x3,x2,x1,x0);
const auto basis=_mm512_set1_epi64(0x0102040810204080ULL);
out[0]=_mm512_gf2p8affine_epi64_epi8(basis,join(x7,x6,x5,x4),0);
out[1]=_mm512_gf2p8affine_epi64_epi8(basis,join(x3,x2,x1,x0),0);
}}
template<bool Basis> struct Maps;
{original}
{transformed}
}}
"""
    print(f"Validated both exact 64x16 expansion maps, both 16x64 feedback maps, "
          f"basis/inverse, full and pruned packet evaluation; coefficient XORs {c0}/{c1}; "
          f"feedback XORs {f0}/{f1}; split-pack feedback XORs {p0}/{p1}; "
          "packet-lane XORs16, inserts15, broadcasts6; "
          "full/pruned evaluator ZMM XORs32/26 (before compiler simplification).", file=sys.stderr)
    return header


def reference_header():
    source = Path(__file__).with_name("gfni_t64_probe.cpp")
    raw = source.read_bytes()
    if sha256(raw).hexdigest() != REFERENCE_SHA256:
        raise ArithmeticError("retained t64 reference source hash changed")
    text = raw.decode().replace("\r\n", "\n")
    marker = "\nstatic void experiment("
    if text.count(marker) != 1:
        raise ArithmeticError("ambiguous retained t64 reference boundary")
    result = "// Retained t64 reference SHA256: " + REFERENCE_SHA256 + "\n#pragma once\n"
    result += text.split(marker)[0] + "\n}\n"
    if source.read_bytes() != raw:
        raise ArithmeticError("retained t64 reference changed while reading")
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, help="generated map header; omitted prints to stdout")
    parser.add_argument("--reference-output", type=Path, help="also extract the authenticated retained reference header")
    args = parser.parse_args()
    header = generate()
    reference = reference_header() if args.reference_output else None
    if args.output:
        args.output.write_text(header, encoding="utf-8", newline="\n")
    else:
        print(header, end="")
    if args.reference_output:
        args.reference_output.write_text(reference, encoding="utf-8", newline="\n")
        print(f"Retained reference SHA256 {REFERENCE_SHA256}", file=sys.stderr)


if __name__ == "__main__":
    main()
