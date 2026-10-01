"""Emit the selected t64/s16 map for an isolated uniform-GL16 experiment.

A is the 64-by-16 matrix declared by t64_s16_selected.json and C=A^T.
Authenticate its rows, columns, rank, spectrum, grouped coordinates and
retained degree-two zeta finisher against generated/SelectedMaps.h.
Only the state-refresh distribution changes in the accompanying probe.
This generator makes no whole-code distance claim.
"""
import ast
from collections import Counter
from hashlib import sha256
import json
from pathlib import Path
import re
import sys

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
MAP_PATH = HERE.parent/'rate_quarter_bch/inner_calibration/maps/t64_s16_selected.json'
SELECTED_PATH = ROOT/'spin/src/kernels/generated/SelectedMaps.h'


def rank(rows):
    pivots = {}
    for word in rows:
        while word:
            bit = word.bit_length()-1
            if bit not in pivots:
                pivots[bit] = word
                break
            word ^= pivots[bit]
    return len(pivots)


def array(section, name):
    match = re.search(r'\b'+name+r'\{([^}]+)\};', section)
    if match is None:
        raise ArithmeticError('missing selected-map array '+name)
    return [int(word.strip(), 0) for word in match[1].split(',')]


def evaluate_expression(expression, variables):
    """Restricted scalar evaluation of the copied XOR-only finish circuit."""
    def visit(node):
        if isinstance(node, ast.Name) and node.id in variables:
            return variables[node.id]
        if (isinstance(node, ast.Subscript) and isinstance(node.value, ast.Name)
                and node.value.id == 'z' and isinstance(node.slice, ast.Constant)
                and type(node.slice.value) is int):
            index = node.slice.value
            if not 0 <= index < 64 or index.bit_count() > 2:
                raise ArithmeticError('selected finisher exceeds pruned degree-two zeta')
            return sum(1 << p for p in range(64) if p & index == index)
        if (isinstance(node, ast.Call) and isinstance(node.func, ast.Name)
                and node.func.id == 'vx' and len(node.args) == 2 and not node.keywords):
            return visit(node.args[0]) ^ visit(node.args[1])
        raise ArithmeticError('unrecognized selected-map finisher expression')
    return visit(ast.parse(expression, mode='eval').body)


def prepare():
    raw_map, raw_selected = MAP_PATH.read_bytes(), SELECTED_PATH.read_bytes()
    record = json.loads(raw_map)
    if record.get('step_bits') != 64 or record.get('state_bits') != 16:
        raise ArithmeticError('selected t64/s16 geometry required')
    rows = [int(word, 16) for word in record['generator_rows_hex']]
    columns = record['columns']
    if (len(rows) != 16 or len(columns) != 64
            or any(not 0 <= row < 1 << 64 for row in rows)
            or any(type(column) is not int or not 0 < column < 1 << 16 for column in columns)
            or rank(rows) != 16 or rank(columns) != 16):
        raise ArithmeticError('rank-16 selected expansion and feedback required')
    reconstructed = [sum(((row >> p) & 1) << j for j, row in enumerate(rows)) for p in range(64)]
    if columns != reconstructed:
        raise ArithmeticError('selected A rows and C=A^T columns disagree')
    images = [0]*(1 << 16)
    for state in range(1, len(images)):
        bit = state & -state
        images[state] = images[state ^ bit] ^ rows[bit.bit_length()-1]
    spectrum = Counter(image.bit_count() for image in images)
    if [spectrum[w] for w in range(65)] != record['a_counts']:
        raise ArithmeticError('fresh selected-map expansion spectrum disagrees')
    text = raw_selected.decode()
    section = text.split('struct Map64S16 {', 1)[1].split('\n};', 1)[0]
    if array(section, 'columns') != columns:
        raise ArithmeticError('JSON and production selected-map columns differ')
    order, grouped = array(section, 'groupOrder'), array(section, 'groupedColumns')
    if sorted(order) != list(range(16)) or grouped != [
            sum(((column >> j) & 1) << i for i, j in enumerate(order)) for column in columns]:
        raise ArithmeticError('selected grouped coordinate permutation differs')
    body = section.split('static inline void finish(const __m128i* z,__m128i* out) {', 1)[1].rsplit('}', 1)[0].strip()
    variables, outputs = {}, {}
    for statement in body.split(';'):
        statement = statement.strip()
        if not statement:
            continue
        assignment = re.fullmatch(r'(?:const auto (v\d+)|out\[(\d+)\])\s*=\s*(.+)', statement)
        if assignment is None:
            raise ArithmeticError('unrecognized selected-map finisher statement')
        value = evaluate_expression(assignment[3], variables)
        if assignment[1]:
            if assignment[1] in variables:
                raise ArithmeticError('duplicate selected-map finisher temporary')
            variables[assignment[1]] = value
        else:
            index = int(assignment[2])
            if index in outputs:
                raise ArithmeticError('duplicate selected-map finisher output')
            outputs[index] = value
    if set(outputs) != set(range(16)) or [outputs[j] for j in range(16)] != rows:
        raise ArithmeticError('selected zeta finisher differs on the exact 64-coordinate basis')
    if MAP_PATH.read_bytes() != raw_map or SELECTED_PATH.read_bytes() != raw_selected:
        raise ArithmeticError('selected source changed during generation')
    canonical = dict(A=list(map(hex, rows)), C=columns)
    digest = sha256(json.dumps(canonical, sort_keys=True, separators=(',', ':')).encode()).hexdigest()
    return rows, columns, order, grouped, body, digest, sha256(raw_map).hexdigest()


def main():
    rows, columns, order, grouped, body, digest, source_digest = prepare()
    print('// Generated by gfni_t64_codegen.py; selected A64, C=A^T, unchanged coordinate order.')
    print('#pragma once\n#include "FusedR4Gfni.h"\nnamespace spin::research::t64 {')
    print('struct Fixed {\nstatic constexpr unsigned T=64,S=16;')
    print('static constexpr const char* name="t64_s16_selected";')
    print(f'static constexpr const char* mapSha256="{digest}";')
    print(f'static constexpr const char* sourceJsonSha256="{source_digest}";')
    print('static constexpr std::array<std::uint32_t,T> columns{'+','.join(map(hex, columns))+'};')
    print('static constexpr auto feedbackColumns=columns; // C=A^T exactly.')
    print('static constexpr std::array<unsigned,S> groupOrder{'+','.join(map(str, order))+'};')
    print('static constexpr std::array<std::uint32_t,T> groupedColumns{'+','.join(map(hex, grouped))+'};')
    print('static constexpr std::uint64_t expansionRows[S]={'+','.join(f'0x{row:016x}ULL' for row in rows)+'};')
    print('static SPIN_FORCEINLINE __m128i vx(__m128i a,__m128i b){return _mm_xor_si128(a,b);}')
    print('static SPIN_FORCEINLINE void finish(const __m128i* z,__m128i* out){\n'+body+'\n}')
    print('static SPIN_FORCEINLINE void emissionTable(const __m128i* state,__m128i table[4][16]){')
    print('alignas(32) __m128i grouped[S];')
    for j in range(16):
        print(f'grouped[{j}]=state[groupOrder[{j}]];')
    print('spin::detail::kernel::tables<S>(grouped,table);\n}')
    print('template<std::size_t R,class Emit> static SPIN_FORCEINLINE void emitPoint(')
    print('const spin::detail::kernel::block* in,__m128i* raw,const __m128i table[4][16],std::size_t base,Emit& emit){')
    print('constexpr unsigned p=T-1-R;const auto value=in[p].mData;raw[p]=value;')
    print('emit(base+p,spin::detail::kernel::block(vx(value,spin::detail::kernel::fixedSum<groupedColumns[p]>(table))));\n}')
    print('template<class Emit,std::size_t... R> static SPIN_FORCEINLINE void emitReverse(')
    print('const spin::detail::kernel::block* in,__m128i* raw,const __m128i table[4][16],std::size_t base,Emit& emit,std::index_sequence<R...>){')
    print('(emitPoint<R>(in,raw,table,base,emit),...);\n}\n};\n}')
    print(f'Exact rank16 A64/C=A^T and retained zeta/grouping validated; map sha256 {digest}', file=sys.stderr)


if __name__ == '__main__':
    main()
