"""Compare imported RM2Sub tables with the selected grid map, using exact integers.

This checks constituent identity and spectra, not the complete encoder or margin.
The C++ correctness tests separately check the optimized kernel against its tables.
"""
import argparse
from collections import Counter
import hashlib
import json
import math
from pathlib import Path
import re


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
IMPLEMENTATION = ROOT / ('constructions/'
    'riffle_parityfanout31x33_bchperm_transpose_bitshuffle_splitstate_preaddmul_'
    'rm2sub_t128_s19/implementation')
TABLE = IMPLEMENTATION / 'libOTe/Tools/RiffleCode/generated/Rm2SubS19Tables.h'
SNAPSHOT = HERE / 'generated/larger_state_inputs_v1'


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def array(text, name, count):
    matches = re.findall(r'\b' + re.escape(name) + r'\s*\{([^{}]*)\}', text)
    if len(matches) != 1:
        raise ValueError(f'Expected one literal array: {name}')
    tokens = [v.strip() for v in matches[0].split(',') if v.strip()]
    values = [int(v, 0) for v in tokens]
    if len(values) != count:
        raise ValueError(f'Wrong array length: {name}')
    return values


def generators(columns, state_bits):
    if any(not 0 < c < (1 << state_bits) for c in columns):
        raise ValueError('Column out of range')
    return [sum(((c >> j) & 1) << i for i, c in enumerate(columns))
            for j in range(state_bits)]


def spectrum(rows):
    pivots = {}
    for word in rows:
        while word:
            p = word.bit_length() - 1
            if p not in pivots:
                pivots[p] = word
                break
            word ^= pivots[p]
    if len(pivots) != len(rows):
        raise ValueError('Dependent generator rows')
    counts = Counter({0: 1})
    word = 0
    for index in range(1, 1 << len(rows)):
        word ^= rows[(index & -index).bit_length() - 1]
        counts[word.bit_count()] += 1
    return dict(sorted(counts.items()))


def dual_spectrum(counts, length, dimension):
    result = {}
    for degree in range(length + 1):
        total = sum(count * sum((-1)**j * math.comb(w, j)
            * math.comb(length-w, degree-j)
            for j in range(max(0, degree-length+w), min(degree, w)+1))
            for w, count in counts.items())
        value, remainder = divmod(total, 1 << dimension)
        if remainder or value < 0:
            raise ValueError('Invalid integral dual spectrum')
        if value:
            result[degree] = value
    if sum(result.values()) != 1 << (length-dimension):
        raise ValueError('Incorrect dual mass')
    return result


def audit(table=TABLE, snapshot=SNAPSHOT):
    manifest_path = snapshot / 'manifest.json'
    manifest = json.loads(manifest_path.read_text())
    names = [f't128_s19_{suffix}.json'
             for suffix in ('selection', 'a_spectrum', 'b_kernel_spectrum')]
    for name in names:
        if sha(snapshot / name) != manifest['snapshot_sha256'][name]:
            raise ValueError(f'Snapshot hash mismatch: {name}')
    selected = json.loads((snapshot / names[0]).read_text())
    if (selected['parameters']['step_bits'], selected['parameters']['state_bits']) != (128, 19):
        raise ValueError('Unexpected selected geometry')
    selected_columns = [int(c, 16) for c in selected['selected']['B_columns_hex']]
    imported_columns = array(table.read_text(), 'Rm2Sub19Columns', 128)
    if len(selected_columns) != 128:
        raise ValueError('Unexpected selected column count')
    imported_rows = generators(imported_columns, 19)
    selected_rows = generators(selected_columns, 19)
    if selected_rows != [int(v, 16) for v in selected['selected']['A_generator_words_hex']]:
        raise ValueError('Selected A and B are not transposes')
    imported_a, selected_a = spectrum(imported_rows), spectrum(selected_rows)
    claimed_a = {r['weight']: r['count'] for r in
                 json.loads((snapshot / names[1]).read_text())['spectrum'] if r['count']}
    if selected_a != claimed_a:
        raise ValueError('Selected spectrum mismatch')
    imported_kernel = dual_spectrum(imported_a, 128, 19)
    selected_kernel = dual_spectrum(selected_a, 128, 19)
    claimed_kernel = {r['total_weight']: r['kernel_words'] for r in
        json.loads((snapshot / names[2]).read_text())['by_total_weight'] if r['kernel_words']}
    if selected_kernel != claimed_kernel:
        raise ValueError('Selected kernel spectrum mismatch')
    mismatches = [i for i, (a, b) in enumerate(zip(imported_columns, selected_columns)) if a != b]
    first = mismatches[0] if mismatches else None
    return dict(
        status='EXACT_CONSTITUENT_COMPARISON_ONLY', configuration='t128_s19',
        ordered_columns_equal=not mismatches, differing_columns=len(mismatches),
        first_mismatch=None if first is None else dict(index=first,
            imported=hex(imported_columns[first]), selected=hex(selected_columns[first])),
        a_spectra_equal=imported_a == selected_a,
        kernel_spectra_equal=imported_kernel == selected_kernel,
        imported_A_spectrum=imported_a, selected_A_spectrum=selected_a,
        imported_minimum_A_weight=min(w for w in imported_a if w),
        selected_minimum_A_weight=min(w for w in selected_a if w),
        imported_minimum_kernel_weight=min(w for w in imported_kernel if w),
        selected_minimum_kernel_weight=min(w for w in selected_kernel if w),
        imported_BA_zero=all(not ((a & b).bit_count() & 1)
                             for a in imported_rows for b in imported_rows),
        full_encoder_equivalence_verified=False, full_margin_bits=None,
        source_sha256={str(p.resolve()): sha(p) for p in
                       [Path(__file__), table, manifest_path, *[snapshot / n for n in names]]},
        limitations=[
            'No equivalence search under coordinate permutations or basis changes.',
            'No audit of the complete outer, route, recurrence, or setup distribution.',
            'Constituent spectra alone do not establish a complete margin.'])


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    report = audit()
    rendered = json.dumps(report, indent=2) + '\n'
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        with args.output.open('x', encoding='utf-8') as stream:
            stream.write(rendered)
    print(rendered, end='')
