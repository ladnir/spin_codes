"""One additional frozen dense-A8 witness; new affine arithmetic throughout."""
from fractions import Fraction as F
from hashlib import sha256
import json
from pathlib import Path

import numpy as np

from hd_two_sum_affine_dense import run, ARITY_CLASSES, CONSTRUCTION


DIRECTORY = Path(__file__).resolve().parent
INITIAL_NAME = 'hd_two_sum_affine_dense_8_0_inherited.json'
# The initial artifact is checked below against a source-derived scope and
# externally recorded byte pin; it is FLOAT proposal evidence, not a proof.
Z = F('8970988817957597/9007199254740992')
LOGP = F('-2266954815/33554432')


def combine(initial, additional):
    expected_masks = {a: list(masks) for a, masks in ARITY_CLASSES}
    banks = []
    for bank in (initial, additional):
        if (bank['construction'] != CONSTRUCTION or bank['cohort_counts'] != [0, 16, 1]
                or bank['source_matching_counts'] != [8, 0] or bank['labels'] != list(range(8, 16))
                or bank['field_order'] != 256 or bank['forbidden_weight'] != 734003
                or bank['total_columns'] != 32 or bank['sparse_pair_assignments'] != 28**16
                or bank['source_colorings'] != 17 or bank['sparse_support_floor'] != 544
                or bank['dense_nonempty_column_floor'] != 17 or bank['sparse_normalizer_exponent'] != 512
                or bank['old_xor_mask_counts_used'] is not False or bank['directed'] is not False):
            raise ValueError('Only the exact affine h1/source8,0 FLOAT scope may be combined')
        rows = {row['active_constituents']: row for row in bank['entries']}
        if len(rows) != len(bank['entries']) or set(rows) != set(expected_masks):
            raise ValueError('The complete six-arity partition is required')
        if any(rows[a]['active_row_masks'] != masks for a, masks in expected_masks.items()):
            raise ValueError('Fixed row-set membership changed')
        banks.append(rows)
    selected = []
    for arity in sorted(expected_masks):
        index = min(range(2), key=lambda i: banks[i][arity]['failure_log2_screen'])
        row = dict(banks[index][arity])
        bank = (initial, additional)[index]
        row.update(witness_label=bank['witness_label'], mark=bank['mark'],
                   log2_potential=bank['log2_potential'])
        selected.append(row)
    return dict(construction=CONSTRUCTION, cohort_counts=[0, 16, 1], source_matching_counts=[8, 0],
        fixed_row_sets=219, deterministic_classes=6, selected=selected,
        failure_log2_screen=float(np.logaddexp2.reduce([row['failure_log2_screen'] for row in selected])),
        setup_dependent_selection=False, common_witness_within_each_arity=True,
        extra_witness_bank_factor=False, old_xor_mask_counts_used=False,
        complete_selected_cohort_and_source=True, directed=False, whole_code_certificate=False)


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        parser.error('Do not overwrite a bounded proposal artifact')
    raw = (DIRECTORY/INITIAL_NAME).read_bytes()
    initial = json.loads(raw)
    additional = run(mark=Z, log_potential=LOGP, witness_label='inherited-old-dense-a8',
        progress=lambda row: print(json.dumps(row, allow_nan=False), flush=True))
    result = dict(initial_filename=INITIAL_NAME, initial_sha256=sha256(raw).hexdigest(),
        additional=additional, combined=combine(initial, additional),
        evaluations=2, fitted_new_witnesses=0, reused_literals_only=True,
        source_sha256=sha256(Path(__file__).read_bytes()).hexdigest(),
        evaluator_sha256=sha256(Path(__file__).with_name('hd_two_sum_affine_dense.py').read_bytes()).hexdigest())
    with args.output.open('x', encoding='utf-8') as target:
        json.dump(result, target, indent=2, allow_nan=False)
        target.write('\n')
    print(json.dumps({key: value for key, value in result.items() if key != 'additional'}, allow_nan=False))
