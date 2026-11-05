"""Bounded low-arity FLOAT screen for new crossed-pair outer ensembles.

One pinned r18,F<=1 regional bank supplies all four state curves.  A
common scalar potential is optimized only for each fixed output slope;
the full regional occupancy range0..18 is checked.  The original input
column floor9 is imposed before fixed-A group powers.  No provider or
native benchmark runs, and no output is a whole-code certificate.
"""
from fractions import Fraction
from hashlib import sha256
import json
from math import log2
from pathlib import Path
from time import perf_counter

import numpy as np
from scipy.optimize import minimize_scalar

from hd_bit_syndrome_scalar_allarity import scalar_envelope
from hd_inverse_column_control import INPUT, INPUT_SHA
from hd_pair_scalar_kernel import pair_block_coefficients, support_counts
from hd_parity_outer import support_counts as parity_counts

DIRECTORY = Path(__file__).resolve().parent
SLOPES = (-10., -9., -8., -7., -6., -5.5, -5., -4.75, -4.70125712,
          -4.5, -4., -3., -2., -1., 0.)


def activity_tail_log(blocks, *, minimum=9):
    """Positive degree convolution in log space; do not subtract low tails."""
    blocks = tuple(np.asarray(block, dtype=float) for block in blocks)
    if (not blocks or any(block.shape != (3,) or np.any(block < 0)
            or not np.all(np.isfinite(block)) or block[0] != 1 for block in blocks)
            or type(minimum) is not int or not 0 <= minimum <= 2*len(blocks)):
        raise ValueError('Positive degree0..2 blocks with constant1 and valid floor required')
    polynomial = np.array((0.,))
    for block in blocks:
        log_block = np.full(3, -np.inf)
        positive = block > 0
        log_block[positive] = np.log2(block[positive])
        following = np.full(len(polynomial)+2, -np.inf)
        for degree in range(3):
            following[degree:degree+len(polynomial)] = np.logaddexp2(
                following[degree:degree+len(polynomial)], polynomial+log_block[degree])
        polynomial = following
    return float(np.logaddexp2.reduce(polynomial[minimum:]))


def load_bank():
    raw = (DIRECTORY/INPUT).read_bytes()
    if sha256(raw).hexdigest() != INPUT_SHA:
        raise ValueError('Pinned regional evidence changed')
    saved = json.loads(raw)
    if (saved['groups'] != 18 or saved['source_class'] != 'Fmax<=1'
            or saved['mark'] != '9927/10000'):
        raise ValueError('Unexpected selected source scope')
    curves = np.array([[float('-inf') if value is None else value for value in row]
                       for row in saved['regional_log2_curves']])
    prefix = log2(int(saved['source_count']))-734003*log2(Fraction(saved['mark']))-log2(255)
    if abs(prefix-saved['common_prefix_log2']) > 1e-8:
        raise ValueError('Source or output-mark normalization changed')
    return saved, curves, prefix


def run(*, max_seconds=60, progress=None):
    start = perf_counter()
    def check_time():
        if perf_counter()-start > max_seconds:
            raise TimeoutError('Bounded pair scalar screen exceeded its budget')
    saved, curves, prefix = load_bank()
    controls = []
    # These H0 arities are complete after multiplicity and class summation.
    # Higher H* cases are only representative fixed-A contributions.
    specifications = [('H0', (0,), 16), ('H0', (0, 1), 56), ('H0', (0, 8), 64),
                      ('Hstar', (0,), 1), ('Hstar', (0, 1), 1),
                      ('Hstar', (0, 1, 2, 3), 1),
                      ('Hstar', (0, 1, 2, 3, 4), 1),
                      ('Hstar', (0, 1, 2, 4, 8), 1)]
    counts = []
    for mixer, source, multiplicity in specifications:
        check_time()
        if mixer == 'H0':
            sparse = parity_counts(sum(1 << i for i in source))
            vector = tuple(sparse.get(mask, 0) for mask in range(65536))
        else:
            vector = support_counts(source, check_time=check_time)
        counts.append(vector)
        if progress:
            progress(dict(stage='exact_support_counts', mixer=mixer, source=source,
                          support_patterns=sum(bool(value) for value in vector)))
    for slope in SLOPES:
        check_time()
        objective = lambda p: 256*scalar_envelope(curves, float(p), slope)['log2_eta']-p
        fit = minimize_scalar(objective, bounds=(-160., 0.), method='bounded',
                              options=dict(xatol=1e-9))
        choices = [(float(fit.fun), float(fit.x)), (objective(0.), 0.), (objective(-160.), -160.)]
        _, potential = min(choices)
        envelope = scalar_envelope(curves, potential, slope)
        common = prefix+256*envelope['log2_eta']+envelope['endpoint_log2']
        rows = []
        for (mixer, source, multiplicity), vector in zip(specifications, counts):
            check_time()
            block = pair_block_coefficients(vector, 2.**slope)['coefficients']
            tail = activity_tail_log(tuple(block.values()))
            group = tail-8*len(source)*log2(255)
            rows.append(dict(mixer=mixer, source=source, labelled_multiplicity=multiplicity,
                             group_log2=group,
                             homogeneous_with_class_multiplicity_log2=common+18*(group+log2(multiplicity))))
        one = common+18*(rows[0]['group_log2']+4)
        two = common+18*float(np.logaddexp2(rows[1]['group_log2']+log2(56),
                                          rows[2]['group_log2']+log2(64)))
        controls.append(dict(log2_t=slope, log2_potential=potential, log2_eta=envelope['log2_eta'],
                             contacts=envelope['contacts'], common_prefix_log2=common,
                             complete_H0_arity1_log2=one, complete_H0_arity2_log2=two,
                             representative_rows=rows))
        if progress:
            progress(dict(stage='scalar_slope', log2_t=slope, H0_a1=one, H0_a2=two))
    allarity_floor = prefix+18*128*log2(253)+256*curves[3, 18]
    pins = {name: sha256((DIRECTORY/name).read_bytes()).hexdigest() for name in
            ('hd_pair_scalar_control.py', 'hd_pair_scalar_kernel.py', 'hd_parity_outer.py',
             'hd_bit_syndrome_scalar_allarity.py')}
    return dict(input_receipt=INPUT, input_sha256=INPUT_SHA, groups=18, source_class='Fmax<=1',
                source_count=saved['source_count'], mark=saved['mark'], output_limit=734003,
                preferred_single_layer_deltas=list(range(8,16)), controls=controls,
                scalar_allarity_positive_branch_floor_log2=float(allarity_floor),
                best_H0_arity1=min(controls,key=lambda row:row['complete_H0_arity1_log2']),
                best_H0_arity2=min(controls,key=lambda row:row['complete_H0_arity2_log2']),
                full_scope_for_H0_arities=[1,2], Hstar_representatives_only=True,
                slopes_prescribed=True, potential_fits=len(SLOPES), regional_provider_runs=0,
                directed=False, whole_code_certificate=False, implementation_measured=False,
                source_pins=pins, elapsed_seconds=perf_counter()-start)


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        parser.error('Refuse to overwrite evidence')
    result = run(progress=lambda item: print(json.dumps(item),flush=True))
    with args.output.open('x',encoding='utf-8') as handle:
        json.dump(result, handle, indent=2, allow_nan=False)
        handle.write('\n')
    print(json.dumps(dict(elapsed_seconds=result['elapsed_seconds'],
        H0_a1=result['best_H0_arity1']['complete_H0_arity1_log2'],
        H0_a2=result['best_H0_arity2']['complete_H0_arity2_log2'],
        allarity_floor=result['scalar_allarity_positive_branch_floor_log2'])),flush=True)
