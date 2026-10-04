"""Bounded floating witness search for the exact retained K18/s22 map.

This preserves the code and changes only proof witnesses. Every reported
score is a proposal, including the union formed with retained endpoints.
Selected improvements require fresh outward replay before certification.
"""
from __future__ import annotations

import argparse
from fractions import Fraction as Q
import hashlib
import json
from math import comb, isfinite, log, log2
import os
from pathlib import Path
import sys
from time import monotonic

for variable in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ[variable] = '1'

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
sys.path.insert(0, str(HERE.parent / 'k16_design'))

from flint import ctx
from scipy.special import logsumexp
import numpy as np
import packet_inner_quadratic_extension as certified
import packet_rs_k20_dense as dense
from packet_rs_length_proposal import ProposalModel, geometry
from packet_outer_geometry_proposal import upper_arrays
import packet_regional_log as logarithmic

K = 262144
MAP_SHA256 = 'e2467b9d5de8403ccd8263f0f55c73df128a8a2dd4c4cc405bc3a8e057fc3581'
RECEIPT_SHA256 = '20cd124704ff83d7ac0e098f5b80ce97e59aabf48f547bd118b1b3b957d5b60b'
RECEIPT = ROOT / 'tmp/rs16-s22-k262144-whole-fresh-p256.json'
DEFAULT_TILTS = ('.081', '.082', '.083', '.084', '.085', '.086')


def endpoint_margin(endpoint):
    mantissa, exponent = endpoint
    if type(mantissa) is not int or mantissa <= 0 or type(exponent) is not int:
        raise ValueError('positive dyadic endpoint required')
    return -log2(mantissa) - exponent


def union_margin(margins):
    values = tuple(margins)
    if not values or not all(isfinite(value) for value in values):
        raise ValueError('nonempty finite floating margins required')
    return -float(logsumexp([-value * log(2) for value in values])) / log(2)


def checked_tilts(values):
    tilts = tuple(Q(value) for value in values)
    if not tilts or len(tilts) > 20 or min(tilts) <= 0 or len(set(tilts)) != len(tilts):
        raise ValueError('one through twenty distinct positive rational tilts required')
    return tilts


def retained(path=RECEIPT):
    raw = path.read_bytes()
    if hashlib.sha256(raw).hexdigest() != RECEIPT_SHA256:
        raise ValueError('exact retained K18 receipt hash required')
    record = json.loads(raw)
    if record['map_record']['map_sha256'] != MAP_SHA256:
        raise ValueError('exact certified22 map required')
    for name, expected in record['source_sha256'].items():
        if hashlib.sha256(Path(name).read_bytes()).hexdigest() != expected:
            raise ValueError(f'retained source pin changed: {name}')
    margins = {int(q): endpoint_margin(value) for q, value in record['occupancy_uppers'].items()}
    if set(margins) != set(range(1, 2049)):
        raise ValueError('complete retained occupancy coverage required')
    return record, margins


class ActivityModel(ProposalModel):
    """One fixed birth-row selection activity; actual packet law stays15/16."""
    def __init__(self, data, record, activity=Q(1, 2), precision=192, candidates=None):
        super().__init__(data, record, precision=precision)
        self.activity = Q(activity)
        self.candidates = {} if candidates is None else candidates
        self.physical_families = {}
        self.row_choices = {}
        if not 0 <= self.activity <= 1:
            raise ValueError('row-selection activity must lie in[0,1]')

    def physical(self, tilt):
        tilt = Q(tilt)
        if tilt <= 0:
            raise ValueError('positive tilt required')
        if tilt not in self.physical_families:
            previous = ctx.prec
            try:
                ctx.prec = self.precision
                if tilt not in self.candidates:
                    self.candidates[tilt] = dense.physical_candidates(self.data, tilt)
                candidates = self.candidates[tilt]
                p = float(self.activity)
                weights = np.array([comb(16, j)*p**j*(1-p)**(16-j) for j in range(17)])
                potential = certified.kernel_t64.kernel_birth_density.continuation(
                    np.tensordot(weights, candidates['arrays'], axes=1))
                family, indices = [], []
                for matrix, options in zip(candidates['local'], candidates['rows']):
                    index = min(range(len(options)), key=lambda k:
                                sum(float(x)*v for x, v in zip(options[k], potential)))
                    value = matrix * certified.kernel_t64.aq(Q(1))
                    for target, mass in enumerate(options[index]):
                        value[0, target] = mass
                    family.append(value)
                    indices.append(index)
                self.physical_families[tilt] = family
                self.row_choices[tilt] = indices
            finally:
                ctx.prec = previous
        return self.physical_families[tilt]

    def local(self, tilt):
        tilt = Q(tilt)
        if tilt not in self.locals:
            previous = ctx.prec
            try:
                ctx.prec = self.precision
                arrays = upper_arrays(certified.kernel_t64.convolve(self.physical(tilt)))
                arrays.setflags(write=False)
                self.locals[tilt] = arrays
            finally:
                ctx.prec = previous
        return self.locals[tilt]


def export_physical(path, model, tilt, pins):
    family = model.physical(tilt)
    physical = model.data['physical_data']
    metadata = {}
    for key, value in physical.items():
        if key == 'birth_class_levels':
            metadata[key] = [int(x) for x in value]
        elif isinstance(value, np.ndarray):
            metadata[key] = dict(shape=list(value.shape), dtype=str(value.dtype))
        elif isinstance(value, (str, int, bool)):
            metadata[key] = value
        elif isinstance(value, (tuple, list)):
            metadata[key] = dict(length=len(value))
    record = dict(schema='k18-certified22-outward-physical-export-1', map_record=model.map_record,
        source_sha256=pins, tilt=str(tilt), local_selection_activity=str(model.activity),
        arithmetic='dyadic upper endpoints from192-bit Arb physical matrices', precision=192,
        birth_row_choices=model.row_choices[tilt], physical_metadata=metadata,
        matrices=[[[list(map(int, matrix[i, j].upper().man_exp()))
                    for j in range(matrix.ncols())] for i in range(matrix.nrows())] for matrix in family],
        whole_code_certificate=False)
    with path.open('x', encoding='utf-8') as handle:
        json.dump(record, handle, indent=2)
        handle.write('\n')
    print(f'Exported physical family:{path}', flush=True)


def estimate(model, qs, tilt):
    """Use the retained scaled evaluator, with its established log fallback."""
    started = monotonic()
    try:
        model.regional(K, tilt, qs[-1])
        margins = {q: -model.score(K, q, tilt) / log(2) for q in qs}
        return margins, 'scaled floating', monotonic() - started
    except FloatingPointError as error:
        failure = str(error)
    previous = ctx.prec
    try:
        ctx.prec = model.precision
        local = certified.kernel_t64.convolve(model.physical(tilt))
    finally:
        ctx.prec = previous
    regional = logarithmic.log_placement(logarithmic.upper_logs(local), qs[-1], epochs=64, windows=32)
    shape = geometry(K)
    margins = {}
    for q in qs:
        mixed = logarithmic.regional_uniform_log(regional, q)
        moment = logarithmic.log_power_matrix(mixed, shape.regions)
        margins[q] = -(log(comb(shape.group_count, q)) + q * model.log_beta
                       + float(tilt) * 52428 + moment) / log(2)
        if not isfinite(margins[q]):
            raise FloatingPointError('nonfinite logarithmic proposal')
    return margins, f'logarithmic floating; scaled failure:{failure}', monotonic() - started


def run(output, *, tilts=DEFAULT_TILTS, minimum=70, maximum=120,
        activities=('1/2',), seconds=300, receipt_path=RECEIPT):
    output = Path(output)
    tilts = checked_tilts(tilts)
    activities = tuple(Q(value) for value in activities)
    if (output.exists() or not 3 <= minimum <= 89 <= maximum <= 384 or seconds <= 0
            or not activities or len(activities) > 8 or len(set(activities)) != len(activities)
            or any(not 0 <= value <= 1 for value in activities)):
        raise ValueError('fresh output, occupancy range containing89, and bounded activities required')
    baseline, old_margins = retained(receipt_path)
    start = monotonic()
    ctx.prec = 192
    print('Preparing exact certified22 census...', flush=True)
    data, map_record = certified.prepare(22, birth_density='capped')
    if map_record['map_sha256'] != MAP_SHA256:
        raise ArithmeticError('prepared map differs from retained certified22 map')
    prepared_seconds = monotonic() - start
    print(f'Prepared exact certified22 in{prepared_seconds:.3f}s', flush=True)
    candidates = {}
    models = {activity: ActivityModel(data, map_record, activity, candidates=candidates)
              for activity in set(activities) | {Q(1, 2)}}
    base_model = models[Q(1, 2)]
    pins = dict(base_model.sources)
    pins[str(Path(__file__).resolve())] = hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    qs = tuple(range(minimum, maximum + 1))
    best, choices, trials = dict(old_margins), {}, []
    schedule = [(Q(1, 2), Q('.085'))]
    schedule.extend((activity, tilt) for activity in activities for tilt in tilts
                    if (activity, tilt) != schedule[0])
    output.parent.mkdir(parents=True, exist_ok=True)
    for index, (activity, tilt) in enumerate(schedule):
        if index and monotonic() - start >= seconds:
            break
        margins, backend, elapsed = estimate(models[activity], qs, tilt)
        if index == 0 and abs(margins[89] - old_margins[89]) > 1e-6:
            raise ArithmeticError(f'q89 baseline mismatch:{margins[89]} vs{old_margins[89]}')
        if tilt == Q('.085'):
            export_physical(output.with_name(output.stem + '-physical-a' + str(activity).replace('/', '_') + '.json'),
                            models[activity], tilt, pins)
        for q, margin in margins.items():
            if margin > best[q]:
                best[q] = margin
                choices[q] = dict(tilt=str(tilt), activity=str(activity))
        trials.append(dict(tilt=str(tilt), local_selection_activity=str(activity),
                           estimated_margin_bits=margins, backend=backend, elapsed_seconds=elapsed))
        for name, expected in pins.items():
            if hashlib.sha256(Path(name).read_bytes()).hexdigest() != expected:
                raise RuntimeError(f'proposal source changed:{name}')
        record = dict(schema='k18-certified22-margin60-proposal-1', proposal_only=True,
            whole_code_certificate=False, has_numerical_upper_endpoints=False,
            K=K, N=2*K, threshold=52428, state_bits=22, physical_t=64,
            zero_initial_state=True, final_flush=False, continuous_state=True,
            packet_activity='15/16', occupancy_values=qs,
            retained_receipt_sha256=RECEIPT_SHA256, map_record=map_record, source_sha256=pins,
            baseline_q89_bits=old_margins[89], baseline_gate_passed=True,
            estimated_combined_union_margin_bits=union_margin(best.values()),
            estimated_searched_union_margin_bits=union_margin(best[q] for q in qs),
            estimated_margin_bits={q: best[q] for q in qs}, selected_witnesses=choices,
            trials=trials, planned_trials=len(schedule), completed_trials=len(trials),
            complete_requested_search=len(trials) == len(schedule),
            seconds_budget=seconds, prepared_seconds=prepared_seconds,
            elapsed_seconds=monotonic()-start,
            scope='Floating proposals only. Combined union uses retained endpoints and floating '
                  'improvements; it is not a probability bound. Every selected improvement '
                  'requires fresh outward replay on this exact map.')
        # The output is a fresh ignored receipt owned by this process.
        with output.open('w' if index else 'x', encoding='utf-8') as handle:
            json.dump(record, handle, indent=2)
            handle.write('\n')
        worst = min(qs, key=best.get)
        print(f'tilt={tilt} activity={activity} q89={margins[89]:.9f} '
              f'combined={record["estimated_combined_union_margin_bits"]:.9f} '
              f'worst=q{worst}:{best[worst]:.6f} trial={elapsed:.2f}s '
              f'total={monotonic()-start:.2f}s', flush=True)
    return record


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--tilts', nargs='+', default=DEFAULT_TILTS)
    parser.add_argument('--activities', nargs='+', default=['1/2'])
    parser.add_argument('--min-q', type=int, default=70)
    parser.add_argument('--max-q', type=int, default=120)
    parser.add_argument('--seconds', type=float, default=300)
    args = parser.parse_args()
    run(args.output, tilts=args.tilts, minimum=args.min_q, maximum=args.max_q,
        activities=args.activities, seconds=args.seconds)


if __name__ == '__main__':
    main()
