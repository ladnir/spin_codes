"""Serial bounded FLOAT proposals for explicitly selected variable-state weights.

This reuses the frozen generic proposal helper and local state provider.
Each completed weight gets an immutable checkpoint. All omitted and pending
weights stay outside scope; even negative rows are not distance certificates.
"""

from fractions import Fraction as F
from functools import partial
from hashlib import sha256
import json
from math import isfinite
from pathlib import Path
from time import monotonic
from types import SimpleNamespace

import hd_field128_candidate_proposal as search
import hd_field128_full_rs_state_screen as state_provider


DEFAULT_WEIGHTS = (1, 2, 4, 8, 16, 32, 64)
DEFAULT_CANDIDATES = ((24, F(2, 5)), (32, F(9, 20)))


def source_pins():
    pins = state_provider.source_pins()
    for module in (search,):
        path = Path(module.__file__).resolve()
        pins[path.name] = sha256(path.read_bytes()).hexdigest()
    path = Path(__file__).resolve()
    pins[path.name] = sha256(path.read_bytes()).hexdigest()
    return dict(sorted(pins.items()))


def checked_plan(candidates, weights):
    candidates = tuple((s, F(d)) for s, d in candidates)
    weights = tuple(weights)
    if (not 1 <= len(candidates) <= 6 or len(set(candidates)) != len(candidates)
            or any(type(s) is not int or s not in state_provider.STATES or not 0 < d < F(1, 2) for s, d in candidates)
            or not 1 <= len(weights) <= 64 or len(set(weights)) != len(weights)
            or any(type(r) is not int or not 1 <= r <= 64 for r in weights)):
        raise ValueError('Use unique states/distances and unique explicitly selected weights1..64')
    return tuple((s, d, r) for s, d in candidates for r in weights)


def initial_witness(r, q=1 << 128):
    """Known field-scaled sparse proposal; no result is inherited from it."""
    alpha = F(1) if r <= 8 else F(2, 5) if r <= 17 else F(1, 5) if r <= 32 else F(1, 10)
    bits = state_provider.old.checked_q(q, construction=True)
    z = (1-F(r, 1200 if r <= 17 else 1400))**(bits//8)
    return float(z), float(alpha*F(8, bits))


def run_ladder(*, candidates=DEFAULT_CANDIDATES, weights=DEFAULT_WEIGHTS, q=1 << 128,
               evaluations=24, max_seconds=3600., max_weight_seconds=180.,
               initial_z=None, initial_alpha=None, checkpoint=None):
    plan = checked_plan(candidates, weights)
    state_provider.old.checked_q(q, construction=True)
    if (type(evaluations) is not int or not 1 <= evaluations <= 1000
            or not isfinite(max_seconds) or not 0 < max_seconds <= 86400
            or not isfinite(max_weight_seconds) or not 0 < max_weight_seconds <= 3600
            or (initial_z is None) != (initial_alpha is None)
            or (initial_z is not None and (not isfinite(initial_z) or not isfinite(initial_alpha)
                                          or not 0 < initial_z < 1 or not 0 < initial_alpha <= 1))):
        raise ValueError('Invalid bounded search resources or initial marks')
    pins, started = source_pins(), monotonic()
    results, failures, attempted = [], [], []

    def snapshot():
        covered = {row['job_index'] for row in results}
        pending = [dict(job_index=i, state=s, distance=str(d), active_groups=r)
                   for i, (s, d, r) in enumerate(plan) if i not in covered]
        return dict(schema='field128-full-rs-state-selected-ladder-float-v1',
                    plan=[dict(job_index=i, state=s, distance=str(d), active_groups=r)
                          for i, (s, d, r) in enumerate(plan)],
                    results=results.copy(), failures=failures.copy(), attempted_job_indices=attempted.copy(),
                    pending_selected_weights=pending, completed_selected_weights=len(results),
                    every_requested_weight_evaluated=not pending,
                    search_stopped_early_job_indices=[row['job_index'] for row in results
                                                     if row.get('resource_limit_reached', False)],
                    field_order=q, field_bits=q.bit_length()-1, max_evaluations=evaluations,
                    max_seconds=max_seconds, max_weight_seconds=max_weight_seconds,
                    elapsed_seconds=monotonic()-started, source_sha256=pins,
                    initial_witness_policy='explicit common marks' if initial_z is not None else 'field-scaled sparse seed per weight',
                    directed=False, whole_code_certificate=False, arithmetic='FLOAT',
                    scope='only the listed completed state/distance/weight triples; no omitted-weight or global-distance claim')

    if checkpoint:
        checkpoint(snapshot())
    for index, (state, distance, r) in enumerate(plan):
        seconds = min(max_weight_seconds, max_seconds-(monotonic()-started))
        if seconds <= 0:
            failures.append(dict(job_index=index, error_type='TimeoutError', message='Total ladder budget exhausted'))
            break
        z, alpha = initial_witness(r, q) if initial_z is None else (initial_z, initial_alpha)
        # One adapter per search; no runtime abstraction is introduced inside numeric kernels.
        provider = SimpleNamespace(weight_screen=partial(state_provider.weight_screen, state=state, q=q))
        attempted.append(index)
        try:
            result = search.proposal(provider, r, z, alpha, distance, evaluations, seconds)
            result.update(job_index=index, source_sha256=pins, initial_z=z, initial_alpha=alpha,
                          requested_state=state, requested_distance=str(distance))
            results.append(result)
        except (TimeoutError, ArithmeticError) as error:
            failures.append(dict(job_index=index, state=state, distance=str(distance), active_groups=r,
                                 error_type=type(error).__name__, message=str(error)))
        if source_pins() != pins:
            raise RuntimeError('A pinned ladder helper changed during execution')
        if checkpoint:
            checkpoint(snapshot())
    if source_pins() != pins:
        raise RuntimeError('A pinned ladder helper changed during execution')
    return snapshot()


def write_new(path, value):
    raw = (json.dumps(value, indent=2, allow_nan=False)+'\n').encode('utf-8')
    with Path(path).open('xb') as stream:
        stream.write(raw)
    return dict(path=str(Path(path).resolve()), sha256=sha256(raw).hexdigest())


def main():
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--candidate', nargs=2, action='append', metavar=('STATE', 'DISTANCE'))
    parser.add_argument('--weight', type=int, action='append')
    parser.add_argument('--field-bits', type=int, default=128)
    parser.add_argument('--max-evaluations', type=int, default=24)
    parser.add_argument('--max-seconds', type=float, default=3600.)
    parser.add_argument('--max-weight-seconds', type=float, default=180.)
    parser.add_argument('--z', type=float)
    parser.add_argument('--alpha', type=float)
    parser.add_argument('--output-dir', type=Path, required=True)
    args = parser.parse_args()
    if args.output_dir.exists():
        parser.error('Output directory must be unused')
    if not 8 <= args.field_bits <= 128 or args.field_bits%8:
        parser.error('Field bits must be a multiple of8in8..128')
    candidates = DEFAULT_CANDIDATES if args.candidate is None else tuple((int(s), F(d)) for s, d in args.candidate)
    weights = args.weight or DEFAULT_WEIGHTS
    # Validate the discrete plan before creating any artifacts.
    checked_plan(candidates, weights)
    args.output_dir.mkdir(exist_ok=False)
    sequence, written = 0, set()

    def save(snapshot):
        nonlocal sequence
        for row in snapshot['results']:
            index = row['job_index']
            if index not in written:
                write_new(args.output_dir/f'weight_{index:03d}.json', row)
                written.add(index)
        write_new(args.output_dir/f'snapshot_{sequence:03d}.json', snapshot)
        sequence += 1
        print(json.dumps(dict(completed_selected_weights=snapshot['completed_selected_weights'],
                              pending_count=len(snapshot['pending_selected_weights']),
                              latest=None if not snapshot['results'] else {key:snapshot['results'][-1][key]
                                  for key in ('requested_state', 'requested_distance', 'active_groups', 'log2_upper_screen')},
                              failures=snapshot['failures'])), flush=True)

    result = run_ladder(candidates=candidates, weights=weights, q=1 << args.field_bits,
                        evaluations=args.max_evaluations, max_seconds=args.max_seconds,
                        max_weight_seconds=args.max_weight_seconds, initial_z=args.z,
                        initial_alpha=args.alpha, checkpoint=save)
    write_new(args.output_dir/'summary.json', result)
    print(json.dumps(dict(output_dir=str(args.output_dir.resolve()),
                          completed_selected_weights=result['completed_selected_weights'],
                          pending_count=len(result['pending_selected_weights']))))


if __name__ == '__main__':
    main()
