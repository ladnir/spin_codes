"""Single original-weight913 directed replay; no saved numerical caps reused."""
import os
os.environ['OPENBLAS_NUM_THREADS'] = '1'

import argparse
from contextlib import contextmanager
from fractions import Fraction as F
from hashlib import sha256
import importlib
import json
from pathlib import Path
import re
import shutil
import signal
import sys
from time import monotonic

HERE = Path(__file__).resolve().parent
NAME = Path(__file__).name
CONTROLLER = 'hd_field128_gamma_portable_band_runner_v25'
CONTROLLER_SHA256 = '5d22c83e3419b1dc2a484b09c9cb2ada17971bd48195987140c9158942d7ef1c'
CONSTRUCTION = 'toeplitz-fullfield-full-rs64-30-gamma16-feedback-alternating-swapped-pairs-k20-rate-half-v1'
WEIGHT, Q, MAX_NODES, MAX_SECONDS = 913, 1 << 128, 4096, 600
MAX_OUTPUT_BYTES, MAX_TABLE_BYTES, MEMORY_BYTES, RESERVE = 32 << 20, 64 << 20, 4 << 30, 20 << 30


def relative_name(value, *, json_file=False, basename=False):
    if (not isinstance(value, str) or '\\' in value or (basename and '/' in value)
            or any(re.fullmatch(r'[A-Za-z0-9_][A-Za-z0-9_.-]*', part) is None or part.endswith('.')
                   for part in value.split('/')) or (json_file and not value.endswith('.json'))):
        raise ValueError('Expected a portable path inside theory')
    (HERE/value).resolve().relative_to(HERE.resolve())
    return value


def digest(value):
    if not isinstance(value, str) or re.fullmatch(r'[0-9a-f]{64}', value) is None:
        raise ValueError('Expected lowercase SHA256')
    return value


def load_proposal(name, expected):
    name = relative_name(name, json_file=True)
    raw = (HERE/name).read_bytes()
    if sha256(raw).hexdigest() != digest(expected):
        raise ValueError('Exact proposal byte pin mismatch')
    proposal = json.loads(raw)
    if (proposal.get('active_groups') != WEIGHT or proposal.get('state_dimension') != 30
            or proposal.get('target_distance') != '2/5' or proposal.get('construction') != CONSTRUCTION
            or proposal.get('field_order', Q) != Q):
        raise ValueError('Expected original913 full-field gamma30/d40 proposal')
    witness = proposal['witness']
    markers = witness['markers']
    if (witness.get('grid') != 128 or witness.get('state_dimension') != 30 or len(markers) != 83
            or proposal.get('marker_count') != 83
            or any(not isinstance(pair, list) or len(pair) != 2 for pair in markers)):
        raise ValueError('Expected exactly83 grid128 positive marker pairs')
    for value in (witness['z'], witness['alpha'], *(x for pair in markers for x in pair)):
        if not isinstance(value, str) or str(F(value)) != value or F(value) <= 0:
            raise ValueError('Expected canonical positive rational literals')
    z, alpha = F(witness['z']), F(witness['alpha'])
    if not z < 1 or not alpha <= 1:
        raise ValueError('Invalid exact z/alpha domain')
    return proposal, (z, alpha, [tuple(map(F, pair)) for pair in markers]), dict(path=name, sha256=expected)


def helper():
    if sha256((HERE/(CONTROLLER+'.py')).read_bytes()).hexdigest() != CONTROLLER_SHA256:
        raise RuntimeError('Reviewed compact controller changed')
    return importlib.import_module(CONTROLLER)


def source_pins(proposal=None, controller=None):
    controller = helper() if controller is None else controller
    pins = controller.source_pins()
    for name, value in ({} if proposal is None else proposal['source_sha256']).items():
        digest(value)
        if name in pins and pins[name] != value:
            raise RuntimeError('Source closure conflict: '+name)
        pins[name] = value
    pins[NAME] = sha256((HERE/NAME).read_bytes()).hexdigest()
    return dict(sorted(pins.items()))


def check_pins(pins, witness_source):
    for name, expected in pins.items():
        if sha256((HERE/name).read_bytes()).hexdigest() != expected:
            raise RuntimeError('Source bytes changed: '+name)
    if sha256((HERE/witness_source['path']).read_bytes()).hexdigest() != witness_source['sha256']:
        raise RuntimeError('Exact proposal changed')


def reserve_guard(planned):
    if shutil.disk_usage(HERE).free-planned < RESERVE:
        raise RuntimeError('Resource stop: planned writes breach20GiB reserve')


@contextmanager
def linux_limits():
    if sys.platform != 'linux':
        raise RuntimeError('Actual replay requires Linux; local tests mock this context')
    import resource
    soft, hard = resource.getrlimit(resource.RLIMIT_AS)
    limit = min([MEMORY_BYTES]+[x for x in (soft, hard) if x != resource.RLIM_INFINITY])
    resource.setrlimit(resource.RLIMIT_AS, (limit, limit))
    previous = signal.getsignal(signal.SIGALRM)
    signal.signal(signal.SIGALRM, lambda unused_signal, unused_frame: os._exit(124))
    signal.setitimer(signal.ITIMER_REAL, MAX_SECONDS)
    try:
        yield limit
    finally:
        signal.setitimer(signal.ITIMER_REAL, 0)
        signal.signal(signal.SIGALRM, previous)


def run(proposal_name, expected_sha256, output_dir, *, cell_failure_bits=58, max_nodes=MAX_NODES):
    if type(cell_failure_bits) is not int or cell_failure_bits not in (54, 58, 64):
        raise ValueError('Choose explicit cell target54,58,or64; aggregate target stays54')
    if type(max_nodes) is not int or not 1 <= max_nodes <= MAX_NODES:
        raise ValueError('Expected at most4096 evaluated nodes')
    destination = HERE/relative_name(output_dir, basename=True)
    proposal, (z, alpha, markers), witness_source = load_proposal(proposal_name, expected_sha256)
    reserve_guard(MAX_OUTPUT_BYTES)
    started, bank = monotonic(), None
    with linux_limits() as memory_limit:
        controller = helper()
        pins = source_pins(proposal, controller)
        check_pins(pins, witness_source)
        destination.mkdir(exist_ok=False)
        outputs = controller.OutputBudget(destination)
        written = []
        def save(name, value, *, phase_check=True):
            if phase_check:
                check_pins(pins, witness_source)
            path = destination/name
            outputs.write_new(path, value)
            written.append((path, sha256(path.read_bytes()).hexdigest()))
        def remaining(cap=MAX_SECONDS):
            seconds = min(cap, MAX_SECONDS-5-(monotonic()-started))
            if seconds <= 0:
                raise TimeoutError('Whole-job budget exhausted')
            return seconds
        try:
            save('driver.json', dict(schema='gamma-r913-directed-driver-v25',
                driver=dict(path=NAME, sha256=pins[NAME]), source_sha256=pins,
                witness_source=witness_source, requested_weight_range=[913, 913],
                numerical_provider=Path(controller.bands.__file__).name,
                partition_controller=CONTROLLER+'.py', controller_sha256=CONTROLLER_SHA256,
                construction=CONSTRUCTION, field_order=Q, state_dimension=30, target_distance='2/5',
                marker_count=83, precision=256, coefficient_bits=80, grid=128,
                cell_failure_bits=cell_failure_bits, aggregate_failure_bits=54,
                max_nodes=max_nodes, whole_job_seconds=MAX_SECONDS, max_table_bytes=MAX_TABLE_BYTES,
                max_output_bytes=MAX_OUTPUT_BYTES, disk_reserve_bytes=RESERVE,
                memory_address_space_bytes=memory_limit, required_normal_process_exit=True,
                fresh_bank_required=True, fresh_bank_rebuilt=False, old_numerical_bounds_reused=False,
                directed=False, whole_code_certificate=False))
            # No progress callback: stdout is only one compact final status row.
            bank = controller.cert.AliasBank(z, alpha, markers, state=30, q=Q, distance=F(2, 5),
                grid=128, precision=256, bits=80, outer_measure='invertible-toeplitz',
                max_seconds=remaining(), max_table_bytes=MAX_TABLE_BYTES)
            bank.check_pins()
            save('bank.json', dict(bank.metadata(), witness_source=witness_source, runner_source_sha256=pins,
                runtime=controller.portable.runtime_metadata(), fresh_bank_rebuilt=True,
                directed=False, whole_code_certificate=False))
            evaluator = controller.bands.Evaluator(bank, 913, 913, max_seconds=remaining(180.),
                witness_source=witness_source)
            result = controller.cover(evaluator, max_nodes=max_nodes, max_seconds=remaining(),
                cell_failure_bits=cell_failure_bits, aggregate_bits=54,
                # Source closure checks surround the entire pinned cover; avoid
                # rehashing 93 files in its checkpoint hot path.
                checkpoint=lambda row: save(f"node_{row['sequence']:06d}.json", row, phase_check=False),
                context_checkpoint=lambda row: save('context.json', row))
            bank.check_pins()
            result.update(execution_driver=dict(path=NAME, sha256=pins[NAME]), execution_source_sha256=pins,
                required_normal_process_exit=True, total_job_seconds=monotonic()-started)
            save('summary.json', result)
            return result
        finally:
            try:
                check_pins(pins, witness_source)
                for path, expected in written:
                    if sha256(path.read_bytes()).hexdigest() != expected:
                        raise RuntimeError('New output bytes changed: '+path.name)
            finally:
                if bank is not None:
                    bank.check_pins()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--proposal', required=True)
    parser.add_argument('--expected-sha256', required=True)
    parser.add_argument('--output-dir', required=True)
    parser.add_argument('--cell-failure-bits', type=int, choices=(54, 58, 64), default=58)
    parser.add_argument('--max-nodes', type=int, default=MAX_NODES)
    args = parser.parse_args()
    result = run(args.proposal, args.expected_sha256, args.output_dir,
                 cell_failure_bits=args.cell_failure_bits, max_nodes=args.max_nodes)
    print(json.dumps({key: result[key] for key in ('status', 'requested_weight_range',
        'evaluated_nodes', 'log2_upper_rational', 'meets_requested_band_target', 'stop_reason')}), flush=True)
    return 0 if (result['complete_requested_band'] and result['directed']
                 and result['meets_requested_band_target']) else 2


if __name__ == '__main__':
    raise SystemExit(main())
