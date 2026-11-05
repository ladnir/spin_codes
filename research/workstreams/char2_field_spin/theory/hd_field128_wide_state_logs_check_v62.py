"""Remote-only all-row check of fast V62 logs against exact V60 entries.

All seven retained marks and all three states are checked. This validates a
FLOAT search evaluator, not directed inequalities or a distance certificate.
"""
import argparse
from fractions import Fraction as F
from hashlib import sha256
import json
from pathlib import Path
import resource
import shutil
from time import monotonic

import numpy as np
import hd_field128_wide_state_logs_v62 as fast

HERE = Path(__file__).resolve().parent
INPUT = 'hd_field128_reversed_dual_marker_v52_20261010.json'
INPUT_SHA = '1ec14fce6aede5f7f45d301fdeae5cb7057670a09f8e244bad52e9bd1b69b3e5'
LIMIT, RESERVE = 32 << 10, 20 << 30


def run(name):
    if HERE != Path('/tmp/spin-field128-proof-52GBmlZD/bundle/theory'):
        raise RuntimeError('Numerical checks are remote-only')
    if Path(name).name != name or not name.endswith('.json'):
        raise ValueError('Expected output basename')
    output = HERE/name
    if output.exists() or shutil.disk_usage(HERE).free-LIMIT < RESERVE:
        raise RuntimeError('Exclusive output and free-space reserve required')
    resource.setrlimit(resource.RLIMIT_AS, (2 << 30, 2 << 30))
    resource.setrlimit(resource.RLIMIT_CPU, (230, 235))
    resource.setrlimit(resource.RLIMIT_FSIZE, (LIMIT, LIMIT))
    raw = (HERE/INPUT).read_bytes()
    if len(raw) > 96 << 10 or sha256(raw).hexdigest() != INPUT_SHA:
        raise RuntimeError('Pinned witness input changed')
    specs = json.loads(raw)['witness_specs']
    if len(specs) != 7:
        raise RuntimeError('Expected all seven retained marks')
    pins, started, rows = fast.source_pins(), monotonic(), []
    for state in fast.STATES:
        for index, spec in enumerate(specs):
            z = F(spec['z'])
            # The serialized witnesses are exact binary64 values; no mark
            # approximation should enter this exact-versus-FLOAT control.
            if F(float(z)) != z:
                raise ArithmeticError('Non-binary64 mark')
            actual = fast.local_logs(float(z), state)
            expected = fast.exact.local_logs(z, state)
            if actual.shape != (65, 2, 2) or expected.shape != actual.shape:
                raise ArithmeticError('Incomplete local domain')
            zeros = np.isneginf(expected)
            if not np.array_equal(zeros, np.isneginf(actual)):
                raise ArithmeticError('Exact zero mask differs')
            if not np.all(np.isfinite(actual[~zeros])):
                raise ArithmeticError('Positive entry became nonfinite')
            error = float(np.max(np.abs(actual[~zeros]-expected[~zeros])))
            if error > 2e-9:
                raise ArithmeticError('FLOAT local-cap disagreement')
            rows.append(dict(state=state, witness=index, max_log2_error=error,
                finite_entries=int(np.count_nonzero(~zeros)),
                exact_log_sha256=sha256(expected.tobytes()).hexdigest(),
                fast_log_sha256=sha256(actual.tobytes()).hexdigest()))
        print(json.dumps(dict(checked_state=state,
                              elapsed_seconds=monotonic()-started)), flush=True)
    fast.check_source_pins(pins)
    if (HERE/INPUT).read_bytes() != raw:
        raise RuntimeError('Witness input changed')
    pins[Path(__file__).name] = sha256(Path(__file__).read_bytes()).hexdigest()
    result = dict(schema='wide-state-logs-check-v62', status='PASS',
        source_sha256=dict(sorted(pins.items())), input=dict(path=INPUT, sha256=INPUT_SHA),
        rows=rows, packets_checked=65*21, scalar_entries_checked=65*21*4,
        all_rows_all_retained_marks_all_states=True,
        max_log2_error=max(row['max_log2_error'] for row in rows),
        directed=False, whole_code_certificate=False, elapsed_seconds=monotonic()-started)
    data = (json.dumps(result, sort_keys=True, separators=(',', ':'), allow_nan=False)+'\n').encode()
    if len(data) > LIMIT or shutil.disk_usage(HERE).free-len(data) < RESERVE:
        raise RuntimeError('Output storage cap')
    with output.open('xb') as stream:
        stream.write(data)
    return dict(path=name, bytes=len(data), sha256=sha256(data).hexdigest(),
        max_log2_error=result['max_log2_error'], elapsed_seconds=result['elapsed_seconds'])


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', required=True)
    print(json.dumps(run(parser.parse_args().output)), flush=True)
