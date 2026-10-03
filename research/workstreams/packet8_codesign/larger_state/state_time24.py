"""Expected source-state residence times in the positive comparison expression."""
import argparse
import hashlib
import json
from pathlib import Path
from time import monotonic

import occupancy24
import numpy as np


def run(args):
    if args.output.exists():
        raise ValueError('fresh output required')
    started = monotonic()
    receipt = json.loads(args.input.read_text(encoding='utf-8'))
    if not receipt['source_pins_verified_at_finish']:
        raise ValueError('complete input required')
    for name, expected in receipt['source_sha256'].items():
        if hashlib.sha256(Path(name).read_bytes()).hexdigest() != expected:
            raise ValueError(f'stale source: {name}')
    local, q, delta = np.asarray(receipt['local_operator_matrices']), receipt['q'], 1e-4
    selectors = dict(zero_source=(slice(None), 0, slice(None)),
        nonzero_source=(slice(None), slice(1, None), slice(None)),
        empty_zero_source=(0, 0, slice(None)),
        occupied_zero_source=(slice(1, None), 0, slice(None)))
    counts = {}
    for label, selector in selectors.items():
        values = []
        for sign in (-1, 1):
            changed = local.copy()
            changed[selector] *= np.exp(sign*delta)
            values.append(occupancy24.moment(changed, q)[0])
        counts[label] = float((values[1]-values[0])/(2*delta))
    if (abs(counts['zero_source']+counts['nonzero_source']-2048) > 2e-4
            or abs(counts['zero_source']-counts['empty_zero_source']-counts['occupied_zero_source']) > 2e-4):
        raise ArithmeticError('state residence checks failed')
    result = dict(schema='actual24-bit-comparison-residence-1', q=q, tilt=receipt['tilt'],
        proposal_only=True, whole_code_certificate=False, actual_failure_distribution=False,
        input_receipt=str(args.input.resolve()), input_receipt_sha256=hashlib.sha256(args.input.read_bytes()).hexdigest(),
        input_source_sha256=receipt['source_sha256'],
        source_sha256={str(path.resolve()): hashlib.sha256(path.read_bytes()).hexdigest()
                      for path in (Path(__file__), Path(occupancy24.__file__))},
        finite_difference_delta=delta, expected_steps=counts, elapsed_seconds=monotonic()-started)
    args.output.write_text(json.dumps(result, indent=2)+'\n', encoding='utf-8')
    print(json.dumps(result, indent=2), flush=True)
    return result


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    run(parser.parse_args())
