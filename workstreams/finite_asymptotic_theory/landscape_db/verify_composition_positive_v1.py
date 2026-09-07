"""Replay selected completed BCH composition witnesses through both kernels."""
import argparse
import json
import math
from pathlib import Path

import numpy as np

import occupation_composition_lazy_v1 as baseline
import occupation_composition_positive_v1 as positive
import study_bch_dominance_v1 as study

HERE = Path(__file__).resolve().parent


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--cover', type=Path, required=True)
    args = parser.parse_args()
    source = args.cover.resolve()
    data = json.loads(source.read_text())
    if data['status'] != 'BINARY64_COMPLETE_SPARSE_INTERVAL':
        raise ValueError('a completed sparse interval is required')
    for name, digest in data['source_sha256'].items():
        if study.sha(study.ROOT/name) != digest:
            raise ValueError('changed source dependency')
    _, spectra, _, _ = study.load_inputs()
    config = data['arguments']; block, exponent = config['block'], config['exponent']
    length = (1 << exponent)//(block//2); cutoff = block*length//10
    maximum = data['occupation_max']; counts = spectra[block]
    models = [(baseline.CompositionBoxes(counts, block, data['bands'], p, maximum),
               positive.CompositionBoxes(counts, block, data['bands'], p, maximum))
              for p in data['probability_banks']]
    tilts = np.unique(np.r_[np.arange(25, 226)/25-math.log(length), np.arange(-30, 6)/5])
    cache = {}; inputs = [source, Path(__file__), Path(positive.__file__), Path(baseline.__file__)]
    checks = []
    for row in data['occupations']:
        q = row['occupation']
        if q not in (5, 64, 128, 256, 512) or row['composition_cover'] is None:
            continue
        boxes = row['composition_cover']['boxes']
        indices = sorted({0, len(boxes)//2, max(range(len(boxes)), key=lambda j: boxes[j]['log_bound'])})
        for index in indices:
            box = boxes[index]; witness = box['witness']; tilt = witness['log_tilt']
            if tilt not in cache:
                hits = np.flatnonzero(tilts == tilt)
                if len(hits) != 1:
                    raise ValueError('unknown recorded tilt')
                path = source.parent/f'tilt_{int(hits[0])}.json'
                cached = json.loads(path.read_text())
                if cached['input_fingerprint'] != data['input_fingerprint'] or cached['log_tilt'] != tilt:
                    raise ValueError('incompatible cached region')
                cache[tilt] = np.array(cached['regions']); inputs.append(path)
            slow, fast = models[witness['probability_bank']]
            arguments = (cache[tilt], box['lower'], box['upper'], q, length, cutoff, math.exp(tilt))
            expected, actual = slow.bound(*arguments), fast.bound(*arguments)
            error = abs(expected-actual)
            if abs(expected-box['log_bound']) > 2e-7 or error > 2e-7:
                raise ArithmeticError('positive composition replay mismatch')
            check = dict(occupation=q, lower=box['lower'], upper=box['upper'], absolute_log_error=error)
            checks.append(check); print(check, flush=True)
    if not checks:
        raise ValueError('no requested composition witnesses in this cover')
    result = dict(status='SELECTED_POSITIVE_COMPOSITION_REPLAY', checks=checks,
                  maximum_log_error=max(c['absolute_log_error'] for c in checks),
                  source_sha256={str(p): study.sha(p) for p in inputs},
                  limitations=['Selected witnesses; the full verifier must replay every new retained bound.',
                               'No outward arithmetic and no performance measurement.'])
    (HERE/'bch_positive_composition_replay_v1.json').write_text(json.dumps(result, indent=2)+'\n')


if __name__ == '__main__':
    main()
