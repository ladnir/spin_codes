"""Replay and clip a saved complete type cover to a higher occupation cutoff."""
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np

import refine_bch_dense_v1 as dense
import study_bch_dominance_v1 as study
import trim_bch_dense_v3 as checks

HERE = Path(__file__).resolve().parent


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input', type=Path, required=True)
    parser.add_argument('--minimum', type=int, required=True)
    args = parser.parse_args()
    raw = args.input.read_bytes(); checkpoint = json.loads(raw)
    settings = checkpoint['arguments']
    block, t, s, exponent = (settings[k] for k in ('block', 'step', 'state', 'exponent'))
    length = (1 << exponent)//(block//2)
    if not settings['minimum'] <= args.minimum <= length:
        raise ValueError('new interval must be contained in the saved cover')
    _, spectra, maps, dependencies = study.load_inputs()
    for name, digest in checkpoint['source_sha256'].items():
        if study.sha(study.ROOT/name) != digest:
            raise ValueError(f'changed checkpoint dependency: {name}')
    original_path = HERE/f'bch_dense_seed_b{block}_t{t}_s{s}_e{exponent}_q{settings["minimum"]}.json'
    original = json.loads(original_path.read_text())['dense']
    dense.check_cover(checkpoint['leaves'], length, settings['minimum'])
    config = maps[t, s]
    ac = {w: n for w, n in enumerate(config['a_counts']) if w and n}
    model = checks.previous.CoverRefiner(spectra[block], block, t, s, ac, config['kernel_counts'], length, original['bands'])
    selected = []
    for index, leaf in enumerate(checkpoint['leaves']):
        old = checks.canonical_witness(model, leaf)
        lo, hi = list(old['lower']), list(old['upper'])
        hi[0] = min(hi[0], length-args.minimum)
        if any(a > b for a, b in zip(lo, hi)) or not dense.lattice_count(lo, hi, length):
            continue
        value = model.direct(lo, hi, old['witness'])
        selected.append(checks.canonical_witness(model, dict(lower=lo, upper=hi, own_log_bound=value, witness=old['witness'])))
        if index % 250 == 0:
            print(f'replayed {index+1}/{len(checkpoint["leaves"])} boxes', flush=True)
    coverage = dense.check_cover(selected, length, args.minimum)
    saved = HERE/f'bch_dense_checkpoint_snapshot_{hashlib.sha256(raw).hexdigest()}.json'
    if not saved.exists():
        saved.write_bytes(raw)
    dependencies.update(checkpoint['source_sha256'])
    for path in (Path(__file__), Path(dense.__file__), Path(checks.__file__), saved, original_path):
        dependencies[path.relative_to(study.ROOT).as_posix()] = study.sha(path)
    total = float(np.logaddexp.reduce([x['own_log_bound'] for x in selected]))
    payload = dict(arguments=dict(block=block, step=t, state=s, exponent=exponent, minimum=args.minimum),
                   source_sha256=dependencies, coverage=coverage,
                   dense=dict(selected_boxes=selected, occupation_min=args.minimum, occupation_max=length,
                              bands=original['bands'], log_union_upper=total))
    path = HERE/f'bch_dense_seed_b{block}_t{t}_s{s}_e{exponent}_q{args.minimum}.json'
    path.write_text(json.dumps(payload, indent=2)+'\n')
    print('clipped full interval margin', -total/np.log(2), coverage, flush=True)


if __name__ == '__main__':
    main()
