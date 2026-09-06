"""Reuse a complete type partition at another inner map, replaying every bound."""
import argparse
import json
import math
from pathlib import Path

import numpy as np

import close_bch_dense_v2 as dense
import study_bch_dominance_v1 as study
import trim_bch_dense_v3 as checks

HERE = Path(__file__).resolve().parent


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input', type=Path, required=True)
    parser.add_argument('--step', type=int, required=True)
    parser.add_argument('--state', type=int, required=True)
    args = parser.parse_args()
    data = json.loads(args.input.read_text())
    block, exponent = data['block_bits'], data['message_exponent']
    length = (1 << exponent)//(block//2)
    if length % args.step:
        raise ValueError('whole epochs required')
    path = HERE/f'bch_dense_seed_b{block}_t{args.step}_s{args.state}_e{exponent}_q{data["occupation_min"]}.json'
    if path.exists():
        raise ValueError('refusing to overwrite an existing seed')
    _, spectra, maps, dependencies = study.load_inputs()
    for name, digest in data['source_sha256'].items():
        if study.sha(study.ROOT/name) != digest:
            raise ValueError(f'changed source cover dependency: {name}')
    leaves = [leaf for root in data['roots'] for leaf in root['leaves']]
    coverage = dense.base.check_cover(leaves, length, data['occupation_min'])
    config = maps[args.step, args.state]
    ac = {w: n for w, n in enumerate(config['a_counts']) if w and n}
    model = dense.CoverRefiner(spectra[block], block, args.step, args.state, ac,
                              config['kernel_counts'], length, data['bands'])
    selected = []
    for index, leaf in enumerate(leaves):
        value = model.direct(leaf['lower'], leaf['upper'], leaf['witness'])
        selected.append(checks.canonical_witness(model, dict(leaf, own_log_bound=value)))
        if index % 250 == 0:
            print(f'replayed target map on {index+1}/{len(leaves)} boxes', flush=True)
    total = float(np.logaddexp.reduce([leaf['own_log_bound'] for leaf in selected]))
    for source in (Path(__file__), Path(dense.__file__), Path(dense.base.__file__), Path(checks.__file__), args.input.resolve()):
        dependencies[source.relative_to(study.ROOT).as_posix()] = study.sha(source)
    result = dict(status='BINARY64_TRANSPORTED_DENSE_COVER', source_sha256=dependencies, coverage=coverage,
                  arguments=dict(block=block, step=args.step, state=args.state, exponent=exponent, minimum=data['occupation_min']),
                  dense=dict(selected_boxes=selected, bands=data['bands'], occupation_min=data['occupation_min'],
                             occupation_max=length, log_union_upper=total))
    path.write_text(json.dumps(result, indent=2)+'\n')
    print('transported dense interval margin', -total/math.log(2), flush=True)


if __name__ == '__main__':
    main()
