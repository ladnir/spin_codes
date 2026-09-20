"""Separate box slack from fixed-type witness slack with direct moment searches."""
import argparse
import json
import math
from pathlib import Path

import numpy as np

import density_dense_witness_v1 as search
import study_bch_dominance_v1 as study

HERE = Path(__file__).resolve().parent


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input', type=Path, required=True)
    parser.add_argument('--limit', type=int, default=3)
    args = parser.parse_args()
    payload = json.loads(args.input.read_text()); settings = payload['arguments']
    b, t, s, e = (settings[k] for k in ('block', 'step', 'state', 'exponent'))
    length = (1 << e)//(b//2)
    _, counts, maps, dependencies = study.load_inputs(); config = maps[t, s]
    seed_path = HERE/f'bch_dense_seed_split_one_b{b}_t{t}_s{s}_e{e}_q{settings["minimum"]}.json'
    bands = json.loads(seed_path.read_text())['dense']['bands']
    ac = {w: n for w, n in enumerate(config['a_counts']) if w and n}
    model = search.DensityCoverRefiner(counts[b], b, t, s, ac, config['kernel_counts'], length, bands)
    table_lookup = model.lookup
    def direct_moment(x, y):
        theta = 1/(1+math.exp(-y))
        raw = search.density.base.epoch_mixture(model.epoch(x), t, [theta])
        return float(search.density.base.terminal_logs(raw, b*length//t)[0])/(b*length)
    def direct_lookup(x, y):
        h = 1e-5
        return (direct_moment(x, y), (direct_moment(x+h, y)-direct_moment(x-h, y))/(2*h),
                (direct_moment(x, y+h)-direct_moment(x, y-h))/(2*h))
    results = []
    for leaf in sorted(payload['leaves'], key=lambda x: x['own_log_bound'], reverse=True)[:args.limit]:
        witness = leaf['witness']
        corners = search.density.base.typed.vertices(leaf['lower'], leaf['upper'], length)
        values = model.value(corners, witness['log_surprisal'], np.array(witness['proposal']),
                             np.array(witness['probabilities']), model.costs(witness['probabilities']))
        vertex = corners[int(np.argmax(values))].tolist()
        point = dict(lower=vertex, upper=vertex, witness=witness, own_log_bound=float(max(values)))
        model.lookup = table_lookup
        refined = model.refine(point)
        model.lookup = direct_lookup
        exact = model.refine(refined)
        replay = model.direct(vertex, vertex, exact['witness'])
        if abs(replay-exact['own_log_bound']) > 2e-6:
            raise ArithmeticError('direct-search point failed replay')
        result = dict(vertex=vertex, box_margin_bits=-leaf['own_log_bound']/math.log(2),
                      point_margin_bits=-point['own_log_bound']/math.log(2),
                      table_search_margin_bits=-refined['own_log_bound']/math.log(2),
                      direct_search_margin_bits=-replay/math.log(2), witness=exact['witness'])
        results.append(result); print(json.dumps(result), flush=True)
    for path in (Path(__file__), Path(search.__file__), Path(search.density.__file__), Path(search.joint.__file__),
                 args.input.resolve(), seed_path):
        dependencies[path.relative_to(study.ROOT).as_posix()] = study.sha(path)
    output = dict(status='FIXED_TYPE_WITNESS_PROBE', arguments=settings, results=results,
                  source_sha256=dependencies, limitations=['Selected type vectors only; no complete interval conclusion.',
                                                          'Nearest binary64 diagnostics.'])
    (HERE/f'bch_dense_point_probe_b{b}_t{t}_s{s}_e{e}.json').write_text(json.dumps(output, indent=2)+'\n')


if __name__ == '__main__':
    main()
