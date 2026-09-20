"""Small sequential BCH occupation probe; generated evidence remains local."""
import argparse
import json
import math
from pathlib import Path

import numpy as np

import occupation_refresh_v1 as refresh
import composition_boxes as boxes
import read_grid_receipts as receipts
import run_occupation_grid as source

HERE = Path(__file__).resolve().parent


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--block', type=int, choices=(64, 128), default=128)
    parser.add_argument('--step', type=int, default=64)
    parser.add_argument('--state', type=int, default=20)
    parser.add_argument('--exponent', type=int, default=20)
    parser.add_argument('--maximum', type=int, default=64)
    parser.add_argument('--nodes', type=int, default=255)
    parser.add_argument('--sparse-only', action='store_true')
    args = parser.parse_args()
    _, observations, _ = receipts.snapshot()
    row = next(r for r in observations.values() if 'BCH' in r['series']
               and tuple(int(r[k]) for k in ('block_bits', 'step_bits', 'state_bits', 'message_exponent'))
               == (args.block, args.step, args.state, args.exponent))
    config = json.loads((source.grid.pilot.ROOT/row['map_source']).read_text())
    counts = source.exact_counts(row)
    ac = {w: n for w, n in enumerate(config['a_counts']) if w and n}
    length = (1 << args.exponent)//(args.block//2)
    cutoff = args.block*length//10
    model = refresh.Epochs(args.step, args.state, ac, config['kernel_counts'])
    # Sparse state-onset tilts and larger tilts for repeated activations.
    tilts = np.unique(np.r_[np.arange(-20, 33, 2)/4-math.log(length), np.arange(-20, 5)/4])
    sparse = {q: refresh.SparseModel(counts, args.block, q, 4, 4) for q in (2, 3, 4)}
    components = {q: np.full(len(m.indices), np.inf) for q, m in sparse.items()}
    groups = [[w for w in counts if 16*w < 3*args.block],
              [w for w in counts if 3*args.block <= 16*w <= 13*args.block],
              [w for w in counts if 16*w > 13*args.block]]
    bands = [g for g in groups if g]
    models = [refresh.CompositionBoxes(counts, args.block, bands,
              [p for g, p in zip(groups, [low, .5, .87]) if g], args.maximum)
              for low in (.22, .28, .35)]
    best = np.full(args.maximum, np.inf)
    for index, tilt in enumerate(tilts):
        lam = math.exp(tilt)
        region = refresh.region_logs(model.at(lam, args.maximum), args.step, length, args.maximum)
        for q, m in sparse.items():
            for shift in np.arange(-2, 11)/2:
                np.minimum(components[q], m.components(region, cutoff, lam, shift), out=components[q])
        for m in models:
            current, matrices = region, []
            for q in range(1, args.maximum+1):
                current = m.envelope.apply(current[:-1], current[1:])
                matrices.append(current[0])
            values = refresh.terminal_logs(np.array(matrices), args.block)+cutoff*lam
            values += np.array([math.log(math.comb(length, q))+q*math.log(len(bands)) for q in range(1, args.maximum+1)])
            np.minimum(best, values, out=best)
        if index % 10 == 0:
            print(f'sparse tilt {index+1}/{len(tilts)}', flush=True)
    for q, m in sparse.items():
        best[q-1] = min(best[q-1], m.aggregate(components[q], length))
    print('sparse margins', {q: -best[q-1]/math.log(2) for q in (2, 3, 4, 8, 16, 32, 64) if q <= args.maximum}, flush=True)
    if args.sparse_only:
        return
    dense = refresh.DenseModel(counts, args.block, args.step, args.state, ac,
                               config['kernel_counts'], length, np.arange(-20, 5)/4, bands=bands)
    result = dense.search(minimum=args.maximum+1, maximum_nodes=args.nodes, target_bits=80)
    print('dense margin', -result['log_union_upper']/math.log(2), flush=True)
    path = HERE/f'bch_occupation_probe_b{args.block}_t{args.step}_s{args.state}_e{args.exponent}.json'
    path.write_text(json.dumps(dict(arguments=vars(args), sparse_log_bounds=best.tolist(), dense=result), indent=2)+'\n')


if __name__ == '__main__':
    main()
