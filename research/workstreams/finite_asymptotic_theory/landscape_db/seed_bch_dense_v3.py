"""Generate a bounded initial type cover for another BCH reference geometry."""
import argparse
import json
import math
from pathlib import Path

import numpy as np

import occupation_refresh_v1 as transfer
import study_bch_dominance_v1 as study

HERE = Path(__file__).resolve().parent


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--block', type=int, choices=(64, 128), default=64)
    parser.add_argument('--step', type=int, default=64)
    parser.add_argument('--state', type=int, default=20)
    parser.add_argument('--exponent', type=int, default=20)
    parser.add_argument('--minimum', type=int, default=257)
    parser.add_argument('--nodes', type=int, default=511)
    args = parser.parse_args()
    _, spectra, maps, dependencies = study.load_inputs()
    block, t, s = args.block, args.step, args.state
    length = (1 << args.exponent)//(block//2); config = maps[t, s]; counts = spectra[block]
    ac = {w: n for w, n in enumerate(config['a_counts']) if w and n}
    bands = [[w for w in counts if w < block]]
    if block in counts:
        bands.append([block])
    model = transfer.DenseModel(counts, block, t, s, ac, config['kernel_counts'], length,
                                np.arange(-24, 5)/4, bands=bands)
    result = model.search(minimum=args.minimum, maximum_nodes=args.nodes, target_bits=60.)
    for path in (Path(__file__), Path(transfer.__file__)):
        dependencies[path.relative_to(study.ROOT).as_posix()] = study.sha(path)
    payload = dict(status='BINARY64_INITIAL_DENSE_COVER', source_sha256=dependencies,
                   arguments=vars(args), dense=result)
    path = HERE/f'bch_dense_seed_bulk_b{block}_t{t}_s{s}_e{args.exponent}_q{args.minimum}.json'
    if path.exists():
        raise ValueError('refusing to overwrite an existing seed')
    path.write_text(json.dumps(payload, indent=2)+'\n')
    print('seed dense margin', -result['log_union_upper']/math.log(2), flush=True)


if __name__ == '__main__':
    main()
