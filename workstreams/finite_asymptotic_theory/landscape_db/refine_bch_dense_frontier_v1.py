"""Refine weak BCH dense bounds where the retained sparse union is useful."""
import argparse
import json
import math
from pathlib import Path
import subprocess
import sys

import numpy as np

import study_bch_dominance_v1 as study
import transport_bch_full_cover_v1 as transport

HERE = Path(__file__).resolve().parent


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--limit', type=int)
    parser.add_argument('--maximum-refinements', type=int, default=300)
    args = parser.parse_args()
    if args.maximum_refinements < 0 or (args.limit is not None and args.limit < 1):
        parser.error('Expected nonnegative refinements and a positive limit')
    selected = []
    for path in sorted(HERE.glob('bch_full_reference_b*_t*_s*_e*.json')):
        data = json.loads(path.read_text())
        if data['status'] != 'VERIFIED_BINARY64_FULL_REFERENCE':
            raise ValueError('Expected verified full references')
        if data['full_margin_bits'] > 0 or data.get('dense_cover') is None:
            continue
        transport.authenticate(data)
        intervals = transport.select_covers(data)
        dense = [(p, d) for p, d in intervals if d['status'] == 'BINARY64_COMPLETE_DENSE_INTERVAL']
        if len(dense) != 1:
            raise ValueError('Expected one selected dense cover')
        dense_path, cover = dense[0]; minimum = cover['occupation_min']
        partial = float(np.logaddexp.reduce([c['log_upper'] for c in data['components'] if c['occupation_max'] < minimum]))
        if partial >= 0:
            continue
        b, t, s, e = (data[k] for k in ('block_bits', 'step_bits', 'state_bits', 'message_exponent'))
        own = HERE/f'bch_dense_v11_b{b}_t{t}_s{s}_e{e}_q{minimum}/cover.json'
        if dense_path.resolve() == own.resolve():
            continue
        selected.append(dict(geometry=f'{b}:{t}:{s}:{e}', block_bits=b, step_bits=t, state_bits=s,
                             message_exponent=e, reference=str(path),
                             preceding_full_margin_bits=data['full_margin_bits'],
                             sparse_union_margin_bits=-partial/math.log(2)))
    selected.sort(key=lambda r: (r['step_bits'] != 64, -r['state_bits'], r['block_bits'], r['message_exponent']))
    if args.limit is not None:
        selected = selected[:args.limit]
    print(f'Selected {len(selected)} weak dense bounds with useful sparse unions', flush=True)
    for row in selected:
        print(row['geometry'], 'sparse margin', row['sparse_union_margin_bits'], flush=True)
    if not selected:
        return
    command = [sys.executable, '-u', str(HERE/'refine_bch_transported_dense_v1.py'),
               '--maximum-refinements', str(args.maximum_refinements)]
    for row in selected:
        command += ['--geometry', row['geometry']]
    subprocess.run(command, cwd=HERE, check=True)
    for row in selected:
        path = Path(row['reference']); current = json.loads(path.read_text())
        transport.authenticate(current)
        row.update(full_margin_bits=current['full_margin_bits'],
                   margin_penalty_bits=current['margin_penalty_bits'],
                   reference_sha256=study.sha(path))
    paths = [Path(__file__), HERE/'refine_bch_transported_dense_v1.py', *[Path(r['reference']) for r in selected]]
    result = dict(status='VERIFIED_DENSE_FRONTIER_REFINEMENT', rows=selected,
                  useful_full_bounds=sum(r['full_margin_bits'] > 0 for r in selected),
                  maximum_refinements=args.maximum_refinements,
                  source_sha256={str(p): study.sha(p) for p in paths},
                  limitations=['Searches dense witnesses at the target map.',
                               'Weak sparse unions require a separate refinement.', 'No outward arithmetic.'])
    (HERE/'bch_dense_frontier_refinement_v1.json').write_text(json.dumps(result, indent=2)+'\n')
    print('New useful full bounds', result['useful_full_bounds'], 'of', len(selected), flush=True)


if __name__ == '__main__':
    main()
