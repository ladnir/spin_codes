"""Bounded r5888 FLOAT search; immutable complete-profile trial checkpoints."""

from hashlib import sha256
import json
from pathlib import Path

import hd_field128_gamma_midpoint_dense_mark_search as search


def main():
    here = Path(__file__).resolve().parent
    runner = Path(__file__).resolve()
    runner_hash = sha256(runner.read_bytes()).hexdigest()
    source = here/'hd_field128_gamma30_d40_r5888_dense_calibration_20261010/adaptive.json'
    source_hash = '804fa8e94946b94403521c8d2a9a4f5644d8175a2f6bdf0cb8b7f6411ccf5c1c'
    if sha256(source.read_bytes()).hexdigest() != source_hash:
        raise RuntimeError('Frozen35-marker literal source changed')
    if sha256(Path(search.__file__).read_bytes()).hexdigest() != '99ca1cde5d29e20eef67f2aa975bec2fc3f8d5e27ac3a621a5685fbe5072ab05':
        raise RuntimeError('Reviewed dense search source changed')
    output = here/'hd_field128_gamma30_d40_r5888_dense_mark_search_20261010'
    output.mkdir(exist_ok=False)
    def check():
        if sha256(runner.read_bytes()).hexdigest() != runner_hash:
            raise RuntimeError('Orchestration source changed')
        if sha256(source.read_bytes()).hexdigest() != source_hash:
            raise RuntimeError('Literal parent changed')
    def checkpoint(row):
        check()
        row['source_sha256'][runner.name] = runner_hash
        pin = search.provider.search.write_new(output/f"evaluation_{row['evaluation_index']:03d}.json", row)
        print(json.dumps(dict(evaluation=row['evaluation_index'], score=row['log2_upper_screen'],
                              largest=row['largest_profile'], elapsed_seconds=row['elapsed_seconds'], checkpoint=pin)), flush=True)
        return pin
    result = search.run_search(source, 5888, max_evaluations=40, max_seconds=120., checkpoint=checkpoint)
    check()
    result['source_sha256'][runner.name] = runner_hash
    pin = search.provider.search.write_new(output/'summary.json', result)
    best = result['best']
    print(json.dumps(dict(summary=pin, completed_trials=result['completed_trials'],
                          best_evaluation=best['evaluation_index'] if best else None,
                          best_score=best['log2_upper_screen'] if best else None,
                          optimizer=result['optimizer'], elapsed_seconds=result['elapsed_seconds'])), flush=True)


if __name__ == '__main__':
    main()
