"""One explicit2048-node retry after the preserved512-node partial band."""

from hashlib import sha256
import gc
import json
from pathlib import Path

import hd_field128_gamma30_d40_r1168_gap_directed_20261010 as original


def run(bank, parent_check):
    here, runner = Path(__file__).resolve().parent, Path(__file__).resolve()
    digest = sha256(runner.read_bytes()).hexdigest()
    root = here/'hd_field128_gamma30_d40_r1168_gap_directed256_20261010'
    previous = root/'band_1145_1191/summary.json'
    previous_hash = 'c105375d56ff3ff09f17a96acbbbd46c3e219561e5329d7dbf87bf015966dd61'
    metadata = json.loads((root/'bank.json').read_bytes())
    def check():
        parent_check()
        if (sha256(runner.read_bytes()).hexdigest() != digest
                or sha256(previous.read_bytes()).hexdigest() != previous_hash):
            raise RuntimeError('Retry source or preserved partial changed')
    check()
    if type(bank) is not original.p.AliasBank or bank.metadata()['witness'] != metadata['witness']:
        raise ValueError('Expected the unchanged actual1168 bank')
    if metadata['calibration_active_groups'] != 1168:
        raise ValueError('Unexpected actual bank calibration')
    output = root/'retry2048'
    output.mkdir(exist_ok=False)
    evaluator = None
    try:
        evaluator = original.wb.Evaluator(bank, 1145, 1191, max_seconds=180.,
            witness_source=metadata['witness_source'], progress=lambda row: print(json.dumps(row), flush=True))
        result = evaluator.run_cover(max_nodes=2048, max_seconds=60., cell_failure_bits=72, aggregate_bits=54,
            context_checkpoint=lambda row: original.write_new(output/'context.json', row),
            checkpoint=lambda row: original.write_new(output/f"node_{row['sequence']:06d}.json", row))
        check()
        result['orchestration_source_sha256'] = {runner.name: digest,
            original.RUNNER.name: original.RUNNER_HASH}
        result['preserved_prior_partial'] = dict(path=str(previous), sha256=previous_hash)
        pin = original.write_new(output/'summary.json', result)
        print(json.dumps(dict(summary=pin, complete=result['complete_requested_band'],
            passed=result['meets_requested_band_target'], log=result['log2_upper_rational'],
            nodes=result['tree_nodes'], stop_reason=result['stop_reason'])), flush=True)
        return result
    finally:
        evaluator = None
        gc.collect()
        check()
