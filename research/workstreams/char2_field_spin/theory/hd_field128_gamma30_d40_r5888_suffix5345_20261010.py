"""Fresh explicit5345..5407 cover; the prior5344 band remains partial."""

from hashlib import sha256
import gc
import json
from pathlib import Path
from time import monotonic

import hd_field128_gamma_midpoint_weight_bands as wb


def run(bank, parent_check):
    here, source = Path(__file__).resolve().parent, Path(__file__).resolve()
    digest, pins = sha256(source.read_bytes()).hexdigest(), wb.source_pins()
    pins[source.name] = digest
    proposal = here/'hd_field128_gamma30_d40_r5888_winner_replay_20261010/proposal.json'
    proposal_hash = '43cea88db107cfc364ab5989978269fa97b17d0068e32ca5d24aa930887e4754'
    def check():
        parent_check()
        wb.check_source_pins(pins)
        if sha256(proposal.read_bytes()).hexdigest() != proposal_hash:
            raise RuntimeError('Exact17 proposal changed')
    check()
    output = here/'hd_field128_gamma30_d40_r5888_winner17_directed256_20261010/band_5345_5407'
    output.mkdir(exist_ok=False)
    evaluator, started = None, monotonic()
    try:
        evaluator = wb.Evaluator(bank, 5345, 5407, max_seconds=180.,
            witness_source=dict(path=str(proposal), sha256=proposal_hash))
        result = evaluator.run_cover(max_nodes=512, max_seconds=60., cell_failure_bits=80, aggregate_bits=54,
            context_checkpoint=lambda row: wb.profiles.write_new(output/'context.json', row),
            checkpoint=lambda row: wb.profiles.write_new(output/f"node_{row['sequence']:06d}.json", row))
        check()
        result['orchestration_source_sha256'] = {source.name: digest}
        result['total_preparation_and_cover_seconds'] = monotonic()-started
        result['previous_partial_band_not_reused'] = True
        pin = wb.profiles.write_new(output/'summary.json', result)
        print(json.dumps(dict(summary=pin, status=result['status'],
            passing=result['meets_requested_band_target'], nodes=result['tree_nodes'],
            log2_upper_rational=result['log2_upper_rational'])), flush=True)
        return result
    finally:
        evaluator = None
        gc.collect()
        check()
