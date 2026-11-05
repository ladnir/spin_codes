"""Four explicit neighboring bands using the retained reviewed r6400 bank."""

from fractions import Fraction as F
from hashlib import sha256
import gc
import json
from pathlib import Path
from time import monotonic

import hd_field128_gamma30_d40_r6400_directed_20261010 as parent


RANGES = ((6304, 6367), (6240, 6303), (6176, 6239), (6432, 6495))
RUNNER = Path(__file__).resolve()
RUNNER_HASH = sha256(RUNNER.read_bytes()).hexdigest()
PARENT_HASH = '7ee4c1f3ddcee34b04382804a95aaf15e29acb7fbe3ea6d3228066872220b069'


def run(bank, parent_check):
    def check():
        parent_check()
        if sha256(RUNNER.read_bytes()).hexdigest() != RUNNER_HASH:
            raise RuntimeError('Extension runner changed')
        if sha256(Path(parent.__file__).read_bytes()).hexdigest() != PARENT_HASH:
            raise RuntimeError('Reviewed parent runner changed')
        bank.check_pins()
    check()
    source = json.loads(parent.PROPOSAL.read_bytes())
    witness = source['witness']
    if (type(bank) is not parent.p.AliasBank or bank.state != 30 or bank.distance != F(2, 5)
            or bank.z != F(witness['z']) or bank.alpha != F(witness['alpha'])
            or bank.markers != tuple(tuple(map(F, row)) for row in witness['markers'])):
        raise ValueError('Retained bank does not match exact reviewed proposal')
    links = []
    for first, last in RANGES:
        output = parent.OUT/f'band_{first:04d}_{last:04d}'
        output.mkdir(exist_ok=False)
        evaluator, started = None, monotonic()
        check()
        try:
            evaluator = parent.wb.Evaluator(bank, first, last, max_seconds=180.,
                witness_source=dict(path=str(parent.PROPOSAL), sha256=parent.PROPOSAL_HASH),
                progress=lambda row: print(json.dumps(row), flush=True))
            result = evaluator.run_cover(max_nodes=512, max_seconds=60.,
                cell_failure_bits=80, aggregate_bits=54,
                context_checkpoint=lambda row: parent.write_new(output/'context.json', row),
                checkpoint=lambda row: parent.write_new(output/f"node_{row['sequence']:06d}.json", row))
            check()
            result['orchestration_source_sha256'] = {
                RUNNER.name: RUNNER_HASH, Path(parent.__file__).name: PARENT_HASH}
            result['total_preparation_and_cover_seconds'] = monotonic()-started
            pin = parent.write_new(output/'summary.json', result)
            link = dict(requested_weight_range=[first, last], summary=pin,
                complete_requested_band=result['complete_requested_band'],
                meets_requested_band_target=result['meets_requested_band_target'],
                log2_upper_rational=result['log2_upper_rational'], tree_nodes=result['tree_nodes'])
            links.append(link)
            print(json.dumps(dict(stage='band-saved', **link)), flush=True)
        except (TimeoutError, MemoryError, OSError) as error:
            check()
            failure = dict(schema='gamma-directed-band-orchestration-failure-v1',
                requested_weight_range=[first, last], complete_requested_band=False,
                log2_upper_rational=None, directed=False, whole_code_certificate=False,
                failure=dict(error_type=type(error).__name__, message=str(error)),
                witness_source=dict(path=str(parent.PROPOSAL), sha256=parent.PROPOSAL_HASH),
                orchestration_source_sha256={RUNNER.name: RUNNER_HASH},
                elapsed_seconds=monotonic()-started)
            try:
                pin = parent.write_new(output/'failure.json', failure)
                print(json.dumps(dict(failure=pin)), flush=True)
            except OSError:
                print(json.dumps(failure), flush=True)
            break
        finally:
            evaluator = None
            gc.collect()
            check()
    check()
    for link in links:
        if sha256(Path(link['summary']['path']).read_bytes()).hexdigest() != link['summary']['sha256']:
            raise RuntimeError('Completed band bytes changed')
    manifest = dict(schema='gamma-directed-explicit-band-links-v1', ranges=[list(x) for x in RANGES],
        bands=links, source_sha256={RUNNER.name: RUNNER_HASH, Path(parent.__file__).name: PARENT_HASH},
        witness_source=dict(path=str(parent.PROPOSAL), sha256=parent.PROPOSAL_HASH),
        nonoverlapping_weight_ranges=True, numerical_weight_interpolation=False,
        whole_code_certificate=False, union_probability_claim=False)
    print(json.dumps(dict(links=parent.write_new(parent.OUT/'extension_links.json', manifest))), flush=True)
    return manifest
