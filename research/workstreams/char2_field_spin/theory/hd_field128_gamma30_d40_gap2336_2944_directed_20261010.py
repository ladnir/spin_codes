"""Fresh reviewed banks and explicit assigned-gap bands, never interpolated bounds."""

import os
os.environ['OPENBLAS_NUM_THREADS'] = '1'

import argparse
from fractions import Fraction as F
from hashlib import sha256
import gc
import json
from pathlib import Path
from shutil import disk_usage
from time import monotonic

import hd_field128_gamma_all_pair_certificate as p
import hd_field128_gamma_midpoint_weight_bands as wb


HERE = Path(__file__).resolve().parent
RUNNER = Path(__file__).resolve()
RUNNER_HASH = sha256(RUNNER.read_bytes()).hexdigest()
CONFIG = {
    2336: (71, '60dbab49585c3a5b4988c2b783fcf0c2c98b23961b923be2490d30c81c194a47',
           ((2304, 2367), (2257, 2303), (2368, 2431))),
    2944: (68, '8a0f5f235bcc87cd83e18d62f428ad67fbf41ab3d3c829bb5fbc2993aa0a0120',
           ((2912, 2975), (2848, 2911), (2976, 3007)))}


def write_new(path, value):
    if disk_usage(HERE).free < 750_000_000:
        raise OSError('Checkpoint guard: less than750MB disk space remains')
    return wb.profiles.write_new(path, value)


def initialize(weight):
    if weight not in CONFIG:
        raise ValueError('Only the two assigned witness centers are allowed')
    marker_count, proposal_hash, ranges = CONFIG[weight]
    proposal_path = HERE/f'hd_field128_gamma30_d40_r{weight}_gap_calibration_20261010/proposal.json'
    raw = proposal_path.read_bytes()
    if sha256(raw).hexdigest() != proposal_hash:
        raise RuntimeError('Reviewed exact proposal changed')
    proposal = json.loads(raw)
    if (proposal['construction'] != p.local.construction(30) or proposal['state_dimension'] != 30
            or proposal['active_groups'] != weight or proposal['target_distance'] != '2/5'
            or proposal['marker_count'] != marker_count or len(proposal['witness']['markers']) != marker_count):
        raise ValueError('Unexpected literal proposal scope')
    pins = wb.source_pins()
    for name, digest in proposal['source_sha256'].items():
        if name in pins and pins[name] != digest:
            raise RuntimeError('Source pin conflict: '+name)
        pins[name] = digest
    pins[RUNNER.name] = RUNNER_HASH
    def check():
        if proposal_path.read_bytes() != raw:
            raise RuntimeError('Exact literal proposal changed')
        wb.check_source_pins(pins)
    check()
    output = HERE/f'hd_field128_gamma30_d40_r{weight}_gap_directed256_20261010'
    output.mkdir(exist_ok=False)
    witness, started = proposal['witness'], monotonic()
    witness_source = dict(path=str(proposal_path), sha256=proposal_hash)
    try:
        bank = p.AliasBank(F(witness['z']), F(witness['alpha']),
            [tuple(map(F, row)) for row in witness['markers']], state=30, distance=F(2, 5),
            grid=128, precision=256, bits=80, max_seconds=180., max_table_bytes=64 << 20,
            progress=lambda row: print(json.dumps(row), flush=True))
        check()
        metadata = dict(bank.metadata(), witness_source=witness_source, calibration_active_groups=weight,
            orchestration_source_sha256={RUNNER.name: RUNNER_HASH},
            proposal_source_sha256=proposal['source_sha256'], old_numerical_bounds_reused=False,
            bank_build_seconds=monotonic()-started)
        pin = write_new(output/'bank.json', metadata)
        print(json.dumps(dict(stage='bank-complete', bank=pin, seconds=metadata['bank_build_seconds'])), flush=True)
    except (TimeoutError, MemoryError, OSError) as error:
        check()
        failure = dict(stage='bank-build', calibration_active_groups=weight, witness_source=witness_source,
            error_type=type(error).__name__, message=str(error), directed=False, whole_code_certificate=False,
            source_sha256=dict(sorted(pins.items())), elapsed_seconds=monotonic()-started)
        try:
            print(json.dumps(dict(failure=write_new(output/'bank_failure.json', failure))), flush=True)
        except OSError:
            print(json.dumps(failure), flush=True)
        raise
    return bank, check, output, witness_source, ranges


def run_bands(bank, check, output, witness_source, ranges):
    if type(bank) is not p.AliasBank:
        raise ValueError('Expected the exact frozen all-pair bank producer')
    links = []
    for first, last in ranges:
        directory = output/f'band_{first:04d}_{last:04d}'
        directory.mkdir(exist_ok=False)
        started, evaluator = monotonic(), None
        check()
        try:
            evaluator = wb.Evaluator(bank, first, last, max_seconds=180., witness_source=witness_source,
                progress=lambda row: print(json.dumps(row), flush=True))
            result = evaluator.run_cover(max_nodes=512, max_seconds=60., cell_failure_bits=72, aggregate_bits=54,
                context_checkpoint=lambda row: write_new(directory/'context.json', row),
                checkpoint=lambda row: write_new(directory/f"node_{row['sequence']:06d}.json", row))
            check()
            result['orchestration_source_sha256'] = {RUNNER.name: RUNNER_HASH}
            result['total_preparation_and_cover_seconds'] = monotonic()-started
            pin = write_new(directory/'summary.json', result)
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
                witness_source=witness_source, orchestration_source_sha256={RUNNER.name: RUNNER_HASH},
                elapsed_seconds=monotonic()-started)
            try:
                pin = write_new(directory/'failure.json', failure)
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
    manifest = dict(schema='gamma-directed-explicit-band-links-v1', ranges=[list(x) for x in ranges],
        bands=links, source_sha256={RUNNER.name: RUNNER_HASH}, witness_source=witness_source,
        nonoverlapping_weight_ranges=True, numerical_weight_interpolation=False,
        whole_code_certificate=False, union_probability_claim=False)
    pin = write_new(output/'band_links.json', manifest)
    print(json.dumps(dict(links=pin)), flush=True)
    return manifest


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--weight', type=int, choices=tuple(CONFIG), required=True)
    args = parser.parse_args()
    bank, check, output, witness_source, ranges = initialize(args.weight)
    band_results = run_bands(bank, check, output, witness_source, ranges)
    print('Assigned bands finished; actual bank retained, all Evaluator tables released.', flush=True)
