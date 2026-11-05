"""Bounded explicit bands using the unchanged directed midpoint checker."""

import gc
from hashlib import sha256
import json
from pathlib import Path
import shutil

import hd_field128_gamma_midpoint_weight_bands as bands


HERE = Path(__file__).resolve().parent


def run(bank, proposal, source_pin, first, last, *, tag):
    if tag not in ('r3936', 'r128', 'r112', 'r1168'):
        raise ValueError('Unexpected explicit run tag')
    bands.checked_weights(4096, first, last)
    if shutil.disk_usage(HERE).free < 750 << 20:
        raise RuntimeError('Less than750MiB free before a new band')
    driver_pin = dict(path=str(Path(__file__).resolve()),
                      sha256=sha256(Path(__file__).read_bytes()).hexdigest())

    def check():
        if sha256(Path(source_pin['path']).read_bytes()).hexdigest() != source_pin['sha256']:
            raise RuntimeError('Positive witness source changed')
        if sha256(Path(__file__).read_bytes()).hexdigest() != driver_pin['sha256']:
            raise RuntimeError('Band driver changed')
        bands.check_source_pins(proposal['source_sha256'])
        bank.check_pins()

    check()
    output = HERE/f'hd_field128_gamma30_d40_weights{first}_{last}_{tag}_band_v21_20261010'
    output.mkdir(exist_ok=False)
    evaluator = None
    try:
        evaluator = bands.Evaluator(bank, first, last, max_seconds=180,
                                    witness_source=source_pin)
        def checkpoint(row):
            if shutil.disk_usage(HERE).free < 750 << 20:
                raise RuntimeError('Disk guard stopped further checkpoints')
            bands.profiles.write_new(output/f"node_{row['sequence']:06d}.json", row)
        result = evaluator.run_cover(max_nodes=512, max_seconds=60,
            cell_failure_bits=72, aggregate_bits=54,
            context_checkpoint=lambda row: bands.profiles.write_new(output/'context.json', row),
            checkpoint=checkpoint)
        check()
        result['explicit_band_driver'] = driver_pin
        pin = bands.profiles.write_new(output/'summary.json', result)
        print(json.dumps(dict(band=[first, last], summary=pin,
            status=result['status'], complete=result['complete_requested_band'],
            target_met=result['meets_requested_band_target'], nodes=result['tree_nodes'],
            log2_upper=result['log2_upper_rational'])), flush=True)
        return result, pin
    finally:
        evaluator = None
        gc.collect()
