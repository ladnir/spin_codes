"""Two bounded gap proposals derived from the frozen weight2496 bank.

Each invocation fits markers once at its declared center and retunes one
common state potential. Alpha stays fixed; log(z) scales with the center.
The complete assigned gap is screened in FLOAT and requires directed replay.
"""
from hashlib import sha256
import json
from pathlib import Path
from time import perf_counter
from unittest.mock import patch

from hd_shifted_outer_density import TOEPLITZ_OUTER_MEASURE
from hd_shifted_potential_band_neighbor import fit_one_central_bank
from hd_shifted_potential_retune import retune_potential
from hd_shifted_run_count_screen import run_log_screen_float
from hd_shifted_swap_certificate import SWAP_CONSTRUCTION


GAPS = {2208: (2154, 2274), 2752: (2684, 2827)}
SOURCE = Path(__file__).with_name('hd_shifted_potential_band_2496_choice0_2154_2827_selected.json')
COUNT_METHOD = 'FLOAT-log-factorial-cyclic-run-counts'


def run(center, output):
    if type(center) is not int or center not in GAPS:
        raise ValueError('Only the two explicitly assigned centers are permitted')
    first, last = GAPS[center]
    output = Path(output)
    if output.exists():
        raise ValueError('Frozen artifacts are never overwritten')
    started = perf_counter()
    emit = lambda row: print(json.dumps(row, allow_nan=False), flush=True)
    raw = SOURCE.read_bytes()
    source = json.loads(raw)
    emit(dict(stage='gap-central-fit-start', center=center, weight_range=[first, last],
              marker_step=32, z_choices=1, outer_measure=TOEPLITZ_OUTER_MEASURE))
    initial = fit_one_central_bank(source, center, marker_step=32,
                                    outer_measure=TOEPLITZ_OUTER_MEASURE)
    frozen_bank = dict(construction=SWAP_CONSTRUCTION, fixed_witness=initial.witness())
    marker_count = len(initial.tables)
    del initial
    emit(dict(stage='gap-central-markers-frozen', center=center, markers=marker_count,
              elapsed_seconds=perf_counter()-started))
    screen, diagnostics = retune_potential(frozen_bank, center,
                                            outer_measure=TOEPLITZ_OUTER_MEASURE)
    # This inherited description predates the explicit common-potential retune.
    screen.source_witness.pop('alpha_and_potential', None)
    screen.source_witness['parameter_rule'] = (
        'z density-scaled; alpha unchanged; markers fitted once; common potential retuned once centrally')
    assert screen.source_witness['alpha'] == source['fixed_witness']['source_exact_witness']['alpha']
    completed = 0
    def progress(row):
        nonlocal completed
        completed += 1
        if completed % 32 == 0 or row['active_groups'] == last:
            emit(dict(stage='gap-fixed-screen', center=center, completed=completed,
                      total=last-first+1, weight=row['active_groups'],
                      log2_upper_screen=row['log2_upper_screen'], seconds=perf_counter()-started))
    with patch('hd_shifted_potential_band_screen.run_log_screen', run_log_screen_float):
        proposal = screen.band(first, last, max_weights=last-first+1,
                                threshold=-60., progress=progress)
    proposal.update(source_artifact_sha256=sha256(raw).hexdigest(),
        source_artifact=SOURCE.name, central_weight=center, potential_retuning=diagnostics,
        marker_step=32, central_marker_fits=marker_count, witness_searches=marker_count,
        z_choices=1, count_method=COUNT_METHOD, requires_fresh_directed_replay=True,
        elapsed_seconds=perf_counter()-started)
    encoded = json.dumps(proposal, allow_nan=False, separators=(',', ':'))
    with output.open('x', encoding='utf-8') as target:
        target.write(encoded+'\n')
    emit(dict(stage='gap-proposal-saved', center=center, output=str(output),
              failed_weights=proposal['failed_weights'],
              union_log2_upper_screen=proposal['union_log2_upper_screen'],
              seconds=perf_counter()-started, directed=False))
    return proposal


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--center', type=int, choices=tuple(GAPS), required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    run(args.center, args.output)
