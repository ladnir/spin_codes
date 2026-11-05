"""Two bounded FLOAT proposal searches; no directed gate or distance claim."""

from fractions import Fraction as F
from hashlib import sha256
import json
from pathlib import Path
from time import monotonic

import hd_field128_gamma_all_pair_search as search
import hd_field128_lift_dense_proposals as dense


def main():
    directory = Path(__file__).resolve().parent
    output = directory/'hd_field128_gamma30_d40_gap_float_search_20261010'
    if output.exists():
        raise ValueError('Research output directory must be new')
    runner = Path(__file__).resolve()
    runner_hash = sha256(runner.read_bytes()).hexdigest()
    pins, dense_pins = search.source_pins(), dense.source_pins()
    source = dense.read_source()
    source_path = Path(source['provenance']['path'])
    source_hash = source['provenance']['sha256']
    started, created = monotonic(), False

    def check_pins():
        if (search.source_pins() != pins or dense.source_pins() != dense_pins
                or sha256(runner.read_bytes()).hexdigest() != runner_hash
                or sha256(source_path.read_bytes()).hexdigest() != source_hash):
            raise RuntimeError('A proposal input or source changed')

    def save(name, value):
        nonlocal created
        check_pins()
        if not created:
            output.mkdir(exist_ok=False)
            created = True
        return dense.write_new(output/name, value)

    for weight, exponent, z, alpha, baseline in (
        (512, 2, '5035991213944493/2305843009213693952',
         '5759504334290649/1152921504606846976', 352.34050668776035),
        (768, 3, '3353526911962953/36893488147419103232',
         '5632108150008701/1152921504606846976', 655.5344117923814),
    ):
        seconds = min(140., 285.-(monotonic()-started))
        if seconds <= 0:
            print(json.dumps(dict(stage='pending-weight', active_groups=weight,
                                  reason='overall285s budget exhausted')), flush=True)
            break
        artifact = dense.make_witness(source, weight, exponent)
        provenance = dict(path=str(runner), sha256=runner_hash,
                          marker_source=source['provenance'], marker_calibration=weight,
                          marker_transform='dense.make_witness from pinned old exact calibrated markers',
                          marks_role='new exact proposal literals from interpolation of successful marks; no bound interpolation',
                          old_numerical_bounds_reused=False)
        last_saved = baseline
        saved = []

        def checkpoint(row):
            nonlocal last_saved
            value = row['log2_upper_screen']
            compact = {key: row[key] for key in ('evaluation_index', 'z', 'alpha',
                       'log2_upper_screen', 'largest_profile', 'regional_alias_diagnostics')}
            compact.update(active_groups=weight, stage='FLOAT-evaluation')
            print(json.dumps(compact, allow_nan=False), flush=True)
            if value < last_saved-1.:
                exported = dict(row, proposal_wrapper_sha256={runner.name: runner_hash},
                                marker_proposal_source_sha256=dense_pins,
                                material_improvement_threshold_bits=1.)
                saved.append(save(f'r{weight}_improved_{row["evaluation_index"]:03d}.json', exported))
                last_saved = value

        result = search.run_search(artifact, provenance, r=weight, state=30, distance=F(2, 5),
                                   z=z, alpha=alpha, grid=128, max_evaluations=40,
                                   max_seconds=seconds, checkpoint=checkpoint)
        check_pins()
        best = result['best']
        result.update(proposal_wrapper_sha256={runner.name: runner_hash},
                      marker_proposal_source_sha256=dense_pins,
                      materially_improved_artifacts=saved, earlier_diagnostic_baseline=baseline,
                      baseline_role='comparison only; not part of any proof',
                      no_directed_gate=True, no_native_benchmark=True)
        summary_pin = None
        if best is not None and best['log2_upper_screen'] < baseline-1.:
            summary_pin = save(f'r{weight}_summary.json', result)
        print(json.dumps(dict(stage='FLOAT-search-complete', active_groups=weight,
                              best_log2_upper_screen=None if best is None else best['log2_upper_screen'],
                              mark_evaluations=result['mark_evaluations'],
                              elapsed_seconds=result['elapsed_seconds'], summary=summary_pin,
                              failure_count=len(result['failures']), directed=False,
                              whole_code_certificate=False), allow_nan=False), flush=True)


if __name__ == '__main__':
    main()
