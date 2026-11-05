"""One bounded FLOAT proposal, explicitly cross-calibrated from r4096/d45."""

from fractions import Fraction as F
from hashlib import sha256
import json
from pathlib import Path
from time import monotonic

import hd_field128_gamma_all_pair_screen as screen


def main():
    directory = Path(__file__).resolve().parent
    source = directory/'hd_field128_gamma30_d45_r4096_allpairs_initial_20261010.json'
    output = directory/'hd_field128_gamma30_d40_r3072_cross_calibrated_float_20261010.json'
    if output.exists():
        raise FileExistsError('Output is immutable')
    raw = source.read_bytes()
    data = json.loads(raw)
    if (data['active_groups'], data['state_dimension'], data['target_distance']) != (4096, 30, '9/20'):
        raise ValueError('Unexpected source calibration')
    witness = data['witness']
    odds = tuple(tuple(float(F(x)) for x in pair) for pair in witness['marker_log_odds_exact'])
    ratio = float(F(witness['log2_potential_ratio_exact']))
    pins = screen.source_pins()
    pins[Path(__file__).name] = sha256(Path(__file__).read_bytes()).hexdigest()
    result = dict(schema='gamma-cross-calibrated-single-FLOAT-proposal-v1',
                  source=dict(path=str(source), sha256=sha256(raw).hexdigest()),
                  source_calibration=dict(active_groups=4096, target_distance='9/20', state_dimension=30),
                  target_calibration=dict(active_groups=3072, target_distance='2/5', state_dimension=30),
                  source_marker_calibration_retained=True, marker_recalibration_performed=False,
                  input_literals={key: witness[key] for key in ('z', 'alpha', 'marker_log_odds_exact', 'log2_potential_ratio_exact', 'grid')},
                  old_numerical_bounds_reused=False, directed=False, whole_code_certificate=False,
                  maximum_banks=1, maximum_evaluations=1, max_seconds=30, potential_optimization=False,
                  source_sha256=pins)
    started, budget = monotonic(), screen.old.Budget(30)
    try:
        bank = screen.AliasBankScreen(float(F(witness['z'])), float(F(witness['alpha'])), odds, ratio,
                    state=30, distance=F(2, 5), grid=witness['grid'], floor_policy='positive-fallback', _budget=budget)
        row = screen.selected_screen(bank, 3072, support_method='runs', optimize_potential=False,
                                     max_evaluations=1, _budget=budget)
        budget.check()
        result.update(status='complete', result=row, witness=bank.witness())
    except (TimeoutError, ArithmeticError) as error:
        result.update(status='incomplete', failure=dict(error_type=type(error).__name__, message=str(error)))
    result['elapsed_seconds'] = monotonic()-started
    after = screen.source_pins()
    after[Path(__file__).name] = sha256(Path(__file__).read_bytes()).hexdigest()
    if after != pins or source.read_bytes() != raw:
        raise RuntimeError('Pinned source changed')
    screen.search.write_new(output, result)
    print(json.dumps(dict(path=str(output), sha256=sha256(output.read_bytes()).hexdigest(),
                          status=result['status'], elapsed_seconds=result['elapsed_seconds'],
                          log2_upper_screen=result.get('result', {}).get('log2_upper_screen'))), flush=True)


if __name__ == '__main__':
    main()
