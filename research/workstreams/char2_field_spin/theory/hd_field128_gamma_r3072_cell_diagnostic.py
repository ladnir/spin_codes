"""Four fresh markers and one regional histogram; not a weight certificate."""

from fractions import Fraction as F
from hashlib import sha256
import json
from math import comb, log2
from pathlib import Path
from time import monotonic

import numpy as np
from scipy.optimize import minimize
from flint import ctx

import hd_field128_gamma_all_pair_screen as screen
import hd_field128_gamma_all_pair_certificate as directed


def run():
    source = Path('hd_field128_gamma30_d40_r3072_fresh_midpoint_float_20261010.json')
    output = Path('hd_field128_gamma_r3072_f1463_four_marker_diagnostic_20261010.json')
    if output.exists():
        raise FileExistsError(output)
    raw = source.read_bytes()
    saved = json.loads(raw)
    z, alpha = F(saved['witness']['z']), F(saved['witness']['alpha'])
    weight, full, single, empty = 3072, 1463, 146, 2487
    assert saved['active_groups'] == weight and saved['state_dimension'] == 30
    histogram = comb(4096, full)*comb(4096-full, single)
    normalizer = log2(histogram)
    pins = {**screen.source_pins(), **directed.source_pins(),
            Path(__file__).name: sha256(Path(__file__).read_bytes()).hexdigest()}
    started, budget = monotonic(), screen.old.Budget(60.)
    coefficients = screen.AliasCoefficients(float(z), float(alpha), state=30, _budget=budget)
    cache = {}
    def numerator(point):
        budget.check()
        key = tuple(map(float, point))
        if key not in cache:
            cache[key] = np.asarray(screen.alias.power(coefficients.pair_logs(key, _budget=budget), 64)).ravel()
        return cache[key]
    fits, marks = [], []
    for entry in range(4):
        seeds = []
        for row in saved['fits']:
            if row['entry'] == entry and row['full'] in (1344, 1520):
                x, y = map(lambda v: float(F(v)), row['marker_logs_exact'])
                sf, ff = row['single'], row['full']
                ef = 4096-weight+ff
                seeds.append((x+log2(single/empty/(sf/ef)), y+log2(full/empty/(ff/ef))))
        best = min(seeds, key=lambda point: numerator(point)[entry]-single*point[0]-full*point[1])
        best_value = float(numerator(best)[entry]-single*best[0]-full*best[1])
        def objective(point):
            nonlocal best, best_value
            value = float(numerator(point)[entry]-single*point[0]-full*point[1])
            if value < best_value:
                best, best_value = tuple(map(float, point)), value
            return value
        fitted = minimize(objective, best, method='Powell', bounds=[(-64., 8.)]*2,
                          options=dict(maxfev=256, xtol=1e-8, ftol=1e-10))
        marker = tuple(F(2.**v) for v in best)
        # Binary64 exp2 outputs are proposed exact positive rationals.
        # They are not asserted equal to an exact exponential of the log seed.
        marks.append(marker)
        fits.append(dict(entry=entry, marker_logs_exact=list(map(lambda v: str(F(v)), best)),
                         marker_positive_rationals=list(map(str, marker)),
                         normalized_positive_cap_FLOAT=best_value-normalizer,
                         calls=int(fitted.nfev), optimizer_success=bool(fitted.success)))
    print(json.dumps(dict(stage='four-marker-fit', fits=fits)), flush=True)
    odds = [tuple(log2(float(value)) for value in marker) for marker in marks]
    float_bank = screen.AliasBankScreen(float(z), float(alpha), odds, 0., state=30, distance=F(2, 5),
                     grid=128, floor_policy='positive-fallback', _coefficients=coefficients, _budget=budget)
    float_rows = []
    for index, (odds_row, table) in enumerate(zip(odds, float_bank.tables)):
        shift = normalizer+single*odds_row[0]+full*odds_row[1]
        alias = table.logs[:, single % 128, full % 128]
        reliable = table.reliable[:, single % 128, full % 128]
        positive = table.positive_numerator
        selected = np.where(reliable, alias, positive)-shift
        float_rows.append(dict(marker=index, reliable=list(map(bool, reliable)),
                          normalized_alias_logs_FLOAT=(alias-shift).tolist(),
                          normalized_positive_logs_FLOAT=(positive-shift).tolist(),
                          normalized_selected_logs_FLOAT=selected.tolist()))
    print(json.dumps(dict(stage='four-marker-FLOAT-cell', rows=float_rows)), flush=True)
    bank = directed.AliasBank(z, alpha, marks, state=30, distance=F(2, 5), grid=128,
                   precision=256, bits=80, max_seconds=120.,
                   progress=lambda row: print(json.dumps(row), flush=True))
    rows = []
    with ctx.workprec(256):
        for index, marker in enumerate(marks):
            denominator = histogram*directed.base.rational(marker[0])**single*directed.base.rational(marker[1])**full
            alias = bank.atlases[index].numerator(((single, single), (full, full)))
            positive = bank.numerators[index]
            selected = bank.numerator(index, ((single, single), (full, full)))
            def encode(values):
                return [str(directed.base.rational_ceiling(directed.tail.log2arb(directed.base.upper(x/denominator))))
                        for x in values]
            rows.append(dict(marker=index, normalized_alias_log2_ceilings=encode(alias),
                             normalized_positive_log2_ceilings=encode(positive),
                             normalized_selected_log2_ceilings=encode(selected)))
        best = [min(F(row['normalized_selected_log2_ceilings'][entry]) for row in rows) for entry in range(4)]
    bank.check_pins()
    current = {**screen.source_pins(), **directed.source_pins(),
               Path(__file__).name: sha256(Path(__file__).read_bytes()).hexdigest()}
    if current != pins or source.read_bytes() != raw:
        raise RuntimeError('Diagnostic source changed')
    result = dict(schema='gamma30-r3072-single-histogram-four-marker-diagnostic-v1',
              active_groups=weight, regional_full=full, regional_single=single, state_dimension=30,
              target_distance='2/5', source=dict(path=str(source.resolve()), sha256=sha256(raw).hexdigest()),
              z=str(z), alpha=str(alpha), exact_new_positive_markers=[list(map(str, row)) for row in marks],
              marker_proposal_rule='binary64 exp2 of fitted log marks, treated as new exact rational witnesses',
              fits=fits, float_marker_rows=float_rows, directed_marker_rows=rows,
              directed_best_regional_log2_ceilings=list(map(str, best)),
              prior_saved_regional_logs_FLOAT=saved['regions'][full],
              bank_metadata=bank.metadata(), source_sha256=pins,
              fitting_budget_seconds=60, directed_bank_budget_seconds=120,
              elapsed_seconds=monotonic()-started, full_profile_sum_evaluated=False,
              regional_histograms_evaluated=1, selected_weight_certificate=False, whole_code_certificate=False,
              scope='One regional histogram with four fresh markers; directed regional caps only, no weight probability bound')
    pin = screen.search.write_new(output, result)
    print(json.dumps(dict(stage='complete', output=str(output), best=list(map(str, best)), artifact=pin)), flush=True)


if __name__ == '__main__':
    run()
