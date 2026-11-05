"""Bounded proposal experiment for one swapped-Toeplitz bank at weight2496.

Fit coefficient markers once at the center, then fit one common state
potential with those matrices fixed. Evaluate every integer weight in the
declared gap using the sharper Toeplitz density and FLOAT count logs.
Saved proposals require fresh directed replay; this driver certifies nothing.
"""
from copy import deepcopy
from fractions import Fraction as F
import json
from math import isfinite
from pathlib import Path
from time import perf_counter

import numpy as np
from scipy.optimize import minimize_scalar

from hd_shifted_outer_density import TOEPLITZ_OUTER_MEASURE
from hd_shifted_potential_band_screen import FixedAtlasScreen
from hd_shifted_run_count_screen import run_log_screen_float


CENTER = 2496
FIRST, LAST = 2154, 2827
BASE_Z = 2.**(-CENTER/1000)
Z_CHOICES = (BASE_Z, BASE_Z*.98, BASE_Z*1.02)
CENTERS = tuple(sorted(set(range(0, CENTER//2+1, 32))
                       | {k for k in (1, 2, 4, 8, 16)}
                       | {CENTER//2-k for k in (1, 2, 4, 8, 16)}))
COUNT_METHOD = 'FLOAT-log-factorial-cyclic-run-counts'


class Bank2496Screen(FixedAtlasScreen):
    def evaluate(self, weight):
        if type(weight) is not int or not 1 <= weight <= 8192:
            raise ValueError('A nonzero integer support weight is required')
        regions = self.regional(weight)
        support_sum = run_log_screen_float(4096, weight, regions, self.log_ratio)
        prefix = self.alpha*self.log_density*weight+self.constant
        return dict(active_groups=weight, regional_histograms=len(regions),
                    log2_support_sum_screen=support_sum, log2_outer_prefix=prefix,
                    log2_upper_screen=prefix+support_sum,
                    covers_all_matching_profiles=True, covers_all_even_odd_splits=True,
                    directed=False, whole_code_certificate=False,
                    count_method=COUNT_METHOD)

    def witness(self):
        return dict(super().witness(), count_method=COUNT_METHOD,
                    requires_fresh_directed_replay=True)


def selected_interval(proposal, *, center=CENTER, threshold=-60.):
    """Select the contiguous passing component containing the fixed center."""
    rows = proposal['results']
    weights = [row['active_groups'] for row in rows]
    if (not rows or weights != list(range(weights[0], weights[-1]+1))
            or center not in weights or not isfinite(threshold)
            or any(not isfinite(row['log2_upper_screen']) for row in rows)):
        raise ValueError('A complete finite integer-weight screen containing the center is required')
    position = center-weights[0]
    if rows[position]['log2_upper_screen'] > threshold:
        return None
    first = last = position
    while first > 0 and rows[first-1]['log2_upper_screen'] <= threshold:
        first -= 1
    while last+1 < len(rows) and rows[last+1]['log2_upper_screen'] <= threshold:
        last += 1
    selected = deepcopy(proposal)
    selected_rows = selected['results'] = selected['results'][first:last+1]
    union = float(np.logaddexp2.reduce([row['log2_upper_screen'] for row in selected_rows]))
    selected.update(requested_weight_range=[weights[first], weights[last]],
                    weights_evaluated=len(selected_rows), max_weights=len(selected_rows),
                    threshold=threshold, failed_weights=[], all_weights_pass_threshold=True,
                    union_log2_upper_screen=union, union_passes_threshold=union <= threshold,
                    selection='contiguous passing component containing the fixed center')
    return selected


def freeze(path, value):
    encoded = json.dumps(value, allow_nan=False, separators=(',', ':'))
    with path.open('x', encoding='utf-8') as target:
        target.write(encoded+'\n')


def run(choice, output, *, threshold=-60.):
    if type(choice) is not int or not 0 <= choice < len(Z_CHOICES):
        raise ValueError('Only the three declared z choices are admissible')
    output = Path(output)
    selected_path = output.with_name(output.stem+'_selected.json')
    if output.exists() or selected_path.exists():
        raise ValueError('Frozen artifacts are never overwritten')
    emit = lambda item: print(json.dumps(item, allow_nan=False), flush=True)
    started = perf_counter()
    z, alpha = Z_CHOICES[choice], .1
    emit(dict(stage='2496-fit-start', choice=choice, z=z, alpha=alpha,
              markers=len(CENTERS), grid=128, outer_measure=TOEPLITZ_OUTER_MEASURE))
    screen = Bank2496Screen.fit_central_once(CENTER, z, alpha, CENTERS, -10.8,
                                            grid=128, outer_measure=TOEPLITZ_OUTER_MEASURE)
    emit(dict(stage='2496-markers-frozen', seconds=perf_counter()-started,
              markers=len(screen.tables)))
    regions = screen.regional(CENTER)
    objective = lambda ratio: run_log_screen_float(4096, CENTER, regions, ratio)
    fitted = minimize_scalar(objective, bounds=(-18., -4.), method='bounded',
                             options={'xatol': 1.e-7, 'maxiter': 64})
    screen.log_ratio = float(fitted.x)
    screen.source_witness.update(potential=str(F(2.**screen.log_ratio)),
        potential_fit='one central scalar optimization, no neighboring refits',
        potential_fit_evaluations=int(fitted.nfev), z_choice=choice,
        z_candidates=[str(F(value)) for value in Z_CHOICES],
        z_rule='binary64 exp2(-2496/1000), optionally multiplied by .98 or1.02',
        central_alpha_rule='fixed binary64 .1')
    emit(dict(stage='2496-potential-frozen', seconds=perf_counter()-started,
              evaluations=int(fitted.nfev), central=screen.evaluate(CENTER)))
    count = 0
    def progress(row):
        nonlocal count
        count += 1
        if count % 32 == 0 or row['active_groups'] == LAST:
            emit(dict(stage='2496-gap-screen', completed=count, total=LAST-FIRST+1,
                      weight=row['active_groups'], value=row['log2_upper_screen'],
                      seconds=perf_counter()-started))
    proposal = screen.band(FIRST, LAST, threshold=threshold,
                            max_weights=LAST-FIRST+1, progress=progress)
    proposal.update(count_method=COUNT_METHOD, requires_fresh_directed_replay=True,
                    central_weight=CENTER, z_choice=choice, z_choice_limit=3,
                    elapsed_seconds=perf_counter()-started)
    freeze(output, proposal)
    selected = selected_interval(proposal, threshold=threshold)
    if selected is not None:
        selected.update(source_screen_file=output.name)
        freeze(selected_path, selected)
    emit(dict(stage='2496-saved', screen=str(output), selected=str(selected_path) if selected else None,
              selected_range=selected['requested_weight_range'] if selected else None,
              selected_union=selected['union_log2_upper_screen'] if selected else None,
              seconds=perf_counter()-started, directed=False))
    return selected


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--choice', type=int, choices=range(3), default=0)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    run(args.choice, args.output)
