"""One FLOAT bank, eight independent complete midpoint-profile weight screens."""

from fractions import Fraction as F
from hashlib import sha256
import json
from math import log2
from pathlib import Path
from time import monotonic

import hd_field128_gamma_midpoint_counts as mode


HERE = Path(__file__).resolve().parent
SOURCE = HERE/'hd_field128_gamma30_d40_r3072_adaptive80_proposal_20261010.json'
OUTPUT = HERE/'hd_field128_gamma30_d40_r3072_neighbor_float_20261010.json'
WEIGHTS = (2944, 3008, 3040, 3104, 3136, 3200, 3264, 3072)
screen = mode.near.provider


def run():
    if OUTPUT.exists():
        raise FileExistsError(OUTPUT)
    raw = SOURCE.read_bytes()
    proposal = json.loads(raw)
    if (proposal['schema'], proposal['active_groups'], proposal['state_dimension'], proposal['target_distance']) != (
            'gamma30-r3072-adaptive80-witness-only-v1', 3072, 30, '2/5'):
        raise ValueError('Unexpected proposal scope')
    pins = mode.source_pins()
    pins[Path(__file__).name] = sha256(Path(__file__).read_bytes()).hexdigest()
    for name, digest in proposal['source_sha256'].items():
        if name in pins and pins[name] != digest:
            raise ValueError('Conflicting input closure: '+name)
        pins[name] = digest
    def check_pins():
        if SOURCE.read_bytes() != raw:
            raise RuntimeError('Proposal changed')
        for name, digest in pins.items():
            if sha256((HERE/name).read_bytes()).hexdigest() != digest:
                raise RuntimeError('Source changed: '+name)
        source = proposal['source']
        if sha256(Path(source['path']).read_bytes()).hexdigest() != source['sha256']:
            raise RuntimeError('Literal-source artifact changed')
    check_pins()
    witness = proposal['witness']
    z, alpha = F(witness['z']), F(witness['alpha'])
    positive = [tuple(map(F, pair)) for pair in witness['exact_positive_marker_proposals']]
    if not (0 < z < 1 and 0 < alpha <= 1 and len(positive) == proposal['marker_count'] == 80
            and all(len(pair) == 2 and min(pair) > 0 for pair in positive)):
        raise ValueError('Invalid exact bank literals')
    odds = [tuple(log2(float(value)) for value in pair) for pair in positive]
    budget, started = screen.old.Budget(180.), monotonic()
    completed, pending, metadata, bank_seconds = [], None, None, None
    status, stage = 'complete', 'one-FLOAT-bank'
    print(json.dumps(dict(stage=stage, weights=WEIGHTS, marker_count=80, max_seconds=180)), flush=True)
    try:
        coeff = screen.AliasCoefficients(float(z), float(alpha), state=30, _budget=budget)
        bank = screen.AliasBankScreen(float(z), float(alpha), odds, 0., state=30, distance=F(2, 5),
                    grid=128, floor_policy='positive-fallback', _coefficients=coeff, _budget=budget)
        bank_seconds, metadata = monotonic()-started, bank.metadata()
        for weight in WEIGHTS:
            budget.check()
            pending, stage = weight, 'complete-midpoint-profile-sum'
            regions = bank.regional(weight)
            result = mode.support_sum(4096, weight, regions, max_singletons=64, budget=budget)
            prefix = bank.constant+bank.alpha*bank.density*weight
            result.update(active_groups=weight, log2_outer_prefix=prefix,
                          log2_upper_screen=prefix+result['log2_support_sum_screen'],
                          regional_alias_diagnostics=dict(bank.regional_diagnostics),
                          regional_logs_sha256=sha256(regions.tobytes()).hexdigest(),
                          complete_selected_weight_screen=True, elapsed_seconds=monotonic()-started)
            completed.append(result)
            pending = None
            print(json.dumps(dict(weight=weight, log2_upper_screen=result['log2_upper_screen'],
                                  largest_profile=result['largest_profile'], elapsed_seconds=result['elapsed_seconds'],
                                  diagnostics=result['regional_alias_diagnostics'])), flush=True)
    except TimeoutError:
        status = 'timeout'
    check_pins()
    output = dict(schema='gamma30-r3072-shared-bank-neighbor-float-v1', status=status,
                  construction=screen.gamma.construction(30), state_dimension=30, target_distance='2/5',
                  calibration_active_groups=3072, requested_weights=list(WEIGHTS), completed_weights=completed,
                  pending_weights=[r for r in WEIGHTS if r not in {row['active_groups'] for row in completed}],
                  incomplete_weight=pending, stopped_stage=stage, bank_metadata=metadata,
                  one_shared_float_bank=True, regional_normalizer_and_prefix_recomputed_per_weight=True,
                  source=dict(path=str(SOURCE), sha256=sha256(raw).hexdigest()),
                  literal_source=proposal['source'], source_sha256=dict(sorted(pins.items())),
                  witness=dict(z=str(z), alpha=str(alpha), exact_positive_marker_proposals=[list(map(str, row)) for row in positive],
                               marker_log_odds_binary64_exact=[[str(F(v)) for v in pair] for pair in odds],
                               grid=128, floor_policy='positive-fallback', maximum_refined_singletons=64),
                  marker_conversion='FLOAT log2 of binary64 conversion of the saved exact positive rationals; no marker refit',
                  max_seconds=180, max_weights=8, float_bank_seconds=bank_seconds, elapsed_seconds=monotonic()-started,
                  all_profiles_retained_for_completed_weights=True, common_potential_used=False,
                  old_numerical_bounds_reused=False, no_weight_interpolation=True, directed=False,
                  selected_weight_certificate=False, whole_code_certificate=False,
                  scope='FLOAT-only candidate calibration at explicitly completed weights; no contiguous-band claim')
    pin = screen.search.write_new(OUTPUT, output)
    print(json.dumps(dict(stage='saved', status=status, artifact=pin)), flush=True)


if __name__ == '__main__':
    run()
