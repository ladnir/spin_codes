"""Complete FLOAT midpoint-weight profile sum for the V46 witness.

Retains all4,198,401 profile pairs at r4096 and uses the original r in
coefficients/prefix. This remains one weight, not a directed full-code cert.
"""
import argparse
from fractions import Fraction as F
from hashlib import sha256
import json
from math import log2
from pathlib import Path
import shutil
from time import monotonic

import hd_field128_reversed_marker_retune_v46 as previous
import hd_field128_gamma_midpoint_counts as counts_module

HERE = Path(__file__).resolve().parent
PREVIOUS_SHA = 'fd6ce9c080b126b104590ab839545f0a795ae6470d7bd9ae28d6e4ea0a25439d'
INPUT_SHA = 'f7314a8fcb5c627a1818db943dc94d94f4f03648fb89f13c6f5348bc654f6dd1'
LIMIT, RESERVE = 32 << 10, 20 << 30


def source_pins():
    pins = previous.source_pins()
    if pins[Path(previous.__file__).name] != PREVIOUS_SHA:
        raise RuntimeError('V46 changed')
    for name, digest in counts_module.source_pins().items():
        if name in pins and pins[name] != digest:
            raise RuntimeError('Conflicting source')
        pins[name] = digest
    pins[Path(__file__).name] = sha256(Path(__file__).read_bytes()).hexdigest()
    return dict(sorted(pins.items()))


def run(source_name, output_name):
    with previous.previous.prior.prior.control.remote_limits():
        return _run(source_name, output_name)


def _run(source_name, output_name):
    import hd_field128_full_rs_gamma_screen as gamma
    import hd_field128_lift_alias_screen as alias
    control = previous.previous.prior.prior.control
    old = control.old
    source, output = old.prior.previous.inside(source_name), old.prior.previous.inside(output_name)
    if output.exists() or source.stat().st_size > LIMIT or shutil.disk_usage(HERE).free-LIMIT < RESERVE:
        raise ValueError('Storage/input guard')
    raw, pins = source.read_bytes(), source_pins()
    if sha256(raw).hexdigest() != INPUT_SHA:
        raise ValueError('Pinned V46 required')
    artifact = json.loads(raw)
    if artifact['active_groups'] != 4096 or artifact['target_distance'] != '9/20':
        raise ValueError('Scope mismatch')
    witness = artifact['best']
    z, alpha = float(F(witness['z'])), float(F(witness['alpha']))
    odds = [tuple(log2(float(F(v))) for v in row) for row in artifact['markers']]
    started, budget = monotonic(), gamma.old.Budget(80.)
    logs = previous.previous.prior.local_logs(z)
    coefficient = old.Coefficients(logs, z, alpha, budget)
    bank = alias.AliasBankScreen(z, alpha, odds, 0., q=1 << 128, distance=F(9, 20),
        grid=128, floor_policy='positive-fallback', max_table_bytes=384 << 20,
        _coefficients=coefficient, _budget=budget)
    regions = bank.regional(4096)
    if regions.shape != (2049, 2, 2):
        raise ArithmeticError('Original-r regional domain')
    counts = counts_module.MidpointCounts(4096, 4096, max_singletons=64, budget=budget)
    selected_count = float(counts.log_counts(2026, 2026, budget=budget))
    prefix = bank.constant+bank.alpha*bank.density*4096
    expected_prefix = old.prior.previous.prefactor(z, alpha, F(9, 20), 4096, bank.density)
    selected = prefix+selected_count+old.prior.profile_transfer(regions[2026], regions[2026])
    if abs(prefix-expected_prefix) > 1.e-7 or abs(selected-witness['score']) > 1.e-4:
        raise ArithmeticError('Prefix or independent selected-profile replay mismatch')
    row = counts_module.support_sum(4096, 4096, regions, counts=counts, budget=budget)
    if row['profile_pairs_evaluated'] != 2049**2 or not row['all_profiles_retained']:
        raise ArithmeticError('Incomplete profiles')
    result = dict(row, schema='reversed-allprofiles-v47-FLOAT', input_sha256=INPUT_SHA,
        source_sha256=pins, active_groups=4096, target_distance='9/20',
        conditional_log2_upper_screen=prefix+row['log2_support_sum_screen'],
        log2_outer_prefix=prefix, selected_control=selected, common_potential=False,
        witness=dict(z=witness['z'], alpha=witness['alpha'], markers=artifact['markers'], grid=128),
        regional_sha256=sha256(regions.tobytes()).hexdigest(),
        local_logs_sha256=sha256(logs.tobytes()).hexdigest(),
        coefficient_sha256=sha256(coefficient.composed.tobytes()).hexdigest(),
        alias_diagnostics=bank.regional_diagnostics,
        original_r_for_coefficients_and_prefix=True, conditional_on_principal_setup=True,
        state32_translated_feedback=True, setup_risk_added=False,
        complete_selected_weight_FLOAT=True, elapsed_seconds=monotonic()-started)
    if source.read_bytes() != raw or source_pins() != pins:
        raise RuntimeError('Source changed')
    data = (json.dumps(old.prior.previous.json_ready(result), sort_keys=True,
                       separators=(',', ':'), allow_nan=False)+'\n').encode()
    if len(data) > LIMIT or shutil.disk_usage(HERE).free-len(data) < RESERVE:
        raise RuntimeError('Storage cap')
    with output.open('xb') as stream:
        stream.write(data)
    return dict(path=output_name, sha256=sha256(data).hexdigest(), bytes=len(data),
                score=result['conditional_log2_upper_screen'], largest_profile=row['largest_profile'],
                profile_pairs_evaluated=row['profile_pairs_evaluated'], elapsed_seconds=result['elapsed_seconds'])


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source', required=True)
    parser.add_argument('--output', required=True)
    args = parser.parse_args()
    print(json.dumps(run(args.source, args.output), allow_nan=False), flush=True)
