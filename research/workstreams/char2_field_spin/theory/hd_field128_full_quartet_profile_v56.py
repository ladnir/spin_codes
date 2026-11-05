"""Full quartet coefficients in a complete r4096 FLOAT profile replay.

Refine only V52 witness0, at its original z/alpha, and retain the six other
complete witnesses. A conditional numerical screen is not a certificate.
"""
import argparse
import ast
from fractions import Fraction as F
from hashlib import sha256
import json
from math import log2
from pathlib import Path
import shutil
from time import monotonic

import numpy as np
import hd_field128_translated_full_quartet_v54 as control
import hd_field128_full_quartet_coefficients_v55 as coefficients
import hd_field128_reversed_return_v44 as caps
import hd_field128_full_rs_gamma_screen as gamma
import hd_field128_lift_alias_screen as alias

HERE = Path(__file__).resolve().parent
CONTROL_SHA = '00d9516a068d9f615caa4c53cb566c62c939ee399f50eb3c5539e75ca090f88c'
COEFFICIENT_SHA = '3dc670779dd085261f682fdb4d201cf90c79304c61ba6a122b651a675026864c'
INPUT54 = 'hd_field128_translated_full_quartet_v54_20261010.json'
INPUT54_SHA = '5545a3712b090fc2392a52392ef85356c3ff6bcbe36da8db45e20ffbf9039653'
LIMIT, RESERVE = 32 << 10, 20 << 30


def source_pins():
    pins = control.source_pins()
    if pins.get(Path(control.__file__).name) != CONTROL_SHA:
        raise RuntimeError('V54 changed')
    pending, seen = [Path(__file__).name], set()
    while pending:
        name = pending.pop()
        if name in seen:
            continue
        seen.add(name)
        raw = (HERE/name).read_bytes()
        digest = sha256(raw).hexdigest()
        if name in pins and pins[name] != digest:
            raise RuntimeError('Conflicting source: '+name)
        pins[name] = digest
        for node in ast.walk(ast.parse(raw)):
            modules = ([v.name for v in node.names] if isinstance(node, ast.Import)
                       else [node.module] if isinstance(node, ast.ImportFrom) and node.module else [])
            for module in modules:
                child = module.split('.')[0]+'.py'
                if (HERE/child).is_file() and child not in seen:
                    pending.append(child)
    if pins.get(Path(coefficients.__file__).name) != COEFFICIENT_SHA:
        raise RuntimeError('V55 changed')
    return dict(sorted(pins.items()))


def run(output_name):
    with control.guard.remote_limits():
        return _run(output_name)


def _run(output_name):
    output = caps.prior.control.old.prior.previous.inside(output_name)
    if output.exists() or shutil.disk_usage(HERE).free-LIMIT < RESERVE:
        raise ValueError('Exclusive output and planned reserve required')
    inputs = {}
    for name, digest in ((control.INPUT52, control.INPUT52_SHA),
                         (control.INPUT53, control.INPUT53_SHA), (INPUT54, INPUT54_SHA)):
        path = HERE/name
        if path.stat().st_size > 96 << 10:
            raise ValueError('Oversized input')
        raw = path.read_bytes()
        if sha256(raw).hexdigest() != digest:
            raise ValueError('Pinned input changed')
        inputs[name] = raw
    pins = source_pins()
    v52, v53, v54 = (json.loads(inputs[name]) for name in
                      (control.INPUT52, control.INPUT53, INPUT54))
    control.validate_inputs(v52, v53)
    if (v54['source_sha256'] != control.source_pins()
            or v54['selected_profile'] != list(control.PROFILE)
            or v54['witness_index'] != 0 or not v54['all_occupancies_retained']
            or v54['ordered_quartets_per_marker'] != 65**4):
        raise ValueError('Full quartet control scope differs')
    started, budget = monotonic(), gamma.old.Budget(165.)
    exchange = control.previous.previous.exchange
    old = caps.prior.control.old
    specs = v52['witness_specs']
    z, alpha = float(F(specs[0]['z'])), float(F(specs[0]['alpha']))
    build_start = monotonic()
    polynomial, construction = coefficients.all_quartet_coefficients(
        caps.local_logs(z), alpha, check=budget.check)
    build_seconds = monotonic()-build_start
    point_errors = []
    for marker in v54['markers']:
        odds = tuple(log2(float(F(v))) for v in marker['marker'])
        actual = coefficients.point_logs(polynomial, odds)
        expected = np.asarray(marker['full_log_four'])
        error = float(np.max(abs(actual-expected)))
        if error > 2.e-8:
            raise ArithmeticError('Independent full65 point control differs')
        point_errors.append(error)
    regions, prefixes, old_hashes = [], [], []
    alias_diagnostics = []
    old_zero, new_zero = None, None
    for index, spec in enumerate(specs):
        budget.check()
        zz, aa = float(F(spec['z'])), float(F(spec['alpha']))
        odds = [tuple(log2(float(F(v))) for v in pair) for pair in spec['markers']]
        pair = old.Coefficients(caps.local_logs(zz), zz, aa, budget)
        bank = alias.AliasBankScreen(zz, aa, odds, 0., q=1 << 128,
            distance=F(9, 20), grid=128, floor_policy='positive-fallback',
            max_table_bytes=384 << 20, _coefficients=pair, _budget=budget)
        if bank.odds != tuple(odds):
            raise ArithmeticError('Marker deduplication changed indices')
        region = bank.regional(4096)
        digest = sha256(region.tobytes()).hexdigest()
        if digest != v52['regional_sha256'][index]:
            raise ArithmeticError('V52 regional control differs')
        old_hashes.append(digest)
        prefixes.append(bank.constant+bank.alpha*bank.density*4096)
        del bank, pair
        if index == 0:
            old_zero = region.copy()
            quartet = np.full_like(region, np.inf)
            for marker in odds:
                table = coefficients.quartet_alias_table(polynomial, marker, check=budget.check)
                proposal = exchange.marker_regional(table, marker)
                np.minimum(quartet, proposal, out=quartet)
                alias_diagnostics.append(table.diagnostic)
            np.minimum(region, quartet, out=region)
            new_zero = region.copy()
        regions.append(region)
    count_module = exchange.previous.previous.counts_module
    counts = count_module.MidpointCounts(4096, 4096, max_singletons=64, budget=budget)
    # Replay the old complete envelope freshly before using the refinement.
    old_regions = [old_zero]+regions[1:]
    replay = exchange.envelope(old_regions, prefixes, counts, budget)
    if replay != v52['final']:
        raise ArithmeticError('Fresh full V52 envelope differs')
    final = exchange.envelope(regions, prefixes, counts, budget)
    if final['score'] > replay['score']+2.e-6:
        raise ArithmeticError('Refinement worsened complete profile bound')
    count = float(counts.log_counts(*control.PROFILE, budget=budget))
    selected_before = prefixes[0]+count+control.four.profile_transfer(*old_zero[list(control.PROFILE)])
    selected_after = prefixes[0]+count+control.four.profile_transfer(*new_zero[list(control.PROFILE)])
    control.close_control(selected_before, v54['scores']['alias_all_markers'], 'V54 alias')
    if source_pins() != pins or any((HERE/n).read_bytes() != raw for n, raw in inputs.items()):
        raise RuntimeError('Source/input changed')
    result = dict(schema='full-quartet-profile-v56-FLOAT', source_sha256=pins,
        inputs=[dict(path=n, sha256=sha256(b).hexdigest()) for n, b in inputs.items()],
        active_groups=4096, target_distance='9/20', witness_index_refined=0,
        z=specs[0]['z'], alpha=specs[0]['alpha'], coefficient_construction=construction,
        coefficient_log_sha256=sha256(polynomial.tobytes()).hexdigest(),
        coefficient_build_seconds=build_seconds, point_control_max_error=max(point_errors),
        point_controls=len(point_errors), old_regional_hashes=old_hashes,
        refined_regional_sha256=sha256(new_zero.tobytes()).hexdigest(),
        baseline=replay, final=final, selected_profile=list(control.PROFILE),
        selected_profile_before=selected_before, selected_profile_after=selected_after,
        complete_profile_gain_bits=replay['score']-final['score'],
        quartet_markers=len(alias_diagnostics),
        maximum_alias_row_error=max(d['normalized_row_sum_error'] for d in alias_diagnostics),
        alias_underflow_terms=sum(d['normalized_input_underflow_terms'] for d in alias_diagnostics),
        all_profiles_retained=True, profiles_evaluated=2049**2,
        conditional_on_principal_setup=True, setup_risk_added=False,
        entry_min_only_same_z_alpha=True, original_r_normalization=True,
        quartet_power=32, counts_outside_alpha=True,
        directed=False, selected_weight_certificate=False, whole_code_certificate=False,
        elapsed_seconds=monotonic()-started)
    receipt = control.write_result(output, result)
    return dict(receipt=receipt, score=final['score'], gain_bits=result['complete_profile_gain_bits'],
                coefficient_build_seconds=build_seconds,
                selected_profile_after=selected_after, elapsed_seconds=result['elapsed_seconds'])


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', required=True)
    args = parser.parse_args()
    print(json.dumps(run(args.output), allow_nan=False), flush=True)
