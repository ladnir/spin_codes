"""Full midpoint-profile FLOAT replay of the exact V57 return refinement.

Refine only the original V52 witness0; retain V56 and all other witnesses.
No numerical score here is a directed or whole-code certificate.
"""
import argparse
from fractions import Fraction as F
from hashlib import sha256
import json
from math import log2
from pathlib import Path
import resource
import shutil
from time import monotonic

import numpy as np
import hd_field128_full_quartet_profile_v56 as previous
import hd_field128_return_tail_moments_v57 as candidate
import hd_field128_return_tail_screen_v58 as diagnostic

HERE = Path(__file__).resolve().parent
INPUT56 = 'hd_field128_full_quartet_profile_v56_20261010.json'
INPUT56_SHA = '571cd9edfd2caee34d0f6cde8f7a79b6fd97df2dd53f19e5b7d41eb7a75578a9'
INPUT58 = 'hd_field128_return_tail_screen_v58_20261010.json'
INPUT58_SHA = 'c97f1281cb485c7da223d051ae6d8bf0fc47dc898302d895a3ed725a84af6a1f'
PREVIOUS_SHA = '2f767032b6d6e46c47b0683cd0ae4f4320b43d4db0ebb19d05424fbc590c7552'
DIAGNOSTIC_SHA = 'a301e2f853f745cbb23dccdf3e31d61ab18890e56b7b9a78ea8467d013509c2a'
LIMIT, RESERVE = 48 << 10, 20 << 30


def source_pins():
    pins = previous.source_pins()
    for name, digest in candidate.source_pins().items():
        if name in pins and pins[name] != digest:
            raise RuntimeError('Conflicting dependency: '+name)
        pins[name] = digest
    for module, expected in ((previous, PREVIOUS_SHA), (diagnostic, DIAGNOSTIC_SHA)):
        path = Path(module.__file__)
        if sha256(path.read_bytes()).hexdigest() != expected:
            raise RuntimeError('Frozen source changed: '+path.name)
        pins[path.name] = expected
    pins[Path(__file__).name] = sha256(Path(__file__).read_bytes()).hexdigest()
    return dict(sorted(pins.items()))


def run(output_name):
    if HERE != Path('/tmp/spin-field128-proof-52GBmlZD/bundle/theory'):
        raise RuntimeError('Remote-only bounded coefficient job')
    if Path(output_name).name != output_name or not output_name.endswith('.json'):
        raise ValueError('Expected a result basename')
    output = HERE/output_name
    if output.exists() or shutil.disk_usage(HERE).free-LIMIT < RESERVE:
        raise ValueError('Exclusive output and planned reserve required')
    resource.setrlimit(resource.RLIMIT_AS, (2 << 30, 2 << 30))
    resource.setrlimit(resource.RLIMIT_CPU, (230, 235))
    resource.setrlimit(resource.RLIMIT_FSIZE, (LIMIT, LIMIT))
    control, coefficients = previous.control, previous.coefficients
    caps, gamma, alias = previous.caps, previous.gamma, previous.alias
    pins = source_pins()
    inputs = {}
    for name, digest in ((control.INPUT52, control.INPUT52_SHA),
                         (INPUT56, INPUT56_SHA), (INPUT58, INPUT58_SHA)):
        path = HERE/name
        if path.stat().st_size > 96 << 10:
            raise ValueError('Oversized input')
        raw = path.read_bytes()
        if sha256(raw).hexdigest() != digest:
            raise ValueError('Pinned result changed')
        inputs[name] = raw
    v52, v56, v58 = (json.loads(inputs[n]) for n in (control.INPUT52, INPUT56, INPUT58))
    if (v56['source_sha256'] != previous.source_pins() or not v56['all_profiles_retained']
            or v56['profiles_evaluated'] != 2049**2 or v56['witness_index_refined'] != 0
            or v58['source_sha256'][candidate.NAME] != diagnostic.PROVIDER_SHA):
        raise ValueError('Predecessor scope differs')
    start, budget = monotonic(), gamma.old.Budget(225.)
    specs = v52['witness_specs']
    exchange = control.previous.previous.exchange
    old = caps.prior.control.old
    regions, prefixes = [], []
    old_poly, old_zero = None, None
    for index, spec in enumerate(specs):
        budget.check()
        z, alpha = float(F(spec['z'])), float(F(spec['alpha']))
        odds = [tuple(log2(float(F(v))) for v in pair) for pair in spec['markers']]
        pair = old.Coefficients(caps.local_logs(z), z, alpha, budget)
        bank = alias.AliasBankScreen(z, alpha, odds, 0., q=1 << 128,
            distance=F(9, 20), grid=128, floor_policy='positive-fallback',
            max_table_bytes=384 << 20, _coefficients=pair, _budget=budget)
        if bank.odds != tuple(odds):
            raise ArithmeticError('Marker sequence changed')
        region = bank.regional(4096)
        if sha256(region.tobytes()).hexdigest() != v52['regional_sha256'][index]:
            raise ArithmeticError('Fresh V52 region mismatch')
        prefixes.append(bank.constant+bank.alpha*bank.density*4096)
        del bank, pair
        if index == 0:
            old_poly, _ = coefficients.all_quartet_coefficients(caps.local_logs(z), alpha,
                                                                 check=budget.check)
            for marker in odds:
                table = coefficients.quartet_alias_table(old_poly, marker, check=budget.check)
                np.minimum(region, exchange.marker_regional(table, marker), out=region)
            if sha256(region.tobytes()).hexdigest() != v56['refined_regional_sha256']:
                raise ArithmeticError('Fresh V56 quartet region mismatch')
            old_zero = region.copy()
        regions.append(region)
    count_module = exchange.previous.previous.counts_module
    counts = count_module.MidpointCounts(4096, 4096, max_singletons=64, budget=budget)
    baseline = exchange.envelope(regions, prefixes, counts, budget)
    if baseline != v56['final']:
        raise ArithmeticError('Fresh V56 complete profile envelope differs')
    z_exact = F(specs[0]['z'])
    z, alpha = float(z_exact), float(F(specs[0]['alpha']))
    if v58['z'] != specs[0]['z'] or v58['alpha'] != specs[0]['alpha']:
        raise ValueError('Local cap diagnostic has a different witness')
    logs, gains = caps.local_logs(z), []
    for j in range(65):
        budget.check()
        old_packet = candidate.provider.exact_packet(j, z_exact)
        new_packet = candidate.exact_packet(j, z_exact)
        for entry in range(4):
            actual = diagnostic.log_fraction(old_packet[entry])
            reference = logs[j].ravel()[entry]
            if (actual is None and reference != -np.inf) or (
                    actual is not None and abs(actual-reference) > 3.e-9):
                raise ArithmeticError('Exact old packet disagrees with FLOAT replay')
            if entry != 2 and new_packet[entry] != old_packet[entry]:
                raise ArithmeticError('Non-return entry changed')
        if new_packet[2] > old_packet[2]:
            raise ArithmeticError('Exact local return became worse')
        raw_return = new_packet[2]*((1 << 4096)-1)
        if diagnostic.rational_digest(raw_return) != v58['occupancies'][j]['selected_exact_sha256']:
            raise ArithmeticError('Independent V58 exact return digest differs')
        new_log = diagnostic.log_fraction(new_packet[2])
        if new_log is not None:
            gains.append(float(logs[j, 1, 0]-new_log))
            logs[j, 1, 0] = min(logs[j, 1, 0], new_log)
    odds = [tuple(log2(float(F(v))) for v in pair) for pair in specs[0]['markers']]
    pair = old.Coefficients(logs, z, alpha, budget)
    bank = alias.AliasBankScreen(z, alpha, odds, 0., q=1 << 128,
        distance=F(9, 20), grid=128, floor_policy='positive-fallback',
        max_table_bytes=384 << 20, _coefficients=pair, _budget=budget)
    if bank.constant+bank.alpha*bank.density*4096 != prefixes[0]:
        raise ArithmeticError('Outer factor or original-r normalization changed')
    np.minimum(regions[0], bank.regional(4096), out=regions[0])
    del bank, pair
    pair_only = exchange.envelope(regions, prefixes, counts, budget)
    polynomial, construction = coefficients.all_quartet_coefficients(logs, alpha, check=budget.check)
    if np.any(polynomial > old_poly+3.e-9):
        raise ArithmeticError('Local refinement increased a quartet coefficient')
    diagnostics = []
    for marker in odds:
        table = coefficients.quartet_alias_table(polynomial, marker, check=budget.check)
        np.minimum(regions[0], exchange.marker_regional(table, marker), out=regions[0])
        diagnostics.append(table.diagnostic)
    final = exchange.envelope(regions, prefixes, counts, budget)
    if final['score'] > pair_only['score']+2.e-6 or pair_only['score'] > baseline['score']+2.e-6:
        raise ArithmeticError('Retained full-profile bound worsened')
    if source_pins() != pins or any((HERE/n).read_bytes() != raw for n, raw in inputs.items()):
        raise RuntimeError('Sources/inputs changed during replay')
    result = dict(schema='return-tail-profile-v59-FLOAT', source_sha256=pins,
        inputs=[dict(path=n, sha256=sha256(raw).hexdigest()) for n, raw in inputs.items()],
        z=specs[0]['z'], alpha=specs[0]['alpha'], witness_index_refined=0,
        target_distance='9/20', active_groups=4096, all_profiles_retained=True,
        profiles_evaluated=2049**2, baseline=baseline, pair_refinement=pair_only, final=final,
        complete_profile_gain_bits=baseline['score']-final['score'],
        maximum_local_return_gain=max(gains), coefficient_construction=construction,
        refined_regional_sha256=sha256(regions[0].tobytes()).hexdigest(),
        coefficient_sha256=sha256(polynomial.tobytes()).hexdigest(),
        maximum_alias_row_error=max(d['normalized_row_sum_error'] for d in diagnostics),
        alias_underflow_terms=sum(d['normalized_input_underflow_terms'] for d in diagnostics),
        entry_min_only_same_z_alpha=True, original_r_normalization=True,
        conditional_on_principal_setup=True, setup_risk_added=False,
        construction_changed=False, directed=False, whole_code_certificate=False,
        elapsed_seconds=monotonic()-start)
    data = (json.dumps(result, separators=(',', ':'), sort_keys=True, allow_nan=False)+'\n').encode()
    if len(data) > LIMIT or shutil.disk_usage(HERE).free-len(data) < RESERVE:
        raise RuntimeError('Output/reserve exceeded')
    with output.open('xb') as stream:
        stream.write(data)
    return dict(path=output_name, bytes=len(data), sha256=sha256(data).hexdigest(),
        baseline=baseline['score'], pair=pair_only['score'], final=final['score'],
        gain_bits=result['complete_profile_gain_bits'], elapsed_seconds=result['elapsed_seconds'])


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', required=True)
    print(json.dumps(run(parser.parse_args().output), allow_nan=False), flush=True)
