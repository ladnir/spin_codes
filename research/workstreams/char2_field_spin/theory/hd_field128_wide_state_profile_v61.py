"""Separate complete midpoint screens for wider-state construction candidates.

Every state uses its own fresh local caps and seven old numerical witnesses.
No bound transfers between constructions. FLOAT proposals, not certificates.
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
import hd_field128_wide_state_v60 as provider
import hd_field128_full_quartet_profile_v56 as framework

HERE = Path(__file__).resolve().parent
SOURCE = 'hd_field128_reversed_dual_marker_v52_20261010.json'
SOURCE_SHA = '1ec14fce6aede5f7f45d301fdeae5cb7057670a09f8e244bad52e9bd1b69b3e5'
LIMIT, RESERVE = 48 << 10, 20 << 30


def source_pins():
    pins = framework.source_pins()
    for name, digest in provider.source_pins().items():
        if name in pins and pins[name] != digest:
            raise RuntimeError('Conflicting source identity')
        pins[name] = digest
    pins[Path(__file__).name] = sha256(Path(__file__).read_bytes()).hexdigest()
    return dict(sorted(pins.items()))


def run(output_name):
    if HERE != Path('/tmp/spin-field128-proof-52GBmlZD/bundle/theory'):
        raise RuntimeError('Remote-only bounded numerical screen')
    if Path(output_name).name != output_name or not output_name.endswith('.json'):
        raise ValueError('Expected result basename')
    output = HERE/output_name
    if output.exists() or shutil.disk_usage(HERE).free-LIMIT < RESERVE:
        raise ValueError('Exclusive output and planned reserve required')
    resource.setrlimit(resource.RLIMIT_AS, (2 << 30, 2 << 30))
    resource.setrlimit(resource.RLIMIT_CPU, (550, 555))
    resource.setrlimit(resource.RLIMIT_FSIZE, (LIMIT, LIMIT))
    path = HERE/SOURCE
    if path.stat().st_size > 96 << 10:
        raise ValueError('Oversized witness input')
    raw = path.read_bytes()
    if sha256(raw).hexdigest() != SOURCE_SHA:
        raise ValueError('Frozen witness input changed')
    specs, pins = json.loads(raw)['witness_specs'], source_pins()
    control, caps, alias, gamma = framework.control, framework.caps, framework.alias, framework.gamma
    exchange, old = control.previous.previous.exchange, caps.prior.control.old
    budget, start, rows = gamma.old.Budget(540.), monotonic(), []
    for state in provider.STATES:
        regions, prefixes, local_hashes = [], [], []
        for index, spec in enumerate(specs):
            budget.check()
            z_exact, alpha_exact = F(spec['z']), F(spec['alpha'])
            z, alpha = float(z_exact), float(alpha_exact)
            logs = np.full((65, 2, 2), -np.inf)
            for j in range(65):
                budget.check()
                packet = provider.exact_packet(j, z_exact, state)
                if len(packet) != 4 or any(value < 0 for value in packet):
                    raise ArithmeticError('Invalid exact packet')
                for entry, value in enumerate(packet):
                    if value:
                        logs[j].ravel()[entry] = log2(value.numerator)-log2(value.denominator)
            if np.any(np.isnan(logs)) or np.any(np.isposinf(logs)):
                raise ArithmeticError('Invalid local logarithm')
            local_hashes.append(sha256(logs.tobytes()).hexdigest())
            odds = [tuple(log2(float(F(v))) for v in pair) for pair in spec['markers']]
            pair = old.Coefficients(logs, z, alpha, budget)
            bank = alias.AliasBankScreen(z, alpha, odds, 0., q=1 << 128,
                distance=F(9, 20), grid=128, floor_policy='positive-fallback',
                max_table_bytes=384 << 20, _coefficients=pair, _budget=budget)
            if bank.odds != tuple(odds):
                raise ArithmeticError('Marker indexing changed')
            regions.append(bank.regional(4096))
            prefixes.append(bank.constant+bank.alpha*bank.density*4096)
            del bank, pair
        counts_module = exchange.previous.previous.counts_module
        counts = counts_module.MidpointCounts(4096, 4096, max_singletons=64, budget=budget)
        envelope = exchange.envelope(regions, prefixes, counts, budget)
        if envelope['profiles_evaluated'] != 2049**2:
            raise ArithmeticError('Incomplete midpoint profile domain')
        row = dict(state=state, final=envelope, local_log_sha256=local_hashes,
            regional_sha256=[sha256(region.tobytes()).hexdigest() for region in regions],
            principal_bad_gamma_upper=str(provider.principal_bad_gamma_upper(state)),
            emission_minimum_distance=65-state)
        rows.append(row)
        print(json.dumps(dict(progress_state=state, score=envelope['score'],
                              elapsed_seconds=monotonic()-start)), flush=True)
        del regions, counts
    if source_pins() != pins or path.read_bytes() != raw:
        raise RuntimeError('Source/input changed during screen')
    result = dict(schema='wide-state-midpoint-v61-FLOAT', source_sha256=pins,
        input=dict(path=SOURCE, sha256=SOURCE_SHA), states=list(provider.STATES), rows=rows,
        target_distance='9/20', active_groups=4096, witnesses=len(specs),
        all_profiles_retained=True, profiles_per_state=2049**2,
        original_r_normalization=True, witness_min_after_own_full_transfer=True,
        cross_construction_minimum=False, old_construction_bounds_reused=False,
        conditional_on_principal_setup=True, setup_risk_added=False,
        directed=False, whole_code_certificate=False, encoder_implemented=False,
        elapsed_seconds=monotonic()-start)
    data = (json.dumps(result, sort_keys=True, separators=(',', ':'), allow_nan=False)+'\n').encode()
    if len(data) > LIMIT or shutil.disk_usage(HERE).free-len(data) < RESERVE:
        raise RuntimeError('Output/reserve budget exceeded')
    with output.open('xb') as stream:
        stream.write(data)
    return dict(path=output_name, bytes=len(data), sha256=sha256(data).hexdigest(),
        scores=[(row['state'], row['final']['score']) for row in rows],
        elapsed_seconds=result['elapsed_seconds'])


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', required=True)
    print(json.dumps(run(parser.parse_args().output), allow_nan=False), flush=True)
