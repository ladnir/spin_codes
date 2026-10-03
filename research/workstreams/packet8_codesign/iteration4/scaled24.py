"""Fresh A-only byte scalings of the literal GF(256)^3 construction.

All numerical gates are floating proposals, not outward certificates.  The
C-character tables are intentionally unchanged when the expansion A changes.
"""
from __future__ import annotations

import os
for _name in ('OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'OMP_NUM_THREADS'):
    os.environ[_name] = '1'

import argparse
from fractions import Fraction
import hashlib
import json
from math import comb, log
from pathlib import Path
import sys
from time import monotonic

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent/'larger_state'))
import maps24
import screen24
sys.path.insert(0, str(HERE.parent/'iteration3'))
import cached_fractional24 as base
import potential_birth_gate as birth


def adjoint_tables(tables, bits):
    result = np.zeros_like(tables)
    for h in range(tables.shape[0]):
        for k in range(3):
            for value in range(1 << bits):
                result[h, k, value] = sum(
                    ((int(tables[h, k, 1 << i]) & value).bit_count() & 1) << i
                    for i in range(bits))
    return result


def make_scaled(scales, *, packet_bits=8, modulus=0x11b, windows=8):
    scales = tuple(scales)
    if len(scales) != windows or any(type(s) is not int or not 0 < s < 1 << packet_bits for s in scales):
        raise ValueError('one nonzero field element per output packet required')
    original = maps24.make_maps(packet_bits=packet_bits, modulus=modulus, windows=windows)
    data = dict(original)
    mul = lambda a, b: maps24.multiply(a, b, bits=packet_bits, modulus=modulus)
    lookup = np.array([[mul(scale, value) for value in range(1 << packet_bits)]
                       for scale in scales], dtype=np.uint8)
    tables = np.array([lookup[h, original['tables'][h]] for h in range(windows)])
    rows = tuple(sum(int(tables[h, i//packet_bits, 1 << (i % packet_bits)]) << (packet_bits*h)
                     for h in range(windows)) for i in range(3*packet_bits))
    identity = dict(packet_bits=packet_bits, windows=windows, state_bits=3*packet_bits,
        modulus=hex(modulus), expansion_rows_hex=list(map(hex, rows)),
        feedback_columns=list(original['columns']))
    record = dict(original['record'])
    record.update(identity)
    record.update(map_sha256=hashlib.sha256(json.dumps(identity, sort_keys=True).encode()).hexdigest(),
        expansion_rank=maps24.rank(rows), CA_columns=[maps24.apply(data['columns'], row) for row in rows],
        output_byte_scales=list(scales), base_map_sha256=original['record']['map_sha256'],
        construction_change='A-only nonzero output-byte scales; C, update, order, outer and route unchanged')
    record['CA_rank'] = maps24.rank(record['CA_columns'])
    if record['expansion_rank'] != 3*packet_bits:
        raise ArithmeticError('scaled expansion must retain full state rank')
    data.update(rows=rows, tables=tables, record=record,
                expansion_transpose_tables=adjoint_tables(tables, packet_bits))
    # This is C^T, not A^T. Do not replace it with the new expansion adjoints.
    if data['columns'] != original['columns'] or not np.array_equal(data['transpose_tables'], original['transpose_tables']):
        raise ArithmeticError('A-only construction unexpectedly changed C')
    return data


def candidate_scales(name):
    if name == 'baseline':
        return (1,)*8
    if name.startswith('bands'):
        return (1,)*4+(int(name[5:], 16),)*4
    if name == 'power3':
        result, value = [], 1
        for _ in range(8):
            result.append(value)
            value = maps24.multiply(value, 3)
        return tuple(result)
    raise ValueError('unknown bounded candidate')


def scale_score(value):
    weights = [a.bit_count()+maps24.multiply(value, a).bit_count() for a in range(1, 256)]
    return dict(scale=value, minimum=min(weights), histogram=np.bincount(weights, minlength=17).tolist())


def sources():
    pins = screen24.sources()
    pins.update(base.fine.source_pins())
    pins.update(birth.source_pins())
    for path in (Path(__file__).resolve(), HERE/'test_scaled24.py', Path(base.__file__).resolve(),
                 Path(base.cs.cg.__file__).resolve(), Path(base.cs.cg.variant_gate.__file__).resolve()):
        pins[str(path)] = hashlib.sha256(path.read_bytes()).hexdigest()
    return pins


def gate(local, theta=Fraction(3, 50), alpha=Fraction(2, 5), q=16):
    old = base.cs.cg.gate.potential_operators(local, 0.)
    potential, change, algebra = birth.rebase_potential(old)
    with np.errstate(under='raise', over='raise', invalid='raise', divide='raise'):
        tuples = base.fine.tuple_products(potential, 4)
        macro = base.fine.fine_operators(tuples, alpha)
        alpha1 = base.fine.fine_operators(tuples, 1.)
    coarse, _ = base.fine.coarse.grouped_operators(potential, 4)
    positive = coarse > 0
    error = float(np.max(np.abs(alpha1[positive]/coarse[positive]-1)))
    if np.any(alpha1[~positive] != 0) or error > 2e-11:
        raise ArithmeticError('alpha1 grouping regression failed')
    ordinary, _ = base.fine.prior.placement(potential, q, epochs=64, windows=8, force_log=True)
    grouped, _ = base.fine.prior.placement(alpha1, q, epochs=16, windows=32, force_log=True)
    mask = np.isfinite(ordinary) & np.isfinite(grouped)
    regional_error = float(np.max(np.abs(ordinary[mask]-grouped[mask])))
    if not np.array_equal(np.isfinite(ordinary), np.isfinite(grouped)) or regional_error > 2e-8:
        raise ArithmeticError('alpha1 regional regression failed')
    fractional, backend = base.fine.prior.placement(macro, q, epochs=16, windows=32, force_log=True)
    moment = base.fine.prior.logarithmic.log_power_matrix(fractional[q], 32)
    beta = base.fine.prior.rs_uniform_envelope.UniformInputEnvelope(16, 8, 8, 2).beta
    exponent = log(comb(512, q))+float(alpha)*(q*(log(beta.numerator)-log(beta.denominator))+float(theta)*13107)+moment
    return dict(q=q, tilt=str(theta), alpha=str(alpha), margin_bits=-exponent/log(2),
        fractional_log_moment=moment, backend=backend, alpha1_macro_relative_error=error,
        alpha1_regional_log_error=regional_error, algebra_checks=algebra,
        change_of_basis=change.tolist(), strict_numpy_underflow_check=True,
        powered_macro_operators=macro.tolist())


def run(args):
    if args.output_dir.exists():
        raise ValueError('fresh output directory required')
    args.output_dir.mkdir(parents=True)
    pins, started = sources(), monotonic()
    baseline_path = HERE.parent/'iteration3/cache24_v2/theta_3_50.json'
    baseline, baseline_saved = base.cs.cg.load_local(baseline_path)
    pins.update(baseline_saved['source_sha256'])
    pins[str(baseline_path)] = hashlib.sha256(baseline_path.read_bytes()).hexdigest()
    for name in args.variants:
        beginning = monotonic()
        data = make_scaled(candidate_scales(name))
        counts = maps24.census(data)
        census_seconds = monotonic()-beginning
        print(f'{name}: fresh census {census_seconds:.2f}s; minweight={np.flatnonzero(counts["spectrum"])[1]}, '
              f'profiles={np.count_nonzero(counts["histogram"])}, CA_rank={data["record"]["CA_rank"]}', flush=True)
        local, diagnostic = screen24.local_operators(data, counts, np.exp(-.06), progress=True)
        regression = None
        if name == 'baseline':
            if data['record']['map_sha256'] != baseline_saved['map_record']['map_sha256']:
                raise ArithmeticError('all-one literal identity regression failed')
            positive = baseline > 0
            regression = dict(bit_identical=bool(np.array_equal(local, baseline)),
                maximum_relative_error=float(np.max(np.abs(local[positive]/baseline[positive]-1))))
            if np.any(local[~positive] != 0) or regression['maximum_relative_error'] > 2e-10:
                raise ArithmeticError('all-one local moment regression failed')
        if not base.cs.cg.variant_gate.checked_pins(pins):
            raise ArithmeticError('source changed during the fresh census')
        cache = dict(schema='packet8-scaled-A-actual24-local-proposal-1', variant=name,
            map_record=data['record'], tilt='3/50', local_operator_matrices=local.tolist(),
            expansion_spectrum={str(i): int(n) for i, n in enumerate(counts['spectrum']) if n},
            expansion_profile_count=int(np.count_nonzero(counts['histogram'])),
            source_sha256=pins, source_pins_verified_at_finish=True, census_seconds=census_seconds,
            local_diagnostics=diagnostic, baseline_regression=regression,
            proposal_only=True, whole_code_certificate=False, has_outward_endpoints=False)
        cache_path = args.output_dir/(name+'_local.json')
        cache_path.write_text(json.dumps(cache, indent=2)+'\n', encoding='utf-8')
        del counts
        result = gate(local)
        receipt = dict(schema='packet8-scaled-A-actual24-G4-proposal-1', variant=name,
            map_record=data['record'], result=result,
            geometry=dict(K=65536, N=131072, groups=512, regions=32, physical_steps_per_region=64,
                          state_bits=24, group_steps=4, cutoff=13107, final_flush=False),
            local_receipt=str(cache_path.resolve()),
            source_sha256=dict(pins, **{str(cache_path.resolve()): hashlib.sha256(cache_path.read_bytes()).hexdigest()}),
            source_pins_verified_at_finish=True, elapsed_seconds=monotonic()-beginning,
            proposal_only=True, whole_code_certificate=False, has_outward_endpoints=False,
            subset_union_is_not_powered=True, message_factor='C(512,q)*beta^(q*alpha)',
            power_order='potential thinning; rebase births; ordered G4 product; entrywise alpha; tuple averaging')
        if not base.cs.cg.variant_gate.checked_pins(receipt['source_sha256']):
            raise ArithmeticError('source changed during the gate')
        (args.output_dir/(name+'_gate.json')).write_text(json.dumps(receipt, indent=2)+'\n', encoding='utf-8')
        print(f'{name}: q16 theta=.06 alpha=.4 G4 = {result["margin_bits"]:.9f} bits; '
              f'elapsed={monotonic()-started:.2f}s', flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--variants', nargs='+', default=['baseline', 'bands17', 'bands2e', 'power3'])
    parser.add_argument('--output-dir', type=Path, required=True)
    run(parser.parse_args())
