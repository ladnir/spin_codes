"""Selected-class directed replay for the measured H* outer.

This fresh replay covers r18, homogeneous arity8, and Fmax<=1 only.
Regional caps and atom/cover weights use Arb; chronological contraction
uses explicitly upward binary64 primitives with exact power-of-two scales.
The final roots, sum and global normalization use Arb again. Saved FLOAT
matrices are never certificate inputs. A successful replay is NOT a
whole-code certificate. The CLI refuses to replace existing receipts.
"""
from fractions import Fraction as F
from hashlib import sha256
import json
from math import nextafter
from pathlib import Path
import platform
from time import perf_counter

import numpy as np
from flint import arb, ctx

import hd_positive_upward as up
from hd_bit_syndrome_atoms_directed import (build_atoms, export_scaled,
    power_upper, primitive_weights, rational, upper_fraction)
from hd_bit_syndrome_columns_directed import counted_columns
from hd_bit_syndrome_gate import cover_gap_geometry, row_geometry
from hd_bit_syndrome_orbit_gate import coefficient_tables, labelled_classes
from hd_bit_syndrome_outer import CONSTRUCTION
from hd_bit_syndrome_regions_directed import regional_bank
from hd_post_refresh_quiet_rows import rowset_gaps

DIRECTORY = Path(__file__).resolve().parent
MARK, GROUPS, ARITY, COLUMNS, LIMIT = F(993, 1000), 18, 8, 16, 734003
YOUNG = (F(1, 16), F(15, 16))
FROZEN_PINS = {
    'hd_bit_syndrome_atoms_directed.py': '8dad6f3a541e91cd6b14b925745b514f0539eb07c7391a7eee4ff098c6b2ceb6',
    'hd_bit_syndrome_regions_directed.py': '431442dbc6c3772a9f818ff624720ecaab0e6208c67601f66cf23ec0613e129b',
    'hd_bit_syndrome_outer.py': 'a2c2d0fee734fad2c4b74da7850589f8e84180798dcf5f79beab30e06939292d',
    'hd_bit_syndrome_kernel.py': 'a27e5557aceaf04273c3e67ce0a01a5d5746c6d01a2ecda8881205ba23541c77',
    'hd_bit_syndrome_gate.py': '3f8fa49998b1aacc00a964e5c6b83611d62510424afeca98875979aea8063505',
    'hd_bit_syndrome_orbit_gate.py': '1aad0e0c7b88397f18f5c81806ba04b4a2aab64015396e740af4155ea7636023',
    'hd_post_refresh_quiet_rows.py': '75e38819d225e9771ade1ed90509b7d4582a2637b24b6522f5b62b0200eb49e7',
}


def checked_sources(reviewed):
    if set(reviewed) != {'hd_positive_upward.py', 'hd_bit_syndrome_columns_directed.py'}:
        raise ValueError('Both newly reviewed arithmetic engines must be pinned explicitly')
    pins = dict(FROZEN_PINS, **reviewed)
    for name, expected in pins.items():
        if len(expected) != 64 or sha256((DIRECTORY/name).read_bytes()).hexdigest() != expected:
            raise ValueError('Reviewed source changed: '+name)
    return pins


def scaled_fraction(value, exponent):
    """Exact interpretation of an already upward binary64 mantissa."""
    value = F(float(value))
    return value*(1 << exponent) if exponent >= 0 else value/(1 << -exponent)


def signed_float_upper(value):
    endpoint = value if value.is_exact() else value.upper()
    q = endpoint.fmpq()
    exact = F(int(q.p), int(q.q))
    answer = float(exact)
    if F(answer) < exact:
        answer = nextafter(answer, float('inf'))
    if F(answer) < exact:
        raise ArithmeticError('Signed binary64 upper export failed')
    return answer


def aggregate(words, *, groups=GROUPS, precision=256):
    """Upper (sum_A J_A^(1/r))^r, preserving every labelled row."""
    if words.mantissas.ndim != 1 or not len(words.mantissas) or groups < 1:
        raise ValueError('A nonempty labelled-word vector and positive group count required')
    with ctx.workprec(precision):
        total = arb(0)
        for mantissa in words.mantissas:
            value = scaled_fraction(mantissa, words.exponent)
            total += power_upper(value, F(1, groups))
        return (total**groups).upper()


def covered_bounds(aggregate_upper, source_counts, *, precision=256):
    """Apply RS exponent1152, projective exponent1 and Chernoff once."""
    with ctx.workprec(precision):
        shared = rational(upper_fraction(aggregate_upper)) * rational(MARK)**(-LIMIT) / arb(255)**1153
        result = {}
        for key in ('F0', 'F1', 'F_le1'):
            bound = (shared*source_counts[key]).upper()
            logarithm = bound.log()/arb(2).log()
            result[key] = dict(bound_upper_dyadic=str(upper_fraction(bound)),
                              log2_upper=signed_float_upper(logarithm),
                              below_2_minus_40=bool(bound < arb(2)**-40))
    return result


def ordered_products(matrices, parents, labels):
    if (len(parents) != len(labels) or not parents or parents[0] != -1 or labels[0] != -1):
        raise ValueError('A rooted chronological gap trie is required')
    result = [up.identity()]
    for index in range(1, len(parents)):
        parent, label = parents[index], labels[index]
        if not 0 <= parent < index or not 0 <= label < len(matrices.mantissas):
            raise ValueError('Invalid gap prefix dependency')
        result.append(up.matmul(result[parent], matrices[label]))
    return up.concatenate([value[None, ...] for value in result])


def run(precision, *, reviewed, cache_path, max_seconds=240, progress=None):
    if precision not in (256, 384) or not 0 < max_seconds <= 600:
        raise ValueError('Precision256/384 and a bounded positive run budget required')
    pins = checked_sources(reviewed)
    environment = up.ieee_checks()
    cache_path = Path(cache_path)
    if cache_path.exists():
        raise ValueError('Do not overwrite a frozen directed cache')
    started = perf_counter()
    stages = []

    def check_time():
        if perf_counter()-started >= max_seconds:
            raise TimeoutError('Directed gate exhausted its explicit budget')

    def stage(name):
        stages.append(dict(stage=name, elapsed_seconds=perf_counter()-started))
        if progress:
            progress(stages[-1])

    bank = regional_bank(precision, check_time=check_time)
    stage('fresh_original_one_state_regions')
    with ctx.workprec(precision):
        curves = np.empty((4, GROUPS+1), dtype=object)
        for entry in range(4):
            for occupancy in range(GROUPS+1):
                # Compare exact dyadic upper endpoints, never ball midpoints.
                caps = [upper_fraction(bank['regions'][occupancy, f][entry])
                        for f in range(min(1, occupancy//2)+1)]
                curves[entry, occupancy] = rational(max(caps)).upper()
    measure = build_atoms(curves, precision=precision)
    stage('all_76_moments_exactly_dominated')
    rows = row_geometry()[0]
    quiet_rows, quiet_codes = rowset_gaps(16, ARITY)
    if rows != quiet_rows or len(rows) != 12870:
        raise AssertionError('The complete fixed labelled row-set order changed')
    classes = labelled_classes()
    covers, cover_values, parents, labels, nonquiet_codes = cover_gap_geometry()
    caps_q, caps_n = coefficient_tables()
    primitive = primitive_weights(measure, cover_values, caps_q, caps_n)
    stage('directed_atom_cover_and_terminal_weights')
    imported = {key: up.from_export(export_scaled(primitive[key])) for key in
                ('edge_masses', 'empty_matrix', 'cover_matrices', 'terminal_quiet', 'terminal_nonquiet')}
    empty = imported['empty_matrix']
    quiet_products = up.concatenate([up.matrix_power(empty, j)[None, ...] for j in range(9)])
    nonquiet_products = ordered_products(imported['cover_matrices'], parents, labels)
    quiet, cost_q = counted_columns(imported['terminal_quiet'], classes, imported['edge_masses'],
        ARITY, quiet_products, quiet_codes, batch_size=128, check_time=check_time)
    stage('quiet_ordered_column_sum')
    nonquiet, cost_n = counted_columns(imported['terminal_nonquiet'], classes, imported['edge_masses'],
        8, nonquiet_products, nonquiet_codes, batch_size=128, check_time=check_time)
    stage('nonquiet_ordered_column_sum')
    with ctx.workprec(precision):
        young = up.from_export(export_scaled(np.array([
            rational(weight)**(1-GROUPS) for weight in YOUNG], dtype=object)))
    combined = up.add(up.multiply(quiet, young[0]), up.multiply(nonquiet, young[1]))
    word_matrices = up.matrix_power(combined, COLUMNS)
    words = up.add(word_matrices[:, 0, 0], word_matrices[:, 0, 1])
    check_time()
    result_aggregate = aggregate(words, precision=precision)
    bounds = covered_bounds(result_aggregate, bank['source_counts'], precision=precision)
    stage('whole_word_group_sum_and_prefix')
    with cache_path.open('xb') as handle:
        np.savez_compressed(handle, quiet=quiet.mantissas, quiet_scale=quiet.exponent,
            nonquiet=nonquiet.mantissas, nonquiet_scale=nonquiet.exponent,
            combined=combined.mantissas, combined_scale=combined.exponent,
            words=words.mantissas, words_scale=words.exponent,
            rowsets=np.asarray(rows, dtype=np.int8), mark=str(MARK), precision=precision)
    return dict(construction=CONSTRUCTION,
        inner='original one-state h4; packet64, state16, four RS[16,4] maps; D2 H16 D1 update',
        route='alternating even/odd adjacent pair matchings; independent pair permutations and fair swaps',
        scope='r18; homogeneous arity8; Fmax<=1; all labelled row sets and boundary paths',
        complete_requested_class=True, whole_code_certificate=False, directed=True,
        message_length=1048576, codeword_length=2097152, record_bytes=16,
        distance='35/100', forbidden_weight=LIMIT, mark=str(MARK), activity_marker='1',
        groups=GROUPS, arity=ARITY, columns=COLUMNS, rowset_count=len(rows),
        source_rowsets_fixed_across_columns=True, source_count_once=True,
        rs_normalizer_exponent=1152, projective_normalizations=1,
        young_weights=list(map(str, YOUNG)), young_before_word_power=True,
        information_cover=True, per_orbit_caps=True, orbits_used_for_counts_only=True,
        all_boundary_state_paths=True, initial_state=0, terminal_vector=[1, 1],
        original_one_state_inner=True, post_refresh_bank_used=False, saved_float_bounds_used=False,
        regional_curves_upper_dyadic=[[str(upper_fraction(x)) for x in row] for row in curves],
        region_metadata=bank['metadata'], source_counts={k: str(bank['source_counts'][k])
            for k in ('F0', 'F1', 'F_le1', 'F_ge2', 'all_sources')},
        exact_atom_measure=measure['exact_dyadic_inputs'], all_76_moments_verified=True,
        aggregate_upper_dyadic=str(upper_fraction(result_aggregate)), covered_bounds=bounds,
        precision_bits=precision, arithmetic_environment=environment,
        runtime_versions=dict(python=platform.python_version(), numpy=np.__version__),
        costs=dict(quiet=cost_q, nonquiet=cost_n), completed_stages=stages,
        elapsed_seconds=perf_counter()-started, source_sha256=sha256(Path(__file__).read_bytes()).hexdigest(),
        reviewed_sources=pins, cache=str(cache_path), cache_sha256=sha256(cache_path.read_bytes()).hexdigest())


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--precision', type=int, choices=(256, 384), required=True)
    parser.add_argument('--upward-sha256', required=True)
    parser.add_argument('--columns-sha256', required=True)
    parser.add_argument('--cache', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists() or args.cache.exists():
        parser.error('Do not overwrite a frozen receipt or cache')
    answer = run(args.precision, reviewed={'hd_positive_upward.py': args.upward_sha256,
        'hd_bit_syndrome_columns_directed.py': args.columns_sha256}, cache_path=args.cache,
        progress=lambda row: print(json.dumps(row), flush=True))
    with args.output.open('x', encoding='utf-8') as handle:
        json.dump(answer, handle, indent=2, allow_nan=False)
        handle.write('\n')
    print(json.dumps(dict(covered_bounds=answer['covered_bounds'], elapsed_seconds=answer['elapsed_seconds']),
                     allow_nan=False), flush=True)
