"""Isolated 17- and 18-coordinate extensions of the selected t64/s16 map.

For s in {17,18}, preserve the selected sixteen expansion vectors. Append
the first s-16 independent quadratic truth tables in the lexicographic
pair order used by packet_inner_s20. This defines a 64-by-s binary map A;
the feedback map is C=A^T. Construction checks rank, C*A=0, distinct nonzero
columns, and full rank for each four-column packet, without a state census.

A physical step emits y=x+A*a and updates a'=M*a+C*x. Each M is fresh
uniform GL(s,2), independent of the input and entering state. The state is
initially zero, persists across steps and regions, and is not flushed.
Two ordered physical steps form the existing 128-bit proof macro.

Only prepare() materializes all 2^s expansion images and regenerates exact
spectra and kernel profiles. No S16/S20 numerical table is reused. Capped
and class birth-density envelopes describe the same ideal update law.
Construction pins identify the base, this adapter, and its construction
helpers; a whole-code replay must additionally pin its full proof engine.
This module claims neither a whole-code distance bound nor a seeded guarantee.
"""
from __future__ import annotations

import hashlib
import json
from pathlib import Path

import packet_inner_s20 as extension


kernel_t64 = extension.kernel_t64
kernel_maps = extension.kernel_maps
s16_maps = extension.s16_maps
BASE_MAP = extension.BASE_MAP
BASE_SHA256 = extension.BASE_SHA256
SUPPORTED_BITS = (17, 18)
STEP_BITS = 64
SCHEMA = 'selected-t64-s16-quadratic-small-extension-1'
aq, up = extension.aq, extension.up


def _bits(bits):
    if type(bits) is not int or bits not in SUPPORTED_BITS:
        raise ValueError('state dimension must be the integer17 or18')
    return bits


def _construction_pins():
    paths = (BASE_MAP, Path(__file__).resolve(), Path(extension.__file__).resolve(),
             Path(s16_maps.__file__).resolve())
    return {str(path): hashlib.sha256(path.read_bytes()).hexdigest() for path in paths}


def _check_pins(record):
    if record.get('source_sha256') != _construction_pins():
        raise RuntimeError('map construction sources changed or required pins are missing')


def construction(bits):
    """Return (rows, columns, record) after structural checks, without enumeration."""
    bits = _bits(bits)
    pins = _construction_pins()
    # This pure constructor performs only bounded row/column algebra. Taking
    # its prefix preserves exactly the declared greedy rank-increasing order.
    full_rows, _, base_record = extension.construction()
    rows = full_rows[:bits]
    columns = extension.columns_from_rows(rows, STEP_BITS)
    packet_ranks = [s16_maps.binary_rank(columns[start:start+4])
                    for start in range(0, STEP_BITS, 4)]
    if (len(rows) != bits or s16_maps.binary_rank(rows) != bits or
            s16_maps.binary_rank(columns) != bits or
            len(set(columns)) != STEP_BITS or not all(columns) or
            any((a & b).bit_count() & 1 for a in rows for b in rows) or
            packet_ranks != [4]*16 or not all(extension._quadratic(row) for row in rows)):
        raise ArithmeticError('small extension failed rank, orthogonality, or packet checks')
    module, helper = Path(__file__).resolve(), Path(extension.__file__).resolve()
    record = dict(schema=SCHEMA, t=STEP_BITS, s=bits, base_state_bits=16,
        construction='preserve all selected16 rows, then append the first rank-increasing '
                     f'quadratic truth tables in lexicographic (i,j) order, i<j<6, until rank{bits}',
        coordinate_order='coordinate x=0..63, variable i=(x>>i)&1; state coordinates follow row order',
        appended_monomials=base_record['appended_monomials'][:bits-16],
        appended_rows_hex=list(map(hex, rows[16:])), expansion_rows_hex=list(map(hex, rows)),
        base_expansion_rows_preserved=True,
        feedback_definition=f'C=A^T for the explicit 64-by-{bits} expansion A',
        feedback_columns=list(columns), expansion_rank=bits, feedback_rank=bits,
        feedback_times_expansion_zero=True, packet_ranks=packet_ranks,
        nonzero_distinct_columns=True, quadratic_expansion=True,
        map_sha256=extension._identity(rows, columns),
        source=dict(path=str(BASE_MAP), sha256=pins[str(BASE_MAP)]),
        source_role=f'selected16 base only; derived{bits} map is specified by rows and construction',
        constructor_source=dict(path=str(module), sha256=pins[str(module)]),
        construction_helper_source=dict(path=str(helper), sha256=pins[str(helper)]),
        source_sha256=pins, distribution='uniform_gl',
        sampling=f'independent uniform GL{bits} for every physical t64 step',
        return_denominator=(1 << bits)-1,
        zero_initial_state=True, final_flush=False,
        state_continuity='retained_between_every_physical_step_and_region',
        macro=dict(physical_steps=2, physical_step_bits=64, macro_step_bits=128,
                   physical_windows=16, macro_windows=32, state_continuity='retained_between_halves',
                   initial_state='zero', flush=False),
        full_state_census=False, whole_code_certificate=False,
        scope=f'Declared fixed t64/s{bits} inner with ideal fresh uniform GL{bits} updates; '
              'no whole-code certificate and no seeded implementation guarantee')
    if pins[str(BASE_MAP)] != BASE_SHA256:
        raise ArithmeticError('the pinned selected16 base source has changed')
    _check_pins(record)
    return rows, columns, record


def prepare(bits, birth_density='capped'):
    """Regenerate the complete 2^s census and return (macro wrapper, map record)."""
    bits = _bits(bits)
    if birth_density not in ('classes', 'capped'):
        raise ValueError('birth density must be classes or capped')
    rows, columns, record = construction(bits)
    images, spectrum, dual = extension.enumerate_census(rows, STEP_BITS)
    minimum = min(weight for weight in spectrum if weight)
    if minimum != 16:
        raise ArithmeticError('the quadratic-prefix extension must have expansion distance16')
    physical = kernel_maps.prepare_maps(images, columns, bits=bits,
        distribution='uniform_gl', birth_density=birth_density)
    data = kernel_t64.wrap(physical)
    record.update(full_state_census=True, state_count=1 << bits,
        expansion_spectrum={str(weight): count for weight, count in spectrum.items()},
        feedback_transpose_spectrum={str(weight): count for weight, count in spectrum.items()},
        feedback_kernel_spectrum=list(map(str, dual)), minimum_expansion_weight=minimum,
        minimum_feedback_kernel_weight=next(weight for weight in range(1, STEP_BITS+1) if dual[weight]),
        birth_density=birth_density, exact_profiles_regenerated=True)
    authenticate(data, record)
    return data, record


def _check_census_record(record, bits):
    """Check completeness metadata without enumerating any expansion image."""
    spectrum = record.get('expansion_spectrum')
    dual_strings = record.get('feedback_kernel_spectrum')
    if (not isinstance(spectrum, dict) or not spectrum or
            any(type(key) is not str or not key.isdecimal() or str(int(key)) != key or
                not 0 <= int(key) <= STEP_BITS or type(count) is not int or count < 1
                for key, count in spectrum.items()) or
            spectrum.get('0') != 1 or sum(spectrum.values()) != 1 << bits or
            record.get('feedback_transpose_spectrum') != spectrum or
            min((int(key) for key in spectrum if key != '0'), default=None) != 16 or
            not isinstance(dual_strings, list) or len(dual_strings) != STEP_BITS+1 or
            any(type(value) is not str or not value.isdecimal() or str(int(value)) != value
                for value in dual_strings)):
        raise ValueError('complete expansion and feedback spectra metadata required')
    dual = tuple(map(int, dual_strings))
    minimum_dual = next((weight for weight in range(1, STEP_BITS+1) if dual[weight]), None)
    if (dual[0] != 1 or sum(dual) != 1 << (STEP_BITS-bits) or
            type(record.get('minimum_feedback_kernel_weight')) is not int or
            record['minimum_feedback_kernel_weight'] != minimum_dual):
        raise ValueError('complete dual spectrum and its nonzero minimum required')


def authenticate(data, record):
    """Bind fresh prepared data to the selected prefix, census scope, and source pins.

    This checks identity and immutable preparation metadata, without rerunning
    a state census. Callers must obtain data from prepare(), not saved tables.
    """
    if not isinstance(record, dict):
        raise ValueError('an explicit small-extension map record is required')
    bits = _bits(record.get('s'))
    rows, columns, expected = construction(bits)
    # Exact JSON equality distinguishes booleans from integer lookalikes in
    # nested scope metadata. The only changed construction flag is the census.
    declared = {key: record.get(key) for key in expected if key != 'full_state_census'}
    wanted = {key: value for key, value in expected.items() if key != 'full_state_census'}
    if (json.dumps(declared, sort_keys=True) != json.dumps(wanted, sort_keys=True) or
            record.get('full_state_census') is not True or
            type(record.get('state_count')) is not int or record['state_count'] != 1 << bits or
            record.get('exact_profiles_regenerated') is not True or
            record.get('birth_density') not in ('classes', 'capped') or
            record.get('minimum_expansion_weight') != 16):
        raise ValueError('current exact selected16 prefix with a fresh full census is required')
    _check_census_record(record, bits)
    physical = kernel_t64.authenticate(data)
    if (data['bits'] != bits or data['physical_step_bits'] != STEP_BITS or
            data['macro_step_bits'] != 128 or data['macro_windows'] != 32 or
            data['birth_density'] != record['birth_density'] or
            data['map_sha256'] != expected['map_sha256'] or
            physical['bits'] != bits or tuple(map(int, physical['columns'])) != columns or
            tuple(int(physical['map_images'][1 << j]) for j in range(bits)) != rows):
        raise ValueError('prepared data does not match the declared t64 small-extension map')
    _check_pins(record)
    return physical
