"""Isolated quadratic-prefix extensions with 19 through 22 state coordinates.

Preserve the selected sixteen expansion rows, then append independent
quadratic truth tables in the lexicographic order defined by packet_inner_s20.
For the requested state dimension s, the resulting map A has shape 64-by-s;
C=A^T is its feedback map. Dimension 22 is the full RM(2,6) space because
there are 1+6+15 independent constant, linear, and quadratic monomials.

A physical step emits x+A*a and updates a to M*a+C*x. Each M is independent
uniform GL(s,2). State starts at zero, continues across steps and regions,
and is not flushed. Two physical steps form the existing 128-bit macro.

construction() performs only rank and structural checks. Only prepare()
materializes all 2^s states and regenerates spectra and exact kernel profiles.
The largest candidate has over four million states; preparation is research
setup work, not an encoder implementation. No distance certificate is claimed.
Construction pins cover all reused construction and metadata-check helpers;
a whole-code replay must also pin its complete numerical proof engine.
"""
from __future__ import annotations

import hashlib
import json
from pathlib import Path

import packet_inner_s20 as extension
import packet_inner_small_extension as metadata


kernel_t64, kernel_maps, s16_maps = extension.kernel_t64, extension.kernel_maps, extension.s16_maps
BASE_MAP, BASE_SHA256 = extension.BASE_MAP, extension.BASE_SHA256
SUPPORTED_BITS, STEP_BITS = (19, 20, 21, 22), 64
SCHEMA = 'selected-t64-s16-quadratic-prefix-19-through-22-1'
aq, up = extension.aq, extension.up


def _bits(bits):
    if type(bits) is not int or bits not in SUPPORTED_BITS:
        raise ValueError('integer state dimension from 19 through 22 required')
    return bits


def _construction_pins():
    paths = (BASE_MAP, Path(__file__).resolve(), Path(extension.__file__).resolve(),
             Path(s16_maps.__file__).resolve(), Path(metadata.__file__).resolve())
    return {str(path): hashlib.sha256(path.read_bytes()).hexdigest() for path in paths}


def _check_pins(record):
    if record.get('source_sha256') != _construction_pins():
        raise RuntimeError('quadratic construction sources changed or required pins are missing')


def construction(bits):
    """Return explicit rows, columns, and a pinned record without a state census."""
    bits = _bits(bits)
    pins = _construction_pins()
    selected_rows = extension.construction()[0][:16]
    rows, additions = extension.greedy_extension(selected_rows, target_bits=bits)
    columns = extension.columns_from_rows(rows, STEP_BITS)
    packet_ranks = [s16_maps.binary_rank(columns[i:i+4]) for i in range(0, STEP_BITS, 4)]
    if (rows[:16] != selected_rows or s16_maps.binary_rank(rows) != bits or
            s16_maps.binary_rank(columns) != bits or len(set(columns)) != STEP_BITS or
            not all(columns) or packet_ranks != [4]*16 or
            not all(extension._quadratic(row) for row in rows) or
            any((a & b).bit_count() & 1 for a in rows for b in rows)):
        raise ArithmeticError('quadratic prefix failed rank, transpose, or packet checks')
    module, helper = Path(__file__).resolve(), Path(extension.__file__).resolve()
    record = dict(schema=SCHEMA, t=STEP_BITS, s=bits, base_state_bits=16,
        construction='preserve selected16 rows, then append rank-increasing quadratic truth tables '
                     f'in lexicographic (i,j) order, i<j<6, until rank {bits}',
        coordinate_order='coordinate x=0..63, variable i=(x>>i)&1; state coordinates follow row order',
        expansion_rows_hex=list(map(hex, rows)), feedback_columns=list(columns),
        appended_monomials=[list(pair) for pair in additions], appended_rows_hex=list(map(hex, rows[16:])),
        base_expansion_rows_preserved=True, expansion_rank=bits, feedback_rank=bits,
        feedback_definition=f'C=A^T for the explicit 64-by-{bits} expansion A',
        feedback_times_expansion_zero=True, packet_ranks=packet_ranks,
        nonzero_distinct_columns=True, quadratic_expansion=True, full_RM_2_6=bits == 22,
        map_sha256=extension._identity(rows, columns),
        source=dict(path=str(BASE_MAP), sha256=pins[str(BASE_MAP)]),
        source_role=f'selected16 base only; explicit rows define the derived {bits}-coordinate map',
        constructor_source=dict(path=str(module), sha256=pins[str(module)]),
        construction_helper_source=dict(path=str(helper), sha256=pins[str(helper)]), source_sha256=pins,
        distribution='uniform_gl', sampling=f'independent uniform GL{bits} for every physical t64 step',
        return_denominator=(1 << bits)-1, zero_initial_state=True, final_flush=False,
        state_continuity='retained_between_every_physical_step_and_region',
        macro=dict(physical_steps=2, physical_step_bits=64, macro_step_bits=128,
                   physical_windows=16, macro_windows=32, state_continuity='retained_between_halves',
                   initial_state='zero', flush=False),
        full_state_census=False, whole_code_certificate=False,
        scope=f'Fixed t64/s{bits} map with ideal independent uniform GL{bits} updates; '
              'no whole-code certificate or seeded guarantee')
    if pins[str(BASE_MAP)] != BASE_SHA256:
        raise ArithmeticError('pinned selected16 base source changed')
    _check_pins(record)
    return rows, columns, record


def prepare(bits, birth_density='capped'):
    """Regenerate all 2^s states and return the physical-step macro plus record."""
    bits = _bits(bits)
    if birth_density not in ('classes', 'capped'):
        raise ValueError('birth density must be classes or capped')
    rows, columns, record = construction(bits)
    images, spectrum, dual = extension.enumerate_census(rows, STEP_BITS)
    minimum = min(weight for weight in spectrum if weight)
    if minimum != 16:
        raise ArithmeticError('quadratic-prefix expansion distance must be 16')
    physical = kernel_maps.prepare_maps(images, columns, bits=bits,
        distribution='uniform_gl', birth_density=birth_density)
    data = kernel_t64.wrap(physical)
    record.update(full_state_census=True, state_count=1 << bits,
        expansion_spectrum={str(weight): count for weight, count in spectrum.items()},
        feedback_transpose_spectrum={str(weight): count for weight, count in spectrum.items()},
        feedback_kernel_spectrum=list(map(str, dual)), minimum_expansion_weight=minimum,
        minimum_feedback_kernel_weight=next(w for w in range(1, STEP_BITS+1) if dual[w]),
        birth_density=birth_density, exact_profiles_regenerated=True)
    authenticate(data, record)
    return data, record


def authenticate(data, record):
    """Check actual prepared basis, columns, census scope, and current source pins."""
    if not isinstance(record, dict):
        raise ValueError('explicit quadratic-extension map record required')
    bits = _bits(record.get('s'))
    rows, columns, expected = construction(bits)
    keys = expected.keys() - {'full_state_census'}
    if (json.dumps({k: record.get(k) for k in keys}, sort_keys=True) !=
            json.dumps({k: expected[k] for k in keys}, sort_keys=True) or
            record.get('full_state_census') is not True or
            type(record.get('state_count')) is not int or record['state_count'] != 1 << bits or
            record.get('exact_profiles_regenerated') is not True or
            record.get('birth_density') not in ('classes', 'capped') or
            record.get('minimum_expansion_weight') != 16):
        raise ValueError('current quadratic prefix and complete fresh census required')
    metadata._check_census_record(record, bits)
    physical = kernel_t64.authenticate(data)
    if (data['bits'] != bits or data['physical_step_bits'] != STEP_BITS or
            data['macro_step_bits'] != 128 or data['macro_windows'] != 32 or
            data['birth_density'] != record['birth_density'] or data['map_sha256'] != expected['map_sha256'] or
            physical['bits'] != bits or tuple(map(int, physical['columns'])) != columns or
            tuple(int(physical['map_images'][1 << j]) for j in range(bits)) != rows):
        raise ValueError('prepared physical maps do not match the quadratic-prefix record')
    _check_pins(record)
    return physical
