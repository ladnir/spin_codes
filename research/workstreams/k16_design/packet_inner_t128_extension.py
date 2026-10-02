"""One physical t128 step with a twenty- or twenty-two-coordinate state.

The base expansion A is the retained t128/s20 map. For s22, preserve its
twenty rows and append the first independent quadratic truth tables in
lexicographic variable-pair order. Set C=A^T. A step emits x+A*a, then
updates a to M*a+C*x, where each M is independent uniform GL(s,2).
State starts at zero, persists across steps and regions, and is not flushed.

Unlike the t64 adapters, this module processes all 32 packet positions in
ONE physical step. Its occupancy operators must not be convolved again.
construction() checks the explicit maps without enumerating their images.
prepare() regenerates the complete census using Python integers, including
the high 64 bits; saved spectra in the base JSON are not accepted as proof.
No whole-code certificate or implementation timing is claimed.
"""
from __future__ import annotations

from collections import Counter
from fractions import Fraction as Q
import hashlib
from itertools import combinations
import json
from pathlib import Path

import packet_q1 as q1


kernel_maps = q1.kernel_t64.kernel_maps
s16_maps = q1.kernel_t64.s16_maps
sparse_kernel = q1.kernel_t64.sparse_kernel
kernel_birth_density = q1.kernel_t64.kernel_birth_density
BASE_MAP = (Path(__file__).resolve().parents[1] /
            'rate_quarter_bch/inner_calibration/maps/t128_s20_nested.json').resolve()
BASE_SHA256 = 'eeca4787b4447607a47d731879f4e3ced1ebfd6ab2266f4b6ecb1244fc108f29'
SUPPORTED_BITS, STEP_BITS, WINDOWS = (20, 22), 128, 32
SCHEMA = 'nested-t128-s20-quadratic-extension-1'
aq, up = kernel_maps.aq, kernel_maps.up


def _bits(bits):
    if type(bits) is not int or bits not in SUPPORTED_BITS:
        raise ValueError('integer state dimension 20 or 22 required')
    return bits


def columns_from_rows(rows, width=STEP_BITS):
    rows = tuple(rows)
    if (type(width) is not int or not 1 <= width <= STEP_BITS or not rows or
            any(type(row) is not int or not 0 <= row < 1 << width for row in rows)):
        raise ValueError('bounded binary rows and an output width from 1 through 128 required')
    return tuple(sum(((row >> x) & 1) << i for i, row in enumerate(rows))
                 for x in range(width))


def _quadratic(row):
    coefficients = [(row >> x) & 1 for x in range(STEP_BITS)]
    for i in range(7):
        for x in range(STEP_BITS):
            if x & (1 << i):
                coefficients[x] ^= coefficients[x ^ (1 << i)]
    return not any(value and mask.bit_count() > 2 for mask, value in enumerate(coefficients))


def _identity(rows, columns):
    payload = dict(bits=len(rows), width=len(columns), expansion_rows=list(map(hex, rows)),
                   feedback_columns=list(columns))
    return hashlib.sha256(json.dumps(payload, sort_keys=True, separators=(',', ':')).encode()).hexdigest()


def _construction_pins():
    paths = (BASE_MAP, Path(__file__).resolve(), Path(s16_maps.__file__).resolve())
    return {str(path): hashlib.sha256(path.read_bytes()).hexdigest() for path in paths}


def _check_pins(record):
    if record.get('source_sha256') != _construction_pins():
        raise RuntimeError('t128 construction sources changed or required pins are missing')


def construction(bits):
    """Return rows, transpose columns, and the pinned one-step definition."""
    bits = _bits(bits)
    pins = _construction_pins()
    if pins[str(BASE_MAP)] != BASE_SHA256:
        raise ArithmeticError('pinned t128/s20 source changed')
    source = json.loads(BASE_MAP.read_bytes())
    base_rows = tuple(int(word, 16) for word in source['generator_rows_hex'])
    if (source.get('step_bits') != STEP_BITS or source.get('state_bits') != 20 or
            len(base_rows) != 20 or s16_maps.binary_rank(base_rows) != 20 or
            list(columns_from_rows(base_rows)) != source.get('columns') or
            not all(_quadratic(row) for row in base_rows)):
        raise ArithmeticError('pinned source must specify the quadratic t128/s20 map')
    rows, additions = list(base_rows), []
    for pair in combinations(range(7), 2):
        if len(rows) == bits:
            break
        i, j = pair
        row = sum(1 << x for x in range(STEP_BITS) if (x >> i) & (x >> j) & 1)
        if s16_maps.binary_rank([*rows, row]) == len(rows)+1:
            rows.append(row)
            additions.append(list(pair))
    rows = tuple(rows)
    columns = columns_from_rows(rows)
    packet_ranks = [s16_maps.binary_rank(columns[x:x+4]) for x in range(0, STEP_BITS, 4)]
    if (len(rows) != bits or s16_maps.binary_rank(rows) != bits or
            s16_maps.binary_rank(columns) != bits or len(set(columns)) != STEP_BITS or
            not all(columns) or packet_ranks != [4]*WINDOWS or
            any((a & b).bit_count() & 1 for a in rows for b in rows)):
        raise ArithmeticError('t128 extension failed rank, packet rank, or CA=0 checks')
    record = dict(schema=SCHEMA, t=STEP_BITS, s=bits, base_state_bits=20,
        construction='preserve nested20 rows, then append independent quadratic truth tables '
                     'in lexicographic (i,j) order, i<j<7',
        coordinate_order='coordinate x=0..127, variable i=(x>>i)&1; state coordinates follow row order',
        expansion_rows_hex=list(map(hex, rows)), feedback_columns=list(columns),
        appended_monomials=additions, base_expansion_rows_preserved=True,
        feedback_definition=f'C=A^T for the explicit 128-by-{bits} expansion A',
        expansion_rank=bits, feedback_rank=bits, feedback_times_expansion_zero=True,
        nonzero_distinct_columns=True, packet_ranks=packet_ranks, quadratic_expansion=True,
        map_sha256=_identity(rows, columns), source=dict(path=str(BASE_MAP), sha256=BASE_SHA256),
        source_sha256=pins, distribution='uniform_gl',
        sampling=f'independent uniform GL{bits} for every physical t128 step',
        return_denominator=(1 << bits)-1, zero_initial_state=True, final_flush=False,
        state_continuity='retained_between_every_physical_step_and_region',
        geometry=dict(physical_steps=1, physical_step_bits=128, physical_windows=32,
                      proof_step_bits=128, proof_windows=32),
        full_state_census=False, whole_code_certificate=False,
        scope='Explicit ideal uniform-GL inner; no whole-code certificate or seeded guarantee')
    _check_pins(record)
    return rows, columns, record


def enumerate_census(rows, width=STEP_BITS):
    """Enumerate exact Python-int images; do not truncate to a machine word."""
    rows = tuple(rows)
    columns_from_rows(rows, width)
    if not 1 <= len(rows) <= 24 or s16_maps.binary_rank(rows) != len(rows):
        raise ValueError('one through twenty-four independent expansion rows required')
    images = tuple(s16_maps.images_from_rows(rows))
    spectrum = dict(sorted(Counter(image.bit_count() for image in images).items()))
    if len(images) != 1 << len(rows) or spectrum.get(0) != 1:
        raise ArithmeticError('full injective expansion census failed')
    dual = tuple(s16_maps.dual_spectrum(spectrum, width, len(rows)))
    return images, spectrum, dual


def prepare(bits, birth_density='capped'):
    """Return freshly censused physical data and its one-step map record."""
    bits = _bits(bits)
    if birth_density not in ('classes', 'capped'):
        raise ValueError('birth density must be classes or capped')
    rows, columns, record = construction(bits)
    images, spectrum, dual = enumerate_census(rows)
    minimum = min(weight for weight in spectrum if weight)
    if minimum != 32:
        raise ArithmeticError('this quadratic extension must have expansion distance 32')
    data = kernel_maps.prepare_maps(images, columns, bits=bits,
        distribution='uniform_gl', birth_density=birth_density)
    record.update(full_state_census=True, state_count=1 << bits,
        expansion_spectrum={str(weight): count for weight, count in spectrum.items()},
        feedback_transpose_spectrum={str(weight): count for weight, count in spectrum.items()},
        feedback_kernel_spectrum=list(map(str, dual)), minimum_expansion_weight=minimum,
        minimum_feedback_kernel_weight=next(w for w in range(1, STEP_BITS+1) if dual[w]),
        birth_density=birth_density, exact_profiles_regenerated=True)
    authenticate(data, record)
    return data, record


def _check_census_record(record, bits):
    spectrum, dual = record.get('expansion_spectrum'), record.get('feedback_kernel_spectrum')
    if (not isinstance(spectrum, dict) or not spectrum or
            any(type(w) is not str or not w.isdecimal() or str(int(w)) != w or
                not 0 <= int(w) <= STEP_BITS or type(n) is not int or n < 1
                for w, n in spectrum.items()) or spectrum.get('0') != 1 or
            sum(spectrum.values()) != 1 << bits or
            record.get('feedback_transpose_spectrum') != spectrum or
            min((int(w) for w in spectrum if w != '0'), default=None) != 32 or
            not isinstance(dual, list) or len(dual) != STEP_BITS+1 or
            any(type(n) is not str or not n.isdecimal() or str(int(n)) != n for n in dual)):
        raise ValueError('complete exact t128 expansion and dual spectrum metadata required')
    counts = tuple(map(int, dual))
    if (counts[0] != 1 or sum(counts) != 1 << (STEP_BITS-bits) or
            type(record.get('minimum_feedback_kernel_weight')) is not int or
            record['minimum_feedback_kernel_weight'] != next((w for w in range(1, STEP_BITS+1) if counts[w]), None)):
        raise ValueError('complete dual census and its nonzero minimum required')


def authenticate(data, record=None):
    """Bind prepared data to this physical map; optionally check census metadata."""
    bits = _bits(data.get('bits'))
    rows, columns, expected = construction(bits)
    kernel_maps.authenticate(data)
    if (data.get('distribution') != 'uniform_gl' or data.get('windows') != WINDOWS or
            data.get('birth_density') not in ('classes', 'capped') or
            data.get('profile_partition', False) or
            data.get('map_sha256') != expected['map_sha256'] or
            tuple(map(int, data['columns'])) != columns or
            tuple(int(data['map_images'][1 << i]) for i in range(bits)) != rows):
        raise ValueError('actual single-step t128 maps and uniform-GL preparation required')
    if record is not None:
        keys = expected.keys() - {'full_state_census'}
        if (not isinstance(record, dict) or
                json.dumps({k: record.get(k) for k in keys}, sort_keys=True) !=
                json.dumps({k: expected[k] for k in keys}, sort_keys=True) or
                record.get('full_state_census') is not True or
                type(record.get('state_count')) is not int or record['state_count'] != 1 << bits or
                record.get('exact_profiles_regenerated') is not True or
                record.get('birth_density') != data['birth_density'] or
                record.get('minimum_expansion_weight') != 32):
            raise ValueError('current t128 construction with a complete fresh census required')
        _check_census_record(record, bits)
        _check_pins(record)
    return data


def local_operators(data, tilt, activity=Q(1, 2)):
    """Return physical operators for occupancies 0..32, without convolution."""
    authenticate(data)
    tilt, activity = Q(tilt), Q(activity)
    if tilt <= 0 or not 0 <= activity <= 1:
        raise ValueError('positive output tilt and valid row-selection activity required')
    z = (-aq(tilt)).exp()
    local = sparse_kernel.outward_at_z(data, z)
    if data['birth_density'] == 'capped':
        local = kernel_birth_density.refine_local(data, local, z, activity)
    if len(local) != WINDOWS+1:
        raise ArithmeticError('complete physical occupancy family 0..32 required')
    return local
