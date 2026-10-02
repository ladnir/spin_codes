"""An isolated twenty-state-coordinate extension of the selected t64 map.

The first sixteen expansion vectors are exactly those in
t64_s16_selected.json. Enumerate pairs (i,j), 0<=i<j<6, in lexicographic
order. The vector for a pair has coordinate x equal to x_i*x_j, where
x=0,...,63 is represented by six bits. Append a vector precisely when it
increases binary rank, stopping at twenty vectors. This procedure defines
the map; it does not use the separate t64_s20_nested candidate.

Let A be the resulting 64-by-20 binary expansion and set C=A^T. A physical
step emits y=x+A*s and updates s'=M*s+C*x. Each M is sampled independently
and uniformly from GL(20,2), independently of the input and entering state.
The state starts at zero and persists between all steps and regions.
There is no final flush. Two ordered physical steps form a 128-bit macro.

construction() verifies ranks, packet ranks, and C*A=0 without a full state
census. prepare() additionally regenerates all 2^20 expansion images, their
weight spectrum, the exact MacWilliams dual spectrum, and the kernel's exact
character and birth-class censuses. No saved spectrum or numerical bound
is accepted. Both APIs identify the base file and this constructor by hash.

The default birth_density='classes' avoids the optional capped-density
attachment. Both choices describe the same ideal GL20 law and differ only
in their envelope representation. Preparation is research setup work, not
an encoder hot path: one million states are materialized, and downstream
capped-density inversion can allocate several large state-by-occupancy arrays.
This module asserts no whole-code distance certificate or seeded guarantee.
"""
from __future__ import annotations

from collections import Counter
from fractions import Fraction as Q
import hashlib
from itertools import combinations
import json
from pathlib import Path

import packet_q1 as q1


kernel_t64 = q1.kernel_t64
kernel_maps = kernel_t64.kernel_maps
s16_maps = kernel_t64.s16_maps
BASE_MAP = kernel_t64.SELECTED_MAP.resolve()
BASE_SHA256 = 'dabe234b3b78d00ae25d1fef29be9f8eebb7dc111495d9ad5fdd0cdbd2d9c121'
STATE_BITS = 20
STEP_BITS = 64
SCHEMA = 'selected-t64-s16-quadratic-extension-s20-1'
aq, up = kernel_t64.aq, kernel_t64.up


def quadratic_monomial(i, j, *, variables=6):
    """Truth table of x_i*x_j; coordinate x is bit x of the returned integer."""
    if (type(variables) is not int or not 2 <= variables <= 6
            or type(i) is not int or type(j) is not int or not 0 <= i < j < variables):
        raise ValueError('ordered variable indices 0<=i<j<variables<=6 required')
    return sum(1 << x for x in range(1 << variables) if (x >> i) & (x >> j) & 1)


def columns_from_rows(rows, width):
    """Transpose the explicit expansion vectors into feedback columns."""
    rows = tuple(rows)
    if (type(width) is not int or not 1 <= width <= 64 or not rows
            or any(type(row) is not int or not 0 <= row < 1 << width for row in rows)):
        raise ValueError('bounded binary rows in a positive output width required')
    return tuple(sum(((row >> coordinate) & 1) << state for state, row in enumerate(rows))
                 for coordinate in range(width))


def greedy_extension(base_rows, target_bits=STATE_BITS, *, variables=6):
    """Append independent quadratic truth tables in lexicographic pair order."""
    if type(variables) is not int or not 2 <= variables <= 6:
        raise ValueError('two through six Boolean variables required')
    rows = list(base_rows)
    if (not rows or type(target_bits) is not int or not len(rows) <= target_bits <= 24
            or any(type(row) is not int or not 0 <= row < 1 << (1 << variables) for row in rows)
            or s16_maps.binary_rank(rows) != len(rows)):
        raise ValueError('independent bounded base rows and a valid target dimension required')
    additions = []
    for pair in combinations(range(variables), 2):
        if len(rows) == target_bits:
            break
        row = quadratic_monomial(*pair, variables=variables)
        if s16_maps.binary_rank([*rows, row]) == len(rows) + 1:
            rows.append(row)
            additions.append(pair)
    if len(rows) != target_bits:
        raise ValueError('quadratic candidates cannot reach the requested dimension')
    return tuple(rows), tuple(additions)


def _quadratic(row):
    coefficients = [(row >> x) & 1 for x in range(STEP_BITS)]
    for bit in (1, 2, 4, 8, 16, 32):
        for x in range(STEP_BITS):
            if x & bit:
                coefficients[x] ^= coefficients[x ^ bit]
    return not any(value and mask.bit_count() > 2 for mask, value in enumerate(coefficients))


def _identity(rows, columns):
    # Match kernel_maps._identity without materializing its full image tuple.
    payload = dict(bits=len(rows), width=len(columns), expansion_rows=list(map(hex, rows)),
                   feedback_columns=list(columns))
    return hashlib.sha256(json.dumps(payload, sort_keys=True, separators=(',', ':')).encode()).hexdigest()


def construction():
    """Return (rows, columns, record) after light exact structural checks."""
    raw = BASE_MAP.read_bytes()
    source_hash = hashlib.sha256(raw).hexdigest()
    if source_hash != BASE_SHA256:
        raise ArithmeticError('the selected sixteen-row base source has changed')
    source = json.loads(raw)
    base_rows = tuple(int(word, 16) for word in source['generator_rows_hex'])
    if (source.get('step_bits') != STEP_BITS or source.get('state_bits') != 16
            or len(base_rows) != 16 or s16_maps.binary_rank(base_rows) != 16
            or not all(_quadratic(row) for row in base_rows)
            or list(columns_from_rows(base_rows, STEP_BITS)) != source.get('columns')):
        raise ArithmeticError('the pinned base does not define the selected quadratic t64/s16 map')
    rows, additions = greedy_extension(base_rows)
    columns = columns_from_rows(rows, STEP_BITS)
    rank = s16_maps.binary_rank(rows)
    packet_ranks = [s16_maps.binary_rank(columns[start:start+4]) for start in range(0, STEP_BITS, 4)]
    if (rank != STATE_BITS or rows[:16] != base_rows
            or s16_maps.binary_rank(columns) != STATE_BITS
            or len(set(columns)) != STEP_BITS or not all(columns)
            or any((a & b).bit_count() & 1 for a in rows for b in rows)
            or packet_ranks != [4] * 16):
        raise ArithmeticError('the derived map failed rank, distinctness, orthogonality, or packet checks')
    module = Path(__file__).resolve()
    module_hash = hashlib.sha256(module.read_bytes()).hexdigest()
    pins = {str(BASE_MAP): source_hash, str(module): module_hash}
    record = dict(schema=SCHEMA, t=STEP_BITS, s=STATE_BITS, base_state_bits=16,
        construction='preserve all selected16 rows, then append the first rank-increasing '
                     'quadratic truth tables in lexicographic (i,j) order, i<j<6, until rank20',
        coordinate_order='coordinate x=0..63, variable i=(x>>i)&1; state coordinates follow row order',
        appended_monomials=[list(pair) for pair in additions],
        appended_rows_hex=list(map(hex, rows[16:])), expansion_rows_hex=list(map(hex, rows)),
        base_expansion_rows_preserved=True, feedback_definition='C=A^T for the explicit 64-by-20 expansion A',
        feedback_columns=list(columns), expansion_rank=rank, feedback_rank=rank,
        feedback_times_expansion_zero=True, packet_ranks=packet_ranks,
        nonzero_distinct_columns=True, quadratic_expansion=True,
        map_sha256=_identity(rows, columns),
        source=dict(path=str(BASE_MAP), sha256=source_hash),
        source_role='selected16 base only; the derived20 map is specified by rows and construction',
        constructor_source=dict(path=str(module), sha256=module_hash), source_sha256=pins,
        distribution='uniform_gl', sampling='independent uniform GL20 for every physical t64 step',
        zero_initial_state=True, final_flush=False,
        state_continuity='retained_between_every_physical_step_and_region',
        macro=dict(physical_steps=2, physical_step_bits=64, macro_step_bits=128,
                   physical_windows=16, macro_windows=32, state_continuity='retained_between_halves',
                   initial_state='zero', flush=False),
        full_state_census=False, whole_code_certificate=False,
        scope='Declared fixed t64/s20 inner with an ideal fresh uniform GL20 law; '
              'no whole-code certificate and no seeded implementation guarantee')
    if BASE_MAP.read_bytes() != raw:
        raise RuntimeError('base source changed during map construction')
    return rows, columns, record


def enumerate_census(rows, width):
    """Regenerate all images and exact primal/dual spectra for independent rows."""
    rows = tuple(rows)
    columns_from_rows(rows, width)  # Validate the declared output width.
    if not 1 <= len(rows) <= 24 or s16_maps.binary_rank(rows) != len(rows):
        raise ValueError('one through twenty-four independent expansion rows required')
    images = tuple(s16_maps.images_from_rows(rows))
    spectrum = dict(sorted(Counter(image.bit_count() for image in images).items()))
    if len(images) != 1 << len(rows) or spectrum.get(0) != 1 or sum(spectrum.values()) != len(images):
        raise ArithmeticError('full injective expansion census failed')
    # Binary rank proves injectivity; no second million-element set is needed.
    dual = tuple(s16_maps.dual_spectrum(spectrum, width, len(rows)))
    return images, spectrum, dual


def _check_pins(record):
    if any(hashlib.sha256(Path(path).read_bytes()).hexdigest() != digest
           for path, digest in record['source_sha256'].items()):
        raise RuntimeError('map construction source changed during fresh preparation')


def prepare(birth_density='classes'):
    """Return fresh (two-step wrapper, map record); materialize the full s20 census."""
    if birth_density not in ('classes', 'capped'):
        raise ValueError('birth density must be classes or capped')
    rows, columns, record = construction()
    images, spectrum, dual = enumerate_census(rows, STEP_BITS)
    minimum = min(weight for weight in spectrum if weight)
    if minimum != 16:
        raise ArithmeticError('the declared extension no longer has expansion distance sixteen')
    physical = kernel_maps.prepare_maps(images, columns, bits=STATE_BITS,
        distribution='uniform_gl', birth_density=birth_density)
    wrapper = kernel_t64.wrap(physical)
    if (wrapper['bits'] != STATE_BITS or wrapper['physical_step_bits'] != STEP_BITS
            or wrapper['map_sha256'] != record['map_sha256']):
        raise ArithmeticError('prepared kernel does not match the declared twenty-coordinate map')
    record.update(full_state_census=True, state_count=1 << STATE_BITS,
        expansion_spectrum={str(weight): count for weight, count in spectrum.items()},
        feedback_transpose_spectrum={str(weight): count for weight, count in spectrum.items()},
        feedback_kernel_spectrum=list(map(str, dual)), minimum_expansion_weight=minimum,
        minimum_feedback_kernel_weight=next(weight for weight in range(1, STEP_BITS+1) if dual[weight]),
        birth_density=birth_density, exact_profiles_regenerated=True)
    _check_pins(record)
    return wrapper, record


def authenticate(wrapper):
    """Require this exact twenty-coordinate extension, not merely any t64 map."""
    physical = kernel_t64.authenticate(wrapper)
    _, _, record = construction()
    if (wrapper['bits'] != STATE_BITS or wrapper['physical_step_bits'] != STEP_BITS
            or wrapper['map_sha256'] != record['map_sha256']):
        raise ValueError('this selected16-to20 extension and uniform GL20 law are required')
    return physical


def local_operators(wrapper, tilt, activity=Q(1, 2)):
    authenticate(wrapper)
    return kernel_t64.local_operators(wrapper, tilt, activity)


def outward(wrapper, probabilities, tilt):
    authenticate(wrapper)
    return kernel_t64.outward(wrapper, probabilities, tilt)


def floating(wrapper, probabilities, tilt):
    authenticate(wrapper)
    return kernel_t64.floating(wrapper, probabilities, tilt)
