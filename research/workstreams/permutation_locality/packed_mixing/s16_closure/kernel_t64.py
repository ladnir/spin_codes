"""Two physical t64 steps inside one explicit t128 proof macro-step.

Each physical step emits y=x+Aq, then updates q'=Mq+Cx.  Its M is fresh
uniform GL(16,2), independent of every input and all previous updates.
The state starts at zero and is retained between halves; there is no flush.

For j active packets among the macro-step's 32 slots, the first half has
a active slots with probability C(16,a)C(16,j-a)/C(32,j).  Conditional on
this split, both subsets and their nonzero labels are independent uniform.
Thus row-vector transfer matrices compose chronologically as L_a L_(j-a).
The IID macro operator is M64(p)^2.  This wrapper declares both geometries;
it neither fabricates 128-column maps nor imports a previous certificate.
"""
from collections import Counter
from fractions import Fraction as Q
import hashlib
import json
from math import comb
from pathlib import Path

from flint import arb, arb_mat
import kernel_maps
import kernel_birth_density
import sparse_kernel
import s16_maps

aq, up = kernel_maps.aq, kernel_maps.up
HERE = Path(__file__).resolve().parent
SELECTED_MAP = HERE.parents[2]/'rate_quarter_bch/inner_calibration/maps/t64_s16_selected.json'
SCHEMA = 'uniform-gl-two-physical-steps-macro-1'


def _source_maps(path):
    """Read only the declared maps; regenerate every algebraic claim."""
    path = Path(path).resolve(); raw = path.read_bytes(); source = json.loads(raw)
    if source.get('step_bits') != 64 or source.get('state_bits') != 16:
        raise ValueError('declared t64/s16 map required')
    rows = [int(word, 16) for word in source['generator_rows_hex']]
    if (len(rows) != 16 or any(not 0 <= row < 1 << 64 for row in rows)
            or s16_maps.binary_rank(rows) != 16):
        raise ValueError('sixteen independent 64-bit expansion rows required')
    columns = [sum(((row >> p)&1) << j for j, row in enumerate(rows)) for p in range(64)]
    if columns != source.get('columns'):
        raise ArithmeticError('declared feedback is not the transpose of the expansion')
    if (s16_maps.binary_rank(columns) != 16 or len(set(columns)) != 64 or not all(columns)
            or any((a & b).bit_count() & 1 for a in rows for b in rows)):
        raise ArithmeticError('selected map rank/distinctness/CA=0 check failed')
    images = s16_maps.images_from_rows(rows)
    spectrum = dict(sorted(Counter(image.bit_count() for image in images).items()))
    kernel = s16_maps.dual_spectrum(spectrum, 64, 16)
    packet_ranks = [s16_maps.binary_rank(columns[i:i+4]) for i in range(0, 64, 4)]
    if len(set(images)) != 1 << 16 or packet_ranks != [4]*16:
        raise ArithmeticError('complete injective expansion and full-rank packets required')
    record = dict(schema='s16-selected-t64-fixed-maps-1', t=64, s=16,
        expansion_rows_hex=list(map(hex, rows)), feedback_columns=columns,
        feedback_definition='C=A^T for the explicitly declared 64-by-16 expansion A',
        expansion_rank=16, feedback_rank=16, feedback_times_expansion_zero=True,
        expansion_spectrum={str(w): n for w, n in spectrum.items()},
        feedback_transpose_spectrum={str(w): n for w, n in spectrum.items()},
        feedback_kernel_spectrum=list(map(str, kernel)), packet_ranks=packet_ranks,
        minimum_expansion_weight=min(w for w in spectrum if w),
        minimum_feedback_kernel_weight=next(w for w in range(1, 65) if kernel[w]),
        source=dict(path=str(path), sha256=hashlib.sha256(raw).hexdigest()),
        whole_code_certificate=False)
    if path.read_bytes() != raw:
        raise ArithmeticError('map source changed during fresh enumeration')
    return images, columns, record


def wrap(physical_data):
    """Make the explicit two-step wrapper; smaller geometries support tests."""
    kernel_maps.authenticate(physical_data)
    W = physical_data['windows']
    if (physical_data['distribution'] != 'uniform_gl' or not 1 <= W <= 16
            or physical_data.get('profile_partition', False)):
        raise ValueError('uniform-GL physical kernel with at most sixteen packets required')
    wrapper = dict(schema=SCHEMA, distribution='two_independent_uniform_gl_physical_steps',
        physical_data=physical_data, bits=physical_data['bits'], windows=2*W,
        physical_windows=W, macro_windows=2*W,
        physical_step_bits=4*W, macro_step_bits=8*W, physical_steps=2,
        state_continuity='retained_between_halves', initial_state='zero', flush=False,
        physical_updates_independent=True, birth_density=physical_data.get('birth_density', 'classes'),
        map_sha256=physical_data['map_sha256'])
    authenticate(wrapper)
    return wrapper


def authenticate(wrapper):
    """Reject malformed geometry and changed physical maps before evaluation."""
    if (not isinstance(wrapper, dict) or wrapper.get('schema') != SCHEMA
            or wrapper.get('distribution') != 'two_independent_uniform_gl_physical_steps'
            or wrapper.get('physical_steps') != 2
            or wrapper.get('state_continuity') != 'retained_between_halves'
            or wrapper.get('initial_state') != 'zero' or wrapper.get('flush') is not False
            or wrapper.get('physical_updates_independent') is not True
            or any(key in wrapper for key in ('columns', 'map_images', 'updates'))):
        raise ValueError('explicit no-reset two-step macro scope required')
    physical = wrapper.get('physical_data')
    if not isinstance(physical, dict):
        raise ValueError('actual physical-map data required')
    kernel_maps.authenticate(physical)
    W = physical['windows']
    if (physical['distribution'] != 'uniform_gl' or not 1 <= W <= 16
            or physical.get('profile_partition', False)
            or wrapper.get('bits') != physical['bits']
            or wrapper.get('physical_windows') != W or wrapper.get('macro_windows') != 2*W
            or wrapper.get('windows') != 2*W or wrapper.get('physical_step_bits') != 4*W
            or wrapper.get('macro_step_bits') != 8*W
            or wrapper.get('birth_density') != physical.get('birth_density', 'classes')
            or wrapper.get('map_sha256') != physical['map_sha256']):
        raise ArithmeticError('macro wrapper no longer matches its physical map or geometry')
    return physical


def prepare(birth_density='capped', *, source=SELECTED_MAP):
    """Fresh selected-map preparation; source spectra/endpoints are ignored."""
    images, columns, record = _source_maps(source)
    physical = kernel_maps.prepare_maps(images, columns, bits=16,
        distribution='uniform_gl', birth_density=birth_density)
    wrapper = wrap(physical)
    record.update(map_sha256=physical['map_sha256'], distribution='uniform_gl',
        sampling='independent uniform GL16 for every physical t64 step',
        macro=dict(physical_steps=2, physical_step_bits=64, macro_step_bits=128,
                   physical_windows=16, macro_windows=32,
                   state_continuity='retained_between_halves', initial_state='zero', flush=False))
    return wrapper, record


def _upper_matrix(value):
    return arb_mat([[up(value[i, j]) for j in range(value.ncols())] for i in range(value.nrows())])


def convolve(left, right=None):
    """Exact hypergeometric split, with left physical step before right."""
    right = left if right is None else right
    if len(left) < 2 or len(left) != len(right):
        raise ValueError('matching complete occupancy families required')
    W = len(left)-1; n = left[0].nrows()
    if (not 1 <= W <= 16 or n < 1
            or any(m.nrows() != n or m.ncols() != n for m in [*left, *right])
            or any(not m[i, j].is_finite() or not m[i, j] >= 0
                   for m in [*left, *right] for i in range(n) for j in range(n))):
        raise ValueError('finite nonnegative matching square physical operators required')
    result = []
    for total in range(2*W+1):
        splits = [(a, Q(comb(W, a)*comb(W, total-a), comb(2*W, total)))
                  for a in range(max(0, total-W), min(W, total)+1)]
        if sum((weight for _, weight in splits), Q(0)) != 1:
            raise ArithmeticError('hypergeometric split does not normalize')
        value = arb_mat(n, n)
        for a, weight in splits:
            value += aq(weight)*(left[a]*right[total-a])
        result.append(_upper_matrix(value))
    return result


def floating(wrapper, probabilities, tilt):
    """Proposal-only IID operator for both physical steps."""
    physical = authenticate(wrapper)
    matrix = kernel_maps.floating(physical, probabilities, tilt)
    return matrix @ matrix


def outward(wrapper, probabilities, tilt):
    """Outward IID operator; global macro-step counts remain unchanged."""
    physical = authenticate(wrapper)
    matrix = kernel_maps.outward(physical, probabilities, tilt)
    return _upper_matrix(matrix*matrix)


def local_operators(wrapper, tilt, activity=Q(1, 2)):
    """Return one macro operator for every occupancy from zero through 32."""
    physical = authenticate(wrapper); tilt, activity = Q(tilt), Q(activity)
    if tilt <= 0 or not 0 <= activity <= 1:
        raise ValueError('positive output tilt and valid row-selection activity required')
    z = (-aq(tilt)).exp()
    local = sparse_kernel.outward_at_z(physical, z)
    if physical.get('birth_density', 'classes') == 'capped':
        local = kernel_birth_density.refine_local(physical, local, z, activity)
    return convolve(local)
