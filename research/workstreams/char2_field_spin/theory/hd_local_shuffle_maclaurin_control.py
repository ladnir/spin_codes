"""Saved-bank FLOAT repair for the independently output-shuffled H* ensemble.

For a fixed labelled input set A, uniform output permutation gives
Kbar_A(t)=sum_w WE_A[w] e_w(t)/choose(16,w).  At atoms in [0,1],
Maclaurin gives Kbar_A(t)<=WE_A(s)<=C_A(s) Z_A(s), s=mean(t).
The second inequality is the scalar translated-fiber bound; it is NOT a
coefficientwise comparison of weight enumerators.  Z_A is the ACTUAL
zero-outside fiber, reconstructed using symmetric-involution duality.

At each common atom-count tuple, this helper takes the pointwise minimum
with the old output-shuffle bound.  For the nonzero-column branch, both
positive polynomials delete their exact constant1 BEFORE evaluation.
Full and nonzero minima are independent valid bounds.  Each fixed-A class
retains BOTH polynomials, then the group moment, ordered word product,
first-nine-active-column bound, and whole-word group Minkowski follow.
No arity/row-set resampling, coefficientwise minimum, or shuffle factorial
is used.  All source/RS/projective normalizers are unchanged.

This is a new FLOAT proposal, not a certificate.  It reuses one pinned
original-one-state common measure and builds no regional provider.
"""
from collections import Counter
from fractions import Fraction
from hashlib import sha256
from itertools import combinations
import json
from math import comb, log2
from pathlib import Path
from time import perf_counter

import numpy as np

from hd_bit_syndrome_scalar_allarity import prepare_geometry, positive_nonempty_coefficients
from hd_inverse_column_control import INPUT, INPUT_SHA
from hd_local_shuffle_control import path_measure, first_active_words, shuffled_coefficients

DIRECTORY = Path(__file__).resolve().parent
HELPER_PINS = {
    'hd_local_shuffle_control.py': '7529573c7a52e039d06c70e41f71796b9e8b85f5d9534c25b7cd55a06c1b09e3',
    'hd_bit_syndrome_scalar_allarity.py': '59d347344a279eec99bf29c2376de831cbd27c286e47fb50229d8b2580cf2d5c',
    'hd_bit_syndrome_flat_geometry.py': '1f535b8d2f3d98cc04ca67cc0af34f75bc6bbfb3fabede8f2717f696138587d1',
    'hd_bit_syndrome_outer.py': 'a2c2d0fee734fad2c4b74da7850589f8e84180798dcf5f79beab30e06939292d',
    'hd_inverse_column_control.py': '40668f4bf673ef0859be429a82119769df169408d74aaaa52dfbb68787b9716b',
    'hd_bit_syndrome_counted_general.py': '980ec6e2d20f7eae25a2819396cbd2b2a5f901560cba4cf5dfeee20267de6283',
    'hd_post_refresh_quiet_rows.py': '75e38819d225e9771ade1ed90509b7d4582a2637b24b6522f5b62b0200eb49e7',
    'hd_two_sum_gf16_state_envelope.py': '5107b945ab18415191c48a3434a354c101d8fc5aeda66ec83a3fd2a5c3498af7',
}


def coefficient_classes(prepared, *, size=16, q=256):
    """Merge only fixed-A classes with BOTH identical positive majorants."""
    result = Counter()
    for item in prepared['classes']:
        a, multiplicity = item['arity'], item['multiplicity']
        if (type(a) is not int or not 1 <= a <= size
                or type(multiplicity) is not int or multiplicity < 1):
            raise ValueError('Exact nonempty arity and labelled multiplicity required')
        outside, fiber = item['counts_by_outputweight'], item['zero_fiber_counts']
        old = shuffled_coefficients(outside, fiber, item['d'], size=size, q=q)
        active = positive_nonempty_coefficients(outside, fiber)
        if len(active) > size+1 or sum(active) != q**a-1:
            raise ValueError('Full-fiber mass or output width disagrees')
        active += (0,) * (size+1-len(active))
        result[a, old, active] += multiplicity
    return result


def scalar_nonempty_logs(coefficients, log_means):
    """Evaluate a positive, exact-zero-constant polynomial; handle s=0."""
    values = tuple(coefficients)
    means = np.asarray(log_means, dtype=float)
    if (not values or values[0] != 0
            or any(type(c) is not int or c < 0 for c in values)
            or means.ndim != 1 or np.any(np.isnan(means)) or np.any(means > 0)):
        raise ValueError('Positive nonempty polynomial and scalar means in[0,1] required')
    answer = np.full(means.shape, -np.inf)
    # Horner in log space: no 0*(-inf), no rounded constant subtraction.
    for coefficient in reversed(values[1:]):
        answer = np.logaddexp2(answer+means, log2(coefficient) if coefficient else -np.inf)
    return answer+means


def kernel_logs(old_coefficients, active_coefficients, elementary):
    """Per-count full/nonzero pointwise minima, before the r-th power."""
    old = tuple(Fraction(c) for c in old_coefficients)
    elementary = np.asarray(elementary, dtype=float)
    if (len(old) < 2 or old[0] != 1 or any(c < 0 for c in old)
            or len(active_coefficients) != len(old) or elementary.ndim != 2
            or elementary.shape[1] != len(old) or np.any(np.isnan(elementary))
            or np.any(np.isposinf(elementary)) or not np.all(elementary[:, 0] == 0)):
        raise ValueError('Compatible elementary tables and exact constant-one majorant required')
    size = len(old)-1
    log_means = elementary[:, 1]-log2(size)
    # path_measure already rejects atoms>1.  This guard also covers direct callers.
    if np.any(log_means > 1e-12):
        raise ValueError('Scalar translated-fiber domination requires mean<=1')
    log_means = np.minimum(log_means, 0.)  # harmless FLOAT roundoff at all-one atoms
    old_logs = np.asarray([log2(c.numerator)-log2(c.denominator) if c else -np.inf for c in old])
    old_active = np.logaddexp2.reduce(elementary[:, 1:]+old_logs[None, 1:], axis=1)
    scalar_active = scalar_nonempty_logs(active_coefficients, log_means)
    old_full, scalar_full = np.logaddexp2(0., old_active), np.logaddexp2(0., scalar_active)
    return dict(full=np.minimum(old_full, scalar_full),
                active=np.minimum(old_active, scalar_active),
                old_full=old_full, old_active=old_active,
                scalar_full=scalar_full, scalar_active=scalar_active)


def column_moments(old_coefficients, active_coefficients, paths, elementary, groups):
    """Integrate the two valid per-count bounds in the same ordered measure."""
    paths = np.asarray(paths, dtype=float)
    if (type(groups) is not int or groups < 1 or paths.shape != (len(elementary), 2, 2)
            or np.any(np.isnan(paths)) or np.any(np.isposinf(paths))):
        raise ValueError('Positive group count and compatible ordered path table required')
    kernels = kernel_logs(old_coefficients, active_coefficients, elementary)
    matrices = {key: np.logaddexp2.reduce(paths+groups*kernels[key][:, None, None], axis=0)
                for key in ('full', 'active')}
    return matrices['full'], matrices['active']


def combine_words(classes, words, *, groups, prefix, row_dimension=8, q=256):
    """Whole-word group norm; labelled multiplicities stay outside each root."""
    if len(classes) != len(words) or type(groups) is not int or groups < 1:
        raise ValueError('One complete word result per exact class required')
    per = {}
    for a in sorted({key[0] for key in classes}):
        terms = [log2(multiplicity)+word/groups-row_dimension*a*log2(q-1)
                 for ((b, _, _), multiplicity), word in zip(classes.items(), words) if b == a]
        per[a] = float(prefix+groups*np.logaddexp2.reduce(terms))
    mixed = float(prefix+groups*np.logaddexp2.reduce([(v-prefix)/groups for v in per.values()]))
    return per, mixed


def run(arities, cache_path, *, max_seconds=240, progress=lambda x: None):
    started = perf_counter()
    arities = tuple(sorted(set(arities)))
    if not arities or any(type(a) is not int or a not in range(1, 17) for a in arities):
        raise ValueError('A nonempty subset of arities1..16 required')
    if not np.isfinite(max_seconds) or max_seconds <= 0:
        raise ValueError('A positive finite time budget is required')
    if Path(cache_path).exists():
        raise FileExistsError('Refuse to overwrite saved evidence')
    for name, digest in HELPER_PINS.items():
        if sha256((DIRECTORY/name).read_bytes()).hexdigest() != digest:
            raise ValueError('Pinned helper changed: '+name)
    raw = (DIRECTORY/INPUT).read_bytes()
    if sha256(raw).hexdigest() != INPUT_SHA:
        raise ValueError('Common input receipt changed')
    saved = json.loads(raw)
    for name, digest in saved['reviewed_sources'].items():
        if sha256((DIRECTORY/name).read_bytes()).hexdigest() != digest:
            raise ValueError('Frozen common dependency changed: '+name)
    if saved['groups'] != 18:
        raise ValueError('This frozen control is scoped to r18')
    def check_time():
        if perf_counter()-started > max_seconds:
            raise TimeoutError('Local-shuffle Maclaurin screen exceeded its budget')
    masks = tuple(sum(1 << i for i in row) for a in arities for row in combinations(range(16), a))
    geometry = prepare_geometry(source_masks=masks, max_masks=len(masks), check_time=check_time)
    classes = coefficient_classes(geometry)
    if sum(classes.values()) != len(masks):
        raise AssertionError('Labelled row-set partition changed')
    progress(dict(stage='exact_geometry', classes=len(classes), elapsed_seconds=perf_counter()-started))
    paths, elementary, empty = path_measure(saved['state_measure'])
    full, active = [], []
    for a, old, scalar_active in classes:
        check_time()
        f, n = column_moments(old, scalar_active, paths, elementary, saved['groups'])
        full.append(f)
        active.append(n)
    full, active = np.stack(full), np.stack(active)
    words = first_active_words(full, active, empty, saved['groups'])
    check_time()
    prefix, r = saved['common_prefix_log2'], saved['groups']
    per, mixed = combine_words(classes, words, groups=r, prefix=prefix)
    with Path(cache_path).open('xb') as handle:
        np.savez_compressed(handle, full=full, active=active, words=words,
            arities=np.array([a for a, _, _ in classes]), multiplicities=np.array(list(classes.values())))
    return dict(construction='Hstar with independent uniform S16 output shuffle per group and column',
        groups=r, source_class=saved['source_class'], source_count=saved['source_count'], mark=saved['mark'],
        arities=list(arities), per_arity_log2=per, mixed_log2=mixed, common_prefix_log2=prefix,
        exact_classes=len(classes), labelled_rowsets=len(masks), all_arities=(arities == tuple(range(1, 17))),
        all_cohorts_within_declared_arities=True, fixed_rowset_across_columns=True,
        active_column_floor=9, uniform_first9_cost_bits=(r-1)*log2(comb(16, 9)),
        bound_selection='per-fixed-A-class per-count pointwise minimum, separately full/nonempty',
        scalar_domination_interval=[0, 1], scalar_domination_coefficientwise=False,
        exact_constant_removed_before_evaluation=True, inverse_input_schedule=False,
        coefficient_classes=[dict(arity=a, multiplicity=m, shuffled_coefficients=[str(c) for c in old],
            scalar_nonempty_coefficients=list(nonempty)) for (a, old, nonempty), m in classes.items()],
        input_sha256=INPUT_SHA, cache=str(Path(cache_path).resolve()),
        cache_sha256=sha256(Path(cache_path).read_bytes()).hexdigest(),
        source_sha256=sha256(Path(__file__).read_bytes()).hexdigest(), helper_sha256=HELPER_PINS,
        elapsed_seconds=perf_counter()-started, directed=False, requires_fresh_directed_replay=True,
        whole_code_certificate=False, new_ensemble=True, native_implemented=False,
        native_benchmarks=0, bad_setup_charge_required=False)


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--arities', type=int, nargs='+', default=list(range(1, 17)))
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--cache', type=Path, required=True)
    parser.add_argument('--max-seconds', type=float, default=240)
    args = parser.parse_args()
    if args.output.exists() or args.cache.exists():
        parser.error('Refuse to overwrite saved evidence')
    answer = run(args.arities, args.cache, max_seconds=args.max_seconds,
                 progress=lambda x: print(json.dumps(x), flush=True))
    with args.output.open('x', encoding='utf-8') as handle:
        json.dump(answer, handle, indent=2, allow_nan=False)
        handle.write('\n')
    print(json.dumps({k: answer[k] for k in ('per_arity_log2', 'mixed_log2', 'elapsed_seconds')}), flush=True)
