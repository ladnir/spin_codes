"""Bounded exact joint counts on selected profiles; all other caps remain.

The selection is explicitly recorded. Automatic count-only ordering is
deterministic; an explicit list may come from a FLOAT transfer screen. The
count identity holds for every selection. Reuse counts across mark searches.
This module's full-weight transfer is FLOAT only, not a certificate.
"""

from hashlib import sha256
from heapq import heappush, heapreplace, nsmallest
from math import log2
from pathlib import Path
from types import MappingProxyType

import numpy as np

import hd_field128_gamma_exact_profile_counts as exact
import hd_field128_gamma_near_face_counts as near


MAX_EXACT_PROFILES = 4096


def source_pins():
    pins = {**near.source_pins(), **exact.source_pins()}
    path = Path(__file__).resolve()
    pins[path.name] = sha256(path.read_bytes()).hexdigest()
    return dict(sorted(pins.items()))


def check_source_pins(pins):
    directory = Path(__file__).resolve().parent
    if any(sha256((directory/name).read_bytes()).hexdigest() != digest for name, digest in pins.items()):
        raise RuntimeError('An exact-count refinement source or proof note changed')


class ExactRefinedCounts:
    """Exact selected counts with a frozen valid upper cap everywhere else.

    Canonical profiles have first>=second. Explicit profiles take priority.
    Automatic candidates have a positive minimum singleton count at most
    exact_singletons. Their order is shortest flow interval, then largest
    existing count cap, then lexicographic coordinates. This ordering is a
    cheap heuristic, not a claim about the unknown exact-count gap.
    """

    def __init__(self, pairs, weight, *, max_singletons=64, exact_singletons=64,
                 max_exact_profiles=4096, selected_profiles=(), budget=None):
        if (type(exact_singletons) is not int or not 0 <= exact_singletons <= 8192
                or type(max_exact_profiles) is not int or not 0 <= max_exact_profiles <= MAX_EXACT_PROFILES):
            raise ValueError('Invalid exact singleton cutoff or exact-profile cap')
        budget = budget or near.provider.old.Budget(60.)
        self.base = near.NearFaceCounts(pairs, weight, max_singletons=max_singletons, budget=budget)
        self.pairs, self.weight, self.full = pairs, weight, self.base.full
        if not isinstance(selected_profiles, (tuple, list)) or len(selected_profiles) > MAX_EXACT_PROFILES:
            raise ValueError('Expected a bounded explicit profile list')
        explicit, seen = [], set()
        for pair in selected_profiles:
            if (not isinstance(pair, (tuple, list)) or len(pair) != 2
                    or any(type(f) is not int or f not in self.full for f in pair)):
                raise ValueError('Explicit profile is outside the feasible matching marginals')
            key = (max(pair), min(pair))
            if key not in seen:
                seen.add(key)
                explicit.append(key)
        if len(explicit) > max_exact_profiles:
            raise ValueError('The explicit unique selection exceeds the exact-profile cap')
        eligible = skipped_sharp = 0
        def candidates():
            nonlocal eligible, skipped_sharp
            for first in self.full:
                singles = weight-2*first
                if singles > exact_singletons:
                    continue
                budget.check()
                if singles == 0:
                    skipped_sharp += first-self.full.start+1
                    continue
                for second in range(self.full.start, first+1):
                    if not second % 128:
                        budget.check()
                    key = first, second
                    if key in seen:
                        continue
                    cap = self.base.count(first, second)
                    if not cap:
                        skipped_sharp += 1
                        continue
                    eligible += 1
                    yield singles, -cap, first, second
        capacity = max_exact_profiles-len(explicit)
        chosen = nsmallest(capacity, candidates()) if capacity else []
        selected = explicit + [(a, b) for unused, unused_cap, a, b in
                               chosen]
        values = {}
        for first, second in selected:
            budget.check()
            value = exact.joint_count(pairs, weight, first, second)
            if not 0 <= value <= self.base.count(first, second):
                raise ArithmeticError('Exact count exceeds the prior proven profile cap')
            values[first, second] = value
        budget.check()
        self.exact_counts = MappingProxyType(values)
        self.selected_profiles = tuple(selected)
        self.explicit_profiles = tuple(explicit)
        self.exact_singletons, self.max_exact_profiles = exact_singletons, max_exact_profiles
        self.eligible_optional_profiles, self.skipped_sharp_profiles = eligible, skipped_sharp
        # One immutable sorted lookup for vector calls, rather than one array
        # mask and allocation per selected profile in every chunk.
        entries = []
        self._stride = pairs+1
        for (a, b), value in values.items():
            logarithm = log2(value) if value else -np.inf
            entries.append((a*self._stride+b, logarithm))
            if a != b:
                entries.append((b*self._stride+a, logarithm))
        entries.sort()
        self._keys = np.frombuffer(np.array([item[0] for item in entries], dtype=np.int64).tobytes(), dtype=np.int64)
        self._logs = np.frombuffer(np.array([item[1] for item in entries], dtype=np.float64).tobytes(), dtype=np.float64)

    def count(self, first, second):
        if type(first) is not int or type(second) is not int:
            raise ValueError('Integer full-pair counts required')
        if first not in self.full or second not in self.full:
            return 0
        key = (first, second) if first >= second else (second, first)
        value = self.exact_counts.get(key)
        return self.base.count(first, second) if value is None else value

    def log_counts(self, first, second, *, budget=None):
        result = self.base.log_counts(first, second, budget=budget)
        if not len(self._keys):
            return result
        first, second = np.broadcast_arrays(np.asarray(first), np.asarray(second))
        keys = first*self._stride+second
        locations = np.searchsorted(self._keys, keys)
        locations = np.minimum(locations, len(self._keys)-1)
        return np.where(self._keys[locations] == keys, self._logs[locations], result)

    def diagnostics(self):
        return dict(method='exact joint count on explicit selected symmetric profiles; frozen near-face cap elsewhere',
                    pairs=self.pairs, active_groups=self.weight,
                    maximum_refined_singletons=self.base.max_singletons,
                    optional_exact_singleton_cutoff=self.exact_singletons,
                    maximum_exact_canonical_profiles=self.max_exact_profiles,
                    explicit_canonical_profiles=[list(pair) for pair in self.explicit_profiles],
                    selected_canonical_profiles=[list(pair) for pair in self.selected_profiles],
                    exact_canonical_profiles=len(self.selected_profiles),
                    exact_ordered_profiles=len(self._keys),
                    eligible_optional_canonical_profiles=self.eligible_optional_profiles,
                    optional_candidates_left_unrefined=self.eligible_optional_profiles-(len(self.selected_profiles)-len(self.explicit_profiles)),
                    skipped_already_exact_optional_profiles=self.skipped_sharp_profiles,
                    optional_candidate_scan_performed=self.max_exact_profiles > len(self.explicit_profiles),
                    selection_order='explicit order first, then increasing min singleton count, decreasing prior integer cap, lexicographic profile',
                    both_matching_orientations_refined=True, all_profiles_retained=True,
                    exact_integer_count_path='hd_field128_gamma_exact_profile_counts.joint_count',
                    indexed_not_rotation_classes=True, frozen_providers_modified=False,
                    automatic_fallback_depends_on_transfer_matrices=False,
                    explicit_selection_may_depend_on_transfer_matrices=True,
                    count_identity_valid_for_every_selection=True)


def get_counts(pairs, weight, counts, budget, options):
    if counts is None:
        return ExactRefinedCounts(pairs, weight, budget=budget, **options)
    if options or type(counts) is not ExactRefinedCounts or (counts.pairs, counts.weight) != (pairs, weight):
        raise ValueError('Expected a matching exact-refinement count object without conflicting options')
    budget.check()
    return counts


def support_sum(pairs, weight, regions, *, repetitions=256, chunk=32,
                counts=None, budget=None, **count_options):
    """FLOAT full-profile transfer; reuse counts across marks at fixed P,r."""
    transfer = near.provider.profiles.transfer
    full = transfer.tail.full_counts(pairs, weight)
    regions = np.asarray(regions)
    if (regions.shape != (len(full), 2, 2) or np.any(np.isnan(regions)) or np.any(np.isposinf(regions))
            or type(repetitions) is not int or repetitions < 2 or repetitions % 2
            or type(chunk) is not int or not 1 <= chunk <= 128):
        raise ValueError('Invalid regional matrices, repetitions, or chunk size')
    budget = budget or near.provider.old.Budget(60.)
    counts = get_counts(pairs, weight, counts, budget, count_options)
    total, largest = -np.inf, None
    for start in range(0, len(full), chunk):
        budget.check()
        logs = counts.log_counts(full[start:start+chunk, None], full[None, :], budget=budget)
        product = transfer.multiply(regions[start:start+chunk, None], regions[None])
        repeated = transfer.power(product, repetitions//2)
        output = np.logaddexp2(repeated[..., 0, 0], repeated[..., 0, 1])
        terms = logs+output
        position = np.unravel_index(np.argmax(terms), terms.shape)
        candidate = float(terms[position])
        if largest is None or candidate > largest['log2_term']:
            even, odd = int(full[start+position[0]]), int(full[position[1]])
            largest = dict(even_full=even, odd_full=odd, runs=weight-even-odd,
                           log2_term=candidate, log2_count_upper=float(logs[position]),
                           log2_transfer=float(output[position]))
        total = float(np.logaddexp2(total, np.logaddexp2.reduce(terms.ravel())))
    budget.check()
    return dict(log2_support_sum_screen=total, largest_profile=largest,
                profile_pairs_evaluated=len(full)**2, common_potential_used=False,
                count_provenance=counts.diagnostics(), all_profiles_retained=True,
                arithmetic='FLOAT transfer with bounded exact joint-count refinement',
                directed=False, whole_code_certificate=False)


def select_profiles(pairs, weight, regions, *, max_profiles=4096, max_singletons=64,
                    eligible_singletons=64, repetitions=256, chunk=32,
                    budget=None, region_source=None):
    """FLOAT proposal for exact refinement, using a bounded top-k heap.

    Rank a canonical pair by the sum of both ordered contributions, with a
    diagonal contribution counted once. Every baseline profile is evaluated
    in chunks. Only eligible nonzero non-face pairs enter the heap. Selection
    affects optional counts only, so its FLOAT ordering need not be certified.
    """
    transfer = near.provider.profiles.transfer
    full = transfer.tail.full_counts(pairs, weight)
    regions = np.asarray(regions, dtype=np.float64)
    if (regions.shape != (len(full), 2, 2) or np.any(np.isnan(regions)) or np.any(np.isposinf(regions))
            or type(max_profiles) is not int or not 0 <= max_profiles <= MAX_EXACT_PROFILES
            or type(eligible_singletons) is not int or not 0 <= eligible_singletons <= 8192
            or type(repetitions) is not int or repetitions < 2 or repetitions % 2
            or type(chunk) is not int or not 1 <= chunk <= 128):
        raise ValueError('Invalid selection matrices, bounds, repetitions, or chunk size')
    budget = budget or near.provider.old.Budget(60.)
    counts = near.NearFaceCounts(pairs, weight, max_singletons=max_singletons, budget=budget)
    heap, eligible = [], 0
    for start in range(0, len(full), chunk):
        budget.check()
        firsts = full[start:start+chunk]
        logs = counts.log_counts(firsts[:, None], full[None, :], budget=budget)
        forward = transfer.power(transfer.multiply(regions[start:start+chunk, None], regions[None]), repetitions//2)
        reverse = transfer.power(transfer.multiply(regions[None], regions[start:start+chunk, None]), repetitions//2)
        forward = np.logaddexp2(forward[..., 0, 0], forward[..., 0, 1])
        reverse = np.logaddexp2(reverse[..., 0, 0], reverse[..., 0, 1])
        for offset, first in enumerate(firsts):
            first = int(first)
            singles = weight-2*first
            if not 0 < singles <= eligible_singletons:
                continue
            for second in range(int(full[0]), first+1):
                column = second-int(full[0])
                if not np.isfinite(logs[offset, column]):
                    continue
                score = float(logs[offset, column] + (forward[offset, column] if first == second else
                              np.logaddexp2(forward[offset, column], reverse[offset, column])))
                if not np.isfinite(score):
                    continue
                eligible += 1
                item = score, -first, -second
                if len(heap) < max_profiles:
                    heappush(heap, item)
                elif max_profiles and item > heap[0]:
                    heapreplace(heap, item)
            budget.check()
    budget.check()
    ranked = sorted(heap, reverse=True)
    return dict(selected_profiles=[[-first, -second] for score, first, second in ranked],
                selected_canonical_profiles=len(ranked), eligible_canonical_profiles=eligible,
                maximum_selected_profiles=max_profiles, maximum_eligible_singletons=eligible_singletons,
                baseline_maximum_refined_singletons=max_singletons,
                selection_scores_FLOAT=[score for score, first, second in ranked],
                selection_rule='largest sum of both ordered baseline contributions; diagonal once; lexicographic tie break',
                baseline_profile_pairs_evaluated=len(full)**2, all_baseline_profiles_evaluated=True,
                selection_does_not_omit_supports=True, exact_singletons_for_explicit_replay=0,
                pairs=pairs, active_groups=weight, repetitions=repetitions, chunk=chunk,
                regional_logs_sha256=sha256(regions.tobytes(order='C')).hexdigest(),
                regional_array_shape=list(regions.shape), region_source=region_source,
                source_sha256=source_pins(), arithmetic='FLOAT selection only',
                directed=False, selected_weight_certificate=False, whole_code_certificate=False)
