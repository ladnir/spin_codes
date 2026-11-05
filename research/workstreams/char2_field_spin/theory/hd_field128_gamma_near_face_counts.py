"""Exact deletion-count refinement and a FLOAT chronological-transfer replay.

For P labeled pairs, weight r, matching full count f, and cyclic run count k,
put s=r-2f. Delete the s mixed pairs, then contract the surviving equal pairs.
The residual word has length P-s, weight f, and between max(0,k-s) and k runs.
Singleton positions, orientations, and the indexed residual word reconstruct
the original word. Thus C(P,s)*2**s times the residual run sum is an upper
count, without a rotation factor. Constant and empty words have zero runs.

Both matching bounds are minimized with the frozen profile count. A cutoff
limits optional refinements only; every profile remains in the support sum.
Integer counts are exact. The matrix replay is FLOAT, not a certificate.
"""

import ast
from fractions import Fraction as F
from hashlib import sha256
import json
from math import comb, log2
from pathlib import Path

import numpy as np

import hd_field128_gamma_all_pair_screen as provider
import hd_field128_gamma_profile_certificate as exact


def near_face_count_upper(pairs, weight, full, runs):
    """One matching's deletion bound, including an empty residual cycle."""
    if any(type(x) is not int for x in (pairs, weight, full, runs)):
        raise ValueError('Integer count arguments required')
    if pairs < 1 or not 0 <= weight <= 2*pairs:
        raise ValueError('Invalid pair count or weight')
    singles = weight-2*full
    if not (0 <= singles <= pairs and 0 <= full <= pairs-singles and 0 <= runs <= weight):
        return 0
    values = exact.tail.circular_run_counts(pairs-singles, full)
    return comb(pairs, singles)*(1 << singles)*sum(values[max(0, runs-singles):runs+1])


class NearFaceCounts:
    """Exact O(1) profile queries after bounded near-face prefix preparation.

    With cutoff c, at most floor(c/2)+1 residual arrays are retained. There
    is no profile-pair table. The public construction uses P=4096; small P
    is supported for exhaustive tests. Timeout never returns a partial sum.
    """

    def __init__(self, pairs, weight, *, max_singletons=32, budget=None):
        if (type(pairs) is not int or not 1 <= pairs <= 4096
                or type(weight) is not int or not 0 <= weight <= 2*pairs
                or type(max_singletons) is not int or not 0 <= max_singletons <= 8192):
            raise ValueError('Invalid count domain or optional singleton cutoff')
        budget = budget or provider.old.Budget(60.)
        self.old = exact.ProfileCounts(pairs, weight, budget)
        self.pairs, self.weight, self.full = pairs, weight, self.old.full
        self.max_singletons, self.prefixes = max_singletons, {}
        self._near_logs = {}
        self._run_logs = np.array([log2(v) if v else -np.inf for v in self.old.runs])
        self._marginal_logs = np.array([log2(self.old.marginals[f]) for f in self.full])
        for f in self.full:
            singles = weight-2*f
            if singles > max_singletons:
                continue
            budget.check()
            values = exact.tail.circular_run_counts(pairs-singles, f)
            prefix, total = [0], 0
            for k, value in enumerate(values):
                if not k % 128:
                    budget.check()
                total += value
                prefix.append(total)
            self.prefixes[f] = (singles, comb(pairs, singles)*(1 << singles), tuple(prefix))
        budget.check()

    def matching_count(self, full, runs):
        """Cached optional bound; None means this matching was not refined."""
        data = self.prefixes.get(full)
        if data is None:
            return None
        singles, factor, prefix = data
        if runs < 0 or runs > self.weight:
            return 0
        lo, hi = min(max(0, runs-singles), len(prefix)-1), min(runs+1, len(prefix)-1)
        return factor*(prefix[hi]-prefix[lo])

    def count(self, first, second):
        if type(first) is not int or type(second) is not int:
            raise ValueError('Integer full-pair counts required')
        if first not in self.full or second not in self.full:
            return 0
        runs = self.weight-first-second
        result = self.old.count(first, second)
        for full in (first, second):
            refined = self.matching_count(full, runs)
            if refined is not None:
                result = min(result, refined)
        return result

    def log_counts(self, first, second, *, budget=None):
        """FLOAT logs of exact integer caps; inputs are feasible full counts."""
        budget = budget or provider.old.Budget(60.)
        first, second = np.broadcast_arrays(np.asarray(first), np.asarray(second))
        if (not np.issubdtype(first.dtype, np.integer) or not np.issubdtype(second.dtype, np.integer)
                or np.any(first < self.full.start) or np.any(first >= self.full.stop)
                or np.any(second < self.full.start) or np.any(second >= self.full.stop)):
            raise ValueError('Complete feasible integer full-count domain required')
        runs = self.weight-first-second
        result = np.array(np.minimum(self._run_logs[runs], np.minimum(
            self._marginal_logs[first-self.full.start], self._marginal_logs[second-self.full.start])), copy=True)
        if self.old.contracted is not None:
            contracted = np.full(self.weight+1, -np.inf)
            contracted[:len(self.old.contracted)] = [log2(v) if v else -np.inf for v in self.old.contracted]
            result = np.where((2*first == self.weight) | (2*second == self.weight), contracted[runs], result)
        for full in self.prefixes:
            mask = (first == full) | (second == full)
            if not np.any(mask):
                continue
            budget.check()
            if full not in self._near_logs:
                values = np.full(self.weight+1, -np.inf)
                for k in range(self.weight+1):
                    if not k % 128:
                        budget.check()
                    value = self.matching_count(full, k)
                    if value:
                        values[k] = log2(value)
                self._near_logs[full] = values
            result[mask] = np.minimum(result[mask], self._near_logs[full][runs[mask]])
        budget.check()
        return result


def support_sum(pairs, weight, regions, *, repetitions=256, chunk=32, max_singletons=32, budget=None):
    """All profiles, unchanged chronological regional matrices, new counts."""
    transfer = provider.profiles.transfer
    full = transfer.tail.full_counts(pairs, weight)
    regions = np.asarray(regions)
    if (regions.shape != (len(full), 2, 2) or np.any(np.isnan(regions)) or np.any(np.isposinf(regions))
            or type(repetitions) is not int or repetitions < 2 or repetitions % 2
            or type(chunk) is not int or not 1 <= chunk <= 128):
        raise ValueError('Invalid regional matrices, repetitions, or chunk size')
    budget = budget or provider.old.Budget(60.)
    counts = NearFaceCounts(pairs, weight, max_singletons=max_singletons, budget=budget)
    total, largest = -np.inf, None
    for start in range(0, len(full), chunk):
        budget.check()
        logs = counts.log_counts(full[start:start+chunk, None], full[None, :], budget=budget)
        product = transfer.multiply(regions[start:start+chunk, None], regions[None])
        repeated = transfer.power(product, repetitions//2)
        output = np.logaddexp2(repeated[..., 0, 0], repeated[..., 0, 1])
        terms = logs+output
        index = np.unravel_index(np.argmax(terms), terms.shape)
        candidate = float(terms[index])
        if largest is None or candidate > largest['log2_term']:
            even, odd = int(full[start+index[0]]), int(full[index[1]])
            largest = dict(even_full=even, odd_full=odd, runs=weight-even-odd,
                           log2_term=candidate, log2_count_upper=float(logs[index]),
                           log2_transfer=float(output[index]))
        total = float(np.logaddexp2(total, np.logaddexp2.reduce(terms.ravel())))
    budget.check()
    return dict(log2_support_sum_screen=total, largest_profile=largest,
                profile_pairs_evaluated=len(full)**2, common_potential_used=False,
                max_singletons=max_singletons, refined_matching_full_counts=list(counts.prefixes),
                support_count='min(frozen profile upper, both cached singleton-deletion bounds)',
                arithmetic='FLOAT transfer with exact integer upper counts',
                directed=False, whole_code_certificate=False)


def source_pins():
    directory, pending, seen = Path(__file__).resolve().parent, [Path(__file__).name], {}
    while pending:
        name = pending.pop()
        if name in seen:
            continue
        raw = (directory/name).read_bytes()
        seen[name] = sha256(raw).hexdigest()
        for node in ast.walk(ast.parse(raw)):
            modules = ([item.name for item in node.names] if isinstance(node, ast.Import)
                       else [node.module] if isinstance(node, ast.ImportFrom) and node.module else [])
            for module in modules:
                target = module.split('.')[0]+'.py'
                if (directory/target).is_file() and target not in seen:
                    pending.append(target)
    return dict(sorted(seen.items()))


def main():
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source', type=Path, required=True)
    parser.add_argument('--r', type=int, required=True)
    parser.add_argument('--state', type=int, choices=provider.gamma.STATES, default=30)
    parser.add_argument('--distance', default='2/5')
    parser.add_argument('--max-singletons', type=int, default=32)
    parser.add_argument('--max-seconds', type=float, default=120.)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        parser.error('Output must be unused')
    raw, pins = args.source.read_bytes(), source_pins()
    provenance = dict(path=str(args.source.resolve()), sha256=sha256(raw).hexdigest())
    proposal = provider.search.prepare_proposal(json.loads(raw), provenance, r=args.r, state=args.state,
                                                distance=F(args.distance))
    witness = proposal['witness']
    budget = provider.old.Budget(args.max_seconds)
    bank = provider.AliasBankScreen(float(F(witness['z'])), float(F(witness['alpha'])),
            tuple(tuple(float(F(x)) for x in pair) for pair in witness['marker_log_odds_exact']),
            float(F(witness['log2_potential_ratio_exact'])), state=args.state,
            distance=F(args.distance), grid=witness['grid'], floor_policy='positive-fallback', _budget=budget)
    regions = bank.regional(args.r)
    baseline = provider.profiles.support_sum(provider.old.PAIRS, args.r, regions, budget=budget)
    result = support_sum(provider.old.PAIRS, args.r, regions, max_singletons=args.max_singletons, budget=budget)
    prefix = bank.alpha*bank.density*args.r+bank.constant
    result.update(bank.metadata(), active_groups=args.r, target_distance=str(F(args.distance)),
                  log2_outer_prefix=prefix, log2_upper_screen=prefix+result['log2_support_sum_screen'],
                  baseline_same_regions=dict(baseline, log2_upper_screen=prefix+baseline['log2_support_sum_screen']),
                  regional_alias_diagnostics=bank.regional_diagnostics, witness=bank.witness(),
                  proposal=proposal, source=provenance, source_sha256=pins,
                  old_numerical_bounds_reused=False, regional_caps_unchanged=True,
                  scope='one selected weight; exact near-face upper counts and all-pair caps; FLOAT only')
    if source_pins() != pins or args.source.read_bytes() != raw:
        raise RuntimeError('Pinned near-face replay source changed')
    provider.search.write_new(args.output, result)
    print(json.dumps({key: result[key] for key in ('active_groups', 'log2_upper_screen', 'largest_profile',
                                                  'directed', 'whole_code_certificate')}))


if __name__ == '__main__':
    main()
