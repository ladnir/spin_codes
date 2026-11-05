"""FLOAT diagnostic for exact-count peeling of ordered source profiles.

The frozen neighbor screen runs once and returns its exact-control report
unchanged. Instrument that same pass to cache exact symmetric source counts
and every ordered conditional chain, without repeating either computation.

For each run count c, sort conditional caps B_i in decreasing order. After
querying k counts, use sum_{i<k}N_i B_i+(N_c-sum_{i<k}N_i)B_k. Remaining
mass is an exact integer; B_k is the largest remaining cap. Stop when the
POST-source tail is at most2^-80/r, or no mass remains. The r tails then
contribute at most2^-80 in exact arithmetic. This implementation uses FLOAT
logs and is diagnostic only; it is not a directed certificate.

Report hypothetical unique exact-count queries needed by peeling, separately
from the complete-control queries actually performed. Opposite terminal
orientations share their count but not their chain. Keep huge chain values
in log form, never Fraction-convert them. Fraction tests are tiny only.
Inherit the unchanged screen CLI and180s/256MiB/20GiB guards. No files,
new parameter search, skipped unresolved coefficients, or automatic retry.
"""

from fractions import Fraction as F
import hashlib
import importlib.util
import json
import math
from pathlib import Path
import sys

sys.dont_write_bytecode = True
SCREEN_PIN = "c7cce6b8763dd8145e525061c82a9c2ad51139f2d0d66d8c870de476bb374870"
TAIL_BITS = 80


def load_screen():
    path = Path(__file__).with_name("cheap_lane_neighbor_screen.py")
    if hashlib.sha256(path.read_bytes()).hexdigest() != SCREEN_PIN:
        raise RuntimeError("frozen neighbor screen hash mismatch")
    spec = importlib.util.spec_from_file_location("peeling_neighbor", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def exact_tests():
    """Every prefix, reordered loose bounds, ties, zeros, and orientations."""
    cases = (
        [(3, F(1, 7), F(1, 7), "a"), (5, F(2, 3), F(2, 3), "b"),
         (2, F(1, 4), F(1, 4), "c")],
        [(3, F(9, 10), F(1), "a"), (5, F(1, 10), F(2), "b"),
         (7, F(1, 3), F(4, 3), "c")],
        [(2, F(1, 4), F(3, 4), "a"), (9, F(1, 2), F(3, 4), "b"),
         (4, F(3, 4), F(3, 4), "c")],
        [(0, F(3), F(4), "zero"), (5, F(1, 3), F(2, 3), "a"),
         (0, F(2), F(2), "zero2")],
        [(6, F(1, 3), F(1, 3), "orientation_pair"),
         (6, F(2, 3), F(2, 3), "orientation_pair"),
         (9, F(1, 5), F(1, 4), "diagonal")],
    )
    prefixes = query_checks = 0
    for entries in cases:
        ordered = sorted(enumerate(entries), key=lambda x: (-x[1][2], x[0]))
        exact = sum((n*actual for n, actual, _, _ in entries), F(0))
        mass = sum(n for n, _, _, _ in entries)
        prefix = F(0); previous = None; seen = {}
        for i in range(len(ordered)+1):
            ceiling = ordered[i][1][2] if i < len(ordered) else F(0)
            envelope = prefix+mass*ceiling
            assert envelope >= exact
            assert previous is None or envelope <= previous
            previous = envelope; prefixes += 1
            if i < len(ordered):
                n, actual, upper, key = ordered[i][1]
                assert 0 <= actual <= upper and 0 <= n <= mass
                if key in seen:
                    assert seen[key] == n
                else:
                    seen[key] = n; query_checks += 1
                prefix += n*upper; mass -= n
        assert mass == 0
        assert prefix == sum((n*upper for n, _, upper, _ in entries), F(0))
    assert len({key for _, _, _, key in cases[-1]}) == 2
    assert sum(n for n, _, _, key in cases[-1] if key == "orientation_pair") == 12
    return {"exact_prefix_bounds": prefixes, "unique_count_key_checks": query_checks,
            "loose_bound_reordering": True, "ties_and_zero_counts": True,
            "shared_orientation_count_not_chain": True}


def install(neighbor, metadata):
    original_join = neighbor.join_weight

    def joined(r, regions, counts, screen, log_source, guard, progress):
        cached_counts = {}
        ordered_runs = [[] for _ in range(r+1)]
        by_identity = {id(region): f for f, region in enumerate(regions)}
        raw_count = counts.profile_count
        raw_chain = screen.paired_chains
        actual_count_queries = 0

        def cached_count(pairs, weight, a, b):
            nonlocal actual_count_queries
            assert pairs == 4096 and weight == r
            key = (min(a, b), max(a, b))
            if key not in cached_counts:
                cached_counts[key] = raw_count(pairs, weight, a, b)
                actual_count_queries += 1
            return cached_counts[key]

        def recorded_chain(left, right):
            a, b = by_identity[id(left)], by_identity[id(right)]
            forward, reverse = raw_chain(left, right)
            assert math.isfinite(forward) and math.isfinite(reverse)
            assert cached_counts[min(a, b), max(a, b)] > 0
            c = r-a-b
            ordered_runs[c].append((forward, a, b))
            if a != b:
                ordered_runs[c].append((reverse, b, a))
            return forward, reverse

        counts.profile_count, screen.paired_chains = cached_count, recorded_chain
        try:
            exact_report = original_join(r, regions, counts, screen, log_source, guard, progress)
        finally:
            counts.profile_count, screen.paired_chains = raw_count, raw_chain
        assert sum(map(len, ordered_runs)) == exact_report["ordered_profiles"]
        zeros = [key for key, number in cached_counts.items() if not number]
        assert zeros == ([(r//2, r//2)] if r % 2 == 0 else [])
        logtwo = math.log(2)
        target = -TAIL_BITS*logtwo-math.log(r)
        queried = set()
        prefix_log = tail_log = peeled_log = exact_replay = -math.inf
        popped = hits = 0
        rows = []
        for c in range(1, r+1):
            guard()
            entries = sorted(ordered_runs[c], key=lambda x: (-x[0], x[1], x[2]))
            initial = counts.run_total(4096, r, c)
            assert initial > 0 and entries
            assert sum(cached_counts[min(a, b), max(a, b)] for _, a, b in entries) == initial
            exact_run = -math.inf
            for chain, a, b in entries:
                exact_run = screen.logadd(exact_run, math.log(cached_counts[min(a, b), max(a, b)])+chain)
            exact_replay = screen.logadd(exact_replay, log_source+exact_run)
            remaining, k, prefix = initial, 0, -math.inf
            while True:
                if not remaining:
                    tail = -math.inf
                else:
                    assert k < len(entries)
                    tail = log_source+math.log(remaining)+entries[k][0]
                if tail <= target:
                    break
                chain, a, b = entries[k]
                key = (min(a, b), max(a, b))
                if key in queried:
                    hits += 1
                else:
                    queried.add(key)
                number = cached_counts[key]
                assert 0 <= number <= remaining
                if number:
                    prefix = screen.logadd(prefix, math.log(number)+chain)
                remaining -= number; k += 1; popped += 1
            weighted_prefix = log_source+prefix
            upper = screen.logadd(weighted_prefix, tail)
            assert upper+1e-7 >= log_source+exact_run
            prefix_log = screen.logadd(prefix_log, weighted_prefix)
            tail_log = screen.logadd(tail_log, tail)
            peeled_log = screen.logadd(peeled_log, upper)
            rows.append({"runs": c, "ordered_profiles": len(entries), "peeled": k,
                "remaining_mass_log2": math.log(remaining)/logtwo if remaining else None,
                "post_source_tail_log2_FLOAT": tail/logtwo if math.isfinite(tail) else None,
                "tail_within_fixed_budget_FLOAT": tail <= target,
                "log2_envelope_FLOAT": upper/logtwo,
                "log2_exact_run_FLOAT": (log_source+exact_run)/logtwo,
                "tail_gap_bits_FLOAT": (upper-log_source-exact_run)/logtwo})
        exact_log = exact_report["log2_exact_profile_sum_FLOAT"]*logtwo
        assert abs(exact_replay-exact_log) < 1e-7
        assert peeled_log+1e-7 >= exact_log
        assert tail_log <= -TAIL_BITS*logtwo+1e-7
        assert abs(screen.logadd(prefix_log, tail_log)-peeled_log) < 1e-7
        assert len(queried) <= actual_count_queries
        largest_tails = sorted(rows, key=lambda x: x["post_source_tail_log2_FLOAT"]
                               if x["post_source_tail_log2_FLOAT"] is not None else -math.inf,
                               reverse=True)[:6]
        largest_peels = sorted(rows, key=lambda x: (x["peeled"], x["runs"]), reverse=True)[:6]
        report = {"weight": r, "directed": False, "negative_is_not_a_certificate": True,
            "log2_peeling_bound_FLOAT": peeled_log/logtwo,
            "log2_original_exact_control_FLOAT": exact_log/logtwo,
            "peeling_loss_bits_FLOAT": (peeled_log-exact_log)/logtwo,
            "post_source_tail_log2_FLOAT": tail_log/logtwo if math.isfinite(tail_log) else None,
            "per_run_post_source_tail_limit_log2": target/logtwo,
            "global_additive_tail_budget_bits": TAIL_BITS,
            "hypothetical_unique_count_queries": len(queried),
            "full_control_unique_count_queries_actually_performed": actual_count_queries,
            "hypothetical_query_fraction": len(queried)/actual_count_queries,
            "peeled_ordered_profiles": popped, "symmetric_count_reuses": hits,
            "additional_exact_count_queries": 0, "additional_finite_chain_calls": 0,
            "zero_count_profiles_checked": len(zeros),
            "all_run_mass_and_exact_control_replays_checked": True,
            "largest_residual_tails": largest_tails, "largest_peel_counts": largest_peels}
        metadata.append(report)
        progress("profile_peeling_complete", **report)
        guard()
        return exact_report

    neighbor.join_weight = joined


def main():
    tests = exact_tests()
    neighbor = load_screen()
    reports = []
    install(neighbor, reports)
    print(json.dumps({"phase": "profile_peeling_wrapper_start", "directed": False,
        "tests": tests, "neighbor_screen_pin": SCREEN_PIN,
        "source_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}), flush=True)
    neighbor.main()
    print(json.dumps({"phase": "profile_peeling_wrapper_complete", "directed": False,
        "whole_code_certificate": False, "negative_is_not_a_certificate": True,
        "original_exact_reports_unmodified": True, "results": reports, "tests": tests,
        "generated_bytes": 0, "neighbor_screen_pin": SCREEN_PIN,
        "source_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}, indent=2))


if __name__ == "__main__":
    main()
