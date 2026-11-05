"""Directed shared route coefficients for requested source counts1..128.

Require explicit rational z and alpha. Apply alpha entrywise to the exact
local rational majorants, then retain the unnormalized Arb matrix polynomial
on singleton/full-pair degrees(a,b) with a+2b<=G. Initially a+b<=32 as well.
Seven ordered positive squarings give128 batches. Nonnegative degrees make
every truncation exact for the retained coefficients.

Check every initial coefficient mass against C32(a,b), where
C_n(a,b)=binom(n,a)binom(n-a,b)2^a. The multinomial Vandermonde identity
gives C4096 after squaring, so divide each required coefficient by that
exact integer once. Arb's arbitrary exponents retain tiny positive entries.

For each requested g, use h0,c=(256-g)/256, and check every allowed matching
profile against exact2^-target_bits. The full source-set binomial remains
outside alpha. Maximize profiles; do not add a second profile union. The
terminal is the positive f<=M*ones envelope, without a physical final reset.

With --run-census, omit the source binomial from each profile bound. Group
profiles by c=g-even_full-odd_full occupied runs and maximize certified
upper endpoints within each c. Multiply by the exact source-set count
(8192/c)binom(g-1,c-1)binom(8192-g-1,c-1), then sum over c. Compare this
entire sum with the exact threshold. Endpoint maxima, counts, and this
sum use exact rationals; source multiplicities are outside alpha.

Default output aggregates all checks by g. --details additionally reports
individual profiles within a512-KiB output cap. No files are written. Root
owns numerical runs; use Python -B with an external time and memory limit.
"""

import argparse
from fractions import Fraction
from functools import lru_cache
import hashlib
import importlib.util
import json
import math
from pathlib import Path
import shutil
import sys
import time

sys.dont_write_bytecode = True
import verify_small_state_exact_route as reference


OUTPUT_LIMIT = 512 * 1024
DEFAULT_GROUPS = (61, 62, 63, 64)


def verify(groups=DEFAULT_GROUPS, z=None, alpha=None, precision=128,
           target_bits=60, seconds=50., memory_mib=1024, details=False,
           run_census=False):
    groups = tuple(sorted(groups))
    if (not groups or len(groups) != len(set(groups))
            or any(not isinstance(g, int) or not 1 <= g <= 128 for g in groups)
            or z is None or alpha is None or precision not in (128, 256)
            or not 1 <= target_bits <= 10000 or not 0 < seconds <= 50
            or not 64 <= memory_mib <= 1024):
        raise ValueError("Require distinct groups1..128, explicit rational z/alpha, "
                         "precision128/256, bounded target, <=50 seconds and<=1024 MiB")
    z, alpha = Fraction(z), Fraction(alpha)
    if not 0 < z < 1 or not 0 < alpha <= 1:
        raise ValueError("Require rational z in(0,1), alpha in(0,1]")
    free = shutil.disk_usage(Path(__file__).resolve().parent).free
    if free < 20*1024**3:
        raise RuntimeError("At least20 GiB free space is required")
    started, next_memory_check = time.monotonic(), 0.

    def guard():
        nonlocal next_memory_check
        now = time.monotonic()
        if now-started >= seconds:
            raise TimeoutError("time-cap")
        if now >= next_memory_check:
            if reference.resident_bytes() > memory_mib*1024**2:
                raise MemoryError("memory-cap")
            next_memory_check = now+.05

    max_groups = max(groups)
    progress = dict(max_groups=max_groups, base_exact_mass_checks=0,
                    completed_pair_slots=0, stages=[], complete_regional_table=False)
    reports, dependencies = [], {}
    status, error, active = "directed-shared-fixed-witness-replay", None, None
    try:
        guard()
        directory = (reference.OLD_THEORY if reference.OLD_THEORY.is_dir()
                     else Path(__file__).resolve().parent)
        sys.path.insert(0, str(directory))
        geometry_path = directory/"cheap_lane_rank_debit_gate.py"
        if hashlib.sha256(geometry_path.read_bytes()).hexdigest() != reference.FROZEN_HASH:
            raise RuntimeError("Frozen striped-geometry hash mismatch")
        spec = importlib.util.spec_from_file_location("small_state_shared_arb_geometry", geometry_path)
        geometry = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(geometry)
        from hd_common_refresh_packets import local_moments, matrix
        from flint import arb, arb_mat, ctx, fmpq
        ctx.prec = precision
        dependencies[geometry_path.name] = reference.FROZEN_HASH

        def rational(value):
            value = Fraction(value)
            return arb(value.numerator)/value.denominator

        exponent = rational(alpha)
        local = local_moments(z, q=256, h=4)
        epsilon = Fraction(1, 256**16)
        guard()

        @lru_cache(maxsize=4845)
        def packet(occupancy):
            guard()
            quiet, birth, ret, live = matrix(occupancy, z, local, q=256, h=4)
            total = quiet+birth
            all_live = epsilon*total+(1-epsilon)*live
            entries = ((quiet, total, epsilon*quiet+(1-epsilon)*ret, all_live)
                       if sum(occupancy) else
                       (Fraction(1), Fraction(0), Fraction(0), all_live))
            if min(entries) < 0 or max(entries) > 1:
                raise ArithmeticError("Exact rational local majorant left[0,1]")
            result = arb_mat(2, 2)
            for index, value in enumerate(entries):
                entry = rational(value)
                if alpha != 1 and value:
                    entry = entry**exponent
                if not entry.is_finite() or not entry >= 0:
                    raise ArithmeticError("Invalid directed local entry")
                result[index//2, index%2] = entry
            return result

        polynomial = {}
        for full in range(min(32, max_groups//2)+1):
            for singles in range(min(32-full, max_groups-2*full)+1):
                guard()
                coefficient, mass = arb_mat(2, 2), 0
                for occupancy, count in geometry.occupancy_counts(singles, full):
                    coefficient += packet(tuple(sorted(occupancy)))*count
                    mass += count
                expected = math.comb(32, singles)*math.comb(32-singles, full)*2**singles
                if mass != expected:
                    raise ArithmeticError("Exact one-batch coefficient mass failed")
                polynomial[(singles, full)] = coefficient
                progress["base_exact_mass_checks"] += 1
        progress["completed_pair_slots"] = 32

        for _ in range(7):
            guard()
            # Degree order stops each right scan at the truncation boundary.
            # Ordered products retain the two-state matrix chronology.
            entries = sorted((a+2*b, a, b, value) for (a, b), value in polynomial.items())
            following = {}
            for degree, a, b, left in entries:
                guard()
                remaining = max_groups-degree
                for right_index, (right_degree, c, d, right) in enumerate(entries):
                    if right_degree > remaining:
                        break
                    if (right_index & 255) == 255:
                        guard()
                    key = (a+c, b+d)
                    value = left*right
                    if key in following:
                        following[key] += value
                    else:
                        following[key] = value
            polynomial = following
            progress["completed_pair_slots"] *= 2
            progress["stages"].append(dict(pair_slots=progress["completed_pair_slots"],
                                            coefficients=len(polynomial)))

        expected_coefficients = sum(max_groups-2*b+1 for b in range(max_groups//2+1))
        if len(polynomial) != expected_coefficients:
            raise ArithmeticError("Incomplete shared triangular coefficient table")
        regions = {}
        for group in groups:
            for full in range(group//2+1):
                guard()
                singles = group-2*full
                key = (singles, full)
                normalizer = math.comb(4096, singles)*math.comb(4096-singles, full)*2**singles
                regions[key] = polynomial[key]/arb(normalizer)
        progress.update(complete_regional_table=True,
                        computed_coefficients=len(polynomial), extracted_regions=len(regions),
                        local_occupancies=packet.cache_info().currsize,
                        normalization="exact C32 checks; multinomial Vandermonde identity; exact C4096 division")

        # Same conservative Gamma<2 receipt as the existing directed verifier.
        if not 2*(16*(reference.LENGTH//64)) < (1 << 32)-1:
            raise ArithmeticError("The Gamma<2 receipt is unavailable")
        beta = arb(256)**160/arb(255)**32
        marker_factor = rational(z)**(-reference.WEIGHT_LIMIT)
        target, log_two = arb(2)**(-target_bits), arb(2).log()
        exact_target = fmpq(1, 1 << target_bits) if run_census else None
        for group in groups:
            active = dict(groups=group, checked_profiles=0)
            guard()
            retained = Fraction(256-group, 256)
            conditional = 2*beta**group*marker_factor/(255*rational(retained))
            powered_conditional = conditional if alpha == 1 else conditional**exponent
            source_sets = math.comb(8192, group)
            prefactor = powered_conditional if run_census else source_sets*powered_conditional
            rows, failures, worst, worst_upper = [], [], None, None
            run_maxima, run_profile_counts = {}, {}
            checked, passed_count = 0, 0
            for even in range(group//2+1):
                for odd in range(group//2+1):
                    if even+odd >= group:
                        continue
                    guard()
                    run = regions[(group-2*even, even)]*regions[(group-2*odd, odd)]
                    for _ in range(7):
                        run = run*run
                    bound = prefactor*(run[0, 0]+run[0, 1])
                    if not bound.is_finite() or not bound >= 0:
                        raise ArithmeticError(f"Invalid profile bound at g={group}, even={even}, odd={odd}")
                    checked += 1
                    active["checked_profiles"] = checked
                    row = dict(even_full=even, odd_full=odd)
                    if run_census:
                        # An Arb upper endpoint is a dyadic rational. Compare
                        # those exact rationals, never overlapping balls.
                        upper = bound.upper().fmpq()
                        runs = group-even-odd
                        if not 1 <= runs <= group:
                            raise ArithmeticError("Invalid occupied-run class")
                        run_profile_counts[runs] = run_profile_counts.get(runs, 0)+1
                        previous = run_maxima.get(runs)
                        if previous is None or upper > previous[0]:
                            run_maxima[runs] = (upper, even, odd)
                        row["cycle_runs"] = runs
                        is_worst = worst_upper is None or upper > worst_upper
                        if details or is_worst:
                            row["log2_conditional_bound_interval"] = str(bound.log()/log_two)
                    else:
                        passed = bool(bound < target)
                        passed_count += passed
                        row["passed"] = passed
                        upper = bound.upper()
                        is_worst = worst_upper is None or upper > worst_upper
                        if details or not passed or is_worst:
                            row["log2_bound_interval"] = str(bound.log()/log_two)
                        if not passed:
                            failures.append(row)
                    if is_worst:
                        worst, worst_upper = row, upper
                    if details:
                        rows.append(row)
            expected_profiles = (group//2+1)**2-(1 if group % 2 == 0 else 0)
            if checked != expected_profiles:
                raise ArithmeticError("Incomplete matching-profile coverage")
            guard()
            if run_census:
                if (set(run_maxima) != set(range(1, group+1))
                        or sum(run_profile_counts.values()) != checked):
                    raise ArithmeticError("Incomplete occupied-run profile coverage")
                total_count, total_upper = 0, fmpq(0)
                run_terms, worst_term, worst_term_upper = [], None, None
                for runs in range(1, group+1):
                    guard()
                    numerator = (8192*math.comb(group-1, runs-1)
                                 * math.comb(8192-group-1, runs-1))
                    count, remainder = divmod(numerator, runs)
                    if remainder or count <= 0:
                        raise ArithmeticError("Invalid exact source-cycle count")
                    endpoint, even, odd = run_maxima[runs]
                    term = count*endpoint
                    total_count += count
                    total_upper += term
                    is_worst_term = worst_term_upper is None or term > worst_term_upper
                    if details or is_worst_term:
                        term_row = dict(runs=runs, source_sets=str(count),
                            checked_profiles=run_profile_counts[runs],
                            even_full=even, odd_full=odd,
                            log2_conditional_upper_interval=str(arb(endpoint).log()/log_two),
                            log2_counted_upper_interval=str(arb(term).log()/log_two))
                        if details:
                            run_terms.append(term_row)
                        if is_worst_term:
                            worst_term, worst_term_upper = term_row, term
                if total_count != source_sets:
                    raise ArithmeticError("Exact source-cycle counts do not partition source sets")
                report = dict(groups=group, checked_profiles=checked, valid_profiles=checked,
                    passed=bool(total_upper < exact_target), complete_profile_cover=True,
                    retained_fraction=str(retained), worst_profile=worst,
                    run_census=True, run_classes=group, complete_run_cover=True,
                    source_sets=str(source_sets), counted_source_sets=str(total_count),
                    exact_source_count_partition=True,
                    cycle_run_union_log2_upper_interval=str(arb(total_upper).log()/log_two),
                    worst_run_term=worst_term,
                    upper_maximum="exact rational comparison of certified Arb upper endpoints",
                    acceptance="exact rational counted sum < exact2^(-target_bits)")
                if details:
                    report["run_terms"] = run_terms
            else:
                report = dict(groups=group, checked_profiles=checked, passed_profiles=passed_count,
                              passed=checked == passed_count, complete_profile_cover=True,
                              retained_fraction=str(retained), worst_profile=worst,
                              failed_profile_count=len(failures), failing_profiles=failures[:8],
                              failing_profiles_truncated=len(failures) > 8)
            if details:
                report["profiles"] = rows
            guard()
            reports.append(report)
            active = None
        # Several retained providers import their dependencies lazily. Record
        # the actually loaded module paths after evaluating the local moments,
        # as the original directed verifier does.
        dependencies.update({name+".py": hashlib.sha256(Path(sys.modules[name].__file__).read_bytes()).hexdigest()
                             for name in reference.DEPENDENCIES})
        guard()
    except (TimeoutError, MemoryError, ArithmeticError, RuntimeError) as failure:
        status = str(failure) if isinstance(failure, (TimeoutError, MemoryError)) else "verification-error"
        error = str(failure)

    completed = [row["groups"] for row in reports]
    complete = len(completed) == len(groups)
    passed = complete and all(row["passed"] for row in reports) and error is None
    dependencies[Path(__file__).name] = hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    dependencies[Path(reference.__file__).name] = hashlib.sha256(Path(reference.__file__).read_bytes()).hexdigest()
    return dict(schema="small-state-shared-route-arb-v1", status=status, error=error,
                passed=passed, requested_groups=list(groups), completed_groups=completed,
                unprocessed_groups=[g for g in groups if g not in completed], incomplete_group=active,
                complete_source_count_cover=complete, reports=reports,
                checked_profiles=sum(row["checked_profiles"] for row in reports),
                precision_bits=precision, target_bits=target_bits, z=str(z), alpha=str(alpha),
                run_census=run_census,
                profile_bound_kind=("conditional powered bound without source multiplicity"
                                    if run_census else "unconditional source-binomial-prefactored bound"),
                acceptance=("exact rational counted sum of certified profile upper endpoints < exact2^(-target_bits)"
                            if run_census else "each complete source-prefactored profile bound < exact2^(-target_bits)"),
                union_verified_bits=target_bits-(len(groups)-1).bit_length() if passed else None,
                source_mode="cutoff", missing_groups=0, gamma_upper="2",
                source_prefactor=("(2*(256^160/255^32)^g*z^(-734003)/(255*(256-g)/256))^alpha"
                                  if run_census else
                                  "binom(8192,g)*(2*(256^160/255^32)^g*z^(-734003)/(255*(256-g)/256))^alpha"),
                coefficient_engine=progress, coefficient_field_bits=32, auxiliary_field_order=256,
                output_length=reference.LENGTH, output_limit=reference.WEIGHT_LIMIT,
                terminal="positive f<=M*ones; no physical final refresh or reset charge",
                profile_accounting=("exact cycle-run counts outside alpha; upper maximum within each run class; sum run classes"
                                    if run_census else "full binomial per profile; maximum, no extra profile sum"),
                seconds=time.monotonic()-started, peak_resident_bytes=reference.resident_bytes(),
                time_cap_seconds=seconds, memory_cap_bytes=memory_mib*1024**2,
                output_cap_bytes=OUTPUT_LIMIT, free_bytes_at_start=free, dependencies=dependencies,
                directed=True, whole_code_certificate=False, generated_file_bytes=0)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--groups", type=int, nargs="+", default=DEFAULT_GROUPS)
    parser.add_argument("--z", type=Fraction, required=True)
    parser.add_argument("--alpha", type=Fraction, required=True)
    parser.add_argument("--precision", type=int, choices=(128, 256), default=128)
    parser.add_argument("--target-bits", type=int, default=60)
    parser.add_argument("--seconds", type=float, default=50.)
    parser.add_argument("--memory-mib", type=int, default=1024)
    parser.add_argument("--details", action="store_true")
    parser.add_argument("--run-census", action="store_true",
                        help="sum exact cycle-run source counts times certified class upper maxima")
    args = parser.parse_args()
    result = verify(args.groups, args.z, args.alpha, args.precision, args.target_bits,
                    args.seconds, args.memory_mib, args.details, args.run_census)
    encoded = json.dumps(result, allow_nan=False, sort_keys=True)
    if len(encoded.encode("utf-8")) > OUTPUT_LIMIT:
        for report in result["reports"]:
            report.pop("profiles", None)
            report.pop("run_terms", None)
        result["details_omitted_for_output_cap"] = True
        encoded = json.dumps(result, allow_nan=False, sort_keys=True)
    if len(encoded.encode("utf-8")) > OUTPUT_LIMIT:
        raise RuntimeError("Compact summary exceeded its output cap")
    print(encoded, flush=True)
    if not result["passed"]:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
