"""Shared positive FLOAT route coefficients for actual source counts1..128.

The regional polynomial uses a for singleton pairs and b for full pairs.
Keep every coefficient with a+2b<=G in one array. A batch has32 pair slots,
so its initial polynomial also requires a+b<=32. Seven ordered squarings
then cover128 batches without recomputing a polynomial for each profile.

For n pair slots, C_n(a,b)=binom(n,a)binom(n-a,b)2^a is exact. Store the
coefficient divided by C_n. Squaring uses the nonnegative weights
C_n(a,b)C_n(c,d)/C_(2n)(a+c,b+d). Their sum at each target is one by the
multinomial Vandermonde identity. Check initial masses with exact integers
and report the floating weight-sum residual at every stage.

This is a FLOAT diagnostic, not directed arithmetic. Positive underflow
raises an error; structural zeros and ordinary addition rounding are allowed.
No FFT, subtraction, persistent cache, or output file is used. Run Python -B;
root owns all numerical jobs. An external timeout and memory limit supplement
the internal guards. The CLI emits one JSON object of at most512 KiB.

Cross-check API:
    regions, receipt = regional_table(32, Fraction(49,50), Fraction(1,5))
    regions[(g-2*f, f)]  # normalized2x2 regional matrix in base-two logs
The returned dictionary contains every a+2b<=G, including the empty case.
"""

import argparse
from fractions import Fraction
import hashlib
import json
import math
import os
from pathlib import Path
import shutil
import sys
import time

sys.dont_write_bytecode = True
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
os.environ.setdefault("OMP_NUM_THREADS", "1")

import numpy as np
import probe_small_state_exact_route as exact


bridge = exact.bridge
DEFAULT_Z = Fraction(24, 25)
DEFAULT_ALPHA = Fraction(1, 10)
OUTPUT_LIMIT = 512 * 1024
WEIGHT_SUM_TOLERANCE = 1e-10


class Budget:
    def __init__(self, seconds, memory_mib):
        if not 0 < seconds <= 50 or not 64 <= memory_mib <= 1024:
            raise ValueError("Require <=50 seconds and64..1024 MiB")
        self.free = shutil.disk_usage(Path(__file__).resolve().parent).free
        if self.free < 20*1024**3:
            raise RuntimeError("At least20 GiB free space is required")
        self.started = time.monotonic()
        self.seconds, self.memory_mib = seconds, memory_mib
        self.next_memory_check = self.started

    def guard(self):
        now = time.monotonic()
        if now-self.started >= self.seconds:
            raise TimeoutError("time-cap")
        if now >= self.next_memory_check:
            if bridge._resident_bytes() > self.memory_mib*1024**2:
                raise MemoryError("memory-cap")
            self.next_memory_check = now+.05

    def receipt(self):
        return dict(seconds=time.monotonic()-self.started,
                    peak_resident_bytes=bridge._resident_bytes(),
                    time_cap_seconds=self.seconds,
                    memory_cap_bytes=self.memory_mib*1024**2,
                    free_bytes_at_start=self.free)


def validate_markers(max_groups, z, alpha):
    z, alpha = Fraction(z), Fraction(alpha)
    if (not isinstance(max_groups, int) or not 1 <= max_groups <= 128
            or not 0 < z < 1 or not 0 < alpha <= 1
            or not 0 < float(z) < 1 or not 0 < float(alpha) <= 1):
        raise ValueError("Require G1..128 and representable rational z in(0,1), alpha in(0,1]")
    return z, alpha


def normalizers(slots, max_groups, guard):
    """Exact combinatorial integers, each converted once to binary64."""
    values = np.zeros((max_groups+1, max_groups//2+1), dtype=float)
    integers = {}
    for full in range(min(slots, max_groups//2)+1):
        guard()
        for singles in range(min(slots-full, max_groups-2*full)+1):
            count = (math.comb(slots, singles)*math.comb(slots-singles, full)
                     * 2**singles)
            value = float(count)
            if not 0 < value < math.inf:
                raise ArithmeticError("Normalizer conversion overflowed or underflowed")
            integers[(singles, full)] = count
            values[singles, full] = value
    return values, integers


def _regional_table(max_groups, z, alpha, guard, progress):
    z, alpha = validate_markers(max_groups, z, alpha)
    guard()
    geometry = exact.load_occupancy_counts()
    local = bridge.local_moments(z, q=bridge.Q0, h=4)
    epsilon = Fraction(1, bridge.Q0**16)
    local_values = {}

    def packet(occupancy):
        occupancy = tuple(sorted(occupancy))
        if occupancy in local_values:
            return local_values[occupancy]
        guard()
        quiet, birth, ret, live = bridge.matrix(occupancy, z, local, q=bridge.Q0, h=4)
        total = quiet+birth
        all_live = epsilon*total+(1-epsilon)*live
        entries = ((quiet, total, epsilon*quiet+(1-epsilon)*ret, all_live)
                   if sum(occupancy) else
                   (Fraction(1), Fraction(0), Fraction(0), all_live))
        if min(entries) < 0 or max(entries) > 1:
            raise ArithmeticError("Exact local majorant left[0,1]")
        logs = np.full(4, -math.inf)
        for index, value in enumerate(entries):
            if value:
                converted = float(value)
                if not 0 < converted < math.inf:
                    raise ArithmeticError("Positive local rational conversion underflowed or overflowed")
                logs[index] = math.log2(converted)
        # The same entrywise power and conversion order as exact.regional.
        with np.errstate(over="raise", invalid="raise", under="raise"):
            result = np.exp2(float(alpha)*logs)
        local_values[occupancy] = result
        return result

    shape = (max_groups+1, max_groups//2+1)
    values = np.zeros((4, *shape), dtype=float)
    counts, integer_counts = normalizers(32, max_groups, guard)
    progress.update(max_groups=max_groups, z=str(z), alpha=str(alpha),
                    base_exact_mass_checks=0, stages=[], completed_pair_slots=0)
    with np.errstate(over="raise", invalid="raise", divide="raise", under="raise"):
        for (singles, full), count in integer_counts.items():
            guard()
            # Domain from normalizers guarantees singles+full<=32, including G>32.
            mass = 0
            for occupancy, multiplicity in geometry(singles, full):
                values[:, singles, full] += multiplicity*packet(occupancy)
                mass += multiplicity
            if mass != count:
                raise ArithmeticError("Exact striped one-batch coefficient mass failed")
            values[:, singles, full] /= counts[singles, full]
            progress["base_exact_mass_checks"] += 1
        progress["completed_pair_slots"] = 32

        # Reuse these small buffers for every shifted right slice. Matrix entries
        # occupy four contiguous planes; no tiny matrix products are allocated.
        weights = np.zeros(shape)
        weighted = np.zeros((4, *shape))
        first, second = np.zeros(shape), np.zeros(shape)
        slots = 32
        for _ in range(7):
            guard()
            following_counts, following_integers = normalizers(2*slots, max_groups, guard)
            following = np.zeros_like(values)
            weight_sums = np.zeros(shape)
            for singles, full in integer_counts:
                guard()
                remaining = max_groups-singles-2*full
                rows, columns = remaining+1, remaining//2+1
                right_counts = counts[:rows, :columns]
                target_counts = following_counts[singles:singles+rows, full:full+columns]
                current_weights = weights[:rows, :columns]
                current_weights.fill(0.)
                valid = (right_counts > 0) & (target_counts > 0)
                # For G<=128 all C_n fit binary64. This positive ordering avoids
                # the much larger unnormalized product C_n(left)*C_n(right).
                np.divide(right_counts, target_counts, out=current_weights, where=valid)
                current_weights *= counts[singles, full]
                weight_sums[singles:singles+rows, full:full+columns] += current_weights
                work = weighted[:, :rows, :columns]
                np.multiply(values[:, :rows, :columns], current_weights, out=work)
                aa, ab, ba, bb = values[:, singles, full]
                target = following[:, singles:singles+rows, full:full+columns]
                one, two = first[:rows, :columns], second[:rows, :columns]

                np.multiply(work[0], aa, out=one)
                np.multiply(work[2], ab, out=two)
                np.add(one, two, out=one)
                np.add(target[0], one, out=target[0])
                np.multiply(work[1], aa, out=one)
                np.multiply(work[3], ab, out=two)
                np.add(one, two, out=one)
                np.add(target[1], one, out=target[1])
                np.multiply(work[0], ba, out=one)
                np.multiply(work[2], bb, out=two)
                np.add(one, two, out=one)
                np.add(target[2], one, out=target[2])
                np.multiply(work[1], ba, out=one)
                np.multiply(work[3], bb, out=two)
                np.add(one, two, out=one)
                np.add(target[3], one, out=target[3])

            permitted = following_counts > 0
            error = float(np.max(np.abs(weight_sums[permitted]-1.)))
            if error > WEIGHT_SUM_TOLERANCE:
                raise ArithmeticError("Floating hypergeometric weights failed normalization diagnostic")
            if np.any(~np.isfinite(following)) or np.any(following < 0):
                raise ArithmeticError("Invalid normalized regional coefficient")
            slots *= 2
            progress["stages"].append(dict(pair_slots=slots,
                coefficient_count=len(following_integers),
                maximum_weight_sum_error=error,
                largest_normalizer_bits=max(value.bit_length() for value in following_integers.values())))
            values, counts, integer_counts = following, following_counts, following_integers
            progress["completed_pair_slots"] = slots
        guard()

    # Only the empty-input polynomial coefficient has structural off-diagonal zeros.
    for (singles, full) in integer_counts:
        entry = values[:, singles, full]
        required = entry[[0, 3]] if singles+full == 0 else entry
        if np.any(required <= 0):
            raise ArithmeticError("A required positive regional entry vanished")
    with np.errstate(divide="ignore"):
        logs = np.log2(values)
    regions = {(singles, full): logs[:, singles, full].reshape(2, 2).copy()
               for singles, full in integer_counts}
    progress.update(regional_coefficients=len(regions), local_occupancies=len(local_values),
                    complete_regional_table=True,
                    normalization="exact C32 masses; multinomial Vandermonde identity; floating stage checks",
                    weight_sum_tolerance=WEIGHT_SUM_TOLERANCE,
                    positive_underflow_policy="reject; structural zeros allowed",
                    matrix_order="left coefficient times right coefficient")
    return regions


def regional_table(max_groups, z=DEFAULT_Z, alpha=DEFAULT_ALPHA,
                   seconds=45., memory_mib=1024):
    """Return {(singles,full): normalized2x2 log matrix}, plus a FLOAT receipt.

    Exceptions propagate, so no partial table can be mistaken for completion.
    """
    budget, progress = Budget(seconds, memory_mib), {}
    regions = _regional_table(max_groups, z, alpha, budget.guard, progress)
    progress.update(budget.receipt(), directed=False, whole_code_certificate=False,
                    generated_file_bytes=0, frozen_geometry_sha256=exact.FROZEN_HASH)
    return regions, progress


def log_matrix_product(left, right):
    """Broadcasted ordered2x2 product in base-two logs."""
    return np.logaddexp2(left[..., :, 0, None]+right[..., None, 0, :],
                        left[..., :, 1, None]+right[..., None, 1, :])


def profile_report(group, regions, z, alpha, guard, run_census, details):
    guard()
    regional = np.asarray([regions[(group-2*full, full)] for full in range(group//2+1)])
    run = log_matrix_product(regional[:, None], regional[None, :])
    # Each matrix is already one even/odd pair; seven squares give128 pairs.
    with np.errstate(under="ignore"):
        # logaddexp2 may discard a tiny summand, but retains its finite logarithm
        # as an input. This is ordinary floating addition, not a vanished entry.
        for _ in range(7):
            guard()
            run = log_matrix_product(run, run)
        operators = np.logaddexp2(run[..., 0, 0], run[..., 0, 1])
    full = np.arange(group//2+1)
    induced = full[:, None]+full[None, :]
    allowed = induced < group
    source = bridge.source_factor(group, "cutoff")
    if source["missing_groups"] != 0 or source["retained_fraction"] != str(Fraction(256-group, 256)):
        raise ArithmeticError("Expected the h0 source correction for g<=128")
    source_choices = math.log2(math.comb(8192, group))
    prefix = source_choices+float(alpha)*(source["log2_prefactor"]-bridge.LIMIT*math.log2(float(z)))
    scores = prefix+operators
    if np.any(~np.isfinite(scores[allowed])):
        raise ArithmeticError("Nonfinite complete profile score")
    candidates = np.argwhere(allowed)
    peak = candidates[int(np.argmax(scores[allowed]))]
    even, odd = map(int, peak)
    result = dict(groups=group, z=str(z), alpha=str(alpha),
                  source_mode="cutoff", missing_groups=0,
                  retained_fraction=source["retained_fraction"],
                  checked_profiles=int(np.count_nonzero(allowed)), complete_profile_cover=True,
                  worst_profile=dict(even_full=even, odd_full=odd,
                                     log2_uncapped_bound=float(scores[even, odd])),
                  log2_uncapped_bound=float(scores[even, odd]))
    if run_census:
        terms, total = [], 0
        for runs in range(1, group+1):
            guard()
            count = (Fraction(8192, runs)*math.comb(group-1, runs-1)
                     * math.comb(8192-group-1, runs-1))
            if count.denominator != 1:
                raise ArithmeticError("Nonintegral cycle-run count")
            mask = allowed & (group-induced == runs)
            if not np.any(mask):
                raise ArithmeticError("Missing source-cycle run class")
            value = float(np.max(scores[mask]))+math.log2(count.numerator)-source_choices
            terms.append(dict(runs=runs, source_sets=str(count.numerator), log2_bound=value))
            total += count.numerator
        if total != math.comb(8192, group):
            raise ArithmeticError("Cycle-run counts failed exact partition")
        result["cycle_run_union_log2"] = bridge.log2_sum(row["log2_bound"] for row in terms)
        result["worst_run_term"] = max(terms, key=lambda row: row["log2_bound"])
    if details:
        result["profiles"] = [dict(even_full=int(e), odd_full=int(o),
                                    log2_uncapped_bound=float(scores[e, o]))
                              for e, o in candidates]
    guard()
    return result


def probe(groups=tuple(range(33, 65)), z=DEFAULT_Z, alpha=DEFAULT_ALPHA,
          seconds=45., memory_mib=1024, run_census=False, details=False):
    groups = tuple(sorted(groups))
    if (not groups or len(set(groups)) != len(groups)
            or any(not isinstance(g, int) or not 1 <= g <= 128 for g in groups)):
        raise ValueError("Require distinct actual source counts1..128")
    z, alpha = validate_markers(max(groups), z, alpha)
    budget, progress, reports = Budget(seconds, memory_mib), {}, []
    status, error = "complete-shared-exact-route-float", None
    try:
        regions = _regional_table(max(groups), z, alpha, budget.guard, progress)
        for group in groups:
            reports.append(profile_report(group, regions, z, alpha, budget.guard, run_census, details))
    except (TimeoutError, MemoryError, ArithmeticError) as failure:
        status = str(failure) if isinstance(failure, (TimeoutError, MemoryError)) else "arithmetic-error"
        error = str(failure)
    completed = [report["groups"] for report in reports]
    complete = len(completed) == len(groups)
    return dict(schema="small-state-shared-route-float-v1", status=status, error=error,
                requested_groups=list(groups), completed_groups=completed,
                unprocessed_groups=[group for group in groups if group not in completed],
                complete_source_count_cover=complete, z=str(z), alpha=str(alpha), reports=reports,
                diagnostic_range_union_log2=bridge.log2_sum(
                    report["log2_uncapped_bound"] for report in reports) if complete else None,
                diagnostic_cycle_run_range_union_log2=bridge.log2_sum(
                    report["cycle_run_union_log2"] for report in reports) if complete and run_census else None,
                coefficient_engine=progress, **budget.receipt(),
                output_cap_bytes=OUTPUT_LIMIT,
                source_accounting=("exact cycle-run counts times maximum profile within each run class; "
                                   "run classes and source counts summed; full-binomial control also reported"
                                   if run_census else
                                   "full binomial per profile; maximum over profiles; source counts summed"),
                local_fractional_grouping="entrywise before route averaging",
                terminal="positive f<=M*ones envelope; no final reset or final refresh charge",
                geometry="striped pairs;128 batches;256 alternating regions",
                frozen_geometry_sha256=exact.FROZEN_HASH,
                dependencies={path.name: hashlib.sha256(path.read_bytes()).hexdigest()
                              for path in (Path(__file__).resolve(), Path(exact.__file__), Path(bridge.__file__))},
                directed=False, whole_code_certificate=False, generated_file_bytes=0)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--groups", type=int, nargs="+", default=list(range(33, 65)))
    parser.add_argument("--z", type=Fraction, default=DEFAULT_Z)
    parser.add_argument("--alpha", type=Fraction, default=DEFAULT_ALPHA)
    parser.add_argument("--seconds", type=float, default=45.)
    parser.add_argument("--memory-mib", type=int, default=1024)
    parser.add_argument("--run-census", action="store_true")
    parser.add_argument("--details", action="store_true")
    args = parser.parse_args()
    result = probe(args.groups, args.z, args.alpha, args.seconds, args.memory_mib,
                   args.run_census, args.details)
    encoded = json.dumps(result, allow_nan=False, sort_keys=True)
    if len(encoded.encode("utf-8")) > OUTPUT_LIMIT:
        for report in result["reports"]:
            report.pop("profiles", None)
        result["details_omitted_for_output_cap"] = True
        encoded = json.dumps(result, allow_nan=False, sort_keys=True)
    if len(encoded.encode("utf-8")) > OUTPUT_LIMIT:
        raise RuntimeError("Compact summary exceeded its output cap")
    print(encoded, flush=True)
    if not result["complete_source_count_cover"]:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
