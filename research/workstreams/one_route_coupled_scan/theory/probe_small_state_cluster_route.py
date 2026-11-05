"""Targeted FLOAT exact-route coefficients for the g4096 one-run class.

Retain singleton degrees0,1,2 and full-pair degree<=2048. Six positive
matrix-polynomial squares cover2048 pair slots. The final square extracts
only degrees(0,2048) and(2,2047), then divides by their exact C4096 masses.
Evaluate both ordered matching profiles(2048,2047) and(2047,2048).

Every coefficient is stored as a base-two logarithm. Positive products use
addition and positive sums use logaddexp2; there is no FFT or subtraction.
Ordinary log-sum rounding can absorb a small summand, but no positive entry
is converted through a possibly underflowing exponential. Structural zeros
are -inf. "Exact" refers to coefficient extraction, NOT directed arithmetic.

Local moments, striped geometry, entrywise alpha, joint source correction,
and the free-terminal envelope are unchanged. Raw empty-input matrices are
the default, matching the dense-anchor control; --pruned selects actual-empty
pruning. A result covers only the one-run class, never all g4096 source sets.

Cross-check API: regional_coefficients(full_pairs, z, alpha, raw=True) returns
{(0,f): log_matrix, (2,f-1): log_matrix} and a receipt, for1<=f<=2048.
The internal extract_coefficients accepts a supplied base log polynomial,
base_slots, and intermediate_squares for smaller or constant-local controls.
No files are written. Root owns all runs; use Python -B and external limits.
"""

import argparse
from fractions import Fraction
import hashlib
import json
import math
import os
from pathlib import Path
import sys

sys.dont_write_bytecode = True
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
os.environ.setdefault("OMP_NUM_THREADS", "1")

import numpy as np
import probe_small_state_exact_route as exact
from probe_small_state_shared_route import Budget


bridge = exact.bridge
OUTPUT_LIMIT = 128 * 1024


def markers(z, alpha):
    if z is None or alpha is None:
        raise ValueError("Explicit rational z and alpha are required")
    z, alpha = Fraction(z), Fraction(alpha)
    if (not 0 < z < 1 or not 0 < alpha <= 1
            or not 0 < float(z) < 1 or not 0 < float(alpha) <= 1):
        raise ValueError("Require representable rational0<z<1 and0<alpha<=1")
    return z, alpha


def check_polynomial(value, slots, full_pairs):
    if value.shape != (4, 3, min(slots, full_pairs)+1):
        raise ArithmeticError("Unexpected bounded-singleton polynomial shape")
    if np.any(np.isnan(value)) or np.any(np.isposinf(value)):
        raise ArithmeticError("Invalid positive log coefficient")
    for singles in range(3):
        allowed = max(0, min(slots-singles, full_pairs)+1)
        if np.any(~np.isneginf(value[:, singles, allowed:])):
            raise ArithmeticError("Coefficient outside the pair-slot domain")


def base_polynomial(full_pairs, z, alpha, raw, guard, progress):
    """Exact striped multiplicities times the unchanged powered local caps."""
    guard()
    occupancy_counts = exact.load_occupancy_counts()
    local = bridge.local_moments(z, q=bridge.Q0, h=4)
    epsilon = Fraction(1, bridge.Q0**16)
    packets = {}
    value = np.full((4, 3, min(32, full_pairs)+1), -math.inf)
    progress.update(base_exact_mass_checks=0, local_occupancies=0)

    def packet(occupancy):
        occupancy = tuple(sorted(occupancy))
        if occupancy not in packets:
            guard()
            quiet, birth, ret, live = bridge.matrix(occupancy, z, local, q=bridge.Q0, h=4)
            total = quiet+birth
            all_live = epsilon*total+(1-epsilon)*live
            if sum(occupancy) or raw:
                entries = (quiet, total, epsilon*quiet+(1-epsilon)*ret, all_live)
            else:
                entries = (Fraction(1), Fraction(0), Fraction(0), all_live)
            if min(entries) < 0 or max(entries) > 1:
                raise ArithmeticError("Exact rational local majorant left[0,1]")
            logs = np.full(4, -math.inf)
            for index, entry in enumerate(entries):
                if entry:
                    # Do not first convert a tiny rational to binary64.
                    logs[index] = float(alpha)*(math.log2(entry.numerator)
                                                - math.log2(entry.denominator))
                    if not math.isfinite(logs[index]):
                        raise ArithmeticError("A positive local logarithm vanished")
            packets[occupancy] = logs
        return packets[occupancy]

    with np.errstate(over="raise", invalid="raise", under="ignore"):
        for singles in range(3):
            for full in range(min(32-singles, full_pairs)+1):
                guard()
                mass = 0
                for occupancy, count in occupancy_counts(singles, full):
                    if count <= 0:
                        raise ArithmeticError("Nonpositive striped multiplicity")
                    term = packet(occupancy)+math.log2(count)
                    np.logaddexp2(value[:, singles, full], term,
                                 out=value[:, singles, full])
                    mass += count
                expected = math.comb(32, singles)*math.comb(32-singles, full)*2**singles
                if mass != expected:
                    raise ArithmeticError("Exact one-batch coefficient mass failed")
                progress["base_exact_mass_checks"] += 1
    progress["local_occupancies"] = len(packets)
    check_polynomial(value, 32, full_pairs)
    return value


def log_convolve(left, right, width, work, guard, progress):
    """Positive scalar convolution, with one reusable shifted-slice buffer."""
    guard()
    result = np.full(width, -math.inf)
    terms = 0
    for index in range(min(len(left), width)):
        if (index & 63) == 0:
            guard()
        if left[index] == -math.inf:
            continue
        count = min(len(right), width-index)
        if not count:
            continue
        np.add(right[:count], left[index], out=work[:count])
        np.logaddexp2(result[index:index+count], work[:count],
                     out=result[index:index+count])
        terms += count
    progress["scalar_convolutions"] += 1
    progress["convolution_pair_terms"] += terms
    return result


def ring_product(left, right, slots, width, work, guard, progress):
    """Multiply scalar polynomials modulo singleton degree3 (six convolutions)."""
    result = np.full((3, width), -math.inf)
    for singles in range(3):
        for first in range(singles+1):
            second = singles-first
            lsize = max(0, min(slots-first+1, left.shape[1]))
            rsize = max(0, min(slots-second+1, right.shape[1]))
            term = log_convolve(left[first, :lsize], right[second, :rsize],
                                width, work, guard, progress)
            np.logaddexp2(result[singles], term, out=result[singles])
    return result


def ring_square(value, slots, width, work, guard, progress):
    """A scalar square needs four convolutions by commutativity."""
    rows = [value[a, :max(0, min(slots-a+1, value.shape[1]))] for a in range(3)]
    result = np.full((3, width), -math.inf)
    result[0] = log_convolve(rows[0], rows[0], width, work, guard, progress)
    result[1] = log_convolve(rows[0], rows[1], width, work, guard, progress)+1.
    result[2] = log_convolve(rows[0], rows[2], width, work, guard, progress)+1.
    term = log_convolve(rows[1], rows[1], width, work, guard, progress)
    np.logaddexp2(result[2], term, out=result[2])
    return result


def matrix_square(value, slots, full_pairs, guard, progress):
    """Ordered matrix square using five positive scalar ring products."""
    width = min(2*slots, full_pairs)+1
    work = np.empty(width)
    aa, ab, ba, bb = value
    diagonal_sum = np.logaddexp2(aa, bb)
    cross = ring_product(ab, ba, slots, width, work, guard, progress)
    following = np.empty((4, 3, width))
    following[0] = ring_square(aa, slots, width, work, guard, progress)
    np.logaddexp2(following[0], cross, out=following[0])
    following[1] = ring_product(ab, diagonal_sum, slots, width, work, guard, progress)
    following[2] = ring_product(ba, diagonal_sum, slots, width, work, guard, progress)
    following[3] = ring_square(bb, slots, width, work, guard, progress)
    np.logaddexp2(following[3], cross, out=following[3])
    check_polynomial(following, 2*slots, full_pairs)
    return following


def target_coefficient(value, slots, singles, full, guard):
    """Extract one coefficient of the final ordered matrix square."""
    result = np.full(4, -math.inf)
    for first in range(singles+1):
        guard()
        second = singles-first
        left_last = min(slots-first, value.shape[2]-1)
        right_last = min(slots-second, value.shape[2]-1)
        low, high = max(0, full-right_last), min(full, left_last)
        if low > high:
            continue
        left = value[:, first, low:high+1]
        right = value[:, second, full-high:full-low+1][:, ::-1]
        for output, x0, y0, x1, y1 in (
                (0, 0, 0, 1, 2), (1, 0, 1, 1, 3),
                (2, 2, 0, 3, 2), (3, 2, 1, 3, 3)):
            terms = np.logaddexp2(left[x0]+right[y0], left[x1]+right[y1])
            result[output] = np.logaddexp2(result[output], np.logaddexp2.reduce(terms))
    if np.any(np.isnan(result)) or np.any(np.isposinf(result)):
        raise ArithmeticError("Invalid final target coefficient")
    return result.reshape(2, 2)


def extract_coefficients(base, full_pairs, guard, progress, *, base_slots=32,
                         intermediate_squares=6):
    """Fixture interface: supplied(4,3,min(base_slots,f)+1) log coefficients.

    Return normalized log matrices for(0,f) and(2,f-1). The supplied base
    must use the same labeled-slot/orientation counting convention as C_n.
    Structural zeros remain allowed, including for constant-matrix controls.
    """
    if (not isinstance(base_slots, int) or not 1 <= base_slots <= 32
            or not isinstance(intermediate_squares, int) or not 0 <= intermediate_squares <= 6
            or not isinstance(full_pairs, int)
            or not 1 <= full_pairs <= min(2048, (base_slots << (intermediate_squares+1))-1)):
        raise ValueError("Invalid bounded coefficient-extraction fixture")
    slots, value = base_slots, np.asarray(base, dtype=float)
    check_polynomial(value, slots, full_pairs)
    progress.update(base_pair_slots=base_slots, stages=[], scalar_convolutions=0, convolution_pair_terms=0,
                    completed_pair_slots=slots, intermediate_squares=intermediate_squares,
                    extracted_coefficients=[], complete_regional_coefficients=False)
    with np.errstate(over="raise", invalid="raise", under="ignore"):
        # logaddexp2 absorption is ordinary addition rounding. No exp2 is used.
        for _ in range(intermediate_squares):
            guard()
            value = matrix_square(value, slots, full_pairs, guard, progress)
            slots *= 2
            progress["completed_pair_slots"] = slots
            progress["stages"].append(dict(pair_slots=slots, full_degree_cap=full_pairs,
                                            coefficient_array_bytes=value.nbytes))
        regions, normalizers = {}, []
        for singles, full in ((0, full_pairs), (2, full_pairs-1)):
            guard()
            coefficient = target_coefficient(value, slots, singles, full, guard)
            count = (math.comb(2*slots, singles)*math.comb(2*slots-singles, full)
                     * 2**singles)
            regions[(singles, full)] = coefficient-math.log2(count)
            normalizers.append(dict(singles=singles, full_pairs=full,
                                    exact_integer_bits=count.bit_length(), log2=math.log2(count)))
            progress["extracted_coefficients"].append([singles, full])
        guard()
    progress.update(completed_pair_slots=2*slots, regional_pair_slots=2*slots,
                    complete_regional_coefficients=True,
                    normalizers=normalizers,
                    normalization="supplied base C_n convention; multinomial Vandermonde identity; exact integer target C_n logs",
                    final_square="two target coefficients only; ordered matrix products",
                    arithmetic="positive log2 convolution; no FFT, subtraction, or exponentiation")
    return regions


def _regional_coefficients(full_pairs, z, alpha, raw, guard, progress):
    if not isinstance(full_pairs, int) or not 1 <= full_pairs <= 2048:
        raise ValueError("Require target full_pairs1..2048")
    base = base_polynomial(full_pairs, z, alpha, raw, guard, progress)
    regions = extract_coefficients(base, full_pairs, guard, progress)
    if any(np.any(~np.isfinite(region)) for region in regions.values()):
        progress["complete_regional_coefficients"] = False
        raise ArithmeticError("A required positive regional entry vanished")
    return regions


def regional_coefficients(full_pairs, z, alpha, *, raw=True, seconds=45., memory_mib=1024):
    """Return the two normalized log matrices and a FLOAT receipt; errors propagate."""
    z, alpha = markers(z, alpha)
    budget, progress = Budget(seconds, memory_mib), {}
    regions = _regional_coefficients(full_pairs, z, alpha, raw, budget.guard, progress)
    progress.update(**budget.receipt(), raw_empty_input=bool(raw), z=str(z), alpha=str(alpha),
                    directed=False, whole_code_certificate=False, generated_file_bytes=0)
    return regions, progress


def probe(z, alpha, *, raw=True, seconds=45., memory_mib=1024):
    z, alpha = markers(z, alpha)
    budget, progress, profiles = Budget(seconds, memory_mib), {}, []
    requested = [[0, 2048], [2, 2047]]
    status, error, result = "complete-targeted-cluster-float", None, {}
    try:
        regions = _regional_coefficients(2048, z, alpha, raw, budget.guard, progress)
        result["regions"] = [dict(singles=a, full_pairs=b, log2_matrix=value.tolist())
                             for (a, b), value in regions.items()]
        source = bridge.source_factor(4096, "joint")
        prefix = float(alpha)*(source["log2_prefactor"]-bridge.LIMIT*math.log2(float(z)))
        source_choices = math.log2(math.comb(8192, 4096))
        for even, odd in ((2048, 2047), (2047, 2048)):
            budget.guard()
            left, right = regions[(4096-2*even, even)], regions[(4096-2*odd, odd)]
            with np.errstate(over="raise", invalid="raise", under="ignore"):
                pair = bridge.log_matrix_product(left, right)
                run = bridge.log_matrix_power(pair, 128)
                operator = float(np.logaddexp2(run[0, 0], run[0, 1]))
            conditional = prefix+operator
            if not math.isfinite(conditional):
                raise ArithmeticError("Nonfinite clustered profile score")
            profiles.append(dict(even_full=even, odd_full=odd, log2_operator=operator,
                                 log2_conditional_bound=conditional,
                                 log2_full_binomial_bound=source_choices+conditional))
        worst = max(profiles, key=lambda row: row["log2_conditional_bound"])
        # A proper half-cycle with one occupied run has8192 possible starts.
        # Both matching orders occur; a maximum charges the class only once.
        result.update(source_factor=source, log2_powered_source_prefactor=prefix,
                      log2_full_source_choices=source_choices, one_run_source_sets=8192,
                      worst_profile=worst,
                      log2_one_run_bound=13.+worst["log2_conditional_bound"],
                      log2_full_binomial_control=worst["log2_full_binomial_bound"])
        budget.guard()
    except (TimeoutError, MemoryError, ArithmeticError, RuntimeError) as failure:
        status = str(failure) if isinstance(failure, (TimeoutError, MemoryError)) else "arithmetic-error"
        error = str(failure)
    completed = progress.get("extracted_coefficients", [])
    complete = error is None and len(completed) == 2 and len(profiles) == 2
    result.update(schema="small-state-cluster-route-float-v1", status=status, error=error,
                  actual_groups=4096, z=str(z), alpha=str(alpha), raw_empty_input=bool(raw),
                  requested_coefficients=requested, completed_coefficients=completed,
                  unprocessed_coefficients=[key for key in requested if key not in completed],
                  requested_profiles=2, checked_profiles=len(profiles), profiles=profiles,
                  unprocessed_profiles=[[e, o] for e, o in ((2048, 2047), (2047, 2048))
                                        if not any(p["even_full"] == e and p["odd_full"] == o for p in profiles)],
                  complete_targeted_profile_cover=complete, complete_source_size_cover=False,
                  coefficient_engine=progress, source_mode="joint",
                  source_accounting="8192 one-run source sets times maximum of both orders; count outside alpha",
                  local_fractional_grouping="entrywise before exact striped coefficient averaging",
                  empty_input_rule="raw" if raw else "actual-empty pruned",
                  terminal="positive f<=M*ones envelope; no final reset or refresh charge",
                  frozen_geometry_sha256=exact.FROZEN_HASH, **budget.receipt(),
                  output_cap_bytes=OUTPUT_LIMIT, directed=False, whole_code_certificate=False,
                  generated_file_bytes=0)
    paths = {Path(__file__).resolve(), Path(exact.__file__).resolve(), Path(bridge.__file__).resolve(),
             Path(sys.modules[Budget.__module__].__file__).resolve()}
    provider = sys.modules.get("hd_common_refresh_packets")
    if provider is not None:
        paths.add(Path(provider.__file__).resolve())
    result["dependencies"] = {path.name: hashlib.sha256(path.read_bytes()).hexdigest()
                              for path in sorted(paths)}
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--z", type=Fraction, required=True)
    parser.add_argument("--alpha", type=Fraction, required=True)
    parser.add_argument("--pruned", action="store_true", help="use actual-empty pruning instead of raw entries")
    parser.add_argument("--seconds", type=float, default=45.)
    parser.add_argument("--memory-mib", type=int, default=1024)
    args = parser.parse_args()
    result = probe(args.z, args.alpha, raw=not args.pruned,
                   seconds=args.seconds, memory_mib=args.memory_mib)
    encoded = json.dumps(result, allow_nan=False, sort_keys=True)
    if len(encoded.encode("utf-8")) > OUTPUT_LIMIT:
        raise RuntimeError("Compact summary exceeded its output cap")
    print(encoded, flush=True)
    if not result["complete_targeted_profile_cover"]:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
