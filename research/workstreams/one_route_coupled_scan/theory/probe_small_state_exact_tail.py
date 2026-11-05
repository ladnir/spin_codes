"""Bounded FLOAT exact regional coefficients near the fully paired endpoint.

For even actual source size g, return every R(s=2j,f=g/2-j),0<=j<=tail.
Retain singleton degrees0..2*tail and total input degree s+2f<=g. Six positive
matrix-polynomial squares take32 pair slots to2048; the final square extracts
only the requested target coefficients in4096 slots. Divide by exact integer
C4096(s,f)=binom(4096,s)binom(4096-s,f)2^s, using its logarithm.

Every scalar coefficient remains in log2 form. Hot convolution loops reuse a
shift buffer and accumulate directly into explicit output planes. Matrix
squaring uses the five-product identity over the commutative scalar polynomial
ring; scalar squares visit only unordered singleton-degree splits. No FFT,
subtractions, or conversions of positive coefficients through exp2 are used.
Ordinary log-sum rounding can absorb a small addend. "Exact" describes the
coefficient being extracted, NOT directed arithmetic or a distance certificate.

The default is g2048,tail8 with raw empty matrices. Smaller g4/g8 controls can
be compared with the shared route engine by the root. No tests run implicitly.
No files are written. Root owns every numerical run; use Python -B.
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
import probe_small_state_cluster_route as cluster


exact, bridge, Budget = cluster.exact, cluster.bridge, cluster.Budget
OUTPUT_LIMIT = 128*1024


def validate(groups, tail):
    if (not isinstance(groups, int) or groups % 2 or not 2 <= groups <= 4096
            or not isinstance(tail, int) or not 0 <= tail <= min(8, groups//2)):
        raise ValueError("Require even groups2..4096 and0<=tail<=min(8,groups/2)")


def last_degree(singles, slots, full_cap):
    """Nonnegative degrees cannot enter a retained target after this cutoff."""
    return min(slots - singles, full_cap - (singles + 1)//2)


def check_polynomial(value, slots, full_cap, max_singles):
    if value.shape != (4, max_singles + 1, min(slots, full_cap) + 1):
        raise ArithmeticError("Unexpected exact-tail coefficient shape")
    if np.any(np.isnan(value)) or np.any(np.isposinf(value)):
        raise ArithmeticError("Invalid positive-log coefficient")
    for singles in range(max_singles + 1):
        allowed = max(0, last_degree(singles, slots, full_cap) + 1)
        if np.any(~np.isneginf(value[:, singles, allowed:])):
            raise ArithmeticError("Coefficient outside the pair-slot or input-degree domain")


def base_polynomial(groups, tail, z, alpha, raw, guard, progress):
    guard()
    full_cap, max_singles = groups//2, 2*tail
    geometry = exact.load_occupancy_counts()
    local = bridge.local_moments(z, q=bridge.Q0, h=4)
    epsilon = Fraction(1, bridge.Q0**16)
    packets = {}
    result = np.full((4, max_singles + 1, min(32, full_cap) + 1), -math.inf)
    progress.update(base_exact_mass_checks=0, local_occupancies=0)

    def packet(occupancy):
        occupancy = tuple(sorted(occupancy))
        if occupancy not in packets:
            guard()
            quiet, birth, ret, live = bridge.matrix(occupancy, z, local, q=bridge.Q0, h=4)
            total = quiet + birth
            all_live = epsilon*total + (1 - epsilon)*live
            entries = ((quiet, total, epsilon*quiet + (1 - epsilon)*ret, all_live)
                       if sum(occupancy) or raw else
                       (Fraction(1), Fraction(0), Fraction(0), all_live))
            if min(entries) < 0 or max(entries) > 1:
                raise ArithmeticError("Exact local majorant left[0,1]")
            logs = np.full(4, -math.inf)
            for index, value in enumerate(entries):
                if value:
                    logs[index] = float(alpha)*(math.log2(value.numerator)
                                                - math.log2(value.denominator))
                    if not math.isfinite(logs[index]):
                        raise ArithmeticError("A positive local logarithm vanished")
            packets[occupancy] = logs
        return packets[occupancy]

    with np.errstate(over="raise", invalid="raise", under="ignore"):
        for singles in range(max_singles + 1):
            for full in range(max(0, last_degree(singles, 32, full_cap) + 1)):
                guard()
                mass = 0
                for occupancy, count in geometry(singles, full):
                    if count <= 0:
                        raise ArithmeticError("Nonpositive striped multiplicity")
                    term = packet(occupancy) + math.log2(count)
                    np.logaddexp2(result[:, singles, full], term,
                                 out=result[:, singles, full])
                    mass += count
                expected = math.comb(32, singles)*math.comb(32 - singles, full)*2**singles
                if mass != expected:
                    raise ArithmeticError("Exact one-batch coefficient mass failed")
                progress["base_exact_mass_checks"] += 1
    progress["local_occupancies"] = len(packets)
    check_polynomial(result, 32, full_cap, max_singles)
    return result


def convolve_into(left, right, output, width, work, guard, progress, twice=False):
    """Accumulate a positive scalar convolution into one existing output plane."""
    guard()
    terms, scale = 0, 1. if twice else 0.
    for index in range(min(len(left), width)):
        if index % 64 == 0:
            guard()
        value = left[index]
        if value == -math.inf:
            continue
        count = min(len(right), width - index)
        if count <= 0:
            continue
        np.add(right[:count], value + scale, out=work[:count])
        np.logaddexp2(output[index:index + count], work[:count],
                     out=output[index:index + count])
        terms += count
    progress["scalar_convolutions"] += 1
    progress["convolution_pair_terms"] += terms


def ring_product_into(output, left, right, slots, full_cap, work, guard, progress):
    for singles in range(output.shape[0]):
        width = max(0, min(output.shape[1], last_degree(singles, 2*slots, full_cap) + 1))
        if not width:
            continue
        for first in range(singles + 1):
            second = singles - first
            left_size = max(0, min(left.shape[1], last_degree(first, slots, full_cap) + 1))
            right_size = max(0, min(right.shape[1], last_degree(second, slots, full_cap) + 1))
            if left_size and right_size:
                convolve_into(left[first, :left_size], right[second, :right_size],
                              output[singles], width, work, guard, progress)


def ring_square_into(output, value, slots, full_cap, work, guard, progress):
    for singles in range(output.shape[0]):
        width = max(0, min(output.shape[1], last_degree(singles, 2*slots, full_cap) + 1))
        if not width:
            continue
        for first in range(singles//2 + 1):
            second = singles - first
            left_size = max(0, min(value.shape[1], last_degree(first, slots, full_cap) + 1))
            right_size = max(0, min(value.shape[1], last_degree(second, slots, full_cap) + 1))
            if left_size and right_size:
                convolve_into(value[first, :left_size], value[second, :right_size],
                              output[singles], width, work, guard, progress, first != second)


def matrix_square(value, slots, full_cap, guard, progress):
    """Five scalar-ring products, preserving the ordered matrix square."""
    width, degrees = min(2*slots, full_cap) + 1, value.shape[1]
    work = np.empty(width)
    following = np.full((4, degrees, width), -math.inf)
    cross = np.full((degrees, width), -math.inf)
    aa, ab, ba, bb = value
    diagonal_sum = np.logaddexp2(aa, bb)
    ring_product_into(cross, ab, ba, slots, full_cap, work, guard, progress)
    ring_square_into(following[0], aa, slots, full_cap, work, guard, progress)
    np.logaddexp2(following[0], cross, out=following[0])
    ring_product_into(following[1], ab, diagonal_sum, slots, full_cap, work, guard, progress)
    ring_product_into(following[2], ba, diagonal_sum, slots, full_cap, work, guard, progress)
    ring_square_into(following[3], bb, slots, full_cap, work, guard, progress)
    np.logaddexp2(following[3], cross, out=following[3])
    check_polynomial(following, 2*slots, full_cap, degrees - 1)
    return following


def extract_coefficients(base, groups, tail, guard, progress, *, base_slots=32,
                         intermediate_squares=6):
    """Supplied-base fixture interface; singleton/full weights use the C_n convention."""
    validate(groups, tail)
    if (not isinstance(base_slots, int) or not 1 <= base_slots <= 32
            or not isinstance(intermediate_squares, int) or not 0 <= intermediate_squares <= 6):
        raise ValueError("Invalid bounded coefficient-extraction fixture")
    final_slots = base_slots << (intermediate_squares + 1)
    full_cap, max_singles = groups//2, 2*tail
    if any(2*j + full_cap - j > final_slots for j in range(tail + 1)):
        raise ValueError("Requested target exceeds the fixture's pair slots")
    slots, value = base_slots, np.asarray(base, dtype=float)
    check_polynomial(value, slots, full_cap, max_singles)
    progress.update(base_pair_slots=base_slots, intermediate_squares=intermediate_squares,
                    stages=[], scalar_convolutions=0, convolution_pair_terms=0,
                    completed_pair_slots=slots, extracted_coefficients=[],
                    complete_regional_coefficients=False)
    with np.errstate(over="raise", invalid="raise", under="ignore"):
        for _ in range(intermediate_squares):
            guard()
            value = matrix_square(value, slots, full_cap, guard, progress)
            slots *= 2
            progress["completed_pair_slots"] = slots
            progress["stages"].append(dict(pair_slots=slots, singleton_degree_cap=max_singles,
                                            full_degree_cap=full_cap,
                                            coefficient_array_bytes=value.nbytes))
        regions, normalizers = {}, []
        for offset in range(tail + 1):
            guard()
            singles, full = 2*offset, full_cap - offset
            coefficient = cluster.target_coefficient(value, slots, singles, full, guard)
            count = (math.comb(2*slots, singles)*math.comb(2*slots - singles, full)
                     * 2**singles)
            regions[full] = coefficient - math.log2(count)
            normalizers.append(dict(singles=singles, full_pairs=full,
                                    exact_integer_bits=count.bit_length(), log2=math.log2(count)))
            progress["extracted_coefficients"].append([singles, full])
        guard()
    progress.update(completed_pair_slots=2*slots, regional_pair_slots=2*slots,
                    complete_regional_coefficients=True, normalizers=normalizers,
                    normalization="supplied base-slot C_n masses; Vandermonde; exact integer target C_n logarithms",
                    final_square="requested target coefficients only; ordered matrix products",
                    arithmetic="positive log2 convolution; no FFT, subtraction, or exp2",
                    allocation="explicit coefficient planes; reused shift buffer; direct output accumulation")
    return regions


def _regional_coefficients(groups, z, alpha, tail, raw, guard, progress):
    validate(groups, tail)
    base = base_polynomial(groups, tail, z, alpha, raw, guard, progress)
    regions = extract_coefficients(base, groups, tail, guard, progress)
    if (len(regions) != tail + 1
            or any(np.any(~np.isfinite(value)) for value in regions.values())):
        progress["complete_regional_coefficients"] = False
        raise ArithmeticError("A requested positive regional entry vanished")
    return regions


def regional_coefficients(groups=2048, z=None, alpha=None, *, tail=8, raw=True,
                          seconds=45., memory_mib=1024):
    """Return{full_pair_count:2x2log matrix},receipt; incomplete calls raise."""
    z, alpha = cluster.markers(z, alpha)
    budget, progress = Budget(seconds, memory_mib), {}
    regions = _regional_coefficients(groups, z, alpha, tail, raw, budget.guard, progress)
    progress.update(**budget.receipt(), actual_groups=groups, tail=tail,
                    raw_empty_input=bool(raw), z=str(z), alpha=str(alpha),
                    directed=False, whole_code_certificate=False, generated_file_bytes=0)
    return regions, progress


def probe(groups=2048, z=None, alpha=None, *, tail=8, raw=True, seconds=45., memory_mib=1024):
    validate(groups, tail)
    z, alpha = cluster.markers(z, alpha)
    budget, progress, regions = Budget(seconds, memory_mib), {}, {}
    status, error = "complete-exact-tail-float", None
    try:
        regions = _regional_coefficients(groups, z, alpha, tail, raw, budget.guard, progress)
    except (TimeoutError, MemoryError, ArithmeticError, RuntimeError) as failure:
        status = str(failure) if isinstance(failure, (TimeoutError, MemoryError)) else "arithmetic-error"
        error = str(failure)
    complete = error is None and len(regions) == tail + 1
    result = dict(schema="small-state-exact-tail-float-v1", status=status, error=error,
                  actual_groups=groups, tail=tail, z=str(z), alpha=str(alpha),
                  raw_empty_input=bool(raw), complete_requested_coefficient_cover=complete,
                  complete_source_size_cover=False,
                  requested_full_pairs=list(range(groups//2 - tail, groups//2 + 1)),
                  completed_full_pairs=sorted(regions),
                  regions=[dict(full_pairs=full, singles=groups - 2*full,
                                log2_matrix=value.tolist()) for full, value in sorted(regions.items())],
                  coefficient_engine=progress, **budget.receipt(),
                  frozen_geometry_sha256=exact.FROZEN_HASH,
                  local_fractional_grouping="entrywise before exact striped coefficient averaging",
                  terminal="positive f<=M*ones envelope; no final reset or refresh charge",
                  output_cap_bytes=OUTPUT_LIMIT, directed=False, whole_code_certificate=False,
                  generated_file_bytes=0)
    paths = {Path(__file__).resolve(), Path(cluster.__file__).resolve(), Path(exact.__file__).resolve(),
             Path(bridge.__file__).resolve(), Path(sys.modules[Budget.__module__].__file__).resolve()}
    provider = sys.modules.get("hd_common_refresh_packets")
    if provider is not None:
        paths.add(Path(provider.__file__).resolve())
    result["dependencies"] = {path.name: hashlib.sha256(path.read_bytes()).hexdigest()
                              for path in sorted(paths)}
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--groups", type=int, default=2048)
    parser.add_argument("--z", type=Fraction, required=True)
    parser.add_argument("--alpha", type=Fraction, required=True)
    parser.add_argument("--tail", type=int, default=8)
    parser.add_argument("--pruned", action="store_true")
    parser.add_argument("--seconds", type=float, default=45.)
    parser.add_argument("--memory-mib", type=int, default=1024)
    args = parser.parse_args()
    result = probe(args.groups, args.z, args.alpha, tail=args.tail, raw=not args.pruned,
                   seconds=args.seconds, memory_mib=args.memory_mib)
    encoded = json.dumps(result, allow_nan=False, sort_keys=True)
    if len(encoded.encode("utf-8")) > OUTPUT_LIMIT:
        raise RuntimeError("Compact summary exceeded its output cap")
    print(encoded, flush=True)
    if not result["complete_requested_coefficient_cover"]:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
