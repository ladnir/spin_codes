"""Bounded FLOAT control for the Boolean full-column source cutoff.

Use the actual striped-pair route and the GF32/GF256 whole-state bridge.
The local state is (column-deficit flag, quiet/live state). If C and D are
the positive clean/dirty 2x2 matrices and M=C+D, fractionalization uses
[[C**alpha,D**alpha],[0,M**alpha]] entrywise, BEFORE route averaging.
Once a column has a deficit, later local tags need not be distinguished.

Four consecutive regions are exactly one RS column: even, odd, even, odd.
Reset the deficit flag at each column, count clean columns, and retain only
counts0..31 across64 columns. A full column supplies4g independent source
constraints;32 full columns force every source to zero. No rank is charged
again: the quiet coefficients already include their uniform-input costs.

The final regular matrix is a positive envelope. For a nonempty final
batch, the next-live column retains its free branch with the SAME tag;
for an empty batch, the diagonal does. No final refresh/reset is required.
The source prefactor is unchanged: binom(8192,g) remains outside alpha,
and [2 beta**g z**(-734003)/(255*(1-g/256))] is raised to alpha.

This is binary64 arithmetic, not a directed certificate. The uncut lifted
control equals the old operator at alpha=1, but can be larger at alpha<1.
No files are written. Run with Python -B; the root owns all numerical runs.
"""

import argparse
from fractions import Fraction
from functools import lru_cache
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
import screen_small_state_bridge as bridge
from probe_small_state_exact_route import (
    FROZEN_HASH, FROZEN_PATH, load_occupancy_counts,
)


Z = Fraction(308022019811679423, 311064346442596352)
ALPHA = Fraction(4, 5)
ZERO = Fraction(0)
ONE = Fraction(1)


def quiet_block_parts(ell, z):
    """Return exact (clean,dirty) sums for reset and free quiet entries."""
    reset = [ZERO, ZERO]
    free = [ZERO, ZERO]
    for requested in range(ell + 1):
        weight = math.comb(ell, requested) * z**(ell-requested) * (1-z)**requested
        rank = requested + min(4, ell-requested)
        reset[int(rank != ell)] += weight / bridge.Q0**rank
        free[int(requested != ell)] += weight / bridge.Q0**requested
    return tuple(reset), tuple(free)


def positive_product_parts(parts):
    """Partition by the first dirty block; never subtract near-equal totals."""
    clean, dirty = ONE, ZERO
    for block_clean, block_dirty in parts:
        dirty = dirty*(block_clean+block_dirty) + clean*block_dirty
        clean *= block_clean
    return clean, dirty


def make_packet_provider(z, alpha, guard):
    """Return cached lifted matrices and exact local validation counters.

    The counters certify identities of the rational partition only. They
    are not floating-point enclosures or a finite-code certificate.
    """
    z, alpha = Fraction(z), Fraction(alpha)
    local = bridge.local_moments(z, q=bridge.Q0, h=4)
    blocks = tuple(quiet_block_parts(ell, z) for ell in range(17))
    epsilon = Fraction(1, bridge.Q0**16)
    statistics = {"exact_packet_partitions": 0, "empty_packets": 0}
    guard()

    @lru_cache(maxsize=4096)
    def packet(occupancy):
        guard()
        occupancy = tuple(occupancy)
        quiet, birth, ret, live = bridge.matrix(occupancy, z, local, q=bridge.Q0, h=4)
        total = quiet + birth
        kc, kd = positive_product_parts(blocks[ell][0] for ell in occupancy)
        tc, td = positive_product_parts(blocks[ell][1] for ell in occupancy)
        if kc+kd != quiet or tc+td != total or min(kc, kd, tc, td, ret, live) < 0:
            raise ArithmeticError("Exact clean/dirty quiet partition failed")
        # These are positive components of the retained provider, including
        # its zero-state mass. No live majorant is subtracted from anything.
        reset = epsilon*quiet + (1-epsilon)*ret
        all_live = epsilon*total + (1-epsilon)*live
        if sum(occupancy):
            clean = (kc, tc, ZERO, ZERO)
            dirty = (kd, td, reset, all_live)
        else:
            # Both forbidden transitions remain zero, even after the flag
            # is already dirty. Empty allowed input has no rank deficit.
            clean = (ONE, ZERO, ZERO, all_live)
            dirty = (ZERO,)*4
            statistics["empty_packets"] += 1
        full = tuple(c+d for c, d in zip(clean, dirty))

        def powered(values):
            result = []
            for value in values:
                converted = float(value)
                if value and converted == 0:
                    raise ArithmeticError("Positive rational local entry underflowed")
                result.append(converted**float(alpha) if value else 0.)
            return np.asarray(result, dtype=float).reshape(2, 2)

        result = np.zeros((4, 4))
        result[:2, :2] = powered(clean)
        result[:2, 2:] = powered(dirty)
        result[2:, 2:] = powered(full)
        if np.any(~np.isfinite(result)) or np.any(result < 0):
            raise ArithmeticError("Invalid lifted local matrix")
        statistics["exact_packet_partitions"] += 1
        return result

    return packet, statistics


def polynomial_product(left, right, singles, full, guard):
    """Ordered positive polynomial product; discard only excessive degrees."""
    guard()
    result = {}
    for (a, b), aa in left.items():
        for (c, d), bb in right.items():
            if a+c <= singles and b+d <= full:
                key = (a+c, b+d)
                value = aa @ bb
                if key in result:
                    result[key] += value
                else:
                    result[key] = value
        guard()
    return result


def regional(packet, singles, full, occupancy_counts, guard):
    """Extract the exact route coefficient; coefficients themselves are floats."""
    base, masses = {}, {}
    for a in range(singles+1):
        for b in range(full+1):
            if a+b > 32:
                continue
            value, mass = np.zeros((4, 4)), 0
            for occupancy, count in occupancy_counts(a, b):
                value += count*packet(tuple(sorted(occupancy)))
                mass += count
            base[(a, b)], masses[(a, b)] = value, mass
        guard()
    value = base
    for _ in range(7):
        value = polynomial_product(value, value, singles, full, guard)
        following = {}
        for (a, b), left in masses.items():
            for (c, d), right in masses.items():
                if a+c <= singles and b+d <= full:
                    key = (a+c, b+d)
                    following[key] = following.get(key, 0)+left*right
        masses = following
        guard()
    normalizer = (math.comb(4096, singles)*math.comb(4096-singles, full)
                  * 2**singles)
    if masses[(singles, full)] != normalizer:
        raise ArithmeticError("Exact striped regional histogram normalization failed")
    result = value[(singles, full)]/normalizer
    if np.any(~np.isfinite(result)) or np.any(result < 0):
        raise ArithmeticError("Invalid floating regional coefficient")
    return result


def log_entries(matrix):
    with np.errstate(divide="ignore"):
        return np.log2(matrix)


def cutoff_log_operator(column, guard, clean_limit=31, columns=64):
    """Positive log-domain recurrence for the truncated clean-column count.

    Setting clean_limit=columns provides an uncut counter control for root
    validation. The default removes exactly histories with32 clean columns.
    """
    if not 0 <= clean_limit <= columns:
        raise ValueError("Invalid clean-column counter limit")
    clean = log_entries(column[:2, :2])
    dirty = log_entries(column[:2, 2:])
    state = np.full((clean_limit+1, 2), -math.inf)
    state[0, 0] = 0.
    for index in range(columns):
        following = np.logaddexp2(state[:, 0, None]+dirty[None, 0, :],
                                 state[:, 1, None]+dirty[None, 1, :])
        if clean_limit:
            increments = np.logaddexp2(state[:-1, 0, None]+clean[None, 0, :],
                                      state[:-1, 1, None]+clean[None, 1, :])
            following[1:] = np.logaddexp2(following[1:], increments)
        state = following
        if index % 8 == 0:
            guard()
    return bridge.log2_sum(state.ravel())


def alpha_one_identity_error(region):
    """Floating validation API: alpha=1 must give A+B=L regionwise."""
    observed, expected = region[:2, :2]+region[:2, 2:], region[2:, 2:]
    if np.any((expected == 0) & (observed != 0)):
        return math.inf
    positive = expected > 0
    return float(np.max(np.abs(observed[positive]-expected[positive])/expected[positive]))


def probe(groups=17, z=Z, alpha=ALPHA, seconds=45., memory_mib=1024,
          profiles=None, details=False):
    if (not isinstance(groups, int) or not 1 <= groups <= 32
            or not 0 < seconds <= 45 or not 64 <= memory_mib <= 1024):
        raise ValueError("Require groups1..32, <=45 seconds and <=1024 MiB")
    z, alpha = Fraction(z), Fraction(alpha)
    if not 0 < z < 1 or not 0 < alpha <= 1:
        raise ValueError("Require rational z in(0,1) and alpha in(0,1]")
    all_profiles = tuple((a, b) for a in range(groups//2+1)
                         for b in range(groups//2+1) if a+b < groups)
    profiles = all_profiles if profiles is None else tuple(dict.fromkeys(profiles))
    if not profiles or any(profile not in all_profiles for profile in profiles):
        raise ValueError("Invalid requested matching profile")
    free = shutil.disk_usage(Path(__file__).resolve().parent).free
    if free < 20*1024**3:
        raise RuntimeError("At least20 GiB free space is required")
    started = time.monotonic()
    last_memory_check = started-1.

    def guard():
        nonlocal last_memory_check
        now = time.monotonic()
        if now-started >= seconds:
            raise TimeoutError("column-cutoff control reached its time cap")
        if now-last_memory_check >= .25:
            if bridge._resident_bytes() > memory_mib*1024**2:
                raise MemoryError("column-cutoff control reached its memory cap")
            last_memory_check = now

    source = bridge.source_factor(groups, "cutoff")
    if source["missing_groups"] != 0 or Fraction(source["retained_fraction"]) != 1-Fraction(groups, 256):
        raise ArithmeticError("Unexpected sparse-source correction")
    prefix = (math.log2(math.comb(8192, groups))
              + float(alpha)*(source["log2_prefactor"]-bridge.LIMIT*math.log2(float(z))))
    regions, records, statistics = {}, [], {}
    max_identity_error = None
    status = "complete-selected-column-cutoff-control"
    try:
        occupancy_counts = load_occupancy_counts()
        packet, statistics = make_packet_provider(z, alpha, guard)
        for full in sorted({value for profile in profiles for value in profile}):
            regions[full] = regional(packet, groups-2*full, full, occupancy_counts, guard)
        if alpha == 1:
            max_identity_error = max(map(alpha_one_identity_error, regions.values()))
            if max_identity_error > 1e-9:
                raise ArithmeticError("Alpha-one regional partition identity failed")
        for even, odd in profiles:
            guard()
            # Region r=4j+a emits component a of RS column j. Its matching
            # parity is a mod2. Each region has128 batches; no flag reset
            # occurs until all four regions of this column are complete.
            pair = regions[even] @ regions[odd]
            column = pair @ pair
            cutoff = cutoff_log_operator(column, guard)
            uncut_matrix = log_entries(column[:2, :2]+column[:2, 2:])
            uncut = bridge.log2_sum(bridge.log_matrix_power(uncut_matrix, 64)[0])
            old_pair = bridge.log_matrix_product(log_entries(regions[even][2:, 2:]),
                                                  log_entries(regions[odd][2:, 2:]))
            old = bridge.log2_sum(bridge.log_matrix_power(old_pair, 128)[0])
            if not all(map(math.isfinite, (cutoff, uncut, old))):
                raise ArithmeticError("Nonfinite floating diagnostic score")
            records.append(dict(even_full=even, odd_full=odd,
                log2_cutoff_operator=cutoff, log2_uncut_lifted_operator=uncut,
                log2_old_untagged_operator=old, log2_uncapped_bound=prefix+cutoff,
                log2_uncut_lifted_bound=prefix+uncut, log2_old_untagged_bound=prefix+old,
                cutoff_gain_bits=uncut-cutoff, fractional_split_cost_bits=uncut-old))
    except (TimeoutError, MemoryError) as error:
        status = str(error)
    complete = len(records) == len(profiles)
    result = dict(schema="small-state-column-cutoff-float-v1", status=status,
        actual_groups=groups, z=str(z), alpha=str(alpha), source_factor=source,
        log2_prefactor_and_source_choices=prefix, origin_log2_penalty=0,
        requested_profiles=len(profiles), evaluated_profiles=len(records),
        complete_selected_profiles=complete,
        complete_profile_cover=complete and set(profiles) == set(all_profiles),
        exact_integer_normalization_checks=len(regions), local_checks=statistics,
        alpha_one_max_regional_relative_error=max_identity_error,
        geometry="128 batches/region; even,odd,even,odd per RS column;64 columns",
        fractional_grouping="entrywise on [[clean,dirty],[0,total]] before route averaging",
        terminal="tag-preserving free branch in final regular matrix; no final refresh charge",
        source_accounting="binom(8192,g) outside alpha; maximum over profiles, no profile union",
        clean_column_counts_retained=[0, 31],
        best_profile=min(records, key=lambda row: row["log2_uncapped_bound"]) if records else None,
        worst_profile=max(records, key=lambda row: row["log2_uncapped_bound"]) if records else None,
        frozen_dependency=FROZEN_PATH.name, frozen_sha256=FROZEN_HASH,
        seconds=time.monotonic()-started, peak_resident_bytes=bridge._resident_bytes(),
        free_bytes_at_start=free, directed=False, whole_code_certificate=False,
        generated_file_bytes=0)
    if details:
        result["profiles"] = records
    return result


def parse_profile(value):
    try:
        profile = tuple(map(int, value.split(",")))
    except ValueError as error:
        raise argparse.ArgumentTypeError("Use even_full,odd_full") from error
    if len(profile) != 2:
        raise argparse.ArgumentTypeError("Use even_full,odd_full")
    return profile


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--groups", type=int, default=17)
    parser.add_argument("--z", type=Fraction, default=Z)
    parser.add_argument("--alpha", type=Fraction, default=ALPHA)
    parser.add_argument("--seconds", type=float, default=45.)
    parser.add_argument("--memory-mib", type=int, default=1024)
    parser.add_argument("--profiles", type=parse_profile, nargs="+")
    parser.add_argument("--details", action="store_true")
    args = parser.parse_args()
    print(json.dumps(probe(args.groups, args.z, args.alpha, args.seconds,
                           args.memory_mib, args.profiles, args.details),
                     allow_nan=False, sort_keys=True), flush=True)


if __name__ == "__main__":
    main()
