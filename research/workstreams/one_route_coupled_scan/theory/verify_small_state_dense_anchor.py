"""Directed complete raw-operator anchor with one rational IID proposal per cell.

Partition every feasible full-pair count at one actual source size into at
most128 cells. A cell's positive rational proposal is fixed from its midpoint
histogram, pseudocount, and tilt. Log histogram probability is concave in the
full-pair count, so its two endpoints supply a common probability lower bound.
Use exact rational LOWER endpoints of Arb probability enclosures for this
divisor. The resulting positive regional matrix bounds every count in the cell.

The striped half-batch law is(p0+p1*(x+y)/2+p2*x*y)^16. Clear its denominator
to obtain an exact integer table H. The two halves have independent proposals.
For sorted occupancies(a,b,c,d), their combined integer numerator is
8*(H_ab*H_cd+H_ac*H_bd+H_ad*H_bc)/product(multiplicity_factorials).
Check divisibility, symmetry, and total integer mass before Arb averaging.
Apply alpha entrywise to exact rational raw local moments BEFORE averaging.

After all cells, take exact rational UPPER endpoints and their entrywise
maximum U. Every actual regional matrix is bounded by U, so U^256 bounds all
even/odd matching profiles, without another profile union. The raw empty
matrix[[1,1],[256^-16,V_empty]] permits inclusion monotonicity at larger source
sizes with the same z and alpha. Sum the joint-Jensen prefactor separately for
every actual source size in the requested interval, with binomial outside alpha.

Accept only after ALL cells and ALL requested source sizes are covered and the
final Arb upper enclosure is strictly below exact2^-target_bits. Partial runs
never pass. No files are written; root owns numerical runs. Use Python -B.
This certifies only the requested source interval, never the whole code.
"""

import argparse
from collections import Counter
from fractions import Fraction
import hashlib
from itertools import combinations_with_replacement
import json
import math
from pathlib import Path
import shutil
import sys
import time

sys.dont_write_bytecode = True
import verify_small_state_exact_route as reference


GROUPS, PAIRS, OUTPUT_LIMIT = 8192, 4096, 128 * 1024


def pair_histogram(groups, full):
    histogram = (PAIRS - groups + full, groups - 2*full, full)
    if min(histogram) < 0 or sum(histogram) != PAIRS:
        raise ArithmeticError("Invalid actual full-pair histogram")
    return histogram


def cells_for_group(groups, maximum):
    low, high = max(0, groups - PAIRS), groups // 2
    width = (high - low + maximum) // maximum
    cells = [(start, min(high, start + width - 1))
             for start in range(low, high + 1, width)]
    if (not 1 <= len(cells) <= maximum or cells[0][0] != low
            or cells[-1][1] != high
            or any(left[1] + 1 != right[0] for left, right in zip(cells, cells[1:]))):
        raise ArithmeticError("Incomplete matching-count cell partition")
    return cells


def proposal(groups, cell, tilt, pseudocount):
    ends = [pair_histogram(groups, full) for full in cell]
    weights = tuple((Fraction(ends[0][i] + ends[1][i], 2) + pseudocount)*tilt**i
                    for i in range(3))
    total = sum(weights)
    probabilities = tuple(weight/total for weight in weights)
    if min(probabilities) <= 0 or sum(probabilities) != 1:
        raise ArithmeticError("Invalid exact rational proposal")
    return probabilities


def integer_occupancy_weights(probabilities, profiles, divisors, guard):
    oriented = (probabilities[0], probabilities[1]/2, probabilities[2])
    denominator = math.lcm(*(value.denominator for value in oriented))
    empty, single, full = (int(value*denominator) for value in oriented)
    if empty + 2*single + full != denominator or min(empty, single, full) <= 0:
        raise ArithmeticError("Failed oriented-pair denominator clearing")
    half = [[1]]
    for _ in range(16):
        guard()
        size = len(half)
        following = [[0]*(size + 1) for _ in range(size + 1)]
        for i in range(size):
            for j in range(size):
                value = half[i][j]
                following[i][j] += empty*value
                following[i + 1][j] += single*value
                following[i][j + 1] += single*value
                following[i + 1][j + 1] += full*value
        half = following
    if (sum(map(sum, half)) != denominator**16
            or any(half[i][j] != half[j][i] for i in range(17) for j in range(17))):
        raise ArithmeticError("Exact half-batch mass or symmetry check failed")
    weights = []
    for index, ((a, b, c, d), divisor) in enumerate(zip(profiles, divisors)):
        if index % 128 == 0:
            guard()
        numerator = 8*(half[a][b]*half[c][d] + half[a][c]*half[b][d]
                       + half[a][d]*half[b][c])
        count, remainder = divmod(numerator, divisor)
        if remainder or count <= 0:
            raise ArithmeticError("Invalid exact sorted occupancy multiplicity")
        weights.append(count)
    mass = denominator**32
    if sum(weights) != mass:
        raise ArithmeticError("Exact complete batch occupancy mass failed")
    return weights, mass


def verify(groups, z, alpha, end_groups=None, max_cells=128,
           tilt=Fraction(3, 2), pseudocount=Fraction(1, 65536), precision=128,
           target_bits=100, seconds=50., memory_mib=1024, details=False):
    end_groups = groups if end_groups is None else end_groups
    z, alpha, tilt, pseudocount = map(Fraction, (z, alpha, tilt, pseudocount))
    if (not isinstance(groups, int) or not isinstance(end_groups, int)
            or not 1 <= groups <= end_groups <= GROUPS
            or not isinstance(max_cells, int) or not 1 <= max_cells <= 128
            or not 0 < z < 1 or not 0 < alpha <= 1 or tilt <= 0 or pseudocount <= 0
            or precision not in (128, 256) or not 1 <= target_bits <= 10000
            or not 0 < seconds <= 50 or not 64 <= memory_mib <= 1024):
        raise ValueError("Require1<=anchor<=end<=8192,1..128 cells, rational0<z<1,"
                         "0<alpha<=1, positive tilt/pseudocount, precision128/256,"
                         "bounded target, <=50 seconds and64..1024 MiB")
    free = shutil.disk_usage(Path(__file__).resolve().parent).free
    if free < 20*1024**3:
        raise RuntimeError("At least20 GiB free space is required")
    started, next_memory = time.monotonic(), 0.

    def guard():
        nonlocal next_memory
        now = time.monotonic()
        if now - started >= seconds:
            raise TimeoutError("time-cap")
        if now >= next_memory:
            if reference.resident_bytes() > memory_mib*1024**2:
                raise MemoryError("memory-cap")
            next_memory = now + .05

    cells = cells_for_group(groups, max_cells)
    status, error = "directed-complete-dense-anchor-replay", None
    completed, completed_sources, reports, result = 0, 0, [], {}
    complete, passed = False, False
    dependencies = {}
    try:
        guard()
        directory = (reference.OLD_THEORY if reference.OLD_THEORY.is_dir()
                     else Path(__file__).resolve().parent)
        sys.path.insert(0, str(directory))
        from hd_common_refresh_packets import local_moments, matrix
        from flint import arb, arb_mat, ctx
        ctx.prec = precision

        def rational(value):
            value = Fraction(value)
            return arb(value.numerator)/value.denominator

        exponent, epsilon = rational(alpha), Fraction(1, 256**16)
        local = local_moments(z, q=256, h=4)
        profiles = tuple(combinations_with_replacement(range(17), 4))
        if len(profiles) != 4845:
            raise ArithmeticError("Incomplete sorted occupancy domain")
        divisors = [math.prod(math.factorial(count) for count in Counter(row).values())
                    for row in profiles]
        rows = []
        for index, occupancy in enumerate(profiles):
            if index % 32 == 0:
                guard()
            quiet, birth, ret, live = matrix(occupancy, z, local, q=256, h=4)
            total = quiet + birth
            entries = (quiet, total, epsilon*quiet + (1 - epsilon)*ret,
                       epsilon*total + (1 - epsilon)*live)
            if min(entries) <= 0 or max(entries) > 1:
                raise ArithmeticError("Exact raw local entry left(0,1]")
            powered = [rational(value) if alpha == 1 else rational(value)**exponent
                       for value in entries]
            if any(not value.is_finite() or not value > 0 for value in powered):
                raise ArithmeticError("Invalid directed raw local entry")
            rows.append(powered)
        local_table = arb_mat(rows)
        guard()
        upper_entries, selected_cells = [None]*4, [None]*4
        log_two = arb(2).log()

        for cell_index, cell in enumerate(cells):
            guard()
            probabilities = proposal(groups, cell, tilt, pseudocount)
            weights, mass = integer_occupancy_weights(probabilities, profiles, divisors, guard)
            average = (arb_mat([weights])*local_table)/arb(mass)
            region = arb_mat([[average[0, 0], average[0, 1]],
                              [average[0, 2], average[0, 3]]])
            for _ in range(7):
                guard()
                region = region*region
            lowers, probability_logs = [], []
            for full in cell:
                guard()
                histogram = pair_histogram(groups, full)
                multinomial = math.comb(PAIRS, histogram[0])*math.comb(
                    PAIRS - histogram[0], histogram[1])
                probability = arb(multinomial)
                for p, n in zip(probabilities, histogram):
                    probability *= rational(p)**n
                if not probability.is_finite() or not probability > 0:
                    raise ArithmeticError("Nonpositive histogram probability enclosure")
                lower = probability.lower().fmpq()
                if lower <= 0:
                    raise ArithmeticError("Nonpositive exact probability lower endpoint")
                lowers.append(lower)
                probability_logs.append(str(probability.log()/log_two))
            divisor = min(lowers)
            region /= arb(divisor)
            for index in range(4):
                value = region[index//2, index%2]
                if not value.is_finite() or not value > 0:
                    raise ArithmeticError("Invalid normalized regional upper matrix")
                upper = value.upper().fmpq()
                if upper_entries[index] is None or upper > upper_entries[index]:
                    upper_entries[index], selected_cells[index] = upper, cell_index
            reports.append(dict(index=cell_index, full_pair_low=cell[0], full_pair_high=cell[1],
                                exact_occupancy_mass_checked=True,
                                proposal=[str(p) for p in probabilities],
                                endpoint_probability_log2_intervals=probability_logs))
            completed += 1
        if completed != len(cells) or any(value is None for value in upper_entries):
            raise ArithmeticError("Incomplete full-pair cell cover")
        guard()
        envelope = arb_mat([[arb(upper_entries[0]), arb(upper_entries[1])],
                            [arb(upper_entries[2]), arb(upper_entries[3])]])
        run = envelope
        for _ in range(8):
            guard()
            run = run*run
        operator = run[0, 0] + run[0, 1]
        if not 2*(16*(reference.LENGTH//64)) < (1 << 32) - 1:
            raise ArithmeticError("The Gamma<2 receipt is unavailable")
        beta = arb(256)**160/arb(255)**32
        kappa = arb(256)**256/beta
        growth = beta*kappa**rational(Fraction(1, 256))
        powered_growth = growth if alpha == 1 else growth**exponent
        source_sum = arb(0)
        term = arb(math.comb(GROUPS, groups))*powered_growth**groups
        for actual in range(groups, end_groups + 1):
            guard()
            source_sum += term
            completed_sources += 1
            if actual < end_groups:
                term *= powered_growth*(GROUPS - actual)
                term /= actual + 1
        common = 2*rational(z)**(-reference.WEIGHT_LIMIT)/255
        source_sum *= common if alpha == 1 else common**exponent
        bound = operator*source_sum
        if not bound.is_finite() or not bound > 0:
            raise ArithmeticError("Invalid final source-interval bound")
        complete = completed_sources == end_groups - groups + 1
        target = arb(2)**(-target_bits)
        passed = complete and bool(bound < target)
        result.update(log2_bound_interval=str(bound.log()/log_two),
                      raw_anchor_operator_log2_interval=str(operator.log()/log_two),
                      interval_source_factor_log2_interval=str(source_sum.log()/log_two),
                      envelope_entry_maximizing_cells=selected_cells,
                      envelope_cells=[dict(index=index, full_pair_low=cells[index][0],
                                           full_pair_high=cells[index][1])
                                      for index in sorted(set(selected_cells))])
        for name in reference.DEPENDENCIES:
            module = sys.modules.get(name)
            if module is None:
                raise RuntimeError("Missing loaded proof dependency: " + name)
            dependencies[name + ".py"] = hashlib.sha256(Path(module.__file__).read_bytes()).hexdigest()
        guard()
    except (TimeoutError, MemoryError, ArithmeticError, RuntimeError) as failure:
        status = str(failure) if isinstance(failure, (TimeoutError, MemoryError)) else "verification-error"
        error = str(failure)
        complete, passed = False, False
    for path in (Path(__file__), Path(reference.__file__)):
        dependencies[path.name] = hashlib.sha256(path.read_bytes()).hexdigest()
    result.update(schema="small-state-dense-anchor-arb-v1", status=status, error=error,
                  passed=passed, complete_profile_cell_cover=complete,
                  complete_source_interval_cover=complete,
                  anchor_groups=groups, end_groups=end_groups,
                  requested_cells=len(cells), completed_cells=completed,
                  requested_source_sizes=end_groups - groups + 1,
                  completed_source_sizes=completed_sources,
                  matching_count_range=[cells[0][0], cells[-1][1]],
                  z=str(z), alpha=str(alpha), tilt=str(tilt), pseudocount=str(pseudocount),
                  precision_bits=precision, target_bits=target_bits,
                  acceptance="complete coverage and Arb bound < exact2^(-target_bits)",
                  source_mode="joint", gamma_upper="2",
                  source_prefactor="sum_g binom(8192,g)*(2*beta^g*kappa^(g/256)*z^(-734003)/255)^alpha",
                  local_fractional_grouping="entrywise before exact striped IID average",
                  cell_denominator="minimum of exact rational LOWER probability endpoints",
                  regional_envelope="entrywise maximum of exact rational UPPER regional endpoints",
                  profile_accounting="all actual regional matrices <= common U; U^256; no profile union",
                  empty_input_rule="raw[[1,1],[256^-16,V_empty]]; omega=1",
                  terminal="positive f<=M*ones; no physical final refresh or reset charge",
                  interval_rule="raw exact-moment inclusion monotonicity; constant markers",
                  seconds=time.monotonic() - started, peak_resident_bytes=reference.resident_bytes(),
                  time_cap_seconds=seconds, memory_cap_bytes=memory_mib*1024**2,
                  output_cap_bytes=OUTPUT_LIMIT, free_bytes_at_start=free,
                  dependencies=dependencies, directed_arithmetic=True, directed=complete,
                  whole_code_certificate=False, generated_file_bytes=0)
    if details:
        result["cells"] = reports
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--groups", type=int, required=True)
    parser.add_argument("--z", type=Fraction, required=True)
    parser.add_argument("--alpha", type=Fraction, required=True)
    parser.add_argument("--end-groups", type=int)
    parser.add_argument("--max-cells", type=int, default=128)
    parser.add_argument("--tilt", type=Fraction, default=Fraction(3, 2))
    parser.add_argument("--pseudocount", type=Fraction, default=Fraction(1, 65536))
    parser.add_argument("--precision", type=int, choices=(128, 256), default=128)
    parser.add_argument("--target-bits", type=int, default=100)
    parser.add_argument("--seconds", type=float, default=50.)
    parser.add_argument("--memory-mib", type=int, default=1024)
    parser.add_argument("--details", action="store_true")
    args = parser.parse_args()
    result = verify(args.groups, args.z, args.alpha, args.end_groups, args.max_cells,
                    args.tilt, args.pseudocount, args.precision, args.target_bits,
                    args.seconds, args.memory_mib, args.details)
    encoded = json.dumps(result, allow_nan=False, sort_keys=True)
    if len(encoded.encode("utf-8")) > OUTPUT_LIMIT:
        result.pop("cells", None)
        result["details_omitted_for_output_cap"] = True
        encoded = json.dumps(result, allow_nan=False, sort_keys=True)
    if len(encoded.encode("utf-8")) > OUTPUT_LIMIT:
        raise RuntimeError("Compact output cap exceeded")
    print(encoded, flush=True)
    if not result["passed"]:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
