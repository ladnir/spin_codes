"""Bounded density-free gate for two through four actual source groups.

The regional average is exact combinatorially: singleton and full pairs
occupy distinct uniformly sampled slots among 4096 pairs, 32 per packet.
The matrix cap ignores the fair orientations because it depends only on
the number of allowed input coordinates. The default FLOAT arithmetic is
diagnostic. --verify replays fixed rational marks with directed Arb bounds
for the selected source-count class, not the whole code.

Keep the pure-quiet path separate throughout every polynomial product.
Nonquiet zero-to-zero excursions must not receive the quiet-column rank
discount. This helper never benchmarks, optimizes, or writes result files.
"""

import argparse
import json
import math
import os
import shutil
import time
from fractions import Fraction

import numpy as np


STATE_COUNT = 17
GROUP_COUNT = 8192
OUTPUT_LENGTH = 1 << 21


def require_free_space():
    """Honor the workstream guard even though this helper writes no data."""
    if shutil.disk_usage(os.path.dirname(os.path.abspath(__file__))).free < 20 * (1 << 30):
        raise RuntimeError("at least 20 GiB free is required before a numerical run")


def packet_matrices(tau, z, max_inputs, activation_d1=False):
    """Caps indexed by the number of distinct allowed input coordinates.

    At most four input columns are independent under the block-diagonal
    RS[16,4] feedback map. Thus a zero-to-zero state transition forces all
    these inputs to vanish. For nonzero incoming state support h, emission
    has at least max(0,13*ceil(h/4)-ell) nonzero coordinates. The MDS refresh
    gives h+h' >=17-4*ell whenever both states are nonzero.
    """
    matrices = []
    for ell in range(max_inputs + 1):
        matrix = np.zeros((STATE_COUNT, STATE_COUNT))
        matrix[0, 0] = 1.0
        for h in range(STATE_COUNT):
            emitted = max(0, 13 * ((h + 3) // 4) - ell)
            for hh in range(STATE_COUNT):
                if h == hh == 0:
                    continue
                multiplicity = math.comb(16, hh)
                if h == 0:
                    if ell and hh >= max(1, 17 - 4 * ell):
                        activation = (min(1.0, ((256 ** ell - 1) // 255)
                                          * 255.0 ** (hh - 16))
                                      if activation_d1 else 1.0)
                        matrix[h, hh] = multiplicity * z * activation
                elif hh == 0:
                    if ell and h <= 4 * ell:
                        matrix[h, hh] = tau ** (-h) * z ** emitted
                elif h + hh >= 17 - 4 * ell:
                    matrix[h, hh] = (multiplicity
                        * tau ** (-min(16 - hh, h)) * z ** emitted)
        matrices.append(matrix)
    return matrices


def pair_product(left, right):
    """Multiply G=quiet*P+remainder, with P=e_0 e_0^T, positively."""
    q1, r1 = left
    q2, r2 = right
    remainder = r1 @ r2
    remainder[:, 0] += q2 * r1[:, 0]
    remainder[0, :] += q1 * r2[0, :]
    return q1 * q2, remainder


def polynomial_product(left, right, singles, full):
    result = {}
    for (a, b), value1 in left.items():
        for (c, d), value2 in right.items():
            key = a + c, b + d
            if key[0] > singles or key[1] > full:
                continue
            q, remainder = pair_product(value1, value2)
            if key in result:
                old_q, old_remainder = result[key]
                old_remainder += remainder
                result[key] = old_q + q, old_remainder
            else:
                result[key] = q, remainder
    return result


def region(singles, full, packets):
    """Return quiet/nonquiet parts of the averaged 128-packet operator."""
    polynomial = {}
    for a in range(singles + 1):
        for b in range(full + 1):
            coefficient = math.comb(32, a) * math.comb(32 - a, b)
            remainder = coefficient * packets[a + 2 * b].copy()
            remainder[0, 0] = 0.0
            polynomial[a, b] = coefficient, remainder
    # 128 is a power of two. All polynomial products remain truncated to
    # the requested histogram and preserve matrix multiplication order.
    for _ in range(7):
        polynomial = polynomial_product(polynomial, polynomial, singles, full)
    quiet, remainder = polynomial[singles, full]
    denominator = (math.comb(4096, singles)
                   * math.comb(4096 - singles, full))
    assert quiet == denominator
    remainder /= denominator
    if not np.all(np.isfinite(remainder)) or np.any(remainder < 0):
        raise ArithmeticError("regional FLOAT range failure")
    return 1, remainder


def profile_gate(groups, even_full, odd_full, regions, tau, z, distance):
    """One matching profile is fixed across all64 outer columns."""
    even, odd = regions[even_full], regions[odd_full]
    column = pair_product(pair_product(pair_product(even, odd), even), odd)
    assert column[0] == 1
    remainder = column[1]
    quiet_weight = tau ** (-4 * groups)
    vector = np.zeros((32, STATE_COUNT))
    vector[0, 0] = 1.0
    logscale = 0.0
    for _ in range(64):
        following = vector @ remainder
        following[1:, 0] += vector[:-1, 0] * quiet_weight
        vector = following
        scale = float(np.max(vector))
        if not math.isfinite(scale) or scale <= 0:
            raise ArithmeticError("column FLOAT range failure")
        vector /= scale
        logscale += math.log2(scale)
    dimensions = 128 * groups
    corrections = np.array([
        1 - tau ** (-(dimensions - 4 * groups * quiet))
        for quiet in range(32)
    ])
    moment = float(np.sum(vector * corrections[:, None]))
    if not math.isfinite(moment) or moment <= 0:
        raise ArithmeticError("terminal FLOAT range failure")
    threshold = math.floor(distance * OUTPUT_LENGTH)
    logmoment = logscale + math.log2(moment)
    log_bound = (dimensions * math.log2(tau) - math.log2(tau - 1)
                 - threshold * math.log2(z) + logmoment)
    return {"even_full": even_full, "odd_full": odd_full,
            "log2_single_support_bound": log_bound,
            "log2_moment": logmoment}


def gate(groups, tau, z, distance, activation_d1=False):
    packets = packet_matrices(tau, z, groups, activation_d1)
    regions = {f: region(groups - 2 * f, f, packets)
               for f in range(groups // 2 + 1)}
    profiles = [profile_gate(groups, a, b, regions, tau, z, distance)
                for a in regions for b in regions if a + b < groups]
    worst = max(profiles, key=lambda item: item["log2_single_support_bound"])
    union_bound = (math.log2(math.comb(GROUP_COUNT, groups))
                   + worst["log2_single_support_bound"])
    return {"method": "density-free-few-groups", "groups": groups,
            "tau": tau, "z": z, "distance": distance,
            "threshold": math.floor(distance * OUTPUT_LENGTH),
            "log2_union_bound": union_bound,
            "worst_profile": [worst["even_full"], worst["odd_full"]],
            "profiles": profiles, "states": STATE_COUNT,
            "quiet_count_maximum": 31,
            "activation_d1": activation_d1,
            "whole_code_certificate": False,
            "status": "FLOAT diagnostic; no certificate"}


def quiet_region_polynomial(even, odd):
    """Column coefficient q counts exactly q completely quiet regions."""
    coefficients = [np.eye(STATE_COUNT)]
    for region_pair in (even, odd, even, odd):
        assert region_pair[0] == 1
        remainder = region_pair[1]
        following = [np.zeros_like(remainder) for _ in range(len(coefficients) + 1)]
        for quiet, coefficient in enumerate(coefficients):
            following[quiet] += coefficient @ remainder
            following[quiet + 1][:, 0] += coefficient[:, 0]
        coefficients = following
    return coefficients


def logsumexp2(values):
    peak = max(values)
    return peak + math.log2(sum(2 ** (x - peak) for x in values))


def top32_gate(groups, tau, z_values, distance, activation_d1=False):
    """Four threshold bounds for the sum of the largest32 quiet counts.

    A path with q_j quiet regions in outer column j has deterministic
    source rank at least r*sum_largest32(q_j). The latter sum equals
    min_v(32*v+sum_j max(q_j-v,0)). A nonzero witness admits a minimizing
    v in0..3; choosing v4 only is possible at full source rank. Partition
    by a minimizer, then enlarge each part before summing its bound.
    """
    threshold = math.floor(distance * OUTPUT_LENGTH)
    profile_values = {(a, b): [] for a in range(groups // 2 + 1)
                      for b in range(groups // 2 + 1) if a + b < groups}
    cache = {}
    for v, z in enumerate(z_values):
        if z not in cache:
            packets = packet_matrices(tau, z, groups, activation_d1)
            cache[z] = {f: region(groups - 2 * f, f, packets)
                        for f in range(groups // 2 + 1)}
        regions = cache[z]
        for (a, b), values in profile_values.items():
            coefficients = quiet_region_polynomial(regions[a], regions[b])
            rank_limit = 128 - 32 * v
            state = np.zeros((rank_limit, STATE_COUNT))
            state[0, 0] = 1.0
            logscale = 0.0
            for _ in range(64):
                following = np.zeros_like(state)
                for q, coefficient in enumerate(coefficients):
                    cost = max(q - v, 0)
                    if cost:
                        following[cost:] += ((state[:-cost] @ coefficient)
                                             * tau ** (-groups * cost))
                    else:
                        following += state @ coefficient
                state = following
                scale = float(np.max(state))
                if not math.isfinite(scale) or scale <= 0:
                    raise ArithmeticError("top32 FLOAT range failure")
                state /= scale
                logscale += math.log2(scale)
            correction = np.array([1 - tau ** (-groups * (rank_limit - cost))
                                   for cost in range(rank_limit)])
            logmoment = logscale + math.log2(float(np.sum(state * correction[:, None])))
            bound = ((128 * groups - 32 * groups * v) * math.log2(tau)
                     - math.log2(tau - 1) - threshold * math.log2(z)
                     + logmoment)
            values.append(bound)
    profiles = [dict(even_full=a, odd_full=b, threshold_bounds=values,
                     log2_single_support_bound=logsumexp2(values))
                for (a, b), values in profile_values.items()]
    worst = max(profiles, key=lambda item: item["log2_single_support_bound"])
    support_count = math.log2(math.comb(GROUP_COUNT, groups))
    return {"method": "density-free-few-groups-top32", "groups": groups,
            "tau": tau, "z_values": z_values, "distance": distance,
            "threshold": threshold,
            "log2_union_bound": support_count + worst["log2_single_support_bound"],
            "threshold_worst_union_bounds": [support_count + max(
                item["threshold_bounds"][v] for item in profiles) for v in range(4)],
            "worst_profile": [worst["even_full"], worst["odd_full"]],
            "profiles": profiles, "states": STATE_COUNT,
            "activation_d1": activation_d1,
            "whole_code_certificate": False,
            "status": "FLOAT diagnostic; no certificate"}


def verify_top32(groups, tau_text, z_texts, distance_text, precision,
                 target_bits, activation_d1):
    """Directed replay of fixed rational marks; no proposal search."""
    from flint import arb, arb_mat, ctx
    require_free_space()
    ctx.prec = precision
    tau_r = Fraction(tau_text)
    z_rs = [Fraction(z) for z in z_texts]
    distance_r = Fraction(distance_text)
    if not (1 <= groups <= 4 and 1 < tau_r <= 65280
            and len(z_rs) == 4 and all(0 < z < 1 for z in z_rs)
            and 0 < distance_r < Fraction(1, 2)
            and 64 <= precision <= 512 and 1 <= target_bits <= 10000):
        raise ValueError("invalid rational replay parameters")

    def value(rational):
        return arb(rational.numerator) / rational.denominator

    tau = value(tau_r)
    n = STATE_COUNT

    def packets(z):
        result = []
        for ell in range(groups + 1):
            matrix = arb_mat(n, n)
            matrix[0, 0] = 1
            for h in range(n):
                emitted = max(0, 13 * ((h + 3) // 4) - ell)
                for hh in range(n):
                    if h == hh == 0:
                        continue
                    multiplicity = math.comb(16, hh)
                    if h == 0:
                        if ell and hh >= max(1, 17 - 4 * ell):
                            activation = (min(Fraction(1), Fraction(
                                (256 ** ell - 1) // 255, 255 ** (16 - hh)))
                                if activation_d1 else Fraction(1))
                            matrix[h, hh] = multiplicity * z * value(activation)
                    elif hh == 0:
                        if ell and h <= 4 * ell:
                            matrix[h, hh] = tau ** (-h) * z ** emitted
                    elif h + hh >= 17 - 4 * ell:
                        matrix[h, hh] = (multiplicity
                            * tau ** (-min(16 - hh, h)) * z ** emitted)
            result.append(matrix)
        return result

    def pair_mul(left, right):
        q1, r1 = left
        q2, r2 = right
        remainder = r1 * r2
        for j in range(n):
            remainder[j, 0] += q2 * r1[j, 0]
            remainder[0, j] += q1 * r2[0, j]
        return q1 * q2, remainder

    def poly_mul(left, right, singles, full):
        result = {}
        for (a, b), val1 in left.items():
            for (c, d), val2 in right.items():
                key = a + c, b + d
                if key[0] > singles or key[1] > full:
                    continue
                q, remainder = pair_mul(val1, val2)
                if key in result:
                    old_q, old_r = result[key]
                    result[key] = old_q + q, old_r + remainder
                else:
                    result[key] = q, remainder
        return result

    def regional(singles, full, packet_caps):
        polynomial = {}
        for a in range(singles + 1):
            for b in range(full + 1):
                coefficient = math.comb(32, a) * math.comb(32 - a, b)
                remainder = coefficient * packet_caps[a + 2 * b]
                remainder[0, 0] = 0
                polynomial[a, b] = coefficient, remainder
        for _ in range(7):
            polynomial = poly_mul(polynomial, polynomial, singles, full)
        quiet, remainder = polynomial[singles, full]
        denominator = (math.comb(4096, singles)
                       * math.comb(4096 - singles, full))
        assert quiet == denominator
        return remainder / denominator

    def column_poly(even, odd):
        identity = arb_mat(n, n)
        for j in range(n):
            identity[j, j] = 1
        coefficients = [identity]
        for remainder in (even, odd, even, odd):
            following = [arb_mat(n, n) for _ in range(len(coefficients) + 1)]
            for quiet, coefficient in enumerate(coefficients):
                following[quiet] += coefficient * remainder
                for j in range(n):
                    following[quiet + 1][j, 0] += coefficient[j, 0]
            coefficients = following
        return coefficients

    profile_terms = {(a, b): [] for a in range(groups // 2 + 1)
                     for b in range(groups // 2 + 1) if a + b < groups}
    cache = {}
    threshold = int(distance_r * OUTPUT_LENGTH)
    for v, z_r in enumerate(z_rs):
        z = value(z_r)
        if z_r not in cache:
            packet_caps = packets(z)
            cache[z_r] = {f: regional(groups - 2 * f, f, packet_caps)
                          for f in range(groups // 2 + 1)}
        regions = cache[z_r]
        rank_limit = 128 - 32 * v
        prefix = (math.comb(GROUP_COUNT, groups) * tau ** (groups * rank_limit)
                  / (tau - 1) * z ** (-threshold))
        for (a, b), terms in profile_terms.items():
            coefficients = column_poly(regions[a], regions[b])
            weighted = {}
            for q, coefficient in enumerate(coefficients):
                cost = max(q - v, 0)
                scaled = coefficient * tau ** (-groups * cost)
                weighted[cost] = weighted.get(cost, arb_mat(n, n)) + scaled
            state = arb_mat(rank_limit, n)
            state[0, 0] = 1
            for _ in range(64):
                following = state * weighted[0]
                for cost, matrix in weighted.items():
                    if cost == 0:
                        continue
                    contribution = state * matrix
                    for c in range(cost, rank_limit):
                        for j in range(n):
                            following[c, j] += contribution[c - cost, j]
                state = following
            moment = sum((1 - tau ** (-groups * (rank_limit - c))) * state[c, j]
                         for c in range(rank_limit) for j in range(n))
            terms.append(prefix * moment)
    limit = arb(2) ** (-target_bits)
    profiles = []
    passed = True
    for (a, b), terms in profile_terms.items():
        bound = sum(terms)
        profile_passed = bool(bound.is_finite() and bound < limit)
        passed = passed and profile_passed
        profiles.append({"even_full": a, "odd_full": b,
                         "passed": profile_passed,
                         "log2_union_bound_interval": str(bound.log() / arb(2).log()),
                         "threshold_log2_intervals": [str(x.log() / arb(2).log())
                                                      for x in terms]})
    return {"method": "density-free-few-groups-top32-directed", "groups": groups,
            "tau": str(tau_r), "z_values": [str(z) for z in z_rs],
            "distance": str(distance_r), "threshold": threshold,
            "precision": precision, "target_bits": target_bits,
            "passed": passed, "profiles": profiles,
            "whole_code_certificate": False,
            "activation_d1": activation_d1,
            "scope": "all nonzero GF128 messages in exactly r of8192 source groups; "
                     "one fixed matching profile per support; ideal setup",
            "states": n}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--groups", type=int, nargs="+", default=[2, 3, 4])
    parser.add_argument("--tau", nargs="+", default=["65280"])
    parser.add_argument("--z", nargs="+", default=["0.999"])
    parser.add_argument("--distance", default="0.35")
    parser.add_argument("--rank-mode", choices=["columns", "top32"], default="columns")
    parser.add_argument("--threshold-z", nargs=4)
    parser.add_argument("--activation-d1", action="store_true")
    parser.add_argument("--verify", action="store_true")
    parser.add_argument("--precision", type=int, default=192)
    parser.add_argument("--target-bits", type=int, default=40)
    args = parser.parse_args()
    if (any(r not in (1, 2, 3, 4) for r in args.groups)
            or any(not 1 < Fraction(t) <= 65280 for t in args.tau)
            or any(not 0 < Fraction(z) < 1 for z in args.z)
            or (args.threshold_z is not None and
                (args.rank_mode != "top32" or len(args.z) != 1
                 or any(not 0 < Fraction(z) < 1 for z in args.threshold_z)))
            or not 0 < Fraction(args.distance) < Fraction(1, 2)):
        parser.error("require r in1..4,1<tau<=65280,0<z<1,0<distance<0.5")
    if len(args.groups) * len(args.tau) * len(args.z) > 24:
        parser.error("at most24 bounded parameter choices")
    if args.verify and (args.rank_mode != "top32" or len(args.groups) != 1
                        or len(args.tau) != 1 or len(args.z) != 1
                        or not 64 <= args.precision <= 512
                        or not 1 <= args.target_bits <= 10000):
        parser.error("replay one top32 group/parameter choice at64..512 bits per invocation")
    require_free_space()
    started = time.monotonic()
    for groups in args.groups:
        for tau in args.tau:
            for z in args.z:
                if time.monotonic() - started > 50:
                    raise SystemExit("bounded50-second invocation budget reached")
                if args.verify:
                    result = verify_top32(groups, tau, args.threshold_z or [z] * 4,
                                          args.distance, args.precision,
                                          args.target_bits, args.activation_d1)
                elif args.rank_mode == "top32":
                    result = top32_gate(groups, float(Fraction(tau)),
                        [float(Fraction(x)) for x in (args.threshold_z or [z] * 4)],
                        float(Fraction(args.distance)), args.activation_d1)
                else:
                    result = gate(groups, float(Fraction(tau)), float(Fraction(z)),
                                  float(Fraction(args.distance)), args.activation_d1)
                print(json.dumps(result,
                                 sort_keys=True), flush=True)
                if args.verify and not result["passed"]:
                    raise SystemExit(1)


if __name__ == "__main__":
    main()
