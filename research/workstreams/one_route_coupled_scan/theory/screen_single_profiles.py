#!/usr/bin/env python3
"""Bounded floating support-profile diagnostic; not a finite certificate.

The default geometry is K=2^20, B=64, packet=4, batch=128, q=2^32.
Only the named, absolute per-batch (h,w,overlap) sequences are covered.
There is no union over other profiles, translations, or excursion placements.
Negative scores concern this relaxed class only; positive scores are failed
bounds, not codewords. The root agent owns execution. Use python -B.

Integer local and selected-time ranks follow SINGLE_DIAGONAL_RANK.md. For a
fixed profile, taking the minimum of the rank/amortized/reset caps dominates
their convex interpolations. SHUFFLED_PACKET_ROUTE.md supplies the actual
packet-power moment and the separate H/S mask groupings. The S grouping uses
the complementary outer-compatibility argument in AFFINE_COMPATIBILITY.md;
its rank is not added to an unrestricted inner rank. A separate S_overlap
bound uses only state-zero pivots on overlap coordinates; compatibility
rows vanish on those columns, so these particular ranks do add.
Regional scores additionally retain the full/touched-region shortening cap.
The old scores remain visible, using only total support size and D(U).

At most 19 profiles, 16384 batches, and 10 fixed fractional markers are used.
The selected-time recurrence uses two prefix-maximum trees: O(n log n) time,
O(n) storage. No state-space search, external dependencies, or files written.
The 30-second default guard aborts an unfinished screen rather than reporting
partial coverage as complete. Intended memory is well below 256 MiB.
"""

import argparse
from bisect import bisect_right
from collections import Counter
from dataclasses import dataclass
from functools import lru_cache
from math import comb, exp2, isfinite, log2
from time import monotonic


B = 64
PACKET = 4
BATCH = 128
FIELD_BITS = 32
LOG_PROJECTIVE = log2((1 << FIELD_BITS) - 1)
LOG_TAU_GRID = (0.25, 0.5, 1.0, 2.0, 4.0, 8.0, 12.0, 16.0, 24.0, 32.0)


@dataclass(frozen=True)
class Profile:
    name: str
    period: tuple
    activation: int = 128
    mode: str = "live"


PROFILES = (
    Profile("partial86_w44_o44", ((86, 44, 44),)),
    Profile("partial86_w44_o22", ((86, 44, 22),)),
    Profile("partial86_w44_o2", ((86, 44, 2),)),
    Profile("partial96_w44_o44", ((96, 44, 44),)),
    Profile("partial65_w44_o44", ((65, 44, 44),)),
    Profile("fullH_w44", ((128, 44, 44),)),
    Profile("above35_partial86_w45", ((86, 45, 45),)),
    Profile("burst3_h1", ((128, 128, 128), (128, 0, 0), (1, 0, 0))),
    Profile("burst3_h32", ((128, 128, 128), (128, 0, 0), (32, 0, 0))),
    Profile("burst3_h64", ((128, 128, 128), (128, 0, 0), (64, 0, 0))),
    Profile("burst3_fullH", ((128, 128, 128), (128, 0, 0), (128, 0, 0))),
    Profile("burst4_h1", ((128, 128, 128), (128, 0, 0), (128, 0, 0), (1, 0, 0))),
    Profile("reset3_a1", ((0, 1, 0), (128, 128, 128), (0, 0, 0)), mode="repeat"),
    Profile("reset3_a32", ((0, 32, 0), (97, 97, 97), (0, 0, 0)), mode="repeat"),
    Profile("reset3_a64", ((0, 64, 0), (65, 65, 65), (0, 0, 0)), mode="repeat"),
    Profile("reset4_a44", ((0, 44, 0), (128, 0, 0), (85, 85, 85), (0, 0, 0)), mode="repeat"),
    Profile("single_impulse", ((128, 0, 0),), activation=1),
    Profile("late_only", (), activation=44, mode="late"),
    Profile("phased_sparse_then_dense", (), mode="phased"),
)


def require(condition, message):
    if not condition:
        raise AssertionError(message)


def guard(deadline):
    if monotonic() > deadline:
        raise TimeoutError("time cap reached; screen incomplete, no coverage claim")


def expand(profile, batches):
    require(len(profile.period) <= 8, "period cap")
    if profile.mode == "late":
        return [(0, 0, 0)] * (batches - 1) + [(0, profile.activation, 0)]
    if profile.mode == "phased":
        # One full activation, then p sparse batches, then m full batches.
        # Choose the largest m with n+127*m+127 <= floor(.35*128*n).
        # This defines one absolute schedule at each accepted n, not a sweep.
        cutoff = (35 * BATCH * batches) // 100
        dense = (cutoff - batches - (BATCH - 1)) // (BATCH - 1)
        require(0 <= dense <= batches - 1, "invalid phased schedule length")
        sparse = batches - 1 - dense
        return [(0, BATCH, 0)] + [(65, 1, 1)] * sparse + [(BATCH, BATCH, BATCH)] * dense
    if profile.mode == "repeat":
        return [profile.period[t % len(profile.period)] for t in range(batches)]
    return [(0, profile.activation, 0)] + [
        profile.period[(t - 1) % len(profile.period)] for t in range(1, batches)]


class PrefixMaximum:
    """Integer Fenwick tree. A query of count0 contains no entries."""

    def __init__(self, count):
        self.data = [-10**12] * (count + 1)

    def add(self, index, value):
        index += 1
        while index < len(self.data):
            self.data[index] = max(self.data[index], value)
            index += index & -index

    def query(self, count):
        result = -10**12
        while count:
            result = max(result, self.data[count])
            count -= count & -count
        return result


def selected_rank(sequence, deadline, source_index=1):
    """Equation(13), using outputs or overlaps; no free-terminal observation."""
    prefix = [0]
    for triple in sequence:
        prefix.append(prefix[-1] + triple[source_index])
    coordinates = sorted(set(prefix[:-1]))
    positions = {value: i for i, value in enumerate(coordinates)}
    count = len(coordinates)
    left, right = PrefixMaximum(count), PrefixMaximum(count)
    left.add(positions[0], 0)
    right.add(count - 1 - positions[0], 0)
    best = 0
    for t in range(1, len(sequence)):
        if t % 2048 == 0:
            guard(deadline)
        capacity = BATCH - sequence[t][0]
        split = bisect_right(coordinates, prefix[t] - capacity)
        # Split min(capacity, prefix[t]-prefix[s]) at the saturation point.
        best = max(left.query(split) + capacity,
                   right.query(count - split) + prefix[t])
        index = positions[prefix[t]]
        left.add(index, best)
        right.add(count - 1 - index, best - prefix[t])
    return best


def check_selected_rank(deadline):
    """Bounded integer regression against the direct quadratic recurrence."""
    seed = 0x51A61E
    checks = 0
    for fixture in range(32):
        sequence = []
        for t in range(1 + fixture % 24):
            seed = (1664525 * seed + 1013904223) & 0xFFFFFFFF
            h = (seed >> 8) % (BATCH + 1)
            seed = (1664525 * seed + 1013904223) & 0xFFFFFFFF
            w = (seed >> 8) % (BATCH + 1)
            if fixture % 4 == 0:
                w = 0 if t % 3 else BATCH
            if fixture % 4 == 1:
                h = 0 if t % 3 else BATCH
            low, high = max(0, h + w - BATCH), min(h, w)
            overlap = low + (seed % (high - low + 1))
            sequence.append((h, w, overlap))
        for source_index in (1, 2):
            prefix = [0]
            for triple in sequence:
                prefix.append(prefix[-1] + triple[source_index])
            direct = [0] * len(sequence)
            for t in range(1, len(sequence)):
                direct[t] = max(direct[s] + min(BATCH - sequence[t][0],
                                               prefix[t] - prefix[s])
                                for s in range(t))
            fast = selected_rank(sequence, deadline, source_index)
            require(fast == direct[-1],
                    f"selected-rank regression fixture={fixture}, source={source_index}")
            checks += 1
    guard(deadline)
    return checks


def rank_caps(sequence, deadline):
    """Uniform integer amplitude caps, plus necessary support exclusions."""
    live = False
    starts = closed = active = excursion_weight = zeros = 0
    previous = (0, 0, 0)
    old_ell = old_shared = 0
    activity = []
    for t, (h, w, overlap) in enumerate(sequence):
        require(0 <= h <= BATCH and 0 <= w <= BATCH, "weight outside batch")
        require(max(0, h + w - BATCH) <= overlap <= min(h, w), "invalid overlap")
        if h == 0:
            if live:
                ph, pw, po = previous
                if not (ph == pw == po):
                    return {"excluded": "reset masks require H=W"}
                if excursion_weight <= BATCH:
                    return {"excluded": "closed excursion has at most b outputs"}
                closed += 1
                live = False
            if w:
                starts += 1
                excursion_weight = 0
                live = True
        elif not live:
            return {"excluded": "live inverse support without activation"}
        if h and old_ell + h <= BATCH:
            return {"excluded": "adjacent MDS support condition"}
        if h and not w:
            zeros += BATCH - h
            if zeros >= BATCH:
                return {"excluded": "zero-output run exceeds polynomial root budget"}
        else:
            zeros = 0
        active += int(live)
        activity.append(live)
        if live:
            excursion_weight += w
        old_shared += min(overlap, max(0, old_ell + h - BATCH))
        old_ell = h + w - overlap
        previous = (h, w, overlap)

    w_total = sum(w for _, w, _ in sequence)
    h_total = sum(h for h, _, _ in sequence)
    local = sum(min(BATCH - sequence[t][0], sequence[t - 1][1])
                for t in range(1, len(sequence)))
    selected = selected_rank(sequence, deadline)
    overlap_rank = selected_rank(sequence, deadline, source_index=2)
    require(selected >= local, "selected-time rank lost the all-times partition")
    require(overlap_rank <= selected, "overlap columns exceed all-output rank bound")
    a = BATCH * (active - starts)
    caps = {
        "ell": min(w_total, old_shared + (old_ell if live else 0)),
        "fresh": w_total - local,
        "selected": w_total - selected,
        "amortized": h_total + (w_total - a + (BATCH if live else 0)) // 2,
        "reset": w_total - BATCH * closed,
    }
    return {"caps": caps, "E": min(caps.values()), "W": w_total, "H": h_total,
            "active": active, "starts": starts, "closed": closed, "free": live,
            "activity": activity, "R_fresh": local, "R_selected": selected,
            "R_overlap": overlap_rank}


def shortening_cap(length, rows, selected, full, touched):
    """Sum of (row occupancy-length/2)+, bounded by its secant on [f,t]."""
    require(length > 0 and length % 2 == 0 and rows > 0, "invalid row geometry")
    require(0 <= full <= touched <= length, "invalid full/touched counts")
    require(full * rows <= selected <= touched * rows, "inconsistent selected total")
    half = length // 2
    if touched <= half:
        return 0
    if full >= half:
        return selected - rows * half
    return (selected - full * rows) * (touched - half) // (touched - full)


def check_shortening_cap(deadline):
    """Exact bounded checks of every f<=i<=t for lengths 4, 8, and 64."""
    triples = aggregate_checks = 0
    for length in (4, 8, 64):
        half = length // 2
        for full in range(length + 1):
            guard(deadline)
            for touched in range(full, length + 1):
                for occupancy in range(full, touched + 1):
                    excess = max(0, occupancy - half)
                    if full < half < touched:
                        require(excess * (touched - full)
                                <= (occupancy - full) * (touched - half),
                                "pointwise secant inequality failed")
                    require(excess <= shortening_cap(length, 1, occupancy, full, touched),
                            "integer single-row cap failed")
                    for occupancies in ((full, occupancy), (occupancy, touched),
                                        (occupancy, occupancy)):
                        actual = sum(max(0, value - half) for value in occupancies)
                        bound = shortening_cap(length, 2, sum(occupancies), full, touched)
                        require(actual <= bound, "integer aggregate secant cap failed")
                        aggregate_checks += 1
                    triples += 1
    guard(deadline)
    return triples, aggregate_checks


def support_geometry(weights):
    """Deterministic cap for any support with fixed per-batch cardinalities."""
    batches = len(weights)
    require(batches % B == 0, "regions must contain whole batches")
    require(all(0 <= value <= BATCH for value in weights), "invalid batch support size")
    per_region = batches // B
    counts = [sum(weights[j * per_region:(j + 1) * per_region]) for j in range(B)]
    touched = sum(value > 0 for value in counts)
    full = sum(value == per_region * BATCH for value in counts)
    n = batches * BATCH
    rows = n // B
    selected = sum(weights)
    cap = shortening_cap(B, rows, selected, full, touched)
    return {"D_cap": cap, "regions": touched, "full_regions": full, "size": selected}


def outer_geometry(activity):
    """The same cap specialized to the actual active-batch union U."""
    geometry = support_geometry([BATCH * int(live) for live in activity])
    return {"D_U": geometry["D_cap"], "regions": geometry["regions"],
            "full_regions": geometry["full_regions"], "U": geometry["size"]}


@lru_cache(maxsize=4096)
def log_comb(n, k):
    require(0 <= k <= n, "invalid binomial")
    return log2(comb(n, k))


def mask_logs(sequence):
    log_h = log_w_given_h = log_s = log_pairs_given_s = 0.0
    s_total = 0
    for (h, w, overlap), count in Counter(sequence).items():
        s = h + w - overlap
        log_h += count * log_comb(BATCH, h)
        log_w_given_h += count * (log_comb(h, overlap) + log_comb(BATCH - h, w - overlap))
        log_s += count * log_comb(BATCH, s)
        log_pairs_given_s += count * (log_comb(s, h) + log_comb(h, overlap))
        s_total += count * s
    require(abs((log_h + log_w_given_h) - (log_s + log_pairs_given_s)) < 1e-6,
            "H/S mask factorizations disagree")
    return log_h, log_w_given_h, log_s, log_pairs_given_s, s_total


def log_sum_exp2(values):
    peak = max(values)
    return peak + log2(sum(exp2(value - peak) for value in values))


def route_log_moment(s, n, log_tau, d_u):
    """Float version of M_p(s,tau), also capped by deterministic dimension."""
    small = min(s, n - s)
    shift = max(0, s - n // 2)
    if small == 0:
        moment = shift * log_tau
    else:
        alpha = small / n
        l0, l1 = log2(1.0 - alpha), log2(alpha)
        terms = [log_comb(B, i) + i * l1 + (B - i) * l0
                 + PACKET * log_tau * max(0, i - B // 2) for i in range(B + 1)]
        moment = shift * log_tau + (n // (PACKET * B)) * log_sum_exp2(terms)
    conditional = min(s // 2, d_u) * log_tau
    return min(moment, conditional), conditional


def score_profile(sequence, ranks, geometry):
    n = len(sequence) * BATCH
    log_h, log_l, log_s, log_ls, s = mask_logs(sequence)
    scores = {name: (float("inf"), None)
              for name in ("H_cond", "H_route", "S_route", "S_overlap", "old_ell",
                           "H_regional", "S_regional")}
    for log_tau in LOG_TAU_GRID:
        fraction = log_tau / FIELD_BITS
        h_moment, h_cond = route_log_moment(ranks["H"], n, log_tau, geometry["D_U"])
        s_moment, _ = route_log_moment(s, n, log_tau, geometry["D_U"])
        projective = -fraction * LOG_PROJECTIVE
        h_base = projective + log_h + fraction * log_l
        values = {
            "H_cond": h_base + (ranks["E"] - ranks["H"]) * log_tau + h_cond,
            "H_route": h_base + (ranks["E"] - ranks["H"]) * log_tau + h_moment,
            "S_route": projective + log_s + fraction * log_ls
                       + (ranks["W"] - s) * log_tau + s_moment,
            "S_overlap": projective + log_s + fraction * log_ls
                         + (ranks["W"] - ranks["R_overlap"] - s) * log_tau + s_moment,
            "old_ell": h_base + (ranks["caps"]["ell"] - ranks["H"]) * log_tau + h_moment,
            "H_regional": h_base + (ranks["E"] - ranks["H"]) * log_tau
                          + min(h_moment, geometry["capH"] * log_tau),
            "S_regional": projective + log_s + fraction * log_ls
                          + (ranks["W"] - ranks["R_overlap"] - s) * log_tau
                          + min(s_moment, geometry["capS"] * log_tau),
        }
        for name, value in values.items():
            require(isfinite(value), "nonfinite diagnostic score")
            if value < scores[name][0]:
                scores[name] = value, log_tau
    return scores


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--batches", type=int, choices=(64, 256, 1024, 16384), default=16384)
    parser.add_argument("--seconds", type=float, default=30.0)
    parser.add_argument("--profile", action="append", choices=tuple(p.name for p in PROFILES))
    args = parser.parse_args()
    require(1 <= args.seconds <= 60, "time cap must be in [1,60] seconds")
    deadline = monotonic() + args.seconds
    chosen = [p for p in PROFILES if not args.profile or p.name in args.profile]
    n = args.batches * BATCH
    require(len(chosen) <= 20 and n % (PACKET * B) == 0, "resource or geometry cap")
    print("FLOAT PROFILE DIAGNOSTIC: NOT A FINITE CERTIFICATE")
    print(f"B={B} p={PACKET} b={BATCH} q=2^{FIELD_BITS} N={n} batches={args.batches}; "
          f"profiles={len(chosen)} markers={len(LOG_TAU_GRID)} seconds_cap={args.seconds:g}")
    print("Absolute profiles only; no translations/placement union. Scores are uncapped log2 bounds.")
    print("H_cond uses D(A)<=min(H/2,D(U)); H_route also uses symmetrized packet moments.")
    print("S_route is a separate outer-compatibility bound; S_overlap adds only disjoint overlap pivots.")
    print("H_regional/S_regional additionally cap shortening using full/touched regions; old scores are retained.")
    print(f"Integer selected-rank regression: {check_selected_rank(deadline)} comparisons passed.")
    triples, aggregate_checks = check_shortening_cap(deadline)
    print(f"Integer regional-cap regression: {triples} f,t,i triples and {aggregate_checks} aggregate checks passed.")
    failures = []
    for profile in chosen:
        guard(deadline)
        sequence = expand(profile, args.batches)
        ranks = rank_caps(sequence, deadline)
        if "excluded" in ranks:
            print(f"{profile.name}: NECESSARY-EMPTY {ranks['excluded']}")
            continue
        geometry = outer_geometry(ranks["activity"])
        if geometry["D_U"] == 0:
            print(f"{profile.name}: OUTER-EMPTY D(U)=0; touched_regions={geometry['regions']}")
            continue
        if ranks["E"] <= 0:
            print(f"{profile.name}: NECESSARY-EMPTY amplitude dimension <=0; caps={ranks['caps']}")
            continue
        h_geometry = support_geometry([h for h, _, _ in sequence])
        s_geometry = support_geometry([h + w - overlap for h, w, overlap in sequence])
        geometry["capH"] = min(geometry["D_U"], ranks["H"] // 2, h_geometry["D_cap"])
        geometry["capS"] = min(geometry["D_U"], s_geometry["size"] // 2, s_geometry["D_cap"])
        scores = score_profile(sequence, ranks, geometry)
        eligible = ranks["W"] <= (35 * n) // 100
        best_name = min(("H_regional", "S_regional"), key=lambda name: scores[name][0])
        best, marker = scores[best_name]
        prefix = rank_caps(sequence[:64], deadline)
        boundary = (f"E/batch64={prefix['E']/64:.4f}" if "E" in prefix else "prefix64_empty")
        print(f"{profile.name}: weight={ranks['W']/n:.6f} eligible35={int(eligible)} "
              f"excursions={ranks['starts']} closed={ranks['closed']} free={int(ranks['free'])} "
              f"U={geometry['U']} D_U_cap={geometry['D_U']}")
        print(f"  capH={geometry['capH']} full/touched={h_geometry['full_regions']}/{h_geometry['regions']} "
              f"capS={geometry['capS']} full/touched={s_geometry['full_regions']}/{s_geometry['regions']}")
        if profile.mode == "phased":
            counts = Counter(sequence)
            print(f"  phases: activation=1 sparse65_1_1={counts[(65, 1, 1)]} "
                  f"dense128_128_128={counts[(BATCH, BATCH, BATCH)]}; terminal is free")
        print(f"  E={ranks['E']} E/batch={ranks['E']/args.batches:.4f} {boundary}; "
              f"caps={ranks['caps']} Rfresh={ranks['R_fresh']} Rselected={ranks['R_selected']} "
              f"Roverlap={ranks['R_overlap']}")
        print("  " + " ".join(f"{name}={value[0]:.3f}@logtau{value[1]:g}" for name, value in scores.items())
              + f" best={best_name}")
        if eligible:
            failures.append((best, profile.name, best_name, marker))
    guard(deadline)
    print("Largest remaining selected-profile bounds (not a search over all supports):")
    for value, name, kind, marker in sorted(failures, reverse=True)[:3]:
        print(f"  {name}: {value:.3f} bits via {kind}, log2(tau)={marker:g}")
    print("Screen complete only for the named profiles and fixed marker grid; no whole-code conclusion.")


if __name__ == "__main__":
    main()
