"""All-route episode bounds with a bounded number of defective batches.

Fixed rejected-uniform SL2 construction; restricted families, not whole-code
closure. Standard library only, no data files, bounded serial checks.
"""
import argparse
import itertools
import math
import time
from functools import lru_cache
from fractions import Fraction

from diagnose_pair_episodes import volume_bits

WIDTH, BATCHES, REGIONS, REGION_BATCHES, BITS = 512, 4096, 256, 16, 32
N = WIDTH * BATCHES


def ceildiv(a, b):
    return (a + b - 1) // b


def logceil_ratio(a, b):
    k = max(0, a.bit_length() - b.bit_length())
    return k + int(a > b << k)


def volume_bound(length, weight, field_bits=BITS):
    """Exact integer majorant using binom(n,z)<=(3n/z)^z, not FLOAT."""
    if weight >= length:
        return field_bits * length
    if (1 << field_bits) * (length - weight + 1) < weight:
        return field_bits * length
    z = min(weight, length - weight)
    entropy = 0 if z == 0 else min(length, z * logceil_ratio(3 * length, z))
    return min(field_bits * length, length + field_bits * weight,
               field_bits * weight + weight.bit_length() + entropy)


@lru_cache(maxsize=256)
def early_region_caps(episodes, partial, initial, region_length=16):
    """Max region cover of early, proper-only intervals for each total length."""
    caps = {0: 0}
    for interval in range(episodes - 1):
        nxt = {}
        for used, regions in caps.items():
            for length in range(1, partial - used + 1):
                # The output/input interval includes its preceding batch.
                extra = (ceildiv(length + 1, region_length) if interval == 0 and initial else
                         (length + 2 * region_length - 1) // region_length)
                new = used + length
                nxt[new] = max(nxt.get(new, -1), regions + extra)
        caps = nxt
    return tuple(sorted(caps.items()))


def profile(h, initial, episodes, live, permille, no_full_reset=False):
    """One interval geometry envelope, before its interval-placement union.

    h is an upper bound on defective incoming batches. initial records
    whether Q_1 is nonzero; episodes counts nonzero incoming-state intervals.
    live counts incoming states Q_1,...,Q_(L-1); Q_L is unconstrained.
    """
    assert initial in (0, 1) and 1 <= episodes <= h + initial
    assert 1 <= live < BATCHES
    partial = min(live, 9 * (h + initial) - episodes)
    # Every full state supplies b random coordinates; every other live
    # state supplies at least two, even if its only live pair is one-sided.
    random_min = WIDTH * (live - partial) + 2 * partial
    free_regions = min(REGIONS, (live + 30 * episodes - 14 * initial) // 16)
    input_regions = min(REGIONS, (live + 31 * episodes - 15 * initial) // 16 + 1)
    if no_full_reset:
        if live <= partial:
            # Do not infer absence of a full state from an upper partial cap.
            # The old input-region bound still covers this small case.
            pass
        else:
            # Some incoming state is full. If it never reaches zero later,
            # the final live interval ends at Q_(L-1). Earlier intervals are
            # entirely proper and have combined length <=partial.
            candidates = [r + ceildiv(live - used + 1, 16)
                          for used, r in early_region_caps(episodes, partial, initial)
                          if used < live]
            if not candidates:
                return None
            input_regions = min(input_regions, max(candidates))
    if input_regions <= 128:
        return None
    assert free_regions > 0
    # Rowwise (f-128)_+ <= f*(r-128)_+/r. Using the integer loss is safe.
    pair_loss = (random_min if free_regions <= 128 else
                 ceildiv(128 * random_min, free_regions))
    # Full states form <=h+1 intervals. If there are h+1, the first starts
    # at Q_1; otherwise there are <=h arbitrary intervals. Both cases give
    # the boundary allowance 30h+16. F contains every coordinate of these
    # full batches, so the outer parity-check rank on F is at least that
    # on their union G. Take the earliest full_batches full states if there
    # are more: truncation preserves the interval bound. This subset step
    # avoids assuming the rounded numerical rank envelope is monotone.
    full_batches = live - partial
    boundary = 30 * h + 16
    full_region_cap = min(REGIONS, (full_batches + boundary) // 16)
    full_mass = WIDTH * full_batches
    full_loss = (full_mass if full_region_cap <= 128 else
                 ceildiv(128 * full_mass, full_region_cap))
    complete_regions = max(0, ceildiv(full_batches - boundary, 16))
    complete_loss = 8192 * min(128, complete_regions)
    inverse_loss = max(pair_loss, full_loss, complete_loss)
    terminal_extra = int(not (no_full_reset and live > partial))
    slots = min(N, WIDTH * (live + episodes + terminal_extra))
    weight = permille * N // 1000
    volume = volume_bound(slots, weight)
    inverse = 1 + volume - BITS * inverse_loss
    # Fix the interval placement BEFORE sampling. All allowed inverse inputs
    # lie in the backward-expanded intervals plus the free terminal batch.
    # The RS subspace supported there has this deterministic dimension cap.
    input_dimension = (0 if input_regions <= 128 else
                       slots * (input_regions - 128) // input_regions)
    # For a fixed message, expose active pair outputs adaptively. Their joint
    # density is <2 times uniform. If there are >=random_min random positions,
    # the first random_min draws already have weight <=W on a bad output.
    # Pad after termination only to define the comparison experiment.
    forward = 1 + volume_bound(random_min, weight) + BITS * (input_dimension - random_min)
    bound = min(inverse, forward)
    return dict(h=h, initial=initial, episodes=episodes, live=live,
                partial_cap=partial, random_min=random_min,
                free_regions=free_regions, input_regions=input_regions,
                slots=slots, pair_loss=pair_loss, full_loss=full_loss,
                complete_loss=complete_loss,
                inverse_loss=inverse_loss, input_dimension=input_dimension,
                forward_upper=forward, inverse_upper=inverse, log2_upper=bound)


def cover(h, permille, no_full_reset=False):
    began = time.monotonic()
    worst = None
    inverse_worst = -(1 << 60)
    cells = excluded = 0
    for initial in (0, 1):
        for episodes in range(1, h + initial + 1):
            for live in range(episodes, BATCHES):
                item = profile(h, initial, episodes, live, permille, no_full_reset)
                cells += 1
                if item is None:
                    excluded += 1
                else:
                    inverse_worst = max(inverse_worst, item['inverse_upper'])
                    if worst is None or item['log2_upper'] > worst['log2_upper']:
                        worst = item
            if time.monotonic() - began > 20:
                raise TimeoutError('20-second coverage cap')
    assert worst is not None
    # <=(L+1)^(2c) interval endpoint choices for each (initial,c).
    # This already sums over lengths; do not union over live again.
    placement_bits = 2 * (h + 1) * BATCHES.bit_length()
    type_bits = (2 * (h + 1) - 1).bit_length()
    return dict(defects_at_most=h, permille=permille, no_full_reset=no_full_reset, cells=cells,
                excluded=excluded, worst=worst,
                placement_bits=placement_bits, type_bits=type_bits,
                family_log2_upper=worst['log2_upper'] + placement_bits + type_bits,
                inverse_only_family_upper=inverse_worst + placement_bits + type_bits,
                all_routes=True, whole_code_certificate=False)


def route_reset_cover():
    """h<=4 at 46%: exact-volume refinement, then a sparse outer route tail.

    All probability arithmetic is integer/rational. This is NOT full history
    coverage: h>=5 is outside the theorem.
    """
    began = time.monotonic()
    h, weight = 4, 460 * N // 1000
    placement_bits = 134
    # I occupies <=130 regions, <=11 of them partial, and <=994816 slots.
    # Each row eligible for a nonzero RS word occupies >=129 allowed regions.
    # Positive elementary-symmetric domination gives the joint-row cap beta^s.
    beta = Fraction(0)
    beta_case = None
    for r in (129, 130):
        for j in range(1, 12):
            remaining = Fraction(994816, 8192) - (r - j)
            if remaining < 0:
                continue
            p = min(Fraction(1), remaining / j)
            k = 129 - (r - j)
            candidate = math.comb(j, k) * p ** k
            if candidate > beta:
                beta, beta_case = candidate, (r, j, k)
    rows = 32
    tail = math.comb(8192, rows) * beta ** rows
    # Integer ceiling of log2 tail, valid for positive fractions <1 too.
    route_bits = tail.numerator.bit_length() - (tail.denominator.bit_length() - 1) + placement_bits
    # Outside the route event there are <=31 eligible rows, each of dimension
    # <=2. A full-state reset requires an active->zero event on every state
    # coordinate; choose its transition index in <=L^b ways.
    dimension = 2 * (rows - 1)
    reset_bits = 1 + BITS * dimension + WIDTH * 12 - BITS * WIDTH
    worst = -(1 << 60)
    refined = rescued = 0
    residual = []
    cells = excluded = 0
    for initial in (0, 1):
        for episodes in range(1, h + initial + 1):
            for live in range(episodes, BATCHES):
                item = profile(h, initial, episodes, live, 460)
                cells += 1
                if item is None:
                    excluded += 1
                    continue
                bound = item['log2_upper']
                if bound + placement_bits > -80:
                    refined += 1
                    bound = min(bound,
                                1 + volume_bits(item['slots'], weight) - BITS * item['inverse_loss'],
                                1 + volume_bits(item['random_min'], weight) +
                                BITS * (item['input_dimension'] - item['random_min']))
                if bound + placement_bits > -80:
                    if item['input_regions'] <= 130 and item['slots'] <= 994816:
                        bound = min(bound, reset_bits)
                        rescued += 1
                    else:
                        residual.append((initial, episodes, live, bound))
                worst = max(worst, bound)
                if time.monotonic() - began > 30:
                    raise TimeoutError('30-second route/reset cover cap')
    no_reset = cover(4, 460, True)['family_log2_upper']
    pieces = (route_bits, worst + placement_bits, no_reset)
    final = max(pieces) + 2  # at most three events
    print(dict(kind='four_defect_route_reset', cells=cells, excluded=excluded,
               exact_volume_refinements=refined, route_reset_rescues=rescued,
               residual=residual, beta_case=beta_case,
               beta=(beta.numerator, beta.denominator), route_failure_upper=route_bits,
               good_route_dimension_cap=dimension, per_placement_reset_upper=reset_bits,
               full_reset_good_route_upper=worst + placement_bits,
               no_full_reset_upper=no_reset, family_log2_upper=final,
               closed=not residual and final < -40, whole_code_certificate=False,
               seconds=time.monotonic() - began), flush=True)


def test_route_reset():
    """Independent exact small-route and adaptive reset probability checks."""
    began = time.monotonic()
    route_checks = 0
    for sizes in itertools.product(range(5), repeat=3):
        masks = [[m for m in range(16) if m.bit_count() == size] for size in sizes]
        outcomes = list(itertools.product(*masks))
        p = tuple(Fraction(size, 4) for size in sizes)
        for threshold in (1, 2, 3):
            elementary = sum((math.prod(p[j] for j in subset)
                              for subset in itertools.combinations(range(3), threshold)), Fraction(0))
            averaged = math.comb(3, threshold) * (sum(p) / 3) ** threshold
            assert elementary <= averaged
            histogram = [0] * 5
            for regions in outcomes:
                eligible = sum(sum((m >> row) & 1 for m in regions) >= threshold for row in range(4))
                histogram[eligible] += 1
            for rows in range(1, 5):
                actual = Fraction(sum(histogram[rows:]), len(outcomes))
                assert actual <= math.comb(4, rows) * elementary ** rows
                route_checks += 1
    # GF4 comparison experiment, two coordinates and two updates; begin full.
    # We count the SUPERSET where each coordinate has some active->zero step,
    # regardless of whether the pair ever resets simultaneously.
    reset_checks = 0
    for inputs in itertools.product(range(16), repeat=2):
        bad = 0
        for noise in itertools.product(range(16), repeat=2):
            state, deaths = 5, 0  # (1,1)
            for x, z in zip(inputs, noise):
                output = z if state else x
                following = state ^ output
                for j in range(2):
                    if (state >> (2*j)) & 3 and not ((following >> (2*j)) & 3):
                        deaths |= 1 << j
                state = following
            bad += deaths == 3
        assert Fraction(bad, 256) <= Fraction(2 ** 2, 4 ** 2)
        reset_checks += 1
    assert time.monotonic() - began < 20
    print(dict(exact_route_tail_checks=route_checks,
               exact_adaptive_reset_checks=reset_checks,
               seconds=time.monotonic() - began), flush=True)


def test_geometry():
    """Exhaustive path DP independent of the numerical bound above."""
    began = time.monotonic()
    d, length = 2, 8
    width = 1 << d
    full = (1 << width) - 1

    def single(mask, t):
        return any(((mask >> j) & 1) != ((mask >> (j ^ (1 << (t % d)))) & 1)
                   for j in range(width))

    # mask, initial, initial_full, h, live_intervals, full_intervals, partial
    states = {(mask, int(bool(mask)), int(mask == full), 0,
               int(bool(mask)), int(mask == full), int(0 < mask < full))
              for mask in range(full + 1)}
    transitions = 0
    for t in range(length - 1):
        nxt = set()
        for mask, initial, initial_full, h, intervals, full_intervals, partial in states:
            for following in range(full + 1):
                bad = int(single(mask, t) or bool(following & ~mask))
                nxt.add((following, initial, initial_full, h + bad,
                         intervals + int(mask == 0 and following != 0),
                         full_intervals + int(mask != full and following == full),
                         partial + int(0 < following < full)))
                transitions += 1
        states = nxt
        assert len(states) < 100000
        if time.monotonic() - began > 20:
            raise TimeoutError('20-second geometry-test cap')
    for mask, initial, initial_full, h, intervals, full_intervals, partial in states:
        h += int(single(mask, length - 1))
        assert intervals <= h + initial
        assert full_intervals <= h + initial_full
        assert partial <= d * (h + initial) - intervals
    # Exhaust intervals and their one-batch backward dilation on a toy route.
    interval_checks = 0
    terminal_checks = 0
    region_length, batches = 4, 16
    for bits in range(1, 1 << (batches - 1)):
        live = {j + 1 for j in range(batches - 1) if bits >> j & 1}
        initial = int(1 in live)
        starts = {t for t in live if t - 1 not in live}
        c, count = len(starts), len(live)
        r = len({t // region_length for t in live})
        complete = sum(all(t in live for t in range(j, j + region_length))
                       for j in range(0, batches, region_length))
        allowed = live | {t - 1 for t in starts} | {batches - 1}
        ri = len({t // region_length for t in allowed})
        assert region_length * r <= count + (2 * region_length - 2) * c - (region_length - 2) * initial
        assert complete >= max(0, ceildiv(count - (2 * region_length - 2) * c, region_length))
        assert ri <= (count + (2 * region_length - 1) * c - (region_length - 1) * initial) // region_length + 1
        assert len(allowed) <= count + c + 1
        if batches - 1 in live:
            last_start = max(starts)
            early_length = sum(t < last_start for t in live)
            cap = max(r0 + ceildiv(count - used + 1, region_length)
                      for used, r0 in early_region_caps(c, early_length, initial, region_length)
                      if used < count)
            assert ri <= cap
            terminal_checks += 1
        interval_checks += 1
        if bits % 4096 == 0 and time.monotonic() - began > 20:
            raise TimeoutError('20-second interval-test cap')
    print(dict(path_states=len(states), path_transitions=transitions,
               interval_geometries=interval_checks,
               terminal_live_geometries=terminal_checks), flush=True)
    comparisons = 0
    for n in range(1, 200):
        for z in range(n + 1):
            k = min(z, n - z)
            entropy = 0 if k == 0 else min(n, k * logceil_ratio(3 * n, k))
            assert math.comb(n, z) <= 1 << entropy
            comparisons += 1
    print(dict(exact_binomial_majorants=comparisons), flush=True)
    volume_checks = 0
    for field_bits in range(1, 5):
        q = 1 << field_bits
        for n in range(1, 32):
            for w in range(n + 1):
                actual = sum(math.comb(n, j) * (q - 1) ** j for j in range(w + 1))
                assert actual <= 1 << volume_bound(n, w, field_bits)
                volume_checks += 1
    print(dict(exact_volume_checks=volume_checks), flush=True)
    # Uniform comparison experiment, not a test of an SL2 ensemble at q=2.
    # Noise is consumed adaptively only when the incoming pair is nonzero.
    stopping_checks = 0
    for inputs in itertools.product(range(4), repeat=4):
        hist = {}
        for noise in itertools.product(range(4), repeat=4):
            state = draws = weight = 0
            for x, z in zip(inputs, noise):
                draws += 2 * int(bool(state))
                output = z if state else x
                weight += output.bit_count()
                state ^= output
            hist[draws, weight] = hist.get((draws, weight), 0) + 1
        for minimum in (2, 4, 6):
            for weight in range(9):
                bad = sum(count for (d0, w0), count in hist.items()
                          if d0 >= minimum and w0 <= weight)
                ball = sum(math.comb(minimum, j) for j in range(minimum + 1)
                           if j <= weight)
                assert bad * (1 << minimum) <= 256 * ball
                stopping_checks += 1
    print(dict(adaptive_uniform_checks=stopping_checks), flush=True)


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--test', action='store_true')
    p.add_argument('--defects', type=int, default=1)
    p.add_argument('--distance', type=int, default=460)
    p.add_argument('--ladder', action='store_true')
    p.add_argument('--no-full-reset', action='store_true')
    p.add_argument('--route-reset', action='store_true')
    p.add_argument('--test-route-reset', action='store_true')
    args = p.parse_args()
    assert 0 <= args.defects <= 32 and 0 <= args.distance <= 500
    if args.test:
        test_geometry()
    elif args.test_route_reset:
        test_route_reset()
    elif args.route_reset:
        route_reset_cover()
    elif args.ladder:
        for h in range(7):
            print(cover(h, args.distance, args.no_full_reset), flush=True)
    else:
        print(cover(args.defects, args.distance, args.no_full_reset), flush=True)


if __name__ == '__main__':
    main()
