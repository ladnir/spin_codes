"""Outward occupancy-one screen for RS16 with eight-bit routed packets.

Only q=1 is covered. The kernel uses exact zero-state births and a valid
uniform-GL envelope for nonzero states. No saved spectrum or numerical
endpoint is accepted. The four coordinate swaps make every packet's
feedback injective. The caller must retain state across steps and regions.
"""
from __future__ import annotations

import argparse
from collections import Counter, defaultdict
from fractions import Fraction as Q
import hashlib
import json
from math import comb
from pathlib import Path
import sys
from time import monotonic

HERE = Path(__file__).resolve().parent
DESIGN = HERE.parents[1] / 'k16_design'
sys.path.insert(0, str(DESIGN))
import packet_q1 as retained
import rs_outer
from flint import arb, arb_mat, ctx

aq, up = retained.kernel_t64.aq, retained.kernel_t64.up


def rank(values):
    pivots = {}
    for value in values:
        while value:
            bit = value.bit_length()-1
            if bit in pivots:
                value ^= pivots[bit]
            else:
                pivots[bit] = value
                break
    return len(pivots)


def images_from_rows(rows):
    result = [0]
    for row in rows:
        result.extend(value ^ row for value in result[:])
    return tuple(result)


def prepare(rows, packet_bits=8):
    """Build exact packet profiles and single-packet birth counts."""
    rows = tuple(rows)
    width = max(row.bit_length() for row in rows)
    if width % packet_bits or rank(rows) != len(rows):
        raise ValueError('full rank and a complete packet geometry required')
    windows, state_bits = width // packet_bits, len(rows)
    columns = tuple(sum(((row >> p) & 1) << i for i, row in enumerate(rows))
                    for p in range(width))
    packet_ranks = [rank(columns[p:p+packet_bits])
                    for p in range(0, width, packet_bits)]
    if packet_ranks != [packet_bits]*windows:
        raise ValueError('q1 kernel requires injective single-packet feedback')
    images = images_from_rows(rows)
    profiles = Counter()
    mask = (1 << packet_bits)-1
    for image in images[1:]:
        profile = [0]*(packet_bits+1)
        for p in range(0, width, packet_bits):
            profile[((image >> p) & mask).bit_count()] += 1
        profiles[tuple(profile)] += 1
    births = Counter()
    for p in range(0, width, packet_bits):
        feedbacks = images_from_rows(columns[p:p+packet_bits])
        for label in range(1, 1 << packet_bits):
            destination = feedbacks[label]
            assert destination != 0
            births[label.bit_count(), images[destination].bit_count()] += 1
    spectrum = Counter(image.bit_count() for image in images)
    assert sum(births.values()) == windows*mask
    return dict(rows=rows, columns=columns, width=width, bits=state_bits,
                windows=windows, packet_bits=packet_bits, packet_ranks=packet_ranks,
                profiles=profiles, births=births, spectrum=spectrum,
                levels=tuple(sorted(set(spectrum)-{0})))


def construction():
    source = retained.kernel_t64.SELECTED_MAP
    original = json.loads(source.read_bytes())
    old_rows = tuple(int(row, 16) for row in original['generator_rows_hex'])
    permutation = list(range(64))
    for left in (7, 23, 39, 55):
        permutation[left], permutation[left+8] = permutation[left+8], permutation[left]
    rows = tuple(sum(((row >> old) & 1) << new
                     for new, old in enumerate(permutation)) for row in old_rows)
    data = prepare(rows)
    assert data['width'] == 64 and data['bits'] == 16
    assert data['spectrum'] == Counter(image.bit_count() for image in images_from_rows(old_rows))
    assert all((a & b).bit_count() % 2 == 0 for a in rows for b in rows)
    record = dict(t=64, s=16, packet_bits=8, windows=8,
                  coordinate_permutation=permutation,
                  expansion_rows_hex=list(map(hex, rows)),
                  feedback_columns=list(data['columns']),
                  feedback_definition='C=A^T after the declared coordinate permutation',
                  packet_ranks=data['packet_ranks'],
                  expansion_spectrum=dict(sorted(data['spectrum'].items())),
                  feedback_times_expansion_zero=True,
                  distribution='independent uniform GL16 at each physical step',
                  zero_initial_state=True, final_flush=False,
                  state_continuity='retained between all physical steps and regions',
                  source_sha256=hashlib.sha256(source.read_bytes()).hexdigest())
    return data, record


def emission(data, profile, z, occupied):
    """Exact tilted emission for a fixed expansion profile and j=0 or 1."""
    total_weight = sum(w*count for w, count in enumerate(profile))
    if occupied == 0:
        return z**total_weight
    if occupied != 1:
        raise ValueError('only zero or one occupied packet is supported')
    b, W = data['packet_bits'], data['windows']
    value = arb(0)
    for w, count in enumerate(profile):
        if count:
            # All b-bit output strings except the original packet are possible.
            local = sum(((comb(b, k)-int(k == w))*z**k
                         for k in range(b+1)), arb(0))
            value += count*z**(total_weight-w)*local
    return value/(W*((1 << b)-1))


def local_operators(data, tilt):
    """Return outward operators for j=0,1; later occupancy is unsupported."""
    z = (-aq(Q(tilt))).exp()
    levels = data['levels']
    n, states = 2+len(levels), (1 << data['bits'])-1
    # Coordinates: zero, uniform-nonzero density, then expansion-weight classes.
    weights = {profile: sum(w*c for w, c in enumerate(profile))
               for profile in data['profiles']}
    result = []
    for occupied in (0, 1):
        rows = [[arb(0)]*n for _ in range(n)]
        if occupied == 0:
            rows[0][0] = arb(1)
        else:
            denominator = data['windows']*((1 << data['packet_bits'])-1)
            for (emitted_weight, birth_weight), count in data['births'].items():
                rows[0][2+levels.index(birth_weight)] += count*z**emitted_weight/denominator
        moments = {profile: emission(data, profile, z, occupied)
                   for profile in data['profiles']}
        mean = up(sum((count*moments[profile] for profile, count in data['profiles'].items()), arb(0))/states)
        bounds = [mean, *[max(up(moments[p]) for p in moments if weights[p] == level)
                         for level in levels]]
        for row, bound in enumerate(bounds, 1):
            rows[row][0] = up(bound/states) if occupied else arb(0)
            rows[row][1] = bound
        result.append(arb_mat([[up(value) for value in row] for row in rows]))
    return result


def run(tilts, precision=256):
    ctx.prec = precision
    start = monotonic()
    data, map_record = construction()
    counts = rs_outer.expected_group_support_counts(n=16, k=8,
                packet_bits=8, packets_per_symbol=2)
    assert sum(counts) == (1 << 128)-1 and len(counts) == 33 and counts[0] == 0
    groups, regions, steps_per_region, threshold = 512, 32, 64, 13107
    best = [arb(1)]*33
    trials = []
    print(f'Prepared {len(data["profiles"])} exact expansion profiles; packet ranks={data["packet_ranks"]}', flush=True)
    for tilt in tilts:
        local = local_operators(data, tilt)
        regional = retained.placement(local, epochs=steps_per_region, windows=8,
                      maximum_groups=1, rounding=retained.rounded)
        moments = retained.support_moments(regional[0], regional[1], regions)
        factor = (aq(Q(tilt))*threshold).exp()
        best = [min(old, up(factor*new)) for old, new in zip(best, moments)]
        terms = [up(groups*aq(count)*value) for count, value in zip(counts, best)]
        upper = up(sum(terms, arb(0)))
        margin = str(-upper.log()/arb(2).log())
        trials.append(dict(tilt=tilt, margin_bits=margin,
                           endpoint=[int(x) for x in upper.man_exp()],
                           dominant_supports=sorted(range(33), key=lambda v: float(terms[v]), reverse=True)[:5]))
        print(f'tilt={tilt} accumulated q1 margin={margin}; dominant={trials[-1]["dominant_supports"]}', flush=True)
    return dict(schema='rs16-packet8-q1-screen-1', K=65536, N=131072,
                groups=groups, regions=regions, steps_per_region=steps_per_region,
                threshold=threshold, precision=precision, map_record=map_record,
                q1_margin_bits=trials[-1]['margin_bits'], trials=trials,
                occupancy_covered=[1], whole_code_certificate=False,
                q1_below_2_minus_40=bool(upper < arb(2)**-40),
                elapsed_seconds=monotonic()-start)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--tilts', nargs='+', default=['.00256', '.00512', '.01024', '.0256', '.0512'])
    parser.add_argument('--precision', type=int, default=256)
    args = parser.parse_args()
    print(json.dumps(run(args.tilts, args.precision), indent=2))


if __name__ == '__main__':
    main()
