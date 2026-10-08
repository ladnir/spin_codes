"""Partial-row shortening in paired regional windows, no native encoder work.

The exact mean lower bound is a theorem over uniform regional permutations.
The optional finite route search is a reproducible witness, not a probability
certificate for the PRNG or the setup ensemble.
"""

from fractions import Fraction as F
from hashlib import sha256
from math import log2

from hd_blockwise_geometry import balanced_profile, failure_lower

WITNESS_HASH = 'c090228e9769e5288d3216dac581ded5e48cf16d521c145468d25b246a94b58c'
WITNESS_LENGTHS = (
    101,82,107,111,102,91,102,94,116,86,112,93,103,102,106,97,
    93,116,99,116,98,90,99,101,83,100,98,93,111,110,91,105,
    84,93,93,102,110,113,105,107,100,101,112,103,106,93,106,113,
    93,110,103,110,98,110,104,102,104,110,99,104,94,102,103,100,
    102,107,94,111,103,94,110,103,93,92,104,106,92,103,102,95,
    103,102,101,98,99,102,108,109,99,111,101,103,104,105,96,100,
    116,108,103,101,112,95,99,120,103,97,103,111,112,105,99,100,
    99,103,103,111,111,111,91,108,111,105,111,97,111,106,103,95,
)


def validate_frozen_witness():
    """Recount without search; physical even regions reverse selection order.

    The first selected W entries therefore form a physical suffix in every
    even region, and a physical prefix in every odd region (zero-based).
    """
    import numpy as np
    rows, regions, width = 16384, 128, 64
    rng = np.random.Generator(np.random.PCG64(0))
    route = np.stack([rng.permutation(rows).astype(np.int32) for _ in range(regions)])
    digest = sha256(route.tobytes()).hexdigest()
    if digest != WITNESS_HASH:
        raise ArithmeticError('The frozen route changed')
    occupancy = np.zeros(rows, dtype=np.int16)
    for c, packets in enumerate(WITNESS_LENGTHS):
        if not 1 <= packets < rows//width or len(np.unique(route[c])) != rows:
            raise ArithmeticError('Invalid window or route permutation')
        occupancy[route[c, :packets*width]] += 1
    histogram = np.bincount(occupancy, minlength=regions+1)
    dimension = sum((hits-64)*int(histogram[hits]) for hits in range(65, regions+1))
    weight = sum(WITNESS_LENGTHS)*width
    if dimension != 720 or weight != 838848 or histogram[128] != 0:
        raise ArithmeticError('The stored partial-row witness did not verify')
    physical = route.copy()
    physical[::2] = route[::2, ::-1]
    return dict(selection_order_sha256=digest, physical_route_sha256=sha256(physical.tobytes()).hexdigest(),
                physical_route_definition='Reverse every even-indexed selection-order region; retain odd regions',
                dimension=dimension, components=64, state8_exit_constraints=512,
                state8_bad_subspace_dimension_lower=dimension-512, low_weight=weight,
                shortened_nonzero_rows=int(histogram[65:].sum()), maximum_hits=int(occupancy.max()),
                fully_contained_rows=int(histogram[128]),
                excess_histogram={str(hits): int(histogram[hits]) for hits in range(65, regions+1) if histogram[hits]},
                setup_probability_claim=False)


def paired_window_lower(state=4, length=128, width=64, distance=F(2, 5)):
    profile = balanced_profile(length, width, distance)
    components = length//2
    probability = failure_lower(profile['mean'], profile['maximum'], components*state)
    return dict(outer=[length, length//2], target_distance=str(distance), state=state,
                components=components, constraints=components*state,
                mean=float(profile['mean']), weight_limit=profile['total_window'],
                probability_lower=str(probability), probability_lower_float=float(probability),
                rejected_at_40_bits=probability > F(1, 1 << 40),
                margin_upper=(-log2(probability.numerator)+log2(probability.denominator)) if probability else None,
                exact=True)


def route_witness(seed=0, steps=256, jumps=(1,), length=128):
    """Bounded deterministic packet-length ascent for one reproducible route.

Selection order is reversed in every even region, giving suffix windows;
odd regions give prefixes. Every consecutive pair is one connected interval.
Lengths stay below R so different pairs never join. The route matrix stores
source-row labels in selection order; reversing even rows gives physical order.
"""
    import numpy as np
    if length not in (128, 256):
        raise ValueError('The bounded search supports outer128 and256 only')
    b, width, parity = length, 64, length//2
    rows, components = (1 << 20)//parity, length//2
    rng = np.random.Generator(np.random.PCG64(seed))
    route = np.stack([rng.permutation(rows).astype(np.int32) for _ in range(b)])
    profile = balanced_profile(b, width, F(2, 5))
    lengths = np.asarray(profile['regional_packet_counts'], dtype=np.int32)
    counts = np.zeros(rows, dtype=np.int16)
    for c in range(b):
        counts[route[c, :int(lengths[c])*width]] += 1
    initial = int(np.maximum(counts-parity, 0).sum())
    history = [initial]
    moves = []
    for _ in range(steps):
        best = (0, -1, -1, 0)
        for jump in jumps:
            if np.any(lengths < jump) or np.any(lengths+jump > rows//width):
                continue
            additions = route[np.arange(b)[:, None], lengths[:, None]*width+np.arange(jump*width)]
            removals = route[np.arange(b)[:, None], (lengths[:, None]-jump)*width+np.arange(jump*width)]
            gains = (counts[additions] >= parity).sum(axis=1)
            losses = (counts[removals] > parity).sum(axis=1)
            gains[lengths+jump >= rows//width] = -10000
            losses[lengths <= jump] = 10000
            for grow in np.argsort(gains)[-8:]:
                adjusted = losses.copy()
                at_threshold = additions[grow][counts[additions[grow]] == parity]
                if len(at_threshold):
                    adjusted += np.isin(removals, at_threshold).sum(axis=1)
                adjusted[grow] = 10000
                shrink = int(np.argmin(adjusted))
                increase = int(gains[grow]-adjusted[shrink])
                if increase > best[0]:
                    best = (increase, int(grow), shrink, jump)
        increase, grow, shrink, jump = best
        if increase <= 0:
            break
        removed = route[shrink, int(lengths[shrink]-jump)*width:int(lengths[shrink])*width]
        added = route[grow, int(lengths[grow])*width:int(lengths[grow]+jump)*width]
        counts[removed] -= 1
        counts[added] += 1
        lengths[shrink] -= jump
        lengths[grow] += jump
        exact = int(np.maximum(counts-parity, 0).sum())
        if exact != history[-1]+increase:
            raise ArithmeticError('Packet move prediction was not exact')
        history.append(exact)
        moves.append((grow, shrink, jump))
    # Independent complete recount verifies the final routed support.
    check = np.zeros(rows, dtype=np.int16)
    for c in range(b):
        if not np.array_equal(np.sort(route[c]), np.arange(rows)):
            raise ArithmeticError('A region was not a permutation')
        check[route[c, :int(lengths[c])*width]] += 1
    if not np.array_equal(check, counts):
        raise ArithmeticError('Incremental support count failed replay')
    physical = route.copy()
    physical[::2] = route[::2, ::-1]
    return dict(seed=seed, bit_generator='NumPy PCG64', selection_order_sha256=sha256(route.tobytes()).hexdigest(),
                physical_route_sha256=sha256(physical.tobytes()).hexdigest(),
                physical_route_definition='Reverse every even-indexed selection-order region; retain odd regions',
                outer=[length, length//2],
                steps=len(moves), jumps=jumps, initial_dimension=initial, final_dimension=history[-1],
                components=components, low_weight=int(lengths.sum())*width,
                min_packets=int(lengths.min()), max_packets=int(lengths.max()),
                lengths=lengths.tolist(), moves=moves,
                state_dimensions_rejected_on_this_route=[s for s in (4, 8, 16, 32) if history[-1] > components*s],
                source_nullity_lower={s: max(0, history[-1]-components*s) for s in (4, 8, 16, 32)},
                setup_failure_probability_not_certified=True)


def main():
    import json
    for state in (4, 8, 16):
        receipt = paired_window_lower(state)
        receipt.pop('probability_lower')
        print(json.dumps(receipt), flush=True)
    print(json.dumps(validate_frozen_witness()), flush=True)


if __name__ == '__main__':
    main()
