"""Compact deterministic fixtures for two realized-route dimension obstructions.

No optimization is run by this verifier. It regenerates each route and counts
window hits independently using concatenation plus bincount. These fixtures
are not a lower bound on the probability of a bad uniformly random route.
"""

from hashlib import sha256
import json

import numpy as np


WITNESSES = (
    {
        "route_seed": 1,
        "numpy_version": "2.4.4",
        "rows": 47680,
        "regions": 32,
        "outer_redundancy": 10,
        "state": 4,
        "alignment": 16,
        "bad_output_limit": 152575,
        "max_dimension_found": 136,
        "excluded_tail_rows": 18,
        "tail_excluded_dimension": 136,
        "total_window_symbols": 152064,
        "route_sha256": "7d3623f6715d405f7409680c4693ae215097a06bcc42c78c08a728a84d1503ed",
        "starts_sha256": "419601156d0c1a28d0d66a1ef12eb7f1a459b61be505c0207f982c0077573e24",
        "widths_sha256": "d06060b8d35c2daf81e3a1488e05c1126179b8fb328be7fcb7eb1324dfc89ea5",
        "starts": [
            15408, 25216, 24624, 23856, 15520, 13664, 4656, 20352,
            28336, 12368, 34752, 5824, 6800, 41696, 11952, 35312,
            15360, 22624, 1232, 30144, 5232, 20816, 20416, 18128,
            3152, 40064, 13728, 30864, 5552, 64, 39504, 432,
        ],
        "widths": [
            5936, 6288, 6416, 4752, 6128, 4640, 4176, 2240,
            3344, 4032, 4688, 3488, 2880, 4336, 5376, 5664,
            4400, 5712, 4160, 6016, 4336, 4608, 3440, 3168,
            5824, 5296, 5488, 5328, 4480, 5936, 5360, 4128,
        ],
        "positive_dimension_row_hits": [
            [792,12], [939,12], [2988,12], [3442,14], [4211,12], [4952,13],
            [5026,11], [6268,15], [6807,12], [8606,12], [8818,12], [9731,14],
            [9991,11], [10161,11], [10166,11], [10601,12], [10988,11], [11728,13],
            [12887,12], [14318,11], [15127,11], [15545,13], [17112,13], [17680,15],
            [19066,11], [19681,12], [20475,14], [20788,12], [21131,12], [21861,12],
            [24347,11], [24414,16], [24517,11], [24699,14], [24913,14], [25055,11],
            [25361,11], [25924,11], [27718,13], [28125,12], [28938,13], [29465,11],
            [31735,11], [31872,11], [32383,11], [32866,14], [33812,11], [34616,14],
            [39054,11], [39108,14], [39624,13], [39938,11], [40828,12], [42076,11],
            [42471,11], [42526,12], [42789,13], [43523,11], [44582,12], [45302,12],
            [45550,13], [45839,12],
        ],
        "high_hit_histogram": {"11":23,"12":19,"13":9,"14":8,"15":2,"16":1},
    },
    {
        "route_seed": 2,
        "numpy_version": "2.4.4",
        "rows": 47680,
        "regions": 32,
        "outer_redundancy": 10,
        "state": 4,
        "alignment": 16,
        "bad_output_limit": 152575,
        "max_dimension_found": 133,
        "excluded_tail_rows": 18,
        "tail_excluded_dimension": 132,
        "total_window_symbols": 152064,
        "route_sha256": "e2e57480572b323d1b505f973b9b544921140cd492daee7eb16d202a4e9e849d",
        "starts_sha256": "d0078283ed6858f8bcea9d7ad88fbbff5fad4735d57c4c3bfea5b2d72488a312",
        "widths_sha256": "abc89444116b81c5d1489a88fc30feed3fe4ea967a685411109a75505fec889c",
        "starts": [
            20128, 13872, 7680, 7088, 25744, 6240, 36512, 20352,
            37360, 41712, 432, 30320, 18416, 36480, 17696, 16544,
            34752, 14160, 18928, 20096, 8128, 27728, 2128, 5520,
            34976, 41696, 18160, 21024, 17376, 10144, 15568, 15392,
        ],
        "widths": [
            4496, 7360, 5744, 3776, 3648, 3808, 3440, 5728,
            3792, 4752, 4736, 6128, 5040, 3920, 5472, 4224,
            6096, 4368, 3328, 2464, 3520, 5072, 3520, 4656,
            6016, 3856, 5968, 4256, 4256, 7936, 6480, 4208,
        ],
        "positive_dimension_row_hits": [
            [47,11], [1643,14], [1924,12], [2125,12], [2214,12], [2554,12],
            [2794,14], [3205,15], [4454,11], [5244,14], [6111,11], [6643,15],
            [6754,13], [8098,13], [8436,14], [8789,12], [10143,12], [10977,12],
            [13632,11], [14051,11], [14426,13], [14620,12], [14951,13], [16746,11],
            [18101,12], [18217,12], [18291,15], [18521,11], [18885,12], [19043,12],
            [20972,11], [22005,11], [22369,14], [23905,13], [25689,12], [27243,12],
            [32366,11], [34970,14], [35184,12], [36602,12], [37018,13], [38051,11],
            [38840,11], [40911,15], [41397,12], [41675,13], [42632,11], [43547,11],
            [43973,16], [45644,13], [46272,11], [46500,14], [47426,15], [47669,11],
        ],
        "high_hit_histogram": {"11":16,"12":17,"13":8,"14":7,"15":5,"16":1},
    },
)


def verify(witness):
    rows, regions = witness["rows"], witness["regions"]
    generator = np.random.Generator(np.random.PCG64(witness["route_seed"]))
    route = np.stack([generator.permutation(rows) for _ in range(regions)]).astype(np.int32)
    route_hash = sha256(route.astype("<i4").tobytes()).hexdigest()
    assert route_hash == witness["route_sha256"]
    starts = np.asarray(witness["starts"], dtype=np.int64)
    widths = np.asarray(witness["widths"], dtype=np.int64)
    assert len(starts) == len(widths) == regions
    assert np.all(starts >= 0) and np.all(widths > 0) and np.all(starts+widths <= rows)
    assert np.all(starts % witness["alignment"] == 0)
    assert np.all(widths % witness["alignment"] == 0)
    assert sha256(starts.astype("<i8").tobytes()).hexdigest() == witness["starts_sha256"]
    assert sha256(widths.astype("<i8").tobytes()).hexdigest() == witness["widths_sha256"]
    selected = np.concatenate([route[j, first:first+length]
                               for j, (first, length) in enumerate(zip(starts, widths))])
    counts = np.bincount(selected, minlength=rows)
    parity = witness["outer_redundancy"]
    raw_dimension = sum(max(int(c)-parity, 0) for c in counts)
    retained = rows-witness["excluded_tail_rows"]
    restricted_dimension = sum(max(int(c)-parity, 0) for c in counts[:retained])
    assert raw_dimension == witness["max_dimension_found"]
    assert restricted_dimension == witness["tail_excluded_dimension"]
    assert int(widths.sum()) == witness["total_window_symbols"] <= witness["bad_output_limit"]
    positive = [[int(i), int(counts[i])] for i in np.flatnonzero(counts > parity)]
    assert positive == witness["positive_dimension_row_hits"]
    histogram = np.bincount(counts, minlength=regions+1)
    high = {str(c): int(histogram[c]) for c in range(parity+1, regions+1) if histogram[c]}
    assert high == witness["high_hit_histogram"]
    exit_constraints = regions*witness["state"]
    assert restricted_dimension > exit_constraints
    return {"route_seed": witness["route_seed"], "verified": True,
            "dimension": raw_dimension, "tail_excluded_dimension": restricted_dimension,
            "exit_constraints_upper": exit_constraints,
            "surviving_dimension_lower": restricted_dimension-exit_constraints,
            "output_weight_upper": int(widths.sum()), "bad_output_limit": witness["bad_output_limit"],
            "positive_dimension_rows": len(positive),
            "distribution_failure_lower_bound_claimed": False}


if __name__ == "__main__":
    for fixture in WITNESSES:
        print(json.dumps(verify(fixture), sort_keys=True))
