"""Separate size replays for the same no-diagonal Toeplitz64/GL3 ensemble.

The fixed K20 engine and witness interface remain unchanged. This wrapper
changes only K, N, regional length, bad-weight threshold, and witnesses.
It also accepts K20 to check equality with that independent fixed replay.
"""

import argparse
from decimal import Decimal
from fractions import Fraction
import json
import math

from gf16_sizes import interval_local_envelope, upper_regions
from rank_continuation_toeplitz_certificate import CONFIG, rational, terminal_power
from subfield_certificate import HIGH, LOW, Interval, integer_margin


FORMAT = "gf16-systematic-toeplitz64-gl3-nodiagonal-size-v1"


def configuration(k_log):
    if k_log not in (16, 18, 20):
        raise ValueError("supported source sizes are K16, K18, and K20")
    k = 1 << k_log
    return dict(CONFIG, k=k, n=2*k, rows=k//32, bad_weight_threshold=(2*k+9)//10-1)


def sparse_upper(config, active, z):
    if not 1 <= active < config["rows"]:
        raise ValueError("sparse class outside supported range")
    zi = rational(z)
    regions = upper_regions(zi, active, config["rows"], config["state"])
    matrix = [Decimal(0)] * 4
    for h, regional in enumerate(regions):
        weight = Decimal(math.comb(active, h) * 15**h)
        matrix = [HIGH.add(a, HIGH.multiply(weight, b)) for a, b in zip(matrix, regional)]
    numerator = HIGH.multiply(Decimal(math.comb(config["rows"], active)), terminal_power(matrix, 64))
    denominator = LOW.multiply(Decimal(15*16**(32*active)), (zi**config["bad_weight_threshold"]).lo)
    return HIGH.divide(numerator, denominator)


def middle_upper(config, active, x, z):
    if not 1 <= active <= config["rows"] or Fraction(x) <= 0:
        raise ValueError("require a valid active class and x > 0")
    xi, zi = rational(x), rational(z)
    local = interval_local_envelope(zi, config["state"])
    matrix = [Decimal(0)] * 4
    for w, entries in enumerate(local):
        weight = math.comb(16, w) * (15*xi)**w * (1+xi)**(16-w)
        matrix = [HIGH.add(a, (weight*b).hi) for a, b in zip(matrix, [v for row in entries for v in row])]
    moment = terminal_power(matrix, config["rows"]*4)
    scalar = rational("1/15") * rational("1/16")**(32*active)
    scalar = scalar / (zi**config["bad_weight_threshold"]) / (xi**(64*active))
    scalar = scalar / Interval.integer(math.comb(config["rows"], active))**63
    return HIGH.multiply(moment, scalar.hi)


def propose(k_log):
    from rank_continuation_random_outer import search, middle_search, middle_screen
    config = configuration(k_log)
    sparse = []
    for active in range(1, 16):
        candidate = search(k_log, active, 3, 32)
        z = Fraction(round(candidate["z"]*2**40), 2**40)
        sparse.append({"active": active, "z": str(z)})
    middle, first = [], 16
    while first <= config["rows"]:
        candidate = middle_search(k_log, first, 3, 32)
        x = Fraction(round(candidate["x"]*10**12), 10**12)
        z = Fraction(round(candidate["z"]*10**12), 10**12)
        if middle_screen(k_log, first, 3, 32, float(x), float(z)) < 110:
            raise ValueError(f"no 110-bit endpoint witness at {first}")
        low, high = first, config["rows"]
        while low < high:
            probe = (low+high+1)//2
            if middle_screen(k_log, probe, 3, 32, float(x), float(z)) >= 110:
                low = probe
            else:
                high = probe-1
        middle.append({"first": first, "last": low, "x": str(x), "z": str(z)})
        first = low+1
    return {"config": config, "sparse": sparse, "middle": middle}


def certify(k_log, witnesses):
    config = configuration(k_log)
    if witnesses["config"] != config:
        raise ValueError("incorrect construction or size")
    previous = [(context.Emax, context.Emin) for context in (HIGH, LOW)]
    try:
        HIGH.Emax = LOW.Emax = 999999999
        HIGH.Emin = LOW.Emin = -999999999
        return _certify(config, witnesses)
    finally:
        for context, (maximum, minimum) in zip((HIGH, LOW), previous):
            context.Emax, context.Emin = maximum, minimum


def _certify(config, witnesses):
    next_active, total = 1, Decimal(0)
    sparse, middle = [], []
    for entry in witnesses["sparse"]:
        if entry["active"] != next_active:
            raise ValueError("sparse coverage gap or overlap")
        bound = sparse_upper(config, next_active, entry["z"])
        total = HIGH.add(total, bound)
        sparse.append(dict(entry, failure_upper=str(bound), integer_bits=integer_margin(bound)))
        next_active += 1
    sparse_total = total
    for entry in witnesses["middle"]:
        if entry["first"] != next_active or not entry["first"] <= entry["last"] <= config["rows"]:
            raise ValueError("middle coverage gap, overlap, or invalid endpoint")
        endpoints = [middle_upper(config, a, entry["x"], entry["z"])
                     for a in (entry["first"], entry["last"])]
        contribution = HIGH.multiply(max(endpoints), Decimal(entry["last"]-entry["first"]+1))
        total = HIGH.add(total, contribution)
        middle.append(dict(entry, endpoint_upper=[str(v) for v in endpoints],
                           failure_upper=str(contribution), integer_bits=integer_margin(contribution)))
        next_active = entry["last"]+1
    if next_active != config["rows"]+1:
        raise ValueError("incomplete coverage")
    return {"format": FORMAT, "config": config, "complete": True,
            "sparse": sparse, "middle": middle, "sparse_upper": str(sparse_total),
            "failure_upper": str(total), "integer_bits": integer_margin(total)}


def verify(record):
    if record.get("format") != FORMAT or record.get("complete") is not True:
        raise ValueError("incorrect format or completeness claim")
    fresh = certify(record["config"]["k"].bit_length()-1, record)
    for field in ("sparse_upper", "failure_upper"):
        if Fraction(record[field]) < Fraction(fresh[field]):
            raise ValueError(f"understated {field}")
    for supplied, actual in [(record, fresh)] + [(a, b) for kind in ("sparse", "middle")
                                                for a, b in zip(record[kind], fresh[kind])]:
        if Fraction(supplied["failure_upper"]) < Fraction(actual["failure_upper"]):
            raise ValueError("understated contribution")
        if (not isinstance(supplied["integer_bits"], int)
                or supplied["integer_bits"] > actual["integer_bits"]
                or supplied["integer_bits"] > integer_margin(Decimal(supplied["failure_upper"]))):
            raise ValueError("overstated or inconsistent margin")
    for supplied, actual in zip(record["middle"], fresh["middle"]):
        if len(supplied["endpoint_upper"]) != 2:
            raise ValueError("missing endpoint bounds")
        for bound, replayed in zip(supplied["endpoint_upper"], actual["endpoint_upper"]):
            if Fraction(bound) < Fraction(replayed):
                raise ValueError("understated endpoint")
        contribution = max(Fraction(v) for v in supplied["endpoint_upper"])*(supplied["last"]-supplied["first"]+1)
        if Fraction(supplied["failure_upper"]) < contribution:
            raise ValueError("interval contribution inconsistent with displayed endpoints")
    if Fraction(record["sparse_upper"]) < sum(Fraction(entry["failure_upper"]) for entry in record["sparse"]):
        raise ValueError("sparse total inconsistent with displayed contributions")
    if Fraction(record["failure_upper"]) < sum(Fraction(entry["failure_upper"])
                                                for kind in ("sparse", "middle") for entry in record[kind]):
        raise ValueError("total inconsistent with displayed contributions")
    return {"verified": True, "k": fresh["config"]["k"], "integer_bits": fresh["integer_bits"],
            "failure_upper": fresh["failure_upper"]}


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--k-log", type=int, choices=(16, 18, 20), default=16)
    parser.add_argument("--propose", action="store_true")
    parser.add_argument("--certify", action="store_true")
    parser.add_argument("--verify")
    args = parser.parse_args()
    if args.verify:
        with open(args.verify, encoding="utf-8") as source:
            print(json.dumps(verify(json.load(source)), indent=2))
    elif args.propose:
        print(json.dumps(propose(args.k_log), indent=2))
    elif args.certify:
        if args.k_log == 20:
            from rank_continuation_toeplitz_witnesses import WITNESSES
            witnesses = WITNESSES
        else:
            from rank_continuation_toeplitz_size_witnesses import WITNESSES
            witnesses = WITNESSES[args.k_log]
        print(json.dumps(certify(args.k_log, witnesses), indent=2))
    else:
        parser.error("select --propose, --certify, or --verify")
