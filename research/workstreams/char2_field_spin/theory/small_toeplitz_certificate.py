"""Standalone directed replay for independent Toeplitz[48,24]/GL3.

The source size is rounded up to a multiple of 24*16. One uniform regional
route, the monomial width-16 inner, and the no-diagonal source law are the
same as in K20_RANK_CONTINUATION.md. Existing engines and witnesses remain
unchanged. Search proposes witnesses; only directed replay certifies them.
"""

import argparse
from decimal import Decimal
from fractions import Fraction
import json
import math

from gf16_sizes import interval_local_envelope, upper_regions
from rank_continuation_toeplitz_certificate import CONFIG, rational, terminal_power
from subfield_certificate import HIGH, LOW, Interval, integer_margin


FORMAT = "gf16-systematic-toeplitz48-gl3-nodiagonal-v1"


def configuration(k_log):
    if k_log not in (16, 18, 20):
        raise ValueError("supported nominal source sizes are K16, K18, and K20")
    k = (((1 << k_log)+383)//384)*384
    return dict(CONFIG, nominal_k=1 << k_log, k=k, n=2*k, outer=[48, 24],
                rows=k//24, regions=48, bad_weight_threshold=(2*k+9)//10-1)


def positive_vector_product(vector, matrix):
    add, mul = HIGH.add, HIGH.multiply
    return [add(mul(vector[0], matrix[0]), mul(vector[1], matrix[2])),
            add(mul(vector[0], matrix[1]), mul(vector[1], matrix[3]))]


def nonempty_systematic_moment(regions, matrix):
    """Positive flag DP retains u != 0 for the one-active-row class."""
    empty, seen = [Decimal(1), Decimal(0)], [Decimal(0), Decimal(0)]
    active_matrix = [HIGH.multiply(Decimal(15), value) for value in regions[1]]
    for _ in range(24):
        continued = positive_vector_product(seen, matrix)
        first = positive_vector_product(empty, active_matrix)
        seen = [HIGH.add(a, b) for a, b in zip(continued, first)]
        empty = positive_vector_product(empty, regions[0])
    for _ in range(24):
        seen = positive_vector_product(seen, matrix)
    return HIGH.add(*seen)


def sparse_upper(config, active, z):
    if not isinstance(active, int) or not 1 <= active < config["rows"]:
        raise ValueError("sparse class outside supported range")
    zi = rational(z)
    if not 0 < zi.lo <= zi.hi < 1:
        raise ValueError("require 0 < z < 1")
    regions = upper_regions(zi, active, config["rows"], 3)
    matrix = [Decimal(0)]*4
    for h, regional in enumerate(regions):
        weight = Decimal(math.comb(active, h)*15**h)
        matrix = [HIGH.add(a, HIGH.multiply(weight, b)) for a, b in zip(matrix, regional)]
    moment = nonempty_systematic_moment(regions, matrix) if active == 1 else terminal_power(matrix, 48)
    numerator = HIGH.multiply(Decimal(math.comb(config["rows"], active)), moment)
    denominator = LOW.multiply(Decimal(15*16**(24*active)), (zi**config["bad_weight_threshold"]).lo)
    return HIGH.divide(numerator, denominator)


def middle_upper(config, active, x, z):
    if not isinstance(active, int) or not 1 <= active <= config["rows"] or Fraction(x) <= 0:
        raise ValueError("require a valid active class and x > 0")
    xi, zi = rational(x), rational(z)
    if not 0 < zi.lo <= zi.hi < 1:
        raise ValueError("require 0 < z < 1")
    local = interval_local_envelope(zi, 3)
    matrix = [Decimal(0)]*4
    for w, entries in enumerate(local):
        weight = math.comb(16, w)*(15*xi)**w*(1+xi)**(16-w)
        matrix = [HIGH.add(a, (weight*b).hi) for a, b in zip(matrix, [v for row in entries for v in row])]
    moment = terminal_power(matrix, config["rows"]*3)
    scalar = rational("1/15")*rational("1/16")**(24*active)
    scalar = scalar/(zi**config["bad_weight_threshold"])/(xi**(48*active))
    scalar = scalar/Interval.integer(math.comb(config["rows"], active))**47
    return HIGH.multiply(moment, scalar.hi)


def propose(k_log):
    from small_toeplitz_screen import search, middle_search, middle_screen
    config, sparse = configuration(k_log), []
    for active in range(1, 19):
        candidate = search(k_log, active, 3, 24, 16)
        sparse.append(dict(active=active, z=str(Fraction(round(candidate["z"]*2**40), 2**40))))
    middle, first = [], 19
    while first <= config["rows"]:
        candidate = middle_search(k_log, first, 3, 24, 16)
        x = Fraction(round(candidate["x"]*10**12), 10**12)
        z = Fraction(round(candidate["z"]*10**12), 10**12)
        if middle_screen(k_log, first, 3, 24, 16, float(x), float(z)) < 90:
            raise ValueError(f"no 90-bit interval witness at {first}")
        low, high = first, config["rows"]
        while low < high:
            probe = (low+high+1)//2
            if middle_screen(k_log, probe, 3, 24, 16, float(x), float(z)) >= 90:
                low = probe
            else:
                high = probe-1
        middle.append(dict(first=first, last=low, x=str(x), z=str(z)))
        first = low+1
    return dict(config=config, sparse=sparse, middle=middle)


def certify(witnesses):
    config = witnesses["config"]
    nominal = config.get("nominal_k", 0)
    if not isinstance(nominal, int) or config != configuration(nominal.bit_length()-1):
        raise ValueError("incorrect construction or rounded size")
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
        first, last = entry["first"], entry["last"]
        if (not isinstance(first, int) or not isinstance(last, int)
                or first != next_active or not first <= last <= config["rows"]):
            raise ValueError("middle coverage gap, overlap, or invalid endpoint")
        endpoints = [middle_upper(config, a, entry["x"], entry["z"]) for a in (first, last)]
        contribution = HIGH.multiply(max(endpoints), Decimal(last-first+1))
        total = HIGH.add(total, contribution)
        middle.append(dict(entry, endpoint_upper=[str(v) for v in endpoints],
                           failure_upper=str(contribution), integer_bits=integer_margin(contribution)))
        next_active = last+1
    if next_active != config["rows"]+1:
        raise ValueError("incomplete coverage")
    return dict(format=FORMAT, config=config, complete=True, sparse=sparse, middle=middle,
                sparse_upper=str(sparse_total), failure_upper=str(total), integer_bits=integer_margin(total))


def verify(record):
    if record.get("format") != FORMAT or record.get("complete") is not True:
        raise ValueError("incorrect format or completeness claim")
    fresh = certify(record)
    for field in ("sparse_upper", "failure_upper"):
        if Fraction(record[field]) < Fraction(fresh[field]):
            raise ValueError(f"understated {field}")
    for supplied, actual in [(record, fresh)]+[(a, b) for kind in ("sparse", "middle")
                                              for a, b in zip(record[kind], fresh[kind])]:
        if Fraction(supplied["failure_upper"]) < Fraction(actual["failure_upper"]):
            raise ValueError("understated contribution")
        if (type(supplied["integer_bits"]) is not int or supplied["integer_bits"] > actual["integer_bits"]
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
    if Fraction(record["sparse_upper"]) < sum(Fraction(e["failure_upper"]) for e in record["sparse"]):
        raise ValueError("sparse total inconsistent with displayed contributions")
    if Fraction(record["failure_upper"]) < sum(Fraction(e["failure_upper"])
                                              for kind in ("sparse", "middle") for e in record[kind]):
        raise ValueError("total inconsistent with displayed contributions")
    return dict(verified=True, k=fresh["config"]["k"], integer_bits=fresh["integer_bits"],
                failure_upper=fresh["failure_upper"])


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--k-log", type=int, choices=(16, 18, 20), default=20)
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
        from small_toeplitz_witnesses import WITNESSES
        print(json.dumps(certify(WITNESSES[args.k_log]), indent=2))
    else:
        parser.error("select --propose, --certify, or --verify")
