"""Directed replay for a new, independent-per-row Toeplitz outer ensemble.

The local GF16 inner transfer is reused from gf16_sizes. This module changes
the outer to systematic [64,32] random Toeplitz blocks and uses 64 regions.
There is no input diagonal: averaging the complete source counting measure
over independent Toeplitz blocks supplies the local independent-value law.
Floating search only proposes rational witnesses; replay uses their exact
fractions and positive directed arithmetic. No native timings occur here.
"""

import argparse
from decimal import Decimal
from fractions import Fraction
import json
import math

from gf16_sizes import interval_local_envelope, upper_regions
from subfield_certificate import HIGH, LOW, Interval, integer_margin


CONFIG = {"k": 2**20, "n": 2**21, "q": 16, "outer": [64, 32],
          "outer_distribution": "systematic-independent-uniform-Toeplitz-per-row",
          "rows": 32768, "regions": 64, "width": 16, "state": 3,
          "route": "independent-uniform-per-region",
          "diagonal": "none",
          "refresh": "independent-uniform-GL3-pre-output",
          "maps": "monomial-Vandermonde-transpose", "initial_state": "zero",
          "final_state": "discarded", "bad_weight_threshold": 209715}


def rational(value):
    value = Fraction(value)
    return Interval.rational(value.numerator, value.denominator)


def positive_product(a, b):
    add, mul = HIGH.add, HIGH.multiply
    return [add(mul(a[0], b[0]), mul(a[1], b[2])),
            add(mul(a[0], b[1]), mul(a[1], b[3])),
            add(mul(a[2], b[0]), mul(a[3], b[2])),
            add(mul(a[2], b[1]), mul(a[3], b[3]))]


def terminal_power(matrix, exponent):
    result = [Decimal(1), Decimal(0), Decimal(0), Decimal(1)]
    while exponent:
        if exponent & 1:
            result = positive_product(result, matrix)
        exponent >>= 1
        if exponent:
            matrix = positive_product(matrix, matrix)
    return HIGH.add(result[0], result[1])


def sparse_upper(active, z):
    if not 1 <= active < CONFIG["rows"]:
        raise ValueError("sparse class outside supported range")
    zi = rational(z)
    regions = upper_regions(zi, active, CONFIG["rows"], CONFIG["state"])
    matrix = [Decimal(0)] * 4
    for h, regional in enumerate(regions):
        weight = Decimal(math.comb(active, h) * 15**h)
        matrix = [HIGH.add(a, HIGH.multiply(weight, b)) for a, b in zip(matrix, regional)]
    numerator = HIGH.multiply(Decimal(math.comb(CONFIG["rows"], active)), terminal_power(matrix, 64))
    denominator = LOW.multiply(Decimal(15 * 16**(32*active)), (zi**209715).lo)
    return HIGH.divide(numerator, denominator)


def middle_upper(active, x, z):
    if not 1 <= active <= CONFIG["rows"] or Fraction(x) <= 0:
        raise ValueError("require a valid active class and x > 0")
    xi, zi = rational(x), rational(z)
    local = interval_local_envelope(zi, CONFIG["state"])
    matrix = [Decimal(0)] * 4
    for w, entries in enumerate(local):
        weight = math.comb(16, w) * (15*xi)**w * (1+xi)**(16-w)
        matrix = [HIGH.add(a, (weight*b).hi) for a, b in zip(matrix, [v for row in entries for v in row])]
    moment = terminal_power(matrix, CONFIG["rows"]*4)
    scalar = rational("1/15")
    scalar = scalar * rational("1/16")**(32*active) / (zi**209715) / (xi**(64*active))
    scalar = scalar / Interval.integer(math.comb(CONFIG["rows"], active))**63
    return HIGH.multiply(moment, scalar.hi)


def propose():
    from rank_continuation_random_outer import search, middle_search, middle_screen
    sparse = []
    for active in range(1, 16):
        candidate = search(20, active, CONFIG["state"], 32)
        z = Fraction(round(candidate["z"] * 2**40), 2**40)
        sparse.append({"active": active, "z": str(z)})
    middle = []
    first = 16
    while first <= CONFIG["rows"]:
        candidate = middle_search(20, first, CONFIG["state"], 32)
        x = Fraction(round(candidate["x"] * 10**12), 10**12)
        z = Fraction(round(candidate["z"] * 10**12), 10**12)
        bits = middle_screen(20, first, CONFIG["state"], 32, float(x), float(z))
        if bits < 110:
            raise ValueError(f"no 110-bit interval witness at {first}: {bits}")
        low, high = first, CONFIG["rows"]
        while low < high:
            probe = (low + high + 1) // 2
            if middle_screen(20, probe, CONFIG["state"], 32, float(x), float(z)) >= 110:
                low = probe
            else:
                high = probe - 1
        middle.append({"first": first, "last": low, "x": str(x), "z": str(z)})
        first = low + 1
    return {"config": CONFIG, "sparse": sparse, "middle": middle}


def certify(witnesses):
    # Large positive matrix powers cancel with small exact rational factors.
    # Widen the exponent range; keep the 80-digit directed precision.
    previous = [(context.Emax, context.Emin) for context in (HIGH, LOW)]
    try:
        HIGH.Emax = LOW.Emax = 999999999
        HIGH.Emin = LOW.Emin = -999999999
        return _certify(witnesses)
    finally:
        for context, (maximum, minimum) in zip((HIGH, LOW), previous):
            context.Emax, context.Emin = maximum, minimum


def _certify(witnesses):
    if witnesses["config"] != CONFIG:
        raise ValueError("incorrect construction")
    next_active, total = 1, Decimal(0)
    sparse, middle = [], []
    for entry in witnesses["sparse"]:
        if entry["active"] != next_active:
            raise ValueError("sparse coverage gap or overlap")
        bound = sparse_upper(next_active, entry["z"])
        total = HIGH.add(total, bound)
        sparse.append(dict(entry, failure_upper=str(bound), integer_bits=integer_margin(bound)))
        next_active += 1
    sparse_total = total
    for entry in witnesses["middle"]:
        if entry["first"] != next_active or not entry["first"] <= entry["last"] <= CONFIG["rows"]:
            raise ValueError("middle coverage gap, overlap, or invalid endpoint")
        endpoints = [middle_upper(a, entry["x"], entry["z"])
                     for a in (entry["first"], entry["last"])]
        contribution = HIGH.multiply(max(endpoints), Decimal(entry["last"] - entry["first"] + 1))
        total = HIGH.add(total, contribution)
        middle.append(dict(entry, endpoint_upper=[str(v) for v in endpoints],
                           failure_upper=str(contribution), integer_bits=integer_margin(contribution)))
        next_active = entry["last"] + 1
    if next_active != CONFIG["rows"] + 1:
        raise ValueError("incomplete coverage")
    return {"format": "gf16-systematic-toeplitz64-gl3-nodiagonal-v1", "config": CONFIG,
            "complete": True, "sparse": sparse, "middle": middle,
            "sparse_upper": str(sparse_total), "failure_upper": str(total),
            "integer_bits": integer_margin(total)}


def verify(record):
    if record.get("format") != "gf16-systematic-toeplitz64-gl3-nodiagonal-v1" or record.get("complete") is not True:
        raise ValueError("incorrect format or completeness claim")
    fresh = certify(record)
    for field in ("sparse_upper", "failure_upper"):
        if Fraction(record[field]) < Fraction(fresh[field]):
            raise ValueError(f"understated {field}")
    if not isinstance(record["integer_bits"], int) or record["integer_bits"] > fresh["integer_bits"]:
        raise ValueError("overstated total margin")
    if record["integer_bits"] > integer_margin(Decimal(record["failure_upper"])):
        raise ValueError("total margin inconsistent with displayed upper bound")
    for kind in ("sparse", "middle"):
        for supplied, actual in zip(record[kind], fresh[kind]):
            if Fraction(supplied["failure_upper"]) < Fraction(actual["failure_upper"]):
                raise ValueError("understated class or interval contribution")
            if not isinstance(supplied["integer_bits"], int) or supplied["integer_bits"] > actual["integer_bits"]:
                raise ValueError("overstated class or interval margin")
            if supplied["integer_bits"] > integer_margin(Decimal(supplied["failure_upper"])):
                raise ValueError("class or interval margin inconsistent with displayed upper bound")
            if kind == "middle":
                if len(supplied["endpoint_upper"]) != 2:
                    raise ValueError("missing endpoint bounds")
                for bound, replayed in zip(supplied["endpoint_upper"], actual["endpoint_upper"]):
                    if Fraction(bound) < Fraction(replayed):
                        raise ValueError("understated endpoint")
                claimed = max(Fraction(v) for v in supplied["endpoint_upper"]) * (supplied["last"] - supplied["first"] + 1)
                if Fraction(supplied["failure_upper"]) < claimed:
                    raise ValueError("interval contribution inconsistent with displayed endpoints")
    if Fraction(record["sparse_upper"]) < sum(Fraction(entry["failure_upper"]) for entry in record["sparse"]):
        raise ValueError("sparse total inconsistent with displayed contributions")
    if Fraction(record["failure_upper"]) < sum(Fraction(entry["failure_upper"])
                                                for kind in ("sparse", "middle") for entry in record[kind]):
        raise ValueError("total inconsistent with displayed contributions")
    return {"verified": True, "integer_bits": fresh["integer_bits"], "failure_upper": fresh["failure_upper"]}


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--propose", action="store_true")
    parser.add_argument("--certify", action="store_true")
    parser.add_argument("--witnesses")
    parser.add_argument("--verify")
    args = parser.parse_args()
    if args.verify:
        with open(args.verify, encoding="utf-8") as source:
            print(json.dumps(verify(json.load(source)), indent=2))
    elif args.propose:
        print(json.dumps(propose(), indent=2))
    elif args.certify:
        if args.witnesses:
            with open(args.witnesses, encoding="utf-8") as source:
                witnesses = json.load(source)
        else:
            from rank_continuation_toeplitz_witnesses import WITNESSES
            witnesses = WITNESSES
        print(json.dumps(certify(witnesses), indent=2))
    else:
        parser.error("select --propose or --certify")
