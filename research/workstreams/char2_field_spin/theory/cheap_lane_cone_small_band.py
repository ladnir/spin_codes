"""Serial directed verification of the fixed source-weight band7..32.

The frozen child verifies all indexed supports and all state paths at one
weight. This wrapper fixes z=1-r/1280 and alpha=1/2 for r<=15, else3/16.
It changes no encoder or child arithmetic. Each printed Arb interval is
parsed as exact decimal rationals. If its upper endpoint is u, then
k=floor(-u) certifies a probability upper bound2^-k. The wrapper adds all
26 dyadic bounds exactly and checks their JOINT40-bit gate.

Children run serially with -B, at most30seconds each and180seconds total.
Every child is hash-checked, uses the same interpreter, and inherits the
20GiB disk guard. Bounded stdout only; no generated files. A failed bound
is a completed diagnostic, not a verifier error. Missing or failed children
cannot produce a band certificate. This is not a whole-code certificate.
"""

import argparse
from fractions import Fraction
import hashlib
import json
import os
from pathlib import Path
import re
import shutil
import subprocess
import sys
import time

sys.dont_write_bytecode = True
PIN = "e32069c849efa8e3abe1d3eee0a252a386e8acb3c161725d9ec83314a7b1678b"
CHILD_BYTES = 12607
MAX_STDOUT = 128*(1 << 10)
TOTAL_SECONDS = 180
CHILD_SECONDS = 30
DECIMAL = r"[+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?"
BALL = re.compile(r"\[\s*("+DECIMAL+r")\s+\+/-\s+("+DECIMAL+r")\s*\]")
ZERO_BALL = re.compile(r"\[\s*\+/-\s+("+DECIMAL+r")\s*\]")
POINT = re.compile(DECIMAL)


def rational_decimal(value):
    if len(value) > 512 or not POINT.fullmatch(value):
        raise ValueError("invalid or oversized finite decimal")
    if "e" in value.lower() and abs(int(value.lower().split("e")[1])) > 4096:
        raise ValueError("decimal exponent exceeds parser budget")
    return Fraction(value)


def upper_endpoint(value):
    if not isinstance(value, str) or len(value) > 1100:
        raise ValueError("invalid directed interval string")
    value = value.strip()
    match = BALL.fullmatch(value)
    if match:
        center, radius = map(rational_decimal, match.groups())
    else:
        match = ZERO_BALL.fullmatch(value)
        if match:
            center, radius = Fraction(0), rational_decimal(match.group(1))
        elif POINT.fullmatch(value):
            center, radius = rational_decimal(value), Fraction(0)
        else:
            raise ValueError("unsupported finite Arb interval format")
    if radius < 0:
        raise ValueError("interval radius must be nonnegative")
    return center+radius


def dyadic(bits):
    if not -65536 <= bits <= 65536:
        raise ValueError("dyadic exponent exceeds bounded certificate budget")
    return Fraction(1, 1 << bits) if bits >= 0 else Fraction(1 << -bits)


def integer_bits(upper):
    return (-upper.numerator)//upper.denominator


def sum_bits(bound):
    if bound <= 0:
        raise ValueError("positive probability upper bound required")
    bits = bound.denominator.bit_length()-bound.numerator.bit_length()
    if bound > dyadic(bits):
        bits -= 1
    assert bound <= dyadic(bits) and bound > dyadic(bits+1)
    return bits


def parser_tests():
    cases = (("[-10.25 +/- 0.125]", Fraction(-81, 8), 10),
             ("[10 +/- 0.01]", Fraction(1001, 100), -11),
             ("[-1.25e2 +/- 1e-20]", Fraction(-125)+Fraction(1, 10**20), 124),
             ("[+/- 0.125]", Fraction(1, 8), -1),
             ("-40", Fraction(-40), 40))
    for text, exact, bits in cases:
        assert upper_endpoint(text) == exact
        assert integer_bits(exact) == bits
        assert -bits >= exact
    assert sum_bits(Fraction(1, 4)) == 2
    assert sum_bits(Fraction(3, 4)) == 0
    assert sum_bits(Fraction(3, 1 << 50)) == 48
    assert dyadic(45)*26 <= dyadic(40)
    assert dyadic(45)*26 > dyadic(41)
    return len(cases)


def extract_result(stdout):
    """The pinned child emits compact phase JSON, then one pretty JSON object."""
    if len(stdout.encode("utf-8")) > MAX_STDOUT:
        raise RuntimeError("child stdout exceeded128KiB")
    decoder = json.JSONDecoder()
    position, objects = 0, []
    while position < len(stdout):
        while position < len(stdout) and stdout[position].isspace():
            position += 1
        if position == len(stdout):
            break
        value, position = decoder.raw_decode(stdout, position)
        if not isinstance(value, dict) or len(objects) >= 10:
            raise RuntimeError("unexpected child JSON stream")
        objects.append(value)
    if not objects or "verifier_completed" not in objects[-1]:
        raise RuntimeError("missing final verifier result")
    if any("phase" not in value for value in objects[:-1]):
        raise RuntimeError("unexpected pre-result child object")
    return objects[-1]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--precision", type=int, choices=(256, 384, 512), default=256)
    args = parser.parse_args()
    started = time.monotonic()
    deadline = started+TOTAL_SECONDS
    path = Path(__file__).with_name("cheap_lane_cone_weight_witness.py").resolve()

    def preflight():
        if time.monotonic() >= deadline:
            raise TimeoutError("180-second serial band budget exhausted")
        if shutil.disk_usage(path.anchor).free < 20*(1 << 30):
            raise RuntimeError("requires20GiB free on source volume")
        contents = path.read_bytes()
        if len(contents) != CHILD_BYTES or hashlib.sha256(contents).hexdigest() != PIN:
            raise RuntimeError("frozen single-weight verifier hash/size mismatch")

    preflight()
    tests = parser_tests()
    environment = os.environ.copy()
    environment["PYTHONDONTWRITEBYTECODE"] = "1"
    environment["PYTHONUNBUFFERED"] = "1"
    results, joint = [], Fraction(0)
    for weight in range(7, 33):
        preflight()
        z = 1-Fraction(weight, 1280)
        alpha = Fraction(1, 2) if weight <= 15 else Fraction(3, 16)
        command = [sys.executable, "-B", str(path), "--weight", str(weight),
                   "--z", str(z), "--alpha", str(alpha),
                   "--precision", str(args.precision)]
        timeout = min(CHILD_SECONDS, deadline-time.monotonic())
        if timeout <= 0:
            raise TimeoutError("serial band budget exhausted before child")
        # The pinned source has bounded prints: <=4phase objects and its final
        # result. No descendant process is launched by that source. run() kills
        # and waits for this child on timeout; another child is never concurrent.
        child = subprocess.run(command, cwd=path.parent, env=environment,
                               stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                               text=True, encoding="utf-8", timeout=timeout)
        if child.returncode != 0:
            raise RuntimeError("child failed at weight%d: %s" %
                               (weight, child.stderr[-2000:]))
        if child.stderr or len(child.stdout.encode("utf-8")) > MAX_STDOUT:
            raise RuntimeError("unexpected stderr or oversized child output")
        result = extract_result(child.stdout)
        parameters = result.get("parameters", {})
        if (result.get("directed") is not True or result.get("verifier_completed") is not True
                or result.get("source_sha256") != PIN
                or result.get("generated_bytes") != 0
                or result.get("precision_bits") != args.precision
                or parameters.get("source_weight") != weight
                or parameters.get("z") != str(z) or parameters.get("alpha") != str(alpha)
                or parameters.get("distance_threshold") != 734003
                or result.get("tests", {}).get("source_count_sum") != "C(8192,r)"
                or result.get("tests", {}).get("per_run_totals") != "passed"
                or result.get("tests", {}).get("terminal_updates") != 32767):
            raise RuntimeError("child result does not match the prescribed witness")
        interval = result["log2_bound_interval"]
        bits = integer_bits(upper_endpoint(interval))
        joint += dyadic(bits)
        results.append({"weight": weight, "z": str(z), "alpha": str(alpha),
                        "integer_bits": bits, "log2_bound_interval": interval})
        print(json.dumps({"phase": "weight_complete", "weight": weight,
                          "integer_bits": bits, "seconds": time.monotonic()-started}), flush=True)
    preflight()
    assert [item["weight"] for item in results] == list(range(7, 33))
    denominator_power = joint.denominator.bit_length()-1
    assert joint.denominator == 1 << denominator_power
    certificate = {"numerator_bits": joint.numerator.bit_length(),
                   "denominator_power_of_two": denominator_power,
                   "integer_bits": sum_bits(joint)}
    if joint.numerator.bit_length() <= 4096:
        certificate["numerator_exact"] = str(joint.numerator)
    print(json.dumps({"scope": "all state paths and all indexed supports with weights7..32",
        "directed": True, "verifier_completed": True, "whole_code_certificate": False,
        "precision_bits": args.precision, "weights": results,
        "joint_probability_upper": certificate,
        "proved_joint_below_2_pow_minus_40": joint < dyadic(40),
        "proved_joint_below_2_pow_minus_100": joint < dyadic(100),
        "rounding": "exact decimal upper endpoint u; sum2^(-floor(-u))",
        "tests": {"exact_interval_parser_cases": tests, "covered_weights": 26,
                  "serial_children": True, "every_child_count_and_terminal_check": "passed"},
        "seconds": time.monotonic()-started, "generated_bytes": 0,
        "child_timeout_seconds": CHILD_SECONDS, "total_timeout_seconds": TOTAL_SECONDS,
        "child_pin": PIN, "source_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest()},
        indent=2), flush=True)


if __name__ == "__main__":
    main()
