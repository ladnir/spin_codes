"""Exact-point tests and replay for the frozen directed cone arithmetic.

Arb.contains(fmpq) may test containment of a rounded enclosing ball, not
just the rational point. At precision256 the frozen toy test rejects the
valid enclosure of4/25. Check exact rational endpoints instead; keep the
entire arithmetic kernel unchanged. This loader replaces only toy_test.
The frozen main retains its48-second,20GiB,stdout-only resource guards.
"""

from collections import defaultdict
from fractions import Fraction as F
import hashlib
import importlib.util
import json
from pathlib import Path
import sys

sys.dont_write_bytecode = True
PIN = "7ccf0cdab823ef15fdd5edd6156d4b49d9a000f2ae28c37d5d43f4b6ea341f01"


def exact_point(ball, fraction, fmpq):
    point = fmpq(fraction.numerator, fraction.denominator)
    assert ball.is_finite()
    assert ball.lower().fmpq() <= point <= ball.upper().fmpq()


def point_tests(self):
    raw = {(0, 0): F(1, 3), (1, 0): F(2, 5), (0, 1): F(3, 7)}
    entry = [self.poly([self.exact(raw[0, 0]), self.exact(raw[1, 0])]),
             self.poly([self.exact(raw[0, 1])])]
    entry += [self.poly([]) for _ in range(self.single-1)]
    actual = self.multiply([[entry]], [[entry]], 2)[0][0]
    expected = defaultdict(F)
    for (f, s), x in raw.items():
        for (ff, ss), y in raw.items():
            if s+ss <= self.single:
                expected[f+ff, s+ss] += x*y
    for (f, s), value in expected.items():
        exact_point(actual[s][f], value, self.fmpq)
    x, y = self.exact(F(1, 3)), self.exact(F(2, 3))
    exact_point(x.min(y), F(1, 3), self.fmpq)
    exact_point(x.max(y), F(2, 3), self.fmpq)
    return len(expected)


def load_verifier():
    path = Path(__file__).with_name("cheap_lane_cone_verify.py")
    if hashlib.sha256(path.read_bytes()).hexdigest() != PIN:
        raise RuntimeError("frozen cone verifier hash mismatch")
    spec = importlib.util.spec_from_file_location("frozen_directed_cone", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)

    class ExactPointVerifier(module.Verifier):
        toy_test = point_tests

    module.Verifier = ExactPointVerifier
    return module


if __name__ == "__main__":
    print(json.dumps({"replay_test_method": "exact rational endpoint comparisons",
        "arithmetic_kernel_sha256": PIN,
        "test_wrapper_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}), flush=True)
    load_verifier().main()
