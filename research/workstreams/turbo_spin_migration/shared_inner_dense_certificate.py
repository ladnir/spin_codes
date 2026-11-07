"""Outward attempt for Turbo's s22 inner and Scalable SPIN's bit route.

This program proves only the positive-occupation exponent inequality.
It recursively partitions the exact rational parameter domain and checks
one rational witness per accepted box with 256-bit Arb arithmetic.
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from fractions import Fraction as F
import json
import math
from pathlib import Path
import sys

from flint import arb, ctx
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[3]))

from research.workstreams.turbo_spin_migration.shared_inner_dense_cover import (
    Box, upper_for_bank,
)
from research.workstreams.turbo_spin_migration.shared_inner_dense_screen import (
    DUAL, M, OUTER, SPECTRUM, WIDTH, optimize, transfer_rate,
)


def num(value):
    value = F(value)
    return arb(value.numerator) / value.denominator


def upper(value):
    mantissa, exponent = map(int, value.upper().man_exp())
    return F(mantissa) * F(2) ** exponent


def kl(x, y):
    x, y = F(x), F(y)
    if x == 0:
        return -(1 - num(y)).log()
    if x == 1:
        return -num(y).log()
    a, b = num(x), num(y)
    return a * (a / b).log() + (1 - a) * ((1 - a) / (1 - b)).log()


def exact_segments():
    payload = json.loads(OUTER.read_text(encoding="utf-8"))
    assert payload["status"] == "proved"
    left = []
    for row in payload["left_segments"]:
        lo, hi = map(F, row["weight"])
        left.append((lo, hi, F(row["slope"]), F(row["intercept"])))
    center_lo = left[-1][1]
    center = [(center_lo, 1 - center_lo, F(0),
               F(payload["central_constant_upper"]))]
    right = [(1 - hi, 1 - lo, -slope, slope + intercept)
             for lo, hi, slope, intercept in reversed(left)]
    result = left + center + right
    assert result[0][0] == F(13, 125)
    assert result[-1][1] == F(112, 125)
    assert all(a[1] == b[0] for a, b in zip(result, result[1:]))
    return result


@dataclass(frozen=True)
class ExactBox:
    alo: F
    ahi: F
    xlo: F
    xhi: F
    slope: F
    intercept: F
    depth: int = 0

    def floating(self):
        return Box(float(self.alo), float(self.ahi), float(self.xlo),
                   float(self.xhi), float(self.slope), float(self.intercept),
                   self.depth)

    def split(self):
        if ((self.ahi - self.alo) / (self.ahi + self.alo)
                > 4 * (self.xhi - self.xlo)):
            mid = (self.alo + self.ahi) / 2
            return (ExactBox(self.alo, mid, self.xlo, self.xhi,
                             self.slope, self.intercept, self.depth + 1),
                    ExactBox(mid, self.ahi, self.xlo, self.xhi,
                             self.slope, self.intercept, self.depth + 1))
        mid = (self.xlo + self.xhi) / 2
        return (ExactBox(self.alo, self.ahi, self.xlo, mid,
                         self.slope, self.intercept, self.depth + 1),
                ExactBox(self.alo, self.ahi, mid, self.xhi,
                         self.slope, self.intercept, self.depth + 1))


def witness_bank():
    bank = []
    for a in (1e-4, .001, .01, .1, .5, 1.):
        for x in (.104, .16, .25, .4, .48, .5, .52, .6, .75, .84, .896):
            row = optimize(a, x, .11)
            bank.append((row["p"], row["y"], row["lam"]))
    for a in (.9, .99, .999, 1.):
        for x in (.485, .49, .494, .498, .502, .506, .51, .515):
            row = optimize(a, x, .11)
            bank.append((row["p"], row["y"], row["lam"]))
    floating = tuple(np.asarray(values) for values in zip(*[
        (*w, transfer_rate(w[0] * w[1], w[2])) for w in bank]))
    rational = [tuple(F(format(v, ".17g")) for v in w) for w in bank]
    return floating, rational


def outward_rate(witness):
    p, y, lam = witness
    beta, z = num(p * y), (-num(lam)).exp()
    g0 = 1 - beta + beta * z
    g1 = beta + (1 - beta) * z
    moments = {w: g0 ** (WIDTH - w) * g1**w for w in SPECTRUM if w}
    fresh = max(upper(v) for v in moments.values())
    uniform = upper(sum((count * moments[w] for w, count in SPECTRUM.items()
                         if w), arb(0)) / M)
    kernel = sum((int(DUAL[j]) * beta**j * (1 - beta)**(WIDTH-j)
                  for j in range(WIDTH + 1)), arb(0))
    zero_zero_raw = sum((int(DUAL[j]) * (beta*z)**j
                         * (1-beta)**(WIDTH-j)
                         for j in range(WIDTH + 1)), arb(0))
    zero_zero = upper(zero_zero_raw)
    zero_fresh = upper(g0**WIDTH - zero_zero_raw)
    active = max(F(0), upper(1 - kernel))
    # The entries below are rational upper bounds on the true matrix.
    matrix = [[zero_zero, max(F(0), zero_fresh), F(0)],
              [min(fresh, active) / M, F(0), fresh],
              [min(uniform, active) / M, F(0), uniform]]
    approx = np.asarray([[float(v) for v in row] for row in matrix])
    values, vectors = np.linalg.eig(approx)
    index = int(np.argmax(values.real))
    vector = np.abs(vectors[:, index].real)
    assert np.all(vector > 0)
    vector /= vector[0]
    rational_vector = [F(format(float(v), ".16g")) for v in vector]
    radius = max(sum(row[j] * rational_vector[j] for j in range(3))
                 / rational_vector[i] for i, row in enumerate(matrix))
    assert radius > 0
    return (num(radius).log() / WIDTH, lam)


def outward_box(box: ExactBox, rate, lam, witness):
    p, y, _ = witness
    values = []
    for alpha in (box.alo, box.ahi):
        for x in (box.xlo, box.xhi):
            value = (num(alpha * (box.slope * x + box.intercept))
                     + kl(alpha, p) + num(alpha) * kl(x, y)
                     + rate + num(F(11, 100) * lam))
            values.append(upper(value))
    return max(values)


def rational_hard_witness(box: ExactBox, witness):
    """Propose rational z, radius, and vector for the Lean hard-box pilot."""
    p, y, lam = witness
    decimal = 10**30
    z = F(round(upper((-num(lam)).exp()) * decimal), decimal)
    beta, zz = num(p * y), num(z)
    g0 = 1 - beta + beta * zz
    g1 = beta + (1 - beta) * zz
    moments = {w: g0 ** (WIDTH - w) * g1**w for w in SPECTRUM if w}
    fresh = max(upper(v) for v in moments.values())
    uniform = upper(sum((count * moments[w] for w, count in SPECTRUM.items()
                         if w), arb(0)) / M)
    kernel = sum((int(DUAL[j]) * beta**j * (1-beta)**(WIDTH-j)
                  for j in range(WIDTH + 1)), arb(0))
    zero_zero_raw = sum((int(DUAL[j]) * (beta*zz)**j
                         * (1-beta)**(WIDTH-j)
                         for j in range(WIDTH + 1)), arb(0))
    zero_zero = upper(zero_zero_raw)
    zero_fresh = upper(g0**WIDTH - zero_zero_raw)
    active = max(F(0), upper(1 - kernel))
    matrix = [[zero_zero, max(F(0), zero_fresh), F(0)],
              [min(fresh, active) / M, F(0), fresh],
              [min(uniform, active) / M, F(0), uniform]]
    approx = np.asarray([[float(v) for v in row] for row in matrix])
    values, vectors = np.linalg.eig(approx)
    index = int(np.argmax(values.real))
    vector = np.abs(vectors[:, index].real)
    vector /= vector[0]
    rational_vector = [F(format(float(v), ".16g")) for v in vector]
    assert all(v > 0 for v in rational_vector)
    radius = max(sum(row[j] * rational_vector[j] for j in range(3))
                 / rational_vector[i] for i, row in enumerate(matrix))
    rho = F(-(-radius.numerator * decimal // radius.denominator), decimal)
    assert rho >= radius > 0
    bound = max(upper(num(alpha * (box.slope * x + box.intercept))
                      + kl(alpha, p) + num(alpha) * kl(x, y)
                      + num(rho).log() / WIDTH
                      - num(F(11, 100)) * num(z).log())
                for alpha in (box.alo, box.ahi)
                for x in (box.xlo, box.xhi))
    return {"z": str(z), "rho": str(rho),
            "vector": [str(v) for v in rational_vector],
            "outward_upper": str(bound),
            "passes_eta_1e_7": bound <= -F(1, 10**7)}


def main():
    if not __debug__:
        raise RuntimeError("Certificate checks require Python assertions")
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--max-boxes", type=int, default=100_000)
    parser.add_argument("--max-depth", type=int, default=36)
    parser.add_argument("--proposal-margin", type=float, default=1e-7)
    parser.add_argument("--precision", type=int, default=256)
    parser.add_argument("--report-hard-box", action="store_true",
                        help="include the tightest accepted box and its proposal")
    parser.add_argument("--sample-accepted", type=int, default=0,
                        help="report this many initial accepted boxes for Lean pilots")
    args = parser.parse_args()
    assert args.precision >= 128
    ctx.prec = args.precision
    floating, witnesses = witness_bank()
    rate_bank = [outward_rate(w) for w in witnesses]
    pending = [ExactBox(F(1, 10_000), F(1), lo, hi, slope, intercept)
               for lo, hi, slope, intercept in exact_segments()]
    processed = accepted = deepest = 0
    largest_upper = None
    hard_box = None
    samples = []
    unresolved = None
    while pending and processed < args.max_boxes:
        box = pending.pop()
        processed += 1
        deepest = max(deepest, box.depth)
        guess, index = upper_for_bank(box.floating(), floating)
        if guess < -args.proposal_margin:
            rate, lam = rate_bank[index]
            bound = outward_box(box, rate, lam, witnesses[index])
            if bound < 0:
                accepted += 1
                if len(samples) < args.sample_accepted:
                    samples.append((box, witnesses[index], index))
                if largest_upper is None or bound > largest_upper:
                    largest_upper = bound
                    hard_box = (box, witnesses[index], index)
                continue
        if box.depth >= args.max_depth:
            midpoint = ((float(box.alo) + float(box.ahi)) / 2,
                        (float(box.xlo) + float(box.xhi)) / 2)
            proposal = optimize(*midpoint, .11)
            witness = tuple(F(format(proposal[key], ".17g"))
                            for key in ("p", "y", "lam"))
            rate, lam = outward_rate(witness)
            bound = outward_box(box, rate, lam, witness)
            if bound < 0:
                accepted += 1
                if len(samples) < args.sample_accepted:
                    samples.append((box, witness, len(witnesses)))
                if largest_upper is None or bound > largest_upper:
                    largest_upper = bound
                    hard_box = (box, witness, len(witnesses))
                witnesses.append(witness)
                rate_bank.append((rate, lam))
                item = (proposal["p"], proposal["y"], proposal["lam"],
                        transfer_rate(proposal["p"] * proposal["y"], proposal["lam"]))
                floating = tuple(np.append(part, value)
                                 for part, value in zip(floating, item))
                continue
            unresolved = {"alpha": list(map(str, (box.alo, box.ahi))),
                          "x": list(map(str, (box.xlo, box.xhi))),
                          "proposal": guess, "midpoint": proposal["exponent"],
                          "outward": str(bound)}
            break
        pending.extend(box.split())
        if processed % 10_000 == 0:
            print(json.dumps({"processed": processed, "accepted": accepted,
                              "pending": len(pending)}), flush=True)
    report = {"status": "proved_dense_exponent" if not pending and
                      unresolved is None else "open",
                      "precision_bits": ctx.prec, "processed": processed,
                      "accepted": accepted, "deepest": deepest,
                      "largest_outward_upper": str(largest_upper),
                       "unresolved": unresolved, "pending": len(pending)}
    if args.report_hard_box and hard_box is not None:
        box, witness, index = hard_box
        report["hard_box"] = {
            "alpha": [str(box.alo), str(box.ahi)],
            "row_density": [str(box.xlo), str(box.xhi)],
            "outer_line": [str(box.slope), str(box.intercept)],
            "witness_index": index,
            "proposal_p_y_lambda": [str(value) for value in witness],
            "rational_pilot": rational_hard_witness(box, witness),
        }
    if samples:
        report["sampled_boxes"] = [{
            "alpha": [str(box.alo), str(box.ahi)],
            "row_density": [str(box.xlo), str(box.xhi)],
            "outer_line": [str(box.slope), str(box.intercept)],
            "witness_index": index,
            "proposal_p_y_lambda": [str(value) for value in witness],
            "rational_pilot": rational_hard_witness(box, witness),
        } for box, witness, index in samples]
    print(json.dumps(report))
    if pending or unresolved is not None:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
