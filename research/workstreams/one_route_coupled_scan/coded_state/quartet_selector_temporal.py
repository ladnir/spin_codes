"""Setup-only temporal selector moment for the unchanged quartet code.

For a fixed eligible row set H, let j be a quartet's H-label count,
P_r the union of quartets with j>=r, and Z_r=sum(j-r)_+.  Let E_r count
coordinate edges joining P_r at adjacent updates.  Every actual profile
with G=F-f and C_out=R-128h obeys

  E_r+8 Z_r >= Ls-512-f+8*((1-r/4)*G-C_out).

Indeed, the selector slack Z_r-((1-r/4)*G-C_out) pays at least one for
each selected quartet outside P_r; such a node touches at most eight
temporal edges.  States outside forced births give at least Ls-512-f
selected-selected edges.  This is simultaneous over all state histories.

Under iid H-slot labels of probability p, the weighted distribution of
|P_r|/4=a is w_a=C(m,a) q0^(m-a) q1^a, where m=width/4 and
q1=sum_{j>=r}Bin(4,p;j)*beta^(8*(j-r)); q0 is the unweighted j<r sum.
Conditional on a, P_r is a uniform 4a-subset.  Thus with

  H_ab=E beta^|U intersect V|
      =sum_j C(4a,j) C(4b,j) (beta-1)^j/C(width,j),

the exact setup moment is e0*(H*diag(w))^updates*1.  The artificial
initial P=empty adds a zero edge weight; the actual terminal is free.
Conditioning every region to exactly h H labels costs pi_h(p)^(-256).
The H-set union costs C(8192,h) once.  No support-history union appears.

The production evaluation is a FLOAT diagnostic, including when positive
matrix powers finish.  It is not a directed-rounding certificate.  Exact
small checks independently enumerate partitions, labels, and intersections.
Only the parent runs jobs.  No files, subprocesses, or witness dumps.
"""

import argparse
from fractions import Fraction
from itertools import combinations, product
from math import comb, isfinite, lgamma, log2
from operator import index
import os
import shutil
import threading
import time

import numpy as np


ST_WIDTH = 512
ST_UPDATES = 4096
ST_ROWS = 8192
ST_REGIONS = 256
ST_LD = np.longdouble
ST_LN2 = np.log(ST_LD(2))


def _st_logchoose(n, k):
    if not 0 <= k <= n:
        return float('-inf')
    if k == 0 or k == n:
        return 0.0
    return (lgamma(n + 1) - lgamma(k + 1) - lgamma(n - k + 1)) / float(ST_LN2)


def _st_normalize(array):
    maximum = np.max(array)
    if not np.isfinite(maximum) or maximum <= 0 or np.any(array <= 0):
        raise ArithmeticError('A positive temporal action became invalid')
    return array / maximum, np.log2(maximum)


def temporal_selector_overlap(log_beta, width=ST_WIDTH):
    """Positive Gram construction of every quartet-population overlap mgf."""
    width = index(width)
    log_beta = ST_LD(log_beta)
    if width < 4 or width % 4 or not np.isfinite(log_beta) or log_beta <= 0:
        raise ValueError('Require a positive beta log and quartet-aligned width')
    delta = np.expm1(log_beta * ST_LN2)
    if not np.isfinite(delta) or delta <= 0:
        raise ArithmeticError('beta-1 is not positive and finite')
    populations = np.arange(0, width + 1, 4, dtype=np.int64)
    features = np.zeros((width // 4 + 1, width + 1), dtype=ST_LD)
    features[:, 0] = 1
    for j in range(width):
        active = populations > j
        ratio = np.sqrt(delta * ST_LD(j + 1) / ST_LD(width - j))
        features[active, j + 1] = (features[active, j]
                                   * (populations[active] - j) / ST_LD(j + 1)
                                   * ratio)
        values = features[active, j + 1]
        if np.any(values <= 0) or not np.all(np.isfinite(values)):
            raise ArithmeticError('Positive Gram feature underflow or overflow')
    overlap = features @ features.T
    if np.any(overlap < 1) or not np.all(np.isfinite(overlap)):
        raise ArithmeticError('Overlap moment is not finite and at least one')
    return _st_normalize(overlap)


def temporal_selector_weights(p, r, log_beta, width=ST_WIDTH):
    """Scaled positive w_a, retaining every population including a=0."""
    width, r = index(width), index(r)
    p, log_beta = ST_LD(p), ST_LD(log_beta)
    if width < 4 or width % 4 or not 1 <= r <= 3:
        raise ValueError('Require quartet width and r in1..3')
    if not np.isfinite(p) or not 0 < p < 1 or not np.isfinite(log_beta) or log_beta <= 0:
        raise ValueError('Require0<p<1 and a positive finite beta log')
    labels = [ST_LD(comb(4, j)) * p ** j * (1 - p) ** (4 - j) for j in range(5)]
    q0 = sum(labels[:r], ST_LD(0))
    q1 = sum((labels[j] * np.exp2(log_beta * (8 * (j - r)))
              for j in range(r, 5)), ST_LD(0))
    if not (np.isfinite(q0) and np.isfinite(q1) and q0 > 0 and q1 > 0):
        raise ArithmeticError('Weighted quartet probabilities are invalid')
    m = width // 4
    logs = np.array([ST_LD(_st_logchoose(m, a)) + (m - a) * np.log2(q0)
                     + a * np.log2(q1) for a in range(m + 1)], dtype=ST_LD)
    scale = np.max(logs)
    weights = np.exp2(logs - scale)
    if np.any(weights <= 0) or not np.all(np.isfinite(weights)):
        raise ArithmeticError('Positive population weight underflow or overflow')
    return weights, scale


def _st_power_row(kernel, kernel_scale, updates):
    """Positive binary powers with scalar log normalizations and free endpoint."""
    updates = index(updates)
    if updates < 0:
        raise ValueError('Nonnegative update count required')
    row = np.zeros(kernel.shape[0], dtype=ST_LD)
    row[0] = 1
    row_scale = ST_LD(0)
    power, power_scale = kernel.copy(), ST_LD(kernel_scale)
    count = updates
    while count:
        if count & 1:
            row, extra = _st_normalize(row @ power)
            row_scale += power_scale + extra
        count >>= 1
        if count:
            power, extra = _st_normalize(power @ power)
            power_scale = 2 * power_scale + extra
    total = np.sum(row)
    if not np.isfinite(total) or total <= 0:
        raise ArithmeticError('Invalid free-endpoint moment')
    return row_scale + np.log2(total)


def temporal_selector_log_moment(p, r, log_beta, width=ST_WIDTH, updates=ST_UPDATES):
    """FLOAT log2 E[beta^(E_r+8Z_r)] before regional quota conditioning."""
    overlap, overlap_scale = temporal_selector_overlap(log_beta, width)
    weights, weight_scale = temporal_selector_weights(p, r, log_beta, width)
    kernel, extra = _st_normalize(overlap * weights[None, :])
    value = _st_power_row(kernel, overlap_scale + weight_scale + extra, updates)
    if not np.isfinite(value):
        raise ArithmeticError('Nonfinite temporal moment')
    return float(value)


def temporal_selector_diagnostic(h=3584, r=2, log_beta=1 / 64, p=None, profile=None):
    """One bounded setup-only FLOAT control; no certificate is claimed."""
    h, r = index(h), index(r)
    log_beta = float(log_beta)
    p = h / ST_ROWS if p is None else float(p)
    if not 1 <= h < ST_ROWS or not 1 <= r <= 3 or not isfinite(log_beta) or log_beta <= 0:
        raise ValueError('Require1<=h<8192,r in1..3,log_beta>0')
    if not isfinite(p) or not 0 < p < 1:
        raise ValueError('Require0<p<1')
    moment = temporal_selector_log_moment(p, r, log_beta)
    choose = _st_logchoose(ST_ROWS, h)
    log_pi = choose + h * log2(p) + (ST_ROWS - h) * log2(1 - p)
    union_and_quota = choose - ST_REGIONS * log_pi
    result = dict(h=h, r=r, p=p, log2_beta=log_beta,
                  setup_log2_moment=moment, log2_H_union=choose,
                  log2_quota_probability=log_pi,
                  external_union_and_quota=union_and_quota,
                  float_tail92_threshold=(moment + union_and_quota + 92) / log_beta,
                  states=ST_WIDTH // 4 + 1, no_state_history_union=True,
                  free_endpoint=True, float_diagnostic_only=True,
                  exact_certificate=False, full_code_certificate=False)
    if profile is not None:
        R, Ls, factual, F = map(index, profile)
        if not 0 <= factual <= Ls <= F <= ST_WIDTH * ST_UPDATES or R < 128 * h:
            raise ValueError('Invalid actual-state profile')
        G = F - factual
        if G % 4:
            raise ValueError('Actual G must be divisible by four')
        z = (4 - r) * (G // 4) - R + 128 * h
        score = Ls - ST_WIDTH - factual + 8 * z
        result.update(R=R, Ls=Ls, factual=factual, F=F, G=G,
                      selector_score=z, temporal_score_required=score,
                      float_log2_tail_upper=moment + union_and_quota - log_beta * score)
    return result


def _st_partitions(width):
    if width == 4:
        return ((15,),)
    if width == 8:
        return tuple((1 + sum(1 << i for i in rest),
                      255 ^ (1 + sum(1 << i for i in rest)))
                     for rest in combinations(range(1, 8), 3))
    raise ValueError('Exact partition toys use width4 or8')


def _st_exact_overlap(width, a, b, beta):
    first = (1 << (4 * a)) - 1
    return sum((beta ** (first & sum(1 << i for i in subset)).bit_count()
                for subset in combinations(range(width), 4 * b)), Fraction(0)) / comb(width, 4 * b)


def _st_exact_mask_weights(width, p, r, beta):
    partitions = _st_partitions(width)
    weights = {}
    for partition in partitions:
        for labels in range(1 << width):
            eligible = labels.bit_count()
            probability = p ** eligible * (1 - p) ** (width - eligible) / len(partitions)
            mask = score = 0
            for quartet in partition:
                j = (quartet & labels).bit_count()
                if j >= r:
                    mask |= quartet
                score += max(0, j - r)
            weights[mask] = weights.get(mask, Fraction(0)) + probability * beta ** (8 * score)
    return weights


def temporal_selector_checks():
    """Independent exact tests; root alone executes this function."""
    overlaps = populations = moments = quota_checks = inequalities = 0
    for width in (4, 8):
        m = width // 4
        for p, r in product((Fraction(1, 2), Fraction(1, 3)), (1, 2, 3)):
            beta = Fraction(3, 2)
            masks = _st_exact_mask_weights(width, p, r, beta)
            labels = [comb(4, j) * p ** j * (1 - p) ** (4 - j) for j in range(5)]
            q0 = sum(labels[:r])
            q1 = sum(labels[j] * beta ** (8 * (j - r)) for j in range(r, 5))
            weights = [comb(m, a) * q0 ** (m - a) * q1 ** a for a in range(m + 1)]
            for mask, value in masks.items():
                a = mask.bit_count() // 4
                assert value == weights[a] / comb(width, 4 * a)
                populations += 1
            matrix = []
            for a in range(m + 1):
                row = []
                for b in range(m + 1):
                    exact = _st_exact_overlap(width, a, b, beta)
                    factorial = sum(Fraction(comb(4 * a, j) * comb(4 * b, j), comb(width, j))
                                    * (beta - 1) ** j for j in range(min(4 * a, 4 * b) + 1))
                    assert exact == factorial
                    row.append(exact * weights[b])
                    overlaps += 1
                matrix.append(row)
            row = [Fraction(int(a == 0)) for a in range(m + 1)]
            assert sum(row) == 1
            for updates in (1, 2):
                row = [sum(row[a] * matrix[a][b] for a in range(m + 1)) for b in range(m + 1)]
                direct = (sum(masks.values()) if updates == 1 else
                          sum(first * second * beta ** (left & right).bit_count()
                              for left, first in masks.items() for right, second in masks.items()))
                assert sum(row) == direct
                floating = temporal_selector_log_moment(float(p), r, log2(float(beta)), width, updates)
                assert abs(floating - log2(float(direct))) < 2e-10
                moments += 1
    # Two width-four updates: condition the eight independent labels on
    # exactly h eligible slots and compare with the iid quota upper bound.
    for h, r in product(range(9), (1, 2, 3)):
        beta = Fraction(3, 2)
        total = Fraction(0)
        for selected in combinations(range(8), h):
            mask = sum(1 << i for i in selected)
            first, second = (mask & 15).bit_count(), (mask >> 4).bit_count()
            edge = 4 if first >= r and second >= r else 0
            score = edge + 8 * (max(0, first - r) + max(0, second - r))
            total += beta ** score
        conditional = total / comb(8, h)
        weights = _st_exact_mask_weights(4, Fraction(1, 2), r, beta)
        unconditional = sum(a * b * beta ** (left & right).bit_count()
                            for left, a in weights.items() for right, b in weights.items())
        assert conditional <= unconditional * 256 / comb(8, h)
        quota_checks += 1
    # A crossed two-layer quartet graph has two edges between each pair
    # of its two left and two right nodes.  Check the deficit lemma itself.
    for labels in product(range(5), repeat=4):
        for selected in range(16):
            bits = [(selected >> i) & 1 for i in range(4)]
            internal = 2 * (bits[0] + bits[1]) * (bits[2] + bits[3])
            for r in (1, 2, 3):
                Z = sum(max(0, j - r) for j in labels)
                required = sum((j - r) * active for j, active in zip(labels, bits))
                deficit = Z - required
                P = [int(j >= r) for j in labels]
                edges = 2 * (P[0] + P[1]) * (P[2] + P[3])
                assert deficit >= sum(active for j, active in zip(labels, bits) if j < r)
                assert internal <= edges + 8 * deficit
                inequalities += 1
    print(dict(exact_overlap_checks=overlaps, exact_uniform_mask_checks=populations,
               exact_and_float_moment_checks=moments, exact_quota_checks=quota_checks,
               exact_temporal_selector_inequalities=inequalities,
               full_code_certificate=False), flush=True)


def temporal_selector_main():
    parser = argparse.ArgumentParser(description=__doc__)
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument('--checks', action='store_true')
    mode.add_argument('--probe', action='store_true')
    parser.add_argument('--h', type=int, default=3584)
    parser.add_argument('--r', type=int, default=2)
    parser.add_argument('--log-beta', type=float, default=1 / 64)
    parser.add_argument('--p', type=float)
    parser.add_argument('--profile', type=int, nargs=4, metavar=('R', 'LS', 'F_BIRTH', 'F'),
                        default=(786432, 786432, 4, 1015808))
    args = parser.parse_args()
    import resource
    if shutil.disk_usage('/tmp').free < 20 * (1 << 30):
        raise RuntimeError('20GiB /tmp free-space guard')
    resource.setrlimit(resource.RLIMIT_AS, (2 * (1 << 30), 2 * (1 << 30)))
    timer = threading.Timer(175, lambda: os._exit(124))
    timer.daemon = True
    timer.start()
    started = time.monotonic()
    try:
        if args.checks:
            temporal_selector_checks()
        else:
            print(temporal_selector_diagnostic(args.h, args.r, args.log_beta, args.p, args.profile), flush=True)
        print(dict(seconds=round(time.monotonic() - started, 3)), flush=True)
    finally:
        timer.cancel()


if __name__ == '__main__':
    temporal_selector_main()
