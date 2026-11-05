"""Accepted-map incidence restriction for the global joint availability count.

In an old-availability-hit quartet let i=|Aprev|>=1 and let Enew=Vnew-Anew.
If Enew is nonempty, the ACTUAL incoming state is nonzero: otherwise every
covered coordinate is a new active birth and therefore belongs to Anew.
Its support size k is at most i.  Superregularity of the accepted A matrix
implies that A*Sprev has at most k-1 zero coordinates.  Every outside-H
input is zero by RS shortening, so on H^c the new state equals A*Sprev.
Every Enew coordinate has zero new state.  Consequently

    |Enew intersect H^c| <= k-1 <= i-1.

When Enew is empty the inequality is automatic.  This is a deterministic
accepted-map constraint, not an additional field probability or run cost.

The new local operator has the same 15 count states as the prior spectral
helper.  Its hit powers also fit the existing joint-subset factorization,
without adding physical state dimensions; the terminal condition Enew=empty
and the general zero operator are unchanged.  All polynomial coefficients are
positive.  For norm probes inject quartet_joint_spectral.py definitions
with its main block removed.  Tests may run standalone.  No local imports.

Optionally retain the output budget from superregular A+delta*I.  In a
hit quartet with j outside-H positions, at least max(0,j-i+1) outputs
are nonzero; each unhit birth also has nonzero output.  Their sum B is at
most W.  output_discount=eta^-1 inserts eta^-B, with external eta^W.
This is a deterministic event marker and does not spend field rank twice.
"""

import argparse
from fractions import Fraction
from itertools import combinations, product
from math import comb, lgamma
import os
import shutil
import threading
import time

import numpy as np


MDS_LD = np.longdouble
MDS_N = 1 << 21
MDS_W = 943718
MDS_STATES = tuple((u, e) for u in range(5) for e in range(5 - u))
MDS_INDEX = {state: i for i, state in enumerate(MDS_STATES)}


def _mds_cast(value):
    if isinstance(value, Fraction):
        return MDS_LD(value.numerator) / MDS_LD(value.denominator)
    return MDS_LD(value)


def joint_mds_local(m, b, eligible, outside, alpha, coverage=1, uncovered=1, exact=False,
                     output_discount=1):
    """Local count kernel with e_out<i; eligible and outside are SEPARATE weights.

    At quota counter z, use eligible=p*z, outside=(1-p)/a_field,
    uncovered=1-p+p*z.  Coverage/class markers multiply the corresponding
    covered weights separately; uncovered retains the raw Bernoulli law.
    Thus the same regional coefficient conditioning remains available.
    """
    cast = Fraction if exact else _mds_cast
    dtype = object if exact else MDS_LD
    m, b, eligible, outside, alpha, coverage, uncovered, output_discount = map(
        cast, (m, b, eligible, outside, alpha, coverage, uncovered, output_discount))
    if min(m, b, eligible, outside, alpha, coverage, uncovered, output_discount) < 0:
        raise ValueError('All local weights must be nonnegative')
    result = np.zeros((15, 15), dtype=dtype)
    birth = m * eligible * coverage * output_discount
    for row, (u, e) in enumerate(MDS_STATES):
        if u:
            if e:
                continue
            for new_a in range(5):
                new_e = 4 - new_a
                state = sum(comb(u, keep) * comb(4 - u, new_a - keep) * alpha ** (u - keep)
                            for keep in range(max(0, new_a - (4 - u)), min(u, new_a) + 1))
                route_weight = sum(
                    comb(new_a, ja) * comb(new_e, je)
                    * eligible ** (4 - ja - je) * outside ** (ja + je)
                    * output_discount ** max(0, ja + je - u + 1)
                    for ja in range(new_a + 1) for je in range(min(new_e, u - 1) + 1))
                result[row, MDS_INDEX[new_a, new_e]] = (
                    coverage ** 4 * b ** new_a * state * route_weight)
        else:
            for new_a in range(5 if e == 0 else 4):
                result[row, MDS_INDEX[new_a, 0]] = (
                    comb(4, new_a) * birth ** new_a * uncovered ** (4 - new_a))
    if not exact and (np.any(result < 0) or not np.all(np.isfinite(result))):
        raise ArithmeticError('Invalid positive MDS local kernel')
    return result


def joint_mds_hit_coefficients(b, eligible, outside, alpha, coverage=1, exact=False,
                               output_discount=1):
    """Coefficients x^i w^f y^(4-f) of one hit quartet, including C(4,i).

    The local count row fixes an old mask.  The additional C(4,i) here is
    needed when averaging partitions by multinomial coefficient extraction.
    """
    local = joint_mds_local(1, b, eligible, outside, alpha, coverage, exact=exact,
                            output_discount=output_discount)
    result = np.zeros((5, 5), dtype=object if exact else MDS_LD)
    for i in range(1, 5):
        for f in range(5):
            result[i, f] = comb(4, i) * local[MDS_INDEX[i, 0], MDS_INDEX[4 - f, f]]
    return result


def joint_mds_hit_powers(groups, hit, exact=False):
    """Positive ragged bivariate powers; power k uses row index u-k.

    H_k(u,f) is the coefficient of x^u w^f y^(4k-f).  For the subset
    transfer, multiply by C(M,k)/multinom(n;u,e,n-u-e), convolve the
    y^(4k-f) factor with G_(M-k,e)(m*eligible*coverage*output_discount*y), and sum f in
    the retained-hole bin.  Uniform subset sampling adds no binomial cost.
    This supplies the same physical (u,e) rows for the existing bin lift.
    This G notation assumes uncovered=1.  At a formal quota counter, G
    must also carry the uncovered-coordinate factor 1-p+p*z.
    """
    dtype = object if exact else MDS_LD
    powers = [np.ones((1, 1), dtype=dtype)]
    for k in range(groups):
        previous = powers[-1]
        after = np.zeros((3 * (k + 1) + 1, 4 * (k + 1) + 1), dtype=dtype)
        for i in range(1, 5):
            for f in range(5):
                if hit[i, f]:
                    after[i - 1:i - 1 + previous.shape[0], f:f + previous.shape[1]] += hit[i, f] * previous
        powers.append(after)
    return powers


def joint_mds_potential_action(width, m, b, eligible, outside, alpha, t, s,
                               coverage=1, uncovered=1, output_discount=1):
    """Exact positive action on t^|A|s^|E| using 129 outer products at width512."""
    if width < 4 or width % 4:
        raise ValueError('Width must be a positive quartet width')
    m, b, eligible, outside, alpha, t, s, coverage, uncovered = map(
        _mds_cast, (m, b, eligible, outside, alpha, t, s, coverage, uncovered))
    groups = width // 4
    local = joint_mds_local(m, b, eligible, outside, alpha, coverage, uncovered,
                            output_discount=output_discount)
    target = np.array([t ** u * s ** e for u, e in MDS_STATES], dtype=MDS_LD)
    action = local @ target
    hit = np.array([MDS_LD(0)] + [comb(4, u) * action[MDS_INDEX[u, 0]] for u in range(1, 5)])
    empty = np.array([comb(4, e) * action[MDS_INDEX[0, e]] for e in range(5)])
    left = np.zeros((width + 1, groups + 1), dtype=MDS_LD)
    right = np.zeros_like(left)
    hit_power = np.ones(1, dtype=MDS_LD)
    empty_power = np.ones(1, dtype=MDS_LD)
    for k in range(groups + 1):
        left[:len(hit_power), k] = comb(groups, k) * hit_power
        right[:len(empty_power), groups - k] = empty_power
        if k < groups:
            hit_power = np.convolve(hit_power, hit)
            empty_power = np.convolve(empty_power, empty)
    coefficients = left @ right.T
    result = np.zeros_like(coefficients)
    for u in range(width + 1):
        denominator = np.array([MDS_LD(comb(width, u) * comb(width - u, e))
                                for e in range(width - u + 1)])
        result[u, :width - u + 1] = coefficients[u, :width - u + 1] / denominator
    if not np.all(np.isfinite(result)):
        raise ArithmeticError('Nonfinite MDS potential action')
    return result


def _joint_mds_physical(m, b, raw_p, outside_factor, alpha, coverage, counter, output_discount=1):
    """Independent 81-configuration matrix, explicitly enumerating 16 H masks."""
    states = tuple(product(range(3), repeat=4))
    masks = [(sum(1 << j for j, value in enumerate(state) if value == 1),
              sum(1 << j for j, value in enumerate(state) if value == 2)) for state in states]
    result = np.full((81, 81), Fraction(0), dtype=object)
    for row, (old_a, old_e) in enumerate(masks):
        for col, (new_a, new_e) in enumerate(masks):
            if old_a:
                if old_e or new_e != (15 ^ new_a):
                    continue
                covered = 15
            else:
                if new_e or (old_e and new_a == 15):
                    continue
                covered = new_a
            for eligible_mask in range(16):
                outside_covered = covered & ~eligible_mask
                if old_a:
                    if (new_e & ~eligible_mask).bit_count() >= old_a.bit_count():
                        continue
                    state_weight = b ** new_a.bit_count() * alpha ** (old_a & ~new_a).bit_count()
                    state_weight *= output_discount ** max(0, outside_covered.bit_count() - old_a.bit_count() + 1)
                else:
                    if outside_covered:
                        continue
                    state_weight = (m * output_discount) ** new_a.bit_count()
                label_weight = ((raw_p * counter) ** eligible_mask.bit_count()
                                * (1 - raw_p) ** (4 - eligible_mask.bit_count()))
                result[row, col] += (label_weight * outside_factor ** outside_covered.bit_count()
                                     * coverage ** covered.bit_count() * state_weight)
    return masks, result


def _joint_mds_field_checks():
    """Exact GF(8) Cauchy-matrix checks of the used superregular lemma only.

    This matrix is a toy algebra witness, not a sample from the construction's
    GF(q^4)-multiplier acceptance distribution.  The new lemma uses A alone.
    """
    def multiply(a, b):
        result = 0
        while b:
            if b & 1:
                result ^= a
            b >>= 1
            a <<= 1
            if a & 8:
                a ^= 0b1011
        return result
    table = [[multiply(a, b) for b in range(8)] for a in range(8)]
    inverse = [0] + [next(b for b in range(1, 8) if table[a][b] == 1) for a in range(1, 8)]
    matrix = [[inverse[x ^ y] for y in range(4, 8)] for x in range(4)]
    minors = 0
    for size in range(1, 5):
        for rows in combinations(range(4), size):
            for cols in combinations(range(4), size):
                work = [[matrix[i][j] for j in cols] for i in rows]
                for pivot in range(size):
                    chosen = next((i for i in range(pivot, size) if work[i][pivot]), None)
                    assert chosen is not None
                    work[pivot], work[chosen] = work[chosen], work[pivot]
                    multiplier = inverse[work[pivot][pivot]]
                    work[pivot] = [table[value][multiplier] for value in work[pivot]]
                    for i in range(pivot + 1, size):
                        multiplier = work[i][pivot]
                        work[i] = [a ^ table[multiplier][b] for a, b in zip(work[i], work[pivot])]
                minors += 1
    vectors = incidences = zero_cases = 0
    maximal_zeros = [0] * 5
    for vector in product(range(8), repeat=4):
        old_support = sum(1 << j for j, value in enumerate(vector) if value)
        if not old_support:
            continue
        output = []
        for row in matrix:
            value = 0
            for coefficient, coordinate in zip(row, vector):
                value ^= table[coefficient][coordinate]
            output.append(value)
        zeros = sum(1 << j for j, value in enumerate(output) if not value)
        support = old_support.bit_count()
        assert zeros.bit_count() <= support - 1
        maximal_zeros[support] = max(maximal_zeros[support], zeros.bit_count())
        for eligible_mask in range(16):
            for new_e in range(16):
                outside_e = new_e & ~eligible_mask
                if outside_e & ~zeros:
                    continue  # Such an E would contradict Qnew=A*Sprev on H^c.
                assert outside_e.bit_count() < support
                incidences += 1
        vectors += 1
    assert maximal_zeros[1:] == [0, 1, 2, 3]
    # With zero incoming state, V=Snew and Snew subset Anew force Enew empty.
    for new_state in range(16):
        for new_available in range(16):
            if new_state & ~new_available:
                continue
            covered = new_state
            assert covered & ~new_available == 0
            zero_cases += 1
    return minors, vectors, incidences, zero_cases


def _joint_mds_pair_checks():
    """Exact GF(256) accepted-pair controls for both state and output bounds."""
    def multiply(a, b):
        result = 0
        while b:
            if b & 1:
                result ^= a
            b >>= 1
            a <<= 1
            if a & 256:
                a ^= 0x11B
        return result
    table = [[multiply(a, b) for b in range(256)] for a in range(256)]
    inverse = [0] + [next(b for b in range(1, 256) if table[a][b] == 1) for a in range(1, 256)]
    selections = [(rows, cols) for size in range(1, 5)
                  for rows in combinations(range(4), size) for cols in combinations(range(4), size)]
    def superregular(matrix):
        for rows, cols in selections:
            size = len(rows)
            work = [[matrix[i][j] for j in cols] for i in rows]
            for pivot in range(size):
                chosen = next((i for i in range(pivot, size) if work[i][pivot]), None)
                if chosen is None:
                    return False
                work[pivot], work[chosen] = work[chosen], work[pivot]
                multiplier = inverse[work[pivot][pivot]]
                work[pivot] = [table[value][multiplier] for value in work[pivot]]
                for i in range(pivot + 1, size):
                    multiplier = work[i][pivot]
                    work[i] = [a ^ table[multiplier][b] for a, b in zip(work[i], work[pivot])]
        return True
    matrix = [[inverse[x ^ y] for y in range(4, 8)] for x in range(4)]
    assert superregular(matrix)
    pair = None
    for delta in range(1, 256):
        candidate = [[matrix[i][j] ^ (delta if i == j else 0) for j in range(4)] for i in range(4)]
        if superregular(candidate):
            pair = candidate
            break
    # At most 80 scalar roots are excluded by the nonconstant minor polynomials.
    assert pair is not None
    vectors = set(product(range(4), repeat=4))
    vectors.remove((0, 0, 0, 0))
    for coordinate in range(4):
        for value in range(1, 256):
            vector = [0, 0, 0, 0]
            vector[coordinate] = value
            vectors.add(tuple(vector))
    output_checks = 0
    for vector in vectors:
        state_support = sum(1 << i for i, value in enumerate(vector) if value)
        def image(matrix):
            result = []
            for row in matrix:
                value = 0
                for coefficient, coordinate in zip(row, vector):
                    value ^= table[coefficient][coordinate]
                result.append(value)
            return result
        state_image = image(matrix)
        output_image = image(pair)
        assert output_image == [value ^ table[delta][coordinate] for value, coordinate in zip(state_image, vector)]
        state_zero = sum(1 << i for i, value in enumerate(state_image) if not value)
        output_nonzero = sum(1 << i for i, value in enumerate(output_image) if value)
        for available in range(16):
            if state_support & ~available:
                continue
            r = available.bit_count()
            for eligible in range(16):
                outside = 15 ^ eligible
                assert (state_zero & outside).bit_count() < r
                assert (output_nonzero & outside).bit_count() >= max(0, outside.bit_count() - r + 1)
                output_checks += 1
    return delta, 2 * len(selections), len(vectors), output_checks


def joint_mds_checks():
    F = Fraction
    m, b, raw_p, outside_factor = F(3, 4), F(5, 4), F(1, 3), F(1, 10)
    entries = terminal_checks = counter_checks = power_checks = norm_checks = 0
    for alpha, coverage, counter, discount in ((F(2, 3), F(4, 5), F(1), F(1)),
                                                (F(2), F(1), F(2), F(1, 2)),
                                                (F(0), F(1), F(1, 2), F(1, 3)),
                                                (F(2, 3), F(0), F(1), F(1))):
        eligible = raw_p * counter
        outside = (1 - raw_p) * outside_factor
        uncovered = 1 - raw_p + raw_p * counter
        local = joint_mds_local(m, b, eligible, outside, alpha, coverage, uncovered, exact=True,
                                output_discount=discount)
        masks, physical = _joint_mds_physical(m, b, raw_p, outside_factor, alpha, coverage, counter, discount)
        for row, (old_a, old_e) in enumerate(masks):
            collapsed = [F(0)] * 15
            for col, (new_a, new_e) in enumerate(masks):
                collapsed[MDS_INDEX[new_a.bit_count(), new_e.bit_count()]] += physical[row, col]
            assert collapsed == list(local[MDS_INDEX[old_a.bit_count(), old_e.bit_count()]])
            entries += 15
        floating = joint_mds_local(m, b, eligible, outside, alpha, coverage, uncovered,
                                   output_discount=discount)
        for i in range(15):
            for j in range(15):
                truth = _mds_cast(local[i, j])
                assert abs(floating[i, j] - truth) <= MDS_LD('2e-14') * max(1, abs(truth))
                counter_checks += 1
        if counter == 1 and discount == 1 and 'joint_spectral_local' in globals():
            old = joint_spectral_local(m, b, raw_p, raw_p + outside, alpha, coverage, exact=True)
            assert all(local[i, j] <= old[i, j] for i in range(15) for j in range(15))
            for i in range(15):
                for new_a in range(5):
                    assert local[i, MDS_INDEX[new_a, 0]] == old[i, MDS_INDEX[new_a, 0]]
                    terminal_checks += 1
        if 'joint_spectral_log_norm' in globals():
            multiplicity = np.array([comb(4, u) * comb(4 - u, e) for u, e in MDS_STATES])
            for logt, logs in ((0, 0), (-1, 1)):
                physical_float = np.array([[_mds_cast(value) for value in row] for row in physical], dtype=float)
                diagonal = np.array([2. ** (logt * a.bit_count() + logs * e.bit_count()) for a, e in masks])
                transformed = physical_float * diagonal[None, :] / diagonal[:, None]
                actual = np.log2(np.linalg.svd(transformed, compute_uv=False)[0])
                expected = joint_spectral_log_norm(floating, logt, logs)
                assert abs(actual - expected) <= MDS_LD('3e-12')
                norm_checks += 1
        hit = joint_mds_hit_coefficients(b, eligible, outside, alpha, coverage, exact=True,
                                        output_discount=discount)
        powers = joint_mds_hit_powers(2, hit, exact=True)
        for u in range(2, 9):
            for f in range(9):
                truth = sum(hit[i, j] * hit[u - i, f - j]
                            for i in range(1, 5) for j in range(5)
                            if 1 <= u - i <= 4 and 0 <= f - j <= 4)
                assert powers[2][u - 2, f] == truth
                power_checks += 1
    # All-outside negative controls exercise each forbidden minor size.
    pure_outside = joint_mds_local(F(1), F(1), F(0), F(1), F(1), exact=True)
    for incoming in range(1, 5):
        assert pure_outside[MDS_INDEX[incoming, 0], MDS_INDEX[4 - incoming, incoming]] == 0
    # Independently check the aggregate-marker substitution on labeled
    # monomials. The old/new A potential telescopes; the local MDS inequality
    # uses old A, so the new-A exponent need not be nonnegative at one step.
    aggregate_checks = 0
    eta = F(2)
    for old_a in range(16):
        for new_a in range(16):
            covered = 15 if old_a else new_a
            holes = covered & ~new_a
            births = 0 if old_a else new_a
            for labels in range(16):
                if births & ~labels:
                    continue
                outside_holes = (holes & ~labels).bit_count()
                if old_a and outside_holes >= old_a.bit_count():
                    continue
                old_exponent = 4 * old_a.bit_count() - covered.bit_count() + births.bit_count() - 4 * outside_holes
                assert old_exponent >= 0
                exponent = 4 * new_a.bit_count() - covered.bit_count() + births.bit_count() - 4 * outside_holes
                # Birth: eta^4. Hit available: eta^3. Hit hole:
                # eta^-1 inside H, eta^-5 outside H.
                substituted = eta ** (4 * births.bit_count())
                if old_a:
                    substituted *= eta ** (3 * new_a.bit_count())
                    substituted *= eta ** (-(holes & labels).bit_count() - 5 * outside_holes)
                assert substituted == eta ** exponent
                assert exponent == old_exponent + 4 * (new_a.bit_count() - old_a.bit_count())
                aggregate_checks += 1
    minors, vectors, incidences, zero_cases = _joint_mds_field_checks()
    delta, pair_minors, pair_vectors, output_checks = _joint_mds_pair_checks()
    print(dict(mds_physical_label_entries=entries, local_numeric_and_quota_checks=counter_checks,
               unchanged_terminal_entries=terminal_checks, exact_hit_power_entries=power_checks,
               count_norm_checks=norm_checks, gf8_superregular_minors=minors,
               gf8_nonzero_state_vectors=vectors, gf8_outside_hole_incidences=incidences,
               zero_incoming_empty_E_checks=zero_cases, aggregate_marker_monomials=aggregate_checks,
               gf256_pair_delta=delta, gf256_pair_superregular_minors=pair_minors,
               gf256_pair_state_vectors=pair_vectors, exact_output_budget_incidence_checks=output_checks,
               full_code_certificate=False), flush=True)


def joint_mds_spectral_value(h=2048, theta=1/16, logz=-.4375, logt=0., logs=0.,
                             coverage=0., outside=0., rankgap=0., tilts=9, outputgap=0.):
    """All-region norm diagnostic with existing coverage/outside/rank-gap caps."""
    if 'joint_spectral_log_norm' not in globals():
        raise RuntimeError('Inject quartet_joint_spectral.py before a norm probe')
    if not 0 < h < 8192 or not 0 < theta <= 1 or logz > 0 or tilts < 1:
        raise ValueError('Invalid field/row/tilt parameters')
    if coverage > 0 or outside > 0 or rankgap < 0 or outputgap < 0:
        raise ValueError('Require coverage/outside<=0 and rankgap/outputgap>=0')
    theta, logz, logt, logs, coverage, outside, rankgap, outputgap = map(
        MDS_LD, (theta, logz, logt, logs, coverage, outside, rankgap, outputgap))
    q = MDS_LD(2) ** 32
    x = np.exp2(logz / theta)
    c = (q + 1) * x
    lam = 1 + c
    d = max(min(MDS_LD(1), lam / (q + 1)), (1 + q * x * x) / lam)
    m, b, alpha = c ** theta, lam ** theta, d ** theta
    m *= np.exp2(rankgap)
    b *= np.exp2(rankgap)
    chi = np.exp2(coverage)
    outside_weight = np.exp2(outside - rankgap - 32 * theta)
    actual = MDS_LD(h) / 8192
    effective_field = 1 / outside_weight
    low = actual / (effective_field - (effective_field - 1) * actual)
    choose = MDS_LD((lgamma(8193) - lgamma(h + 1) - lgamma(8193 - h)) / float(np.log(MDS_LD(2))))
    best = None
    for fraction in np.linspace(0, 1, tilts):
        odds = (low / (1 - low)) ** (1 - fraction) * (actual / (1 - actual)) ** fraction
        p = odds / (1 + odds)
        logpi = choose + h * np.log2(p) + (8192 - h) * np.log2(1 - p)
        local = joint_mds_local(m, b, p, (1 - p) * outside_weight, alpha, chi,
                                output_discount=np.exp2(-outputgap))
        lognorm = joint_spectral_log_norm(local, logt, logs)
        logregion = 2048 * lognorm - logpi
        if best is None or logregion < best[0]:
            best = (logregion, p, lognorm, logpi)
    terminal = 256 * np.logaddexp2(MDS_LD(0), -2 * logt)
    external = (theta - 32 * theta - MDS_W * logz + choose - 32 * theta * 128 * h
                - (1009255 - 128 * h) * outside - 1736703 * coverage
                + (65536 - 128 * h) * rankgap + MDS_W * outputgap)
    return float(external + 256 * best[0] + terminal), dict(
        h=h, theta=float(theta), logz=float(logz), logt=float(logt), logs=float(logs),
        logcoverage=float(coverage), logoutside=float(outside), logrankgap=float(rankgap),
        outputgap=float(outputgap),
        p=float(best[1]), log2_local_norm=float(best[2]), log2_quota_probability=float(best[3]),
        log2_terminal=float(terminal), local_dimension=15, tilts=tilts,
        all_regions_only=True, occupied_129_retained=False, certificate=False)


def joint_mds_main():
    parser = argparse.ArgumentParser(description=__doc__)
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument('--checks', action='store_true')
    mode.add_argument('--probe', action='store_true')
    parser.add_argument('--h', type=int, default=2048)
    parser.add_argument('--theta', type=float, default=1/16)
    parser.add_argument('--logz', type=float, default=-.4375)
    parser.add_argument('--logt', type=float, default=0.)
    parser.add_argument('--logs', type=float, default=0.)
    parser.add_argument('--coverage', type=float, default=0.)
    parser.add_argument('--outside', type=float, default=0.)
    parser.add_argument('--rankgap', type=float, default=0.)
    parser.add_argument('--outputgap', type=float, default=0.)
    parser.add_argument('--tilts', type=int, default=9)
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
            joint_mds_checks()
        else:
            bound, detail = joint_mds_spectral_value(args.h, args.theta, args.logz, args.logt, args.logs,
                                                     args.coverage, args.outside, args.rankgap, args.tilts,
                                                     args.outputgap)
            print(dict(log2_bound=bound, seconds=round(time.monotonic() - started, 3), **detail), flush=True)
    finally:
        timer.cancel()


if __name__ == '__main__':
    joint_mds_main()
