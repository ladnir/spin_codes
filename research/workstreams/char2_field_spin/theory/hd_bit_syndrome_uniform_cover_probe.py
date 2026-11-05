"""Bounded exact algebra probe for uniform information-set covers.

No regional provider or distance contraction is evaluated.  A successful
probe returns explicit balanced bases for the selected labelled A only.
Failure to find a template is inconclusive; a counterexample is reported
only after checking16*rank(H[S,A]) < |A|*|S| exactly.
"""
from itertools import combinations
import random
import time

from gf8_packet_maps import MUL
import hd_bit_syndrome_outer as old
import hd_bit_syndrome_star_outer as star
from hd_bit_syndrome_flat_geometry import affine_flats

FULL = (1 << 16) - 1


def labels(mask):
    while mask:
        bit = mask & -mask
        yield bit.bit_length() - 1
        mask ^= bit


def checked_mask(mask):
    if type(mask) is not int or not 0 <= mask <= FULL:
        raise ValueError('Expected a sixteen-coordinate mask')
    return mask


def binary_plus_rank_one_rank(rows, u, v):
    """Exact GF256 rank of binary P plus u*v^T, with P rows as bitmasks.

    Binary elimination transforms u using only XORs.  A field-valued row
    reduction tests v's membership in row(P) and, when both memberships
    hold, obtains v*x for any P*x=u.  The remaining decision is one scalar.
    """
    rows, u, v = tuple(rows), tuple(u), tuple(v)
    if len(rows) != len(u) or any(x < 0 or x >= (1 << len(v)) for x in rows):
        raise ValueError('Binary dimensions do not agree')
    if any(type(x) is not int or not 0 <= x < 256 for x in u + v):
        raise ValueError('AES-field bytes required')
    basis = {}
    u_outside = False
    for row, value in zip(rows, u):
        while row:
            pivot = (row & -row).bit_length() - 1
            if pivot not in basis:
                basis[pivot] = (row, value)
                break
            row ^= basis[pivot][0]
            value ^= basis[pivot][1]
        if not row and value:
            u_outside = True
    residual = list(v)
    scalar = 0
    for pivot in sorted(basis):
        row, value = basis[pivot]
        coefficient = residual[pivot]
        scalar ^= MUL[coefficient][value]
        for column in labels(row):
            residual[column] ^= coefficient
    v_outside = any(residual)
    r = len(basis)
    if u_outside and v_outside:
        return r + 1
    if u_outside or v_outside:
        return r
    return r - int(scalar == 1)


class RankOracle:
    """Actual H[S,A] ranks, INCLUDING its diagonal on overlapping sets."""
    def __init__(self, construction='current'):
        if construction not in ('current', 'star-coupled'):
            raise ValueError('Unknown fixed construction')
        self.construction = construction
        self.matrix = old.matrix() if construction == 'current' else star.matrix('coupled')
        self.u = old.LABELS
        binary = tuple(tuple(self.matrix[i][j] ^ MUL[self.u[i]][self.u[j]]
                             for j in range(16)) for i in range(16))
        if any(x not in (0, 1) for row in binary for x in row):
            raise AssertionError('The binary-plus-rank-one decomposition changed')
        self.rows = tuple(sum(x << j for j, x in enumerate(row)) for row in binary)
        self.queries = 0

    def rank(self, output_mask, source_mask):
        checked_mask(output_mask)
        checked_mask(source_mask)
        self.queries += 1
        outputs = tuple(labels(output_mask))
        # Keep the original sixteen column labels; zero outside A.
        return binary_plus_rank_one_rank(
            tuple(self.rows[i] & source_mask for i in outputs),
            tuple(self.u[i] for i in outputs),
            tuple(self.u[j] if source_mask >> j & 1 else 0 for j in range(16)))


def translate(mask, shift):
    return sum(1 << (i ^ shift) for i in labels(mask))


def verify_balanced_bases(oracle, source, bases):
    """Check a positive witness without relying on the search procedure."""
    checked_mask(source)
    a = source.bit_count()
    bases = tuple(bases)
    if not a or not bases:
        return False
    if any(checked_mask(b).bit_count() != a or oracle.rank(b, source) != a for b in bases):
        return False
    degrees = [sum((b >> i) & 1 for b in bases) for i in range(16)]
    return len(set(degrees)) == 1 and degrees[0] * 16 == len(bases) * a


def find_translated_cover(oracle, source, candidates):
    """Try fixed balanced templates; no search failure is a counterexample."""
    a = source.bit_count()
    for candidate in candidates:
        if candidate.bit_count() != a:
            raise ValueError('Template has the wrong cardinality')
        bases = tuple(sorted(set(translate(candidate, s) for s in range(16))))
        if all(oracle.rank(b, source) == a for b in bases):
            if not verify_balanced_bases(oracle, source, bases):
                raise AssertionError('A proposed positive cover is not balanced')
            return bases
    return None


def density_counterexample(oracle, source, candidates):
    a = source.bit_count()
    for output in candidates:
        r = oracle.rank(output, source)
        if 16 * r < a * output.bit_count():
            return dict(source_mask=source, output_mask=output, arity=a,
                        output_size=output.bit_count(), rank=r,
                        deficit=a * output.bit_count() - 16 * r)
    return None


def probe(*, max_seconds=8, samples_per_arity=12, templates=32):
    """Bounded selected-row-set probe, not a complete all-A certificate."""
    if not 0 < max_seconds <= 10 or not 1 <= samples_per_arity <= 32 or not 1 <= templates <= 64:
        raise ValueError('Probe budget exceeds the explicit local-algebra guard')
    started = time.monotonic()
    rng = random.Random(20261010)
    flat_masks = set(affine_flats())
    candidates_s = tuple(sorted(flat_masks | {FULL ^ s for s in flat_masks}))
    cases = []
    for a in range(3, 9):
        values = {sum(1 << i for i in range(a)),
                  sum(1 << (2 * i) for i in range(min(a, 8)))}
        while len(values) < samples_per_arity:
            values.add(sum(1 << i for i in rng.sample(range(16), a)))
        cases.extend(sorted(values))
    results = {}
    for construction in ('current', 'star-coupled'):
        oracle = RankOracle(construction)
        successes, unresolved = [], []
        counter = None
        for source in cases:
            if time.monotonic() - started >= max_seconds:
                break
            a = source.bit_count()
            first = {source, FULL ^ source} if a == 8 else {source}
            while len(first) < templates:
                first.add(sum(1 << i for i in rng.sample(range(16), a)))
            cover = find_translated_cover(oracle, source, sorted(first))
            if cover is not None:
                successes.append(dict(source_mask=source, arity=a, bases=cover))
                continue
            counter = density_counterexample(oracle, source, candidates_s + (source, FULL ^ source))
            if counter is not None:
                break
            unresolved.append(source)
        results[construction] = dict(successes=successes, unresolved=unresolved,
                                     counterexample=counter, rank_queries=oracle.queries,
                                     requested_cases=len(cases),
                                     processed_cases=len(successes) + len(unresolved) + int(counter is not None))
    return dict(schema='bounded-exact-uniform-cover-probe-v1', results=results,
                elapsed_seconds=time.monotonic() - started,
                whole_code_certificate=False, all_row_sets_covered=False,
                regional_provider_evaluated=False)
