"""Exact input-weight/output-support counts for the double mixer, a<=6.

For x in GF(q)^a, v -> sum(v_i*x_i) is GF4-linear. Its canonical
relation kernel determines every zero input/output coordinate. A kernel
of codimension d has prod_{h<d}(q-4^h) preimages. We enumerate only
positive-mass kernels, without merging original source labels.

Packed RREF rows are stored by their pivot in bounded uint16 arrays.
Membership tests operate on whole arrays, not expanded kernel bitsets.
At a6 there are only64 possible dense Hnew row forms (lambda or lambda+1
in each position), plus six original input forms. The controller can
reuse those tests across all sixteen scheduled columns.

All histogram additions use uint64. Total mass is at most256^6=2^48;
there is no floating-point arithmetic in the enumerator. It supplies
raw counts, not RS/source normalization, a regional bound, or a certificate.
"""
from hashlib import sha256
from itertools import combinations, product
from math import comb
from pathlib import Path
from time import perf_counter

import numpy as np

import hd_lambda_double_transvection as outer
import hd_quaternary_outer_relations as base
from hd_pair_scheduled_outer import column_permutations


DIRECTORY = Path(__file__).resolve().parent
SCHEDULE = 'input-inverse-output-xor'
PINS = {
    'hd_lambda_double_transvection.py': '973175e25051621470bd19f8c6932806e4a830b908a8a1eb04186e4b616430bd',
    'hd_quaternary_outer_relations.py': 'd8be4e7a1ab77b98c76695f13c1bc9a9cf867df2aae6408ac3e08c6fc12346c0',
    'hd_pair_scheduled_outer.py': '64dc50bf18171d371ae16f4484087e9d42d89cc7c9b55a7576be55b1251e4741',
}
DECODE = {value: i for i, value in enumerate(base.AES_EMBEDDING)}


def dependencies():
    for name, digest in PINS.items():
        if sha256((DIRECTORY/name).read_bytes()).hexdigest() != digest:
            raise ValueError('Reviewed dependency changed: '+name)
    return {**outer.checked_dependencies(), **PINS,
            Path(__file__).name: sha256(Path(__file__).read_bytes()).hexdigest()}


def _parameters(arity, q):
    if type(arity) is not int or not 1 <= arity <= 6 or type(q) is not int or q not in (4, 256):
        raise ValueError('Only arity1..6 over GF4 or GF256 is supported')


def _charts(arity, q):
    """Canonical RREF charts, with their exact positive raw mass."""
    _parameters(arity, q)
    for dimension in range(arity+1):
        count = base.injection_count(arity-dimension, q)
        if not count:
            continue
        for pivots in combinations(range(arity), dimension):
            free = tuple((i, j) for i, pivot in enumerate(pivots)
                         for j in range(pivot+1, arity) if j not in pivots)
            yield pivots, free, count


def relation_arrays(arity, q=256, *, check_time=None):
    """One array column per distinct positive-mass GF4 relation kernel."""
    charts = tuple(_charts(arity, q))
    size = sum(4**len(free) for _, free, _ in charts)
    basis = np.zeros((arity, size), dtype=np.uint16)
    masses = np.empty(size, dtype=np.uint64)
    cursor = 0
    for pivots, free, count in charts:
        if check_time: check_time()
        chart_size = 4**len(free)
        masses[cursor:cursor+chart_size] = count
        for entries in product(range(4), repeat=len(free)):
            rows = [1 << (2*pivot) for pivot in pivots]
            for (i, j), value in zip(free, entries):
                rows[i] |= value << (2*j)
            for pivot, row in zip(pivots, rows):
                basis[pivot, cursor] = row
            cursor += 1
        if check_time: check_time()
    if cursor != size or int(masses.sum(dtype=np.uint64)) != q**arity:
        raise AssertionError('Canonical relation kernels lost exact input mass')
    scale2 = np.array([sum(base.multiply_gf4(2, row >> (2*i) & 3) << (2*i)
                              for i in range(arity)) for row in range(4**arity)], dtype=np.uint16)
    return basis, scale2[basis], masses


def zero_membership(form, basis, twice):
    """Exact membership of one packed form in each canonical kernel."""
    arity, size = basis.shape
    if type(form) is not int or not 0 <= form < 4**arity or twice.shape != basis.shape:
        raise ValueError('A packed form and matching canonical row arrays are required')
    residue = np.full(size, form, dtype=np.uint16)
    for pivot in range(arity):
        coefficient = (residue >> (2*pivot)) & 3
        residue ^= np.where(coefficient & 1, basis[pivot], np.uint16(0))
        residue ^= np.where(coefficient & 2, twice[pivot], np.uint16(0))
    return residue == 0


def scheduled_forms(source, column):
    source = tuple(source)
    if (not 1 <= len(source) <= 6 or len(set(source)) != len(source)
            or any(type(i) is not int or not 0 <= i < 16 for i in source)):
        raise ValueError('One to six distinct ordered original rows required')
    before, after = column_permutations(SCHEDULE, column)
    matrix, forms = outer.matrix(), [0]*16
    for i, row in enumerate(matrix):
        forms[after[i]] = sum(DECODE[row[before[j]]] << (2*k) for k, j in enumerate(source))
    return tuple(forms)


def input_weights(basis, twice):
    weights = np.zeros(basis.shape[1], dtype=np.uint8)
    for i in range(basis.shape[0]):
        weights += ~zero_membership(1 << (2*i), basis, twice)
    return weights


def exact_columns(source, *, q=256, columns=tuple(range(16)), check_time=None):
    """List of exact {(original input weight,physical output mask): count}."""
    source, columns = tuple(source), tuple(columns)
    _parameters(len(source), q)
    if (not columns or len(set(columns)) != len(columns)
            or any(type(j) is not int or not 0 <= j < 16 for j in columns)):
        raise ValueError('One to sixteen distinct scheduled columns required')
    started = perf_counter()
    all_forms = [scheduled_forms(source, j) for j in columns]
    basis, twice, masses = relation_arrays(len(source), q, check_time=check_time)
    weights = input_weights(basis, twice)
    # Cache at most64 dense-form membership arrays (about35MiB at a6).
    membership = {}
    for form in sorted(set(form for row in all_forms for form in row)):
        if check_time: check_time()
        membership[form] = zero_membership(form, basis, twice)
    result = []
    for j, forms in zip(columns, all_forms):
        if check_time: check_time()
        masks = np.zeros(len(masses), dtype=np.uint16)
        for i, form in enumerate(forms):
            masks |= (~membership[form]).astype(np.uint16) << i
        index = weights.astype(np.uint32)*65536+masks
        histogram = np.zeros((len(source)+1)*65536, dtype=np.uint64)
        np.add.at(histogram, index, masses)
        if int(histogram.sum(dtype=np.uint64)) != q**len(source):
            raise AssertionError('Exact joint histogram lost input mass')
        for weight, subtotal in enumerate(histogram.reshape(len(source)+1, 65536).sum(axis=1, dtype=np.uint64)):
            if int(subtotal) != comb(len(source), weight)*(q-1)**weight:
                raise AssertionError('Original input-weight count changed')
        counts = {divmod(int(i), 65536): int(histogram[i]) for i in np.flatnonzero(histogram)}
        if counts.get((0, 0)) != 1 or sum(c for (w, _), c in counts.items() if not w) != 1:
            raise AssertionError('The zero input atom changed')
        result.append(counts)
    return result, dict(original_source=source, columns=columns, field_order=q,
        positive_mass_relation_kernels=len(masses), distinct_output_forms=len(membership),
        exact_mass_per_column=q**len(source), original_input_weights=True,
        original_source_fixed_across_columns=True, integer_histograms=True,
        elapsed_seconds=perf_counter()-started)
