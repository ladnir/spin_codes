"""Exact column catalogs and bounded FLOAT shared-ratio norms for a<=6.

The construction is the ORIGINAL Hlambda (AES lambda0xbc), followed by
the fixed P/Q schedule and staggered column shear.  Only the generic
GF4 relation enumerator is reused from hd_lambda_double_exact6; its old
double-mixer matrix is NOT used here.  Actual column forms come from
hd_lambda_column_shear_ranks.column_rows.

For a original input coordinates, all Hlambda row forms have coefficients
in {lambda,lambda+1}.  There are2^a such forms, at most64.  A GF4 relation
kernel determines their zero bits and the original input weight.  Merge
kernels only when both agree.  Each kernel has its exact GF256 injection
mass; zero-mass kernels are omitted.  No256^a amplitude enumeration occurs.

The catalog is common to all source sets and columns, but a relation
kernel is NOT fixed across columns.  Each call gathers the actual sixteen
physical forms for the same original labelled A.  Exact single-column
histograms retain Python-int counts before FLOAT conversion.  Their raw
pair mass may exceed uint64, so products use sums of individual log counts.

The shared-ratio lemma, atom domain[0,1], and Lr Minkowski coefficients
are unchanged from the frozen helper.  The only extensions are bounded
histogram/input-weight guards and a reusable exact dense-form catalog.
pair_coefficients already returns the branch-root polynomial: do not
take another root.  The caller applies entrywise r-powers and eight
chronological matrix-polynomial factors, with tail >=9*a*r.  RS/source
normalization is not included.  No catalog or gate is built on import.
"""
from dataclasses import dataclass
from hashlib import sha256
import json
from math import comb, log2
from pathlib import Path

import numpy as np

import hd_lambda_double_exact6 as relation
import hd_lambda_column_shear_ranks as shear
import hd_lambda_shear_shared_ratio as ratio


DIRECTORY = Path(__file__).resolve().parent
MAX_ARITY = 6
MAX_HISTOGRAM = 4096
MAX_PAIR_PRODUCTS = 1_000_000
MAX_CHUNK = 8192
PINS = {
    'hd_lambda_double_exact6.py': 'c6065ca7f6b9605aec096283d8ff2db94630692f99a4d2deed9e20b30d4dc960',
    'hd_lambda_column_shear_ranks.py': '4067013201cadef412ae647cc81e394e90b49d9403dfdd4713a739a191a96fa8',
    'hd_lambda_shear_shared_ratio.py': '5b60c22fa6e120c9e4faee276de0f029dc52121c21dc53c2d7c12eafdac9b363',
}


class BudgetExceeded(ValueError):
    """Exact preparation succeeded but the prescribed contraction is too large."""


def dependencies():
    # Old double-matrix pins identify the inherited enumerator source;
    # they do not bind the construction evaluated by column_rows.
    result = {**relation.dependencies(), **ratio.dependencies(), **shear.dependencies()}
    for name, digest in PINS.items():
        if sha256((DIRECTORY/name).read_bytes()).hexdigest() != digest:
            raise ValueError('Frozen dependency changed: '+name)
        result[name] = digest
    result[Path(__file__).name] = sha256(Path(__file__).read_bytes()).hexdigest()
    return result


def _parameters(arity, q):
    if (type(arity) is not int or not 1 <= arity <= MAX_ARITY
            or type(q) is not int or q not in (4, 256)):
        raise ValueError('Only arities1..6 over GF4/GF256 are supported')


def universal_forms(arity):
    _parameters(arity, 256)
    return tuple(sum((2 | (bits >> i & 1)) << (2*i) for i in range(arity))
                 for bits in range(1 << arity))


def catalog_plan(arity, q=256):
    """Exact chart/resource counts only; no relation arrays are allocated."""
    _parameters(arity, q)
    charts = tuple(relation._charts(arity, q))
    kernels = sum(4**len(free) for _, free, _ in charts)
    mass = sum(4**len(free)*count for _, free, count in charts)
    if mass != q**arity:
        raise ArithmeticError('The relation-chart partition lost raw mass')
    return dict(arity=arity, field_order=q, positive_mass_kernels=kernels,
                dense_forms=1 << arity, exact_raw_mass=mass,
                basis_and_twice_bytes=4*arity*kernels,
                mass_zero_weight_bytes=17*kernels)


def _readonly(values, dtype):
    result = np.array(values, dtype=dtype, copy=True)
    result.setflags(write=False)
    return result


@dataclass(frozen=True)
class DenseColumnCatalog:
    arity: int
    q: int
    zero_masks: np.ndarray
    input_weights: np.ndarray
    counts: np.ndarray
    positive_mass_kernels: int
    source_pins: tuple

    def validate(self):
        _parameters(self.arity, self.q)
        n = len(self.counts)
        if (not n or self.counts.dtype != np.uint64 or self.zero_masks.dtype != np.uint64
                or self.input_weights.dtype != np.uint8
                or self.zero_masks.shape != (n,) or self.input_weights.shape != (n,)
                or np.any(self.counts == 0) or np.any(self.input_weights > self.arity)
                or self.positive_mass_kernels < n):
            raise ValueError('Malformed exact dense-column catalog')
        full = np.uint64((1 << (1 << self.arity))-1)
        if np.any(self.zero_masks > full):
            raise ValueError('Zero bits exceed the dense-form dictionary')
        order = np.lexsort((self.zero_masks, self.input_weights))
        if not np.array_equal(order, np.arange(n)):
            raise ValueError('Expected input-weight/zero-mask order')
        if n > 1 and np.any((self.zero_masks[1:] == self.zero_masks[:-1])
                & (self.input_weights[1:] == self.input_weights[:-1])):
            raise ValueError('Duplicate catalog category')
        if sum(map(int, self.counts)) != self.q**self.arity:
            raise ValueError('Lost exact column mass')
        for weight in range(self.arity+1):
            if sum(map(int, self.counts[self.input_weights == weight])) != comb(self.arity, weight)*(self.q-1)**weight:
                raise ValueError('Lost original-input-weight mass')
        zero = self.input_weights == 0
        if (np.count_nonzero(zero) != 1 or int(self.counts[zero][0]) != 1
                or int(self.zero_masks[zero][0]) != int(full)
                or np.count_nonzero(self.zero_masks == full) != 1):
            raise ValueError('Dense forms must annihilate only the unique zero input')
        return self


def build_catalog(arity, q=256, *, check_time=None):
    """Root-owned exact preparation; production a5/a6 builds are not tests."""
    _parameters(arity, q)
    pins = dependencies()
    basis, twice, masses = relation.relation_arrays(arity, q, check_time=check_time)
    weights = relation.input_weights(basis, twice)
    zeros = np.zeros(len(masses), dtype=np.uint64)
    for index, form in enumerate(universal_forms(arity)):
        if check_time is not None:
            check_time()
        zeros |= relation.zero_membership(form, basis, twice).astype(np.uint64) << np.uint64(index)
    order = np.lexsort((zeros, weights))
    zeros, weights, masses = zeros[order], weights[order], masses[order]
    starts = np.r_[0, 1+np.flatnonzero((zeros[1:] != zeros[:-1]) | (weights[1:] != weights[:-1]))]
    result = DenseColumnCatalog(arity, q, _readonly(zeros[starts], np.uint64),
        _readonly(weights[starts], np.uint8),
        _readonly(np.add.reduceat(masses, starts, dtype=np.uint64), np.uint64),
        len(masses), tuple(sorted(pins.items()))).validate()
    if dependencies() != pins:
        raise ValueError('Sources changed during catalog preparation')
    return result


def column_histogram(source, column, catalog, *, schedule=ratio.SCHEDULE):
    """Exact (original input weight,physical16-bit mask,Python-int count)."""
    source = tuple(source)
    if not isinstance(catalog, DenseColumnCatalog) or len(source) != catalog.arity:
        raise ValueError('Original source arity and catalog must agree')
    forms = shear.column_rows(source, column, schedule=schedule)
    lookup = {form: index for index, form in enumerate(universal_forms(catalog.arity))}
    if any(form not in lookup for form in forms):
        raise ValueError('Actual Hlambda output escaped the dense GF4 dictionary')
    masks = np.zeros(len(catalog.counts), dtype=np.uint16)
    for row, form in enumerate(forms):
        zero = (catalog.zero_masks >> np.uint64(lookup[form])) & np.uint64(1)
        masks |= (1-zero).astype(np.uint16) << row
    keys = catalog.input_weights.astype(np.uint32)*65536+masks
    order = np.argsort(keys)
    keys, masses = keys[order], catalog.counts[order]
    starts = np.r_[0, 1+np.flatnonzero(keys[1:] != keys[:-1])]
    counts = np.add.reduceat(masses, starts, dtype=np.uint64)
    result = tuple((int(key)//65536, int(key)%65536, int(count))
                   for key, count in zip(keys[starts], counts))
    if (sum(row[2] for row in result) != catalog.q**catalog.arity
            or result[0] != (0, 0, 1) or sum(row[1] == 0 for row in result) != 1):
        raise ArithmeticError('Actual column map lost injectivity or raw mass')
    for weight in range(catalog.arity+1):
        if sum(n for w, _, n in result if w == weight) != comb(catalog.arity, weight)*(catalog.q-1)**weight:
            raise ArithmeticError('Original input-weight histogram changed')
    if len(result) > MAX_HISTOGRAM:
        raise BudgetExceeded('Exact histogram exceeds the declared4096-row contraction budget')
    return result


def _histogram(values, rows):
    values = tuple(values)
    if (not 1 <= len(values) <= MAX_HISTOGRAM
            or any(len(row) != 3 or any(type(x) is not int for x in row)
                   or not 0 <= row[0] <= MAX_ARITY or not 0 <= row[1] < 1 << rows
                   or row[2] <= 0 for row in values)
            or len({row[:2] for row in values}) != len(values)):
        raise ValueError('Distinct exact positive histograms with weights0..6 and <=4096 rows required')
    return values


def pair_coefficients(first_histogram, second_histogram, prepared, *, q=256,
                      chunk_size=8192, max_products=MAX_PAIR_PRODUCTS, check_time=None):
    """Already-rooted log polynomial, shape(2a+1,2,2), with degree at most12."""
    if (not isinstance(prepared, ratio.RatioKernel) or type(chunk_size) is not int
            or not 1 <= chunk_size <= MAX_CHUNK or type(max_products) is not int
            or not 1 <= max_products <= MAX_PAIR_PRODUCTS):
        raise ValueError('Prepared kernel, chunk<=8192 and positive product budget<=1M required')
    first = _histogram(first_histogram, prepared.rows)
    second = _histogram(second_histogram, prepared.rows)
    products = len(first)*len(second)
    if products > max_products:
        raise BudgetExceeded(f'Exact histogram product {products} exceeds budget {max_products}')
    wa = np.array([row[0] for row in first], dtype=np.int8)
    wb = np.array([row[0] for row in second], dtype=np.int8)
    sa = np.array([row[1] for row in first], dtype=np.uint16)
    sb = np.array([row[1] for row in second], dtype=np.uint16)
    na = np.array([log2(row[2]) for row in first])
    nb = np.array([log2(row[2]) for row in second])
    maximum = int(wa.max()+wb.max())
    result = np.full((maximum+1, 2, 2), -np.inf)
    width = len(second)
    for start in range(0, products, chunk_size):
        if check_time is not None:
            check_time()
        indices = np.arange(start, min(start+chunk_size, products))
        left, right = indices//width, indices%width
        terms = ratio.mixture_log_norms(sa[left], sb[right], prepared, q=q)
        terms += (na[left]+nb[right])[:, None, None]
        weights = wa[left]+wb[right]
        for degree in range(maximum+1):
            selected = terms[weights == degree]
            if len(selected):
                result[degree] = np.logaddexp2(result[degree], np.logaddexp2.reduce(selected, axis=0))
    return result


def pair_polynomial(source, column, catalog, prepared, *, schedule=ratio.SCHEDULE,
                    chunk_size=8192, max_products=MAX_PAIR_PRODUCTS, check_time=None):
    if (type(column) is not int or column not in range(0, 16, 2)
            or not isinstance(prepared, ratio.RatioKernel) or prepared.rows != 16):
        raise ValueError('Even chronological column0..14 and16-row prepared kernel required')
    return pair_coefficients(column_histogram(source, column, catalog, schedule=schedule),
        column_histogram(source, column+1, catalog, schedule=schedule), prepared, q=catalog.q,
        chunk_size=chunk_size, max_products=max_products, check_time=check_time)


def save_catalog(catalog, path):
    catalog.validate()
    if dict(catalog.source_pins) != dependencies():
        raise ValueError('Catalog source pins differ from current frozen dependencies')
    metadata = dict(arity=catalog.arity, q=catalog.q, complete=True,
        positive_mass_kernels=catalog.positive_mass_kernels,
        merged_categories=len(catalog.counts), source_pins=dict(catalog.source_pins),
        universal_forms=list(universal_forms(catalog.arity)),
        category_order='input-weight-zero-mask', exact_raw_mass=catalog.q**catalog.arity,
        original_input_weights=True, construction='original-Hlambda-bc-not-double-mixer',
        whole_code_certificate=False)
    path = Path(path)
    with path.open('xb') as handle:
        np.savez_compressed(handle, zero_masks=catalog.zero_masks, input_weights=catalog.input_weights,
                            counts=catalog.counts, metadata=np.array(json.dumps(metadata, sort_keys=True)))
    return dict(metadata, artifact=str(path), sha256=sha256(path.read_bytes()).hexdigest())


def load_catalog(path, expected_sha256):
    path = Path(path)
    if sha256(path.read_bytes()).hexdigest() != expected_sha256:
        raise ValueError('Pinned catalog hash differs')
    with np.load(path, allow_pickle=False) as data:
        if set(data.files) != {'zero_masks', 'input_weights', 'counts', 'metadata'}:
            raise ValueError('Unexpected catalog schema')
        metadata = json.loads(str(data['metadata']))
        if (metadata.get('complete') is not True or metadata.get('source_pins') != dependencies()
                or metadata.get('universal_forms') != list(universal_forms(metadata['arity']))
                or metadata.get('category_order') != 'input-weight-zero-mask'
                or metadata.get('construction') != 'original-Hlambda-bc-not-double-mixer'):
            raise ValueError('Catalog scope or provenance differs')
        result = DenseColumnCatalog(metadata['arity'], metadata['q'], _readonly(data['zero_masks'], np.uint64),
            _readonly(data['input_weights'], np.uint8), _readonly(data['counts'], np.uint64),
            metadata['positive_mass_kernels'], tuple(sorted(metadata['source_pins'].items()))).validate()
    if (metadata.get('merged_categories') != len(result.counts)
            or metadata.get('exact_raw_mass') != result.q**result.arity):
        raise ValueError('Catalog metadata and exact arrays disagree')
    return result


def metadata():
    return dict(ratio.metadata(), supported_arities=list(range(1, 7)),
        catalog_exact_integer=True, catalog_matrix_independent=True,
        actual_matrix='Hlambda-bc', old_double_matrix_used=False,
        maximum_dense_forms=64, maximum_histogram=MAX_HISTOGRAM,
        maximum_pair_products=MAX_PAIR_PRODUCTS, maximum_chunk=MAX_CHUNK,
        maximum_pair_input_weight=12, raw_pair_counts_fit_uint64_claimed=False)
