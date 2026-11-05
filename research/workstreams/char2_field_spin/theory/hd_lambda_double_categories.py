"""Universal exact dense-form categories for the double mixer, arity<=6.

Repeated RREF enumeration is unnecessary for different mapped source sets.
Every matrix entry is lambda or lambda+1. In ordered input coordinates, its
row form is one of2^a forms lambda*1+b, b in F2^a. A GF4 relation kernel
determines which forms vanish and the total original input weight. Merge
kernels with the same pair (input_weight, dense_zero_pattern), adding their
exact raw injection counts. Any selected output matrix then gathers pattern
bits and forms a joint input-weight/output-mask histogram without elimination.

The catalog is independent of the mapped source set. Merging preserves total
input weight, NOT the identities of active input rows. A pre-Q cache may sort
mapped source labels because permuting inputs preserves that weight and raw
ambient counting measure. Each scheduled column still maps the SAME original
A through P_j, and applies its actual Q_j to every output support.

Arrays are uint64 counts/patterns and uint8 weights; total count<=256^6=2^48.
No float, RS normalization, regional calculation, or proof gate is included.
The optional CLI builds ONE catalog only. Its completeness is the complete
positive-mass relation-kernel enumeration, not an all-source proof claim.
"""
from collections import OrderedDict
from dataclasses import dataclass
from hashlib import sha256
import json
from math import comb
from pathlib import Path
from time import perf_counter

import numpy as np

import hd_lambda_double_exact6 as exact


DIRECTORY = Path(__file__).resolve().parent
EXACT_SHA = 'c6065ca7f6b9605aec096283d8ff2db94630692f99a4d2deed9e20b30d4dc960'


def dependencies():
    if sha256((DIRECTORY/'hd_lambda_double_exact6.py').read_bytes()).hexdigest() != EXACT_SHA:
        raise ValueError('Frozen exact relation enumerator changed')
    return {**exact.dependencies(), Path(__file__).name: sha256(Path(__file__).read_bytes()).hexdigest()}


def dense_forms(arity):
    exact._parameters(arity, 256)
    return tuple(sum((2 | (pattern >> i & 1)) << (2*i) for i in range(arity))
                 for pattern in range(1 << arity))


def dense_form_index(form, arity):
    exact._parameters(arity, 256)
    if (type(form) is not int or not 0 <= form < 4**arity
            or any((form >> (2*i) & 3) not in (2, 3) for i in range(arity))):
        raise ValueError('Every coefficient must be GF4 symbol2 or3, not an AES byte')
    return sum((form >> (2*i) & 1) << i for i in range(arity))


def _readonly(values, dtype):
    result = np.array(values, dtype=dtype, copy=True)
    result.setflags(write=False)
    return result


@dataclass(frozen=True)
class DenseCategories:
    arity: int
    q: int
    zero_patterns: np.ndarray
    input_weights: np.ndarray
    counts: np.ndarray
    positive_mass_kernels: int
    source_pins: tuple

    def validate(self):
        exact._parameters(self.arity, self.q)
        n = len(self.counts)
        if (not n or self.zero_patterns.shape != (n,) or self.input_weights.shape != (n,)
                or self.counts.dtype != np.uint64 or self.zero_patterns.dtype != np.uint64
                or self.input_weights.dtype != np.uint8 or np.any(self.counts == 0)
                or np.any(self.input_weights > self.arity)
                or self.positive_mass_kernels < n):
            raise ValueError('Inconsistent exact dense-category arrays')
        order = np.lexsort((self.zero_patterns, self.input_weights))
        if not np.array_equal(order, np.arange(n)):
            raise ValueError('Catalog must be ordered by input weight then zero pattern')
        if n > 1 and np.any((self.input_weights[1:] == self.input_weights[:-1])
                             & (self.zero_patterns[1:] == self.zero_patterns[:-1])):
            raise ValueError('Duplicate merged category')
        full = (1 << (1 << self.arity))-1
        if self.arity < 6 and np.any(self.zero_patterns > full):
            raise ValueError('Zero-pattern bits exceed the declared dense forms')
        # Python summation also detects a corrupted archive with uint64 wraparound.
        if sum(map(int, self.counts)) != self.q**self.arity:
            raise ValueError('Catalog lost exact ambient mass')
        for weight in range(self.arity+1):
            subtotal = sum(map(int, self.counts[self.input_weights == weight]))
            if subtotal != comb(self.arity, weight)*(self.q-1)**weight:
                raise ValueError('Catalog changed original input-weight counts')
        zero = self.input_weights == 0
        if np.count_nonzero(zero) != 1 or int(self.counts[zero][0]) != 1 or int(self.zero_patterns[zero][0]) != full:
            raise ValueError('The zero input must be a single all-forms-zero category')
        if np.count_nonzero(self.zero_patterns == np.uint64(full)) != 1:
            raise ValueError('Dense forms span the input space; only zero annihilates all of them')
        return self


def build_catalog(arity=6, q=256, *, check_time=None):
    """Enumerate once; merge exactly; no source-set inventory or proof bound."""
    pins = dependencies()
    basis, twice, masses = exact.relation_arrays(arity, q, check_time=check_time)
    weights = exact.input_weights(basis, twice)
    patterns = np.zeros(len(masses), dtype=np.uint64)
    for index, form in enumerate(dense_forms(arity)):
        if check_time: check_time()
        zero = exact.zero_membership(form, basis, twice)
        patterns |= zero.astype(np.uint64) << np.uint64(index)
    order = np.lexsort((patterns, weights))
    patterns, weights, masses = patterns[order], weights[order], masses[order]
    begins = np.r_[0, 1+np.flatnonzero((patterns[1:] != patterns[:-1]) | (weights[1:] != weights[:-1]))]
    merged = np.add.reduceat(masses, begins, dtype=np.uint64)
    catalog = DenseCategories(arity, q, _readonly(patterns[begins], np.uint64),
        _readonly(weights[begins], np.uint8), _readonly(merged, np.uint64), len(masses), tuple(sorted(pins.items())))
    catalog.validate()
    if dependencies() != pins:
        raise ValueError('Sources changed during catalog construction')
    return catalog


def joint_arrays(row_forms, catalog):
    """Sparse flat indices weight*2^rows+mask and exact masses, both read-only."""
    forms = tuple(row_forms)
    if not 1 <= len(forms) <= 16:
        raise ValueError('One to sixteen dense output forms required')
    indexes = tuple(dense_form_index(form, catalog.arity) for form in forms)
    masks = np.zeros(len(catalog.counts), dtype=np.uint16)
    for row, index in enumerate(indexes):
        active = 1 ^ ((catalog.zero_patterns >> np.uint64(index)) & np.uint64(1))
        masks |= active.astype(np.uint16) << row
    stride = 1 << len(forms)
    addresses = catalog.input_weights.astype(np.uint32)*stride+masks
    histogram = np.zeros((catalog.arity+1)*stride, dtype=np.uint64)
    np.add.at(histogram, addresses, catalog.counts)
    if int(histogram.sum(dtype=np.uint64)) != catalog.q**catalog.arity:
        raise AssertionError('Gathering dense forms lost raw mass')
    occupied = np.flatnonzero(histogram)
    return _readonly(occupied, np.uint32), _readonly(histogram[occupied], np.uint64)


def joint_counts(row_forms, catalog):
    """Fresh dictionary of (total original input weight,output mask)->raw count."""
    forms = tuple(row_forms)
    indexes, counts = joint_arrays(forms, catalog)
    stride = 1 << len(forms)
    return {divmod(int(index), stride): int(count) for index, count in zip(indexes, counts)}


class MappedHistogramCache:
    """Bounded pre-Q LRU; reuse changes execution order only, never source scope."""
    def __init__(self, catalog, *, max_entries=512, max_terms=1 << 20):
        catalog.validate()
        if (type(max_entries) is not int or not 0 <= max_entries <= 8008
                or type(max_terms) is not int or not 0 <= max_terms <= 1 << 24):
            raise ValueError('Explicit bounded histogram-cache capacity required')
        self.catalog = catalog
        self.max_entries, self.max_terms = max_entries, max_terms
        self._cache, self._terms = OrderedDict(), 0
        self.hits = self.misses = 0
        self._matrix = exact.outer.matrix()

    def _source(self, source):
        source = tuple(source)
        if (len(source) != self.catalog.arity or len(set(source)) != len(source)
                or any(type(i) is not int or not 0 <= i < 16 for i in source)):
            raise ValueError('Exactly arity distinct source labels0..15 required')
        return source

    def mapped_arrays(self, mapped_source):
        # Sorting is valid ONLY for total input weight, not labelled input support.
        source = tuple(sorted(self._source(mapped_source)))
        if source in self._cache:
            self.hits += 1
            result = self._cache.pop(source)
            self._cache[source] = result
            return result
        self.misses += 1
        forms = tuple(sum(exact.DECODE[row[label]] << (2*i) for i, label in enumerate(source))
                      for row in self._matrix)
        result = joint_arrays(forms, self.catalog)
        size = len(result[0])
        if self.max_entries and size <= self.max_terms:
            while self._cache and (len(self._cache) >= self.max_entries or self._terms+size > self.max_terms):
                _, old = self._cache.popitem(last=False)
                self._terms -= len(old[0])
            self._cache[source] = result
            self._terms += size
        return result

    def mapped_counts(self, mapped_source):
        indexes, counts = self.mapped_arrays(mapped_source)
        return {divmod(int(index), 65536): int(count) for index, count in zip(indexes, counts)}

    def scheduled_columns(self, original_source, *, columns=tuple(range(16)), check_time=None):
        source, columns = self._source(original_source), tuple(columns)
        if (not columns or len(set(columns)) != len(columns)
                or any(type(j) is not int or not 0 <= j < 16 for j in columns)):
            raise ValueError('One to sixteen distinct original scheduled columns required')
        result = []
        for column in columns:
            if check_time: check_time()
            before, after = exact.column_permutations(exact.SCHEDULE, column)
            indexes, counts = self.mapped_arrays(tuple(before[i] for i in source))
            masks, weights = indexes & 65535, indexes >> 16
            physical = np.zeros(len(indexes), dtype=np.uint16)
            for i, destination in enumerate(after):
                physical |= (((masks >> i) & 1) << destination).astype(np.uint16)
            keys = weights*65536+physical
            order = np.argsort(keys)
            result.append({divmod(int(keys[k]), 65536): int(counts[k]) for k in order})
        return result, dict(original_source=source, columns=columns, field_order=self.catalog.q,
            exact_mass_per_column=self.catalog.q**self.catalog.arity,
            merged_dense_categories=len(self.catalog.counts),
            positive_mass_relation_kernels=self.catalog.positive_mass_kernels,
            original_source_fixed_across_columns=True, total_input_weight_preserved=True,
            individual_input_row_identity_preserved=False, mapped_cache_entries=len(self._cache),
            mapped_cache_terms=self._terms, cache_hits=self.hits, cache_misses=self.misses,
            cache_key_sorting_is_not_a_source_orbit_quotient=True)


def save_catalog(catalog, path):
    """Exclusive NPZ export with exact arrays; caller records returned file hash."""
    catalog.validate()
    path = Path(path)
    metadata = dict(arity=catalog.arity, q=catalog.q,
        positive_mass_relation_kernels=catalog.positive_mass_kernels,
        merged_categories=len(catalog.counts), raw_mass=catalog.q**catalog.arity,
        dense_form_symbols=[2, 3], aes_embedding=exact.base.AES_EMBEDDING,
        category_order='input-weight-then-uint64-dense-zero-pattern',
        input_weight_is_total_not_labelled=True, complete=True,
        source_pins=dict(catalog.source_pins), regional_provider_runs=0,
        whole_code_certificate=False)
    with path.open('xb') as handle:
        np.savez_compressed(handle, zero_patterns=catalog.zero_patterns,
            input_weights=catalog.input_weights, counts=catalog.counts,
            metadata=np.array(json.dumps(metadata, sort_keys=True)))
    return dict(metadata, artifact=str(path), sha256=sha256(path.read_bytes()).hexdigest())


def load_catalog(path, expected_sha256):
    """Consume a pinned catalog, never silently accept an unbound count table."""
    path = Path(path)
    if sha256(path.read_bytes()).hexdigest() != expected_sha256:
        raise ValueError('Catalog artifact hash mismatch')
    with np.load(path, allow_pickle=False) as archive:
        if set(archive.files) != {'zero_patterns', 'input_weights', 'counts', 'metadata'}:
            raise ValueError('Unexpected catalog archive schema')
        metadata = json.loads(str(archive['metadata']))
        if (metadata.get('complete') is not True or metadata.get('source_pins') != dependencies()
                or metadata.get('category_order') != 'input-weight-then-uint64-dense-zero-pattern'
                or metadata.get('input_weight_is_total_not_labelled') is not True):
            raise ValueError('Catalog scope or source bindings differ')
        catalog = DenseCategories(metadata['arity'], metadata['q'],
            _readonly(archive['zero_patterns'], np.uint64), _readonly(archive['input_weights'], np.uint8),
            _readonly(archive['counts'], np.uint64), metadata['positive_mass_relation_kernels'],
            tuple(sorted(metadata['source_pins'].items())))
    if (metadata.get('merged_categories') != len(catalog.counts)
            or metadata.get('raw_mass') != catalog.q**catalog.arity
            or metadata.get('dense_form_symbols') != [2, 3]
            or metadata.get('aes_embedding') != list(exact.base.AES_EMBEDDING)):
        raise ValueError('Catalog metadata disagrees with its exact arrays')
    return catalog.validate()


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--arity', type=int, default=6)
    parser.add_argument('--q', type=int, choices=(4, 256), default=256)
    parser.add_argument('--max-seconds', type=float, default=120)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists(): parser.error('Refuse to overwrite a catalog')
    if not 0 < args.max_seconds <= 120: parser.error('Use a wall budget at most120 seconds')
    started = perf_counter()
    def check_time():
        if perf_counter()-started >= args.max_seconds: raise TimeoutError('Catalog construction exhausted its budget')
    catalog = build_catalog(args.arity, args.q, check_time=check_time)
    report = save_catalog(catalog, args.output)
    report['elapsed_seconds'] = perf_counter()-started
    print(json.dumps(report, sort_keys=True), flush=True)
