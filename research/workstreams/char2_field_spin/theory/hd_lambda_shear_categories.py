"""Exact joint support categories for at most three SHEAR source rows.

Each pair has two ordered input blocks of size a. Its output forms are
left-only dense, right-only dense, or dense on both blocks. Dense means
every nonzero block coefficient is GF4 symbol2 or3 (AES0xbc or0xbd).
For a3 this gives8+8+64=80 forms, not just64 full-input dense forms.

Enumerate the positive-mass GF4 relation kernels on2a inputs once. Merge
only kernels with equal total original input weight and equal zero pattern
on ALL universal forms. Store the80 bits as uint64+uint16, never truncate
the high16 bits. Raw injection masses are exact uint64, total256^(2a)<=2^48.
Every actual32-row paired map then gathers its form bits and merges the
joint output mask with the input weight. No independent-column law is used.

The relation kernel is local to one input pair; no kernel is held fixed
across pairs. The caller instead keeps the original source row set fixed,
uses its actual P/Q schedule in every pair, and restores the RS normalization.
This helper contains no regional provider, norm, or proof gate.
"""
from dataclasses import dataclass
from hashlib import sha256
import json
from math import comb
from pathlib import Path
from time import perf_counter

import numpy as np

import hd_lambda_double_exact6 as relation
import hd_lambda_column_shear_ranks as shear


DIRECTORY = Path(__file__).resolve().parent
PINS = {
    'hd_lambda_double_exact6.py': 'c6065ca7f6b9605aec096283d8ff2db94630692f99a4d2deed9e20b30d4dc960',
    'hd_lambda_column_shear_ranks.py': '4067013201cadef412ae647cc81e394e90b49d9403dfdd4713a739a191a96fa8',
}


def dependencies():
    for name, expected in PINS.items():
        if sha256((DIRECTORY/name).read_bytes()).hexdigest() != expected:
            raise ValueError('Frozen dependency changed: '+name)
    # relation_arrays is matrix-independent; its old double-matrix pins are
    # provenance only, NOT the construction used by paired_rows below.
    return {**relation.dependencies(), **shear.dependencies(), **PINS,
            Path(__file__).name: sha256(Path(__file__).read_bytes()).hexdigest()}


def _arity(a):
    if type(a) is not int or not 1 <= a <= 3:
        raise ValueError('Only source arities1..3 are in scope')


def universal_forms(a):
    _arity(a)
    block = tuple(sum((2 | (b >> i & 1)) << (2*i) for i in range(a)) for b in range(1 << a))
    return (block+tuple(form << (2*a) for form in block)
            +tuple(block[b & ((1 << a)-1)] | (block[b >> a] << (2*a)) for b in range(1 << (2*a))))


def _readonly(values, dtype):
    result = np.array(values, dtype=dtype, copy=True)
    result.setflags(write=False)
    return result


@dataclass(frozen=True)
class ShearCategories:
    arity: int
    q: int
    zero_low: np.ndarray
    zero_high: np.ndarray
    input_weights: np.ndarray
    counts: np.ndarray
    positive_mass_kernels: int
    source_pins: tuple

    def validate(self):
        _arity(self.arity)
        relation._parameters(2*self.arity, self.q)
        n, forms = len(self.counts), len(universal_forms(self.arity))
        if (not n or any(array.shape != (n,) for array in (self.zero_low, self.zero_high, self.input_weights))
                or self.zero_low.dtype != np.uint64 or self.zero_high.dtype != np.uint16
                or self.input_weights.dtype != np.uint8 or self.counts.dtype != np.uint64
                or np.any(self.counts == 0) or np.any(self.input_weights > 2*self.arity)
                or self.positive_mass_kernels < n):
            raise ValueError('Invalid exact category arrays')
        full_low, full_high = (1 << min(forms, 64))-1, (1 << max(0, forms-64))-1
        if np.any(self.zero_low > np.uint64(full_low)) or np.any(self.zero_high > full_high):
            raise ValueError('Zero pattern exceeds the universal-form bit width')
        order = np.lexsort((self.zero_low, self.zero_high, self.input_weights))
        if not np.array_equal(order, np.arange(n)):
            raise ValueError('Expected input-weight/high-zero/low-zero order')
        if n > 1 and np.any((self.input_weights[1:] == self.input_weights[:-1])
                & (self.zero_high[1:] == self.zero_high[:-1]) & (self.zero_low[1:] == self.zero_low[:-1])):
            raise ValueError('Duplicate merged category')
        if sum(map(int, self.counts)) != self.q**(2*self.arity):
            raise ValueError('Lost exact raw pair mass')
        for weight in range(2*self.arity+1):
            if sum(map(int, self.counts[self.input_weights == weight])) != comb(2*self.arity, weight)*(self.q-1)**weight:
                raise ValueError('Changed original input-weight mass')
        zero = self.input_weights == 0
        if (np.count_nonzero(zero) != 1 or int(self.counts[zero][0]) != 1
                or int(self.zero_low[zero][0]) != full_low or int(self.zero_high[zero][0]) != full_high
                or np.count_nonzero((self.zero_low == np.uint64(full_low)) & (self.zero_high == full_high)) != 1):
            raise ValueError('Only the unique zero input annihilates every universal form')
        return self


def build_catalog(arity=3, q=256, *, check_time=None):
    _arity(arity)
    pins = dependencies()
    basis, twice, masses = relation.relation_arrays(2*arity, q, check_time=check_time)
    weights = relation.input_weights(basis, twice)
    low = np.zeros(len(masses), dtype=np.uint64)
    high = np.zeros(len(masses), dtype=np.uint16)
    for index, form in enumerate(universal_forms(arity)):
        if check_time: check_time()
        zero = relation.zero_membership(form, basis, twice)
        if index < 64:
            low |= zero.astype(np.uint64) << np.uint64(index)
        else:
            high |= zero.astype(np.uint16) << (index-64)
    order = np.lexsort((low, high, weights))
    low, high, weights, masses = low[order], high[order], weights[order], masses[order]
    starts = np.r_[0, 1+np.flatnonzero((low[1:] != low[:-1]) | (high[1:] != high[:-1]) | (weights[1:] != weights[:-1]))]
    counts = np.add.reduceat(masses, starts, dtype=np.uint64)
    result = ShearCategories(arity, q, _readonly(low[starts], np.uint64), _readonly(high[starts], np.uint16),
        _readonly(weights[starts], np.uint8), _readonly(counts, np.uint64), len(masses), tuple(sorted(pins.items())))
    result.validate()
    if dependencies() != pins: raise ValueError('Sources changed while constructing catalog')
    return result


def joint_arrays(row_forms, catalog):
    """Sparse exact (weight,32-bit output mask,count); no2^32 allocation."""
    forms = tuple(row_forms)
    lookup = {form: i for i, form in enumerate(universal_forms(catalog.arity))}
    if (not 1 <= len(forms) <= 32 or any(type(form) is not int or form not in lookup for form in forms)):
        raise ValueError('Every output form must belong to the exact three-family dictionary')
    masks = np.zeros(len(catalog.counts), dtype=np.uint32)
    for row, form in enumerate(forms):
        bit = lookup[form]
        zero = ((catalog.zero_low >> np.uint64(bit)) & np.uint64(1)) if bit < 64 else ((catalog.zero_high >> (bit-64)) & 1)
        masks |= (1-zero).astype(np.uint32) << np.uint32(row)
    keys = (catalog.input_weights.astype(np.uint64) << np.uint64(32)) | masks.astype(np.uint64)
    order = np.argsort(keys)
    keys, masses = keys[order], catalog.counts[order]
    starts = np.r_[0, 1+np.flatnonzero(keys[1:] != keys[:-1])]
    counts = np.add.reduceat(masses, starts, dtype=np.uint64)
    keys = keys[starts]
    if int(counts.sum(dtype=np.uint64)) != catalog.q**(2*catalog.arity):
        raise ArithmeticError('Joint histogram lost raw pair mass')
    return (_readonly(keys >> np.uint64(32), np.uint8), _readonly(keys & np.uint64(0xffffffff), np.uint32),
            _readonly(counts, np.uint64))


def pair_arrays(source, column, catalog):
    source = tuple(source)
    if len(source) != catalog.arity or type(column) is not int or column not in range(0, 16, 2):
        raise ValueError('Catalog arity and an even chronological column are required')
    left, right = shear.paired_rows(source, column, 1, schedule=shear.SCHEDULE)
    weights, masks, counts = joint_arrays(left+right, catalog)
    if int(counts[0]) != 1 or int(weights[0]) != 0 or int(masks[0]) != 0 or np.count_nonzero(masks == 0) != 1:
        raise ArithmeticError('The injective pair map lost its unique zero input')
    return weights, masks, counts


def pair_counts(source, column, catalog):
    weights, masks, counts = pair_arrays(source, column, catalog)
    return {(int(w), int(mask)): int(n) for w, mask, n in zip(weights, masks, counts)}


def save_catalog(catalog, path):
    catalog.validate()
    metadata = dict(arity=catalog.arity, q=catalog.q, complete=True,
        coefficient_field_order=4, aes_embedding=relation.base.AES_EMBEDDING,
        universal_forms=list(universal_forms(catalog.arity)),
        positive_mass_kernels=catalog.positive_mass_kernels, merged_categories=len(catalog.counts),
        exact_raw_mass=catalog.q**(2*catalog.arity), category_order='input-weight-high-zero-low-zero',
        original_input_weight_is_total=True, common_cross_pair_relation_kernel_claimed=False,
        source_pins=dict(catalog.source_pins), whole_code_certificate=False)
    path = Path(path)
    with path.open('xb') as handle:
        np.savez_compressed(handle, zero_low=catalog.zero_low, zero_high=catalog.zero_high,
            input_weights=catalog.input_weights, counts=catalog.counts, metadata=np.array(json.dumps(metadata, sort_keys=True)))
    return dict(metadata, artifact=str(path), sha256=sha256(path.read_bytes()).hexdigest())


def load_catalog(path, expected_sha256):
    path = Path(path)
    if sha256(path.read_bytes()).hexdigest() != expected_sha256: raise ValueError('Pinned catalog hash differs')
    with np.load(path, allow_pickle=False) as data:
        if set(data.files) != {'zero_low', 'zero_high', 'input_weights', 'counts', 'metadata'}:
            raise ValueError('Unexpected catalog archive schema')
        metadata = json.loads(str(data['metadata']))
        if (metadata.get('complete') is not True or metadata.get('source_pins') != dependencies()
                or metadata.get('universal_forms') != list(universal_forms(metadata['arity']))
                or metadata.get('category_order') != 'input-weight-high-zero-low-zero'):
            raise ValueError('Catalog provenance or universal-form definition differs')
        catalog = ShearCategories(metadata['arity'], metadata['q'], _readonly(data['zero_low'], np.uint64),
            _readonly(data['zero_high'], np.uint16), _readonly(data['input_weights'], np.uint8),
            _readonly(data['counts'], np.uint64), metadata['positive_mass_kernels'], tuple(sorted(metadata['source_pins'].items())))
    if metadata.get('merged_categories') != len(catalog.counts) or metadata.get('exact_raw_mass') != catalog.q**(2*catalog.arity):
        raise ValueError('Catalog metadata and exact arrays disagree')
    return catalog.validate()


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--arity', type=int, choices=(1, 2, 3), default=3)
    parser.add_argument('--q', type=int, choices=(4, 256), default=256)
    parser.add_argument('--max-seconds', type=float, default=120)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists(): parser.error('Refuse to overwrite catalog')
    if not 0 < args.max_seconds <= 120: parser.error('Use at most120 seconds')
    started = perf_counter()
    def check_time():
        if perf_counter()-started >= args.max_seconds: raise TimeoutError('Catalog budget exhausted')
    result = save_catalog(build_catalog(args.arity, args.q, check_time=check_time), args.output)
    result['elapsed_seconds'] = perf_counter()-started
    print(json.dumps(result), flush=True)
