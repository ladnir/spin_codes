"""Immutable marker subsets for the frozen generic directed band driver only.

Each retained marker independently supplies a valid entrywise upper bound.
Taking the minimum over fewer markers can only weaken that bound. No atlas,
coefficient, normalization, or construction is recomputed or relabelled here.
Content hashes detect mutation; they do not authenticate the bank producer.
"""

import ast
from dataclasses import dataclass
from hashlib import sha256
import json
from pathlib import Path
from types import MappingProxyType

import numpy as np

import hd_field128_gamma_all_pair_certificate as provider
import hd_field128_lift_alias_cover as cover


def source_pins():
    directory, pending, pins = Path(__file__).resolve().parent, [Path(__file__).name], {}
    while pending:
        name = pending.pop()
        if name in pins:
            continue
        raw = (directory/name).read_bytes()
        pins[name] = sha256(raw).hexdigest()
        for node in ast.walk(ast.parse(raw)):
            modules = ([item.name for item in node.names] if isinstance(node, ast.Import)
                       else [node.module] if isinstance(node, ast.ImportFrom) and node.module else [])
            for module in modules:
                target = module.split('.')[0]+'.py'
                if (directory/target).is_file() and target not in pins:
                    pending.append(target)
    return dict(sorted(pins.items()))


def _json(value):
    return json.dumps(value, sort_keys=True, separators=(',', ':'), allow_nan=False)


def _data_digest(markers, atlases, numerators, grid, scalars):
    if not len(markers) == len(atlases) == len(numerators):
        raise RuntimeError('Marker, atlas, and fallback rows do not match')
    rows = []
    for atlas, fallback in zip(atlases, numerators):
        array = atlas.upper
        if (type(atlas) is not provider.interval.RegionalAliases
                or not isinstance(array, np.ndarray) or array.dtype != np.dtype('float64')
                or array.shape != (4, grid, grid) or not array.flags.c_contiguous
                or len(atlas.factors) != 4 or len(fallback) != 4):
            raise RuntimeError('Expected the actual all-pair interval-DFT atlas layout')
        rows.append((sha256(memoryview(array).cast('B')).hexdigest(),
                     tuple(value.repr() for value in atlas.factors),
                     tuple(value.repr() for value in fallback)))
    encoded = (tuple(tuple(map(str, row)) for row in markers), grid,
               tuple(value.repr() for value in scalars), rows)
    return sha256(_json(encoded).encode()).hexdigest()


def _parent_signature(parent):
    # Exact Arb repr includes the full dyadic midpoint/radius, unlike str().
    data = _data_digest(parent.markers, parent.atlases, parent.numerators, parent.grid,
                        (parent.constant, parent.slope, parent.density_log))
    metadata = provider.AliasBank.metadata(parent)
    return data, sha256(_json(metadata).encode()).hexdigest()


@dataclass(frozen=True, slots=True, init=False)
class MarkerSubsetBank:
    """A duck-typed band interface, deliberately not an AliasBank subclass.

    The selected arrays are copied once into immutable bytes-backed storage.
    The parent is never modified. All parent atlas/fallback values and witness
    metadata are checked at entry/exit, outside the per-weight extraction loop.
    Use band_receipt below for checks on exceptional exits as well.
    """

    _parent: object
    indices: tuple
    markers: tuple
    atlases: tuple
    numerators: tuple
    precision: int
    grid: int
    constant: object
    slope: object
    pins: object
    _parent_signature: tuple
    _subset_digest: str

    def __init__(self, parent, indices):
        if type(parent) is not provider.AliasBank:
            raise TypeError('Expected the actual gamma all-pair AliasBank, not an adapter or subclass')
        indices = tuple(indices)
        if (not indices or any(type(index) is not int or not 0 <= index < len(parent.markers)
                               for index in indices) or len(set(indices)) != len(indices)):
            raise ValueError('Expected a nonempty ordered subset of distinct valid marker indices')
        provider.AliasBank.check_pins(parent)
        pins, signature = source_pins(), _parent_signature(parent)
        selected = tuple(parent.atlases[index] for index in indices)
        # Bytes-backed arrays cannot be made writable with setflags(). Arb
        # scalar objects have value semantics; arithmetic returns new objects.
        atlases = tuple(provider.interval.RegionalAliases(
            np.frombuffer(atlas.upper.tobytes(), dtype=np.float64).reshape(atlas.upper.shape),
            tuple(atlas.factors)) for atlas in selected)
        values = dict(_parent=parent, indices=indices,
                      markers=tuple(parent.markers[index] for index in indices),
                      atlases=atlases,
                      numerators=tuple(tuple(parent.numerators[index]) for index in indices),
                      precision=parent.precision, grid=parent.grid,
                      constant=parent.constant, slope=parent.slope,
                      pins=MappingProxyType(pins), _parent_signature=signature)
        for name, value in values.items():
            object.__setattr__(self, name, value)
        digest = _data_digest(self.markers, self.atlases, self.numerators, self.grid,
                              (self.constant, self.slope, parent.density_log))
        object.__setattr__(self, '_subset_digest', digest)
        self.check_pins()

    def check_pins(self):
        provider.AliasBank.check_pins(self._parent)
        directory = Path(__file__).resolve().parent
        if any(sha256((directory/name).read_bytes()).hexdigest() != digest
               for name, digest in self.pins.items()):
            raise RuntimeError('A pinned marker-subset or generic band source changed')
        if _parent_signature(self._parent) != self._parent_signature:
            raise RuntimeError('Parent bank witness, atlas, fallback, or scalar data changed')
        if (not self.indices or len(set(self.indices)) != len(self.indices)
                or any(type(index) is not int or not 0 <= index < len(self._parent.markers)
                       for index in self.indices)
                or self.markers != tuple(self._parent.markers[index] for index in self.indices)
                or self.precision != self._parent.precision or self.grid != self._parent.grid):
            raise RuntimeError('Immutable marker-subset selection or numerical context changed')
        if _data_digest(self.markers, self.atlases, self.numerators, self.grid,
                        (self.constant, self.slope, self._parent.density_log)) != self._subset_digest:
            raise RuntimeError('Immutable marker-subset data changed')

    def metadata(self):
        self.check_pins()
        # Detach every nested container: editing a returned receipt must not
        # alter either the parent metadata or a later receipt from this view.
        result = json.loads(_json(provider.AliasBank.metadata(self._parent)))
        parent_witness = result['witness']
        result['witness'] = dict(parent_witness, markers=[list(map(str, row)) for row in self.markers])
        result['regional_diagnostics'] = [result['regional_diagnostics'][index] for index in self.indices]
        result.update(marker_subset=True, subset_marker_indices=list(self.indices),
                      subset_marker_count=len(self.markers), parent_marker_count=len(self._parent.markers),
                      parent_witness=parent_witness,
                      parent_bank_source_sha256=dict(self._parent.pins),
                      parent_bank_data_sha256=self._parent_signature[0],
                      parent_metadata_sha256=self._parent_signature[1],
                      subset_bank_data_sha256=self._subset_digest,
                      source_sha256=dict(self.pins), bank_source_sha256=dict(self.pins),
                      marker_subset_source_sha256=dict(self.pins),
                      bank_pin_scope='subset wrapper, generic band driver, and complete parent closure',
                      subset_rule='entrywise minimum over selected valid marker upper bounds only',
                      supported_consumer='hd_field128_lift_alias_cover.band_receipt',
                      parent_mutation_checks='exact witness/diagnostics, all atlas and fallback values, scalar constants',
                      bank_content_hash_scope='mutation detection only; not producer authentication',
                      coefficient_recomputation=False, parent_modified=False)
        return result


def band_receipt(bank, first, last, ratios, **options):
    """Call the frozen generic driver with full checks before and after it."""
    if type(bank) is not MarkerSubsetBank:
        raise TypeError('Expected a MarkerSubsetBank')
    bank.check_pins()
    try:
        result = cover.band_receipt(bank, first, last, ratios, **options)
    finally:
        bank.check_pins()
    return dict(result, subset_source_verified_at_entry=True, subset_source_verified_at_exit=True)
