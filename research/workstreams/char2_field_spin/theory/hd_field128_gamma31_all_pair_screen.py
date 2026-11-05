"""Isolated FLOAT all-pair adapter for deterministic gamma/state31.

State30 is an explicit replay control. Generic frozen local formulas already
accept31; this adapter does not alter any public whitelist or frozen module.
It reuses the complete all-pair coefficient and generic alias computations.
No result is a directed certificate, encoder implementation, or timing claim.
"""

import ast
from fractions import Fraction as F
from functools import lru_cache
from hashlib import sha256
from math import isfinite
from pathlib import Path

import numpy as np

import hd_field128_full_rs_gamma_screen as gamma
import hd_field128_gamma_all_pair_screen as frozen


old, alias = frozen.old, frozen.alias
FIELD_ORDER, STATES = frozen.FIELD_ORDER, (30, 31)
NOTE = Path(__file__).resolve().parent.parent/'HD_FIELD128_GAMMA31.md'


def construction(state=31, q=FIELD_ORDER):
    if type(state) is not int or state not in STATES or q != FIELD_ORDER:
        raise ValueError('Use GF(2^128), state31, or the explicit state30 replay control')
    return f'toeplitz-fullfield-full-rs64-{state}-gamma16-feedback-alternating-swapped-pairs-k20-rate-half-v1'


@lru_cache(2)
def fixed_maps(state=31):
    construction(state)
    a = gamma.frozen.lagrange(range(64), range(128, 128+state))
    h = gamma.frozen.lagrange(range(state, 2*state), range(state))
    return a, tuple(zip(*a)), h


def local_logs(z, q=FIELD_ORDER, *, state=31):
    construction(state, q)
    return gamma.local_logs(z, q, 64, state)


def exact_packet(j, z, q=FIELD_ORDER, *, state=31):
    construction(state, q)
    return gamma.exact_packet(j, z, q, 64, state)


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


class AliasCoefficients:
    """All65 local occupancies and all561² ordered histogram pairs."""

    def __init__(self, z, alpha, *, q=FIELD_ORDER, state=31, max_seconds=60., _budget=None):
        construction(state, q)
        if not isfinite(z) or not isfinite(alpha) or not 0 < z < 1 or not 0 < alpha <= 1:
            raise ValueError('Expected finite admissible z and alpha')
        budget = _budget or old.Budget(max_seconds)
        self.z, self.alpha, self.q, self.state = float(z), float(alpha), q, state
        self.joint, self.gamma_degree, self.all_pairs = True, 16, True
        self.composed = frozen.all_pair_coefficients(local_logs(self.z, q, state=state),
                                                      self.alpha, _budget=budget)
        self.single = np.full((4, 33, 33), frozen.NEG)
        self.packet = self
        budget.check()

    pair_logs = frozen.AliasCoefficients.pair_logs


class AliasBankScreen(alias.AliasBankScreen):
    """The generic alias interface, with explicit state31/control metadata."""

    def __init__(self, z, alpha, marker_log_odds, log_ratio, *, state=31, q=FIELD_ORDER,
                 distance=F(2, 5), max_seconds=60., _coefficients=None, _budget=None, **kwargs):
        name = construction(state, q)
        budget = _budget or old.Budget(max_seconds)
        coefficients = (_coefficients if _coefficients is not None else
                        AliasCoefficients(z, alpha, q=q, state=state, _budget=budget))
        if (type(coefficients) is not AliasCoefficients
                or (coefficients.z, coefficients.alpha, coefficients.q, coefficients.state)
                != (float(z), float(alpha), q, state)):
            raise ValueError('Coefficient producer must match this adapter, state, field, z, and alpha')
        self.state, self.construction = state, name
        self.pins, self.note_sha256 = source_pins(), sha256(NOTE.read_bytes()).hexdigest()
        super().__init__(z, alpha, marker_log_odds, log_ratio, q=q, distance=distance,
                         _coefficients=coefficients, _budget=budget, **kwargs)
        self.check_pins()

    def check_pins(self):
        directory = Path(__file__).resolve().parent
        if (any(sha256((directory/name).read_bytes()).hexdigest() != digest
                for name, digest in self.pins.items())
                or sha256(NOTE.read_bytes()).hexdigest() != self.note_sha256):
            raise RuntimeError('Pinned state31 adapter source or proof note changed')

    def metadata(self):
        return dict(construction=self.construction, state_dimension=self.state, inner_width=64,
                    emission_minimum=65-self.state, joint_rank_return_cap=True,
                    pair_swaps=True, state_persists_across_regions=True, full_field_random_coefficients=True,
                    feedback='gamma A transpose', gamma_extension_degree=16, gamma_base_field_order=256,
                    gamma_choice='fixed degree16 element; sparse-tower Y is permitted',
                    principal_nonsingularity='symmetric unit-diagonal determinant square factor has degree<=15',
                    fixed_nodes='GF256 additive subspace0..63', state_nodes=f'128..{127+self.state}',
                    refresh=f'D2 H_{self.state} D1; Lagrange0..{self.state-1} to{self.state}..{2*self.state-1}',
                    refresh_diagonal_fullfield_products_per_packet=2*self.state,
                    state_denominator=str(self.q**self.state-1),
                    refresh_density=str(F(self.q**self.state-1, (self.q-1)**self.state)),
                    CA_zero=True, initial_state='zero', final_state='unrestricted',
                    coefficient_method='all-pair-before-alpha-normalized-FFT-alias-FLOAT',
                    raw_occupancy_pair_products=65**2, raw_histogram_pair_products=561**2,
                    composed_histogram_terms=2145, single_components='identically zero',
                    maximum_packet_occupancy=64, maximum_paired_occupancy=128,
                    pair_slots_per_packet=32, pair_slots_per_region=4096,
                    packets_per_region=128, composed_pairs_per_region=64, regions=256,
                    all_pair_before_alpha=True, counts_outside_alpha=True, no_coefficient_truncation=True,
                    explicit_state30_control=self.state == 30, native_implementation=False,
                    source_sha256=dict(self.pins), proof_note_sha256=self.note_sha256,
                    old_numerical_bounds_reused=False, directed=False, whole_code_certificate=False,
                    arithmetic='FLOAT', scope='selected supports only; new state31 or explicit state30 control; no timing claim')

    def evaluate(self, weight, **options):
        return dict(super().evaluate(weight, **options), **self.metadata())

    def witness(self):
        return dict(super().witness(), **self.metadata())


def selected_screen(bank, weight, *, support_method='exact-profiles', optimize_potential=False,
                    max_evaluations=40, _budget=None):
    """Reuse the frozen complete profile/chronological or runs objective."""
    if type(bank) is not AliasBankScreen:
        raise TypeError('Expected this isolated state31/control bank')
    bank.check_pins()
    try:
        result = frozen.selected_screen(bank, weight, support_method=support_method,
                                        optimize_potential=optimize_potential,
                                        max_evaluations=max_evaluations, _budget=_budget)
    finally:
        bank.check_pins()
    return result
