"""Fresh outer authentication and explicitly scoped S16 dense models.

The saved S19 file supplies a rational comparison-mixture proposal only.
Every shell of that proposal is checked against newly authenticated outer
counts. No S19 inner operator or numerical endpoint is imported.
"""
import copy
from fractions import Fraction as Q
import hashlib
import json
from pathlib import Path
import sys

HERE = Path(__file__).resolve().parent
PACKED = HERE.parent
PERMUTATION = PACKED.parent
for path in (PERMUTATION, PERMUTATION/'gf16_packets', PACKED):
    sys.path.insert(0, str(path))

from flint import arb, ctx
import cell_search
from outer_hill_intersection import authenticated_bch_cdf
from monotone import transport_shells
from local_models import full_block
import shared_mixture
import scalar_cover as sc
import variance_partition
import kernel_maps
import adapters_regional as regional


def fingerprint(value):
    # Spectra use integer keys before saving and string keys after JSON load.
    # Normalize first, so this identity survives a saved receipt round trip.
    canonical = json.loads(json.dumps(value))
    return hashlib.sha256(json.dumps(canonical, sort_keys=True, separators=(',', ':')).encode()).hexdigest()


class Model(sc.Model):
    """Use the unchanged outer comparison with separately bound local maps."""

    def __init__(self, mixture, data, threshold=209715, minimum_groups=33,
                 variance_bins=16):
        super().__init__(shared_mixture.as_components(mixture), data, threshold,
            minimum_groups, Q(3, 16), inner=kernel_maps, variance_shuffle=True,
            variance_bins=variance_bins, regional_count=False)
        # The base class only recognizes its legacy inner by module name.
        # Never spoof that name or route S16 data into its actual-map helpers.
        self.regional_count = True
        self.proposal_stop_bits = 54

    def outward(self, cell, witness):
        if 'regional_count_parts' in witness:
            return regional.outward(self, cell, witness)[0]
        return super().outward(cell, witness)

    def candidates(self, cell, scales=(Q(1), Q(1, 2), Q(1, 4)), callback=None):
        base = self.propose_with(cell, self.tilt)
        best = base
        if callback:
            callback(dict(method='iid', score=float(base[0])), best)
        if not 0 < cell[0] <= cell[1] < 1:
            return best
        variance = variance_partition.propose(self, cell, base)
        if variance[0] < best[0]:
            best = variance
        if callback:
            callback(dict(method='variance', score=float(variance[0])), best)
        if cell[1]-cell[0] > Q(1, 1024):
            return best
        witness = copy.deepcopy(variance[1])
        witness.update(regional_tilted_atom=True, regional_fine_tilts=True,
            regional_tilted_variance=True, regional_direct_counts=True,
            regional_exact_zero=True)
        if self.data.get('distribution') != 'uniform_gl':
            witness.update(regional_lazy_density_through=6,
                regional_joint_return_through=3,
                regional_feedback_classes_from=3,
                regional_feedback_classes_through=32,
                regional_feedback_uniform_classes=True,
                regional_feedback_uniform_replace=True)
        lam = Q(witness['parameters'][0])
        for scale in scales:
            if Q(scale) <= 0:
                raise ValueError('positive output-tilt scale required')
            trial = copy.deepcopy(witness)
            trial['parameters'][0] = str(lam*Q(scale))
            score, checked = regional.propose(self, cell, trial)
            if score < best[0]:
                best = (score, checked)
            if callback:
                callback(dict(method='regional', scale=str(scale),
                              output_tilt=trial['parameters'][0], score=float(score)), best)
        return best

    def proposal(self, cell):
        return self.candidates(cell)


def fresh_model(source, *, feedback='bch16', seed=0, refresh='uniform',
                updates=8, precision=192, variance_bins=16, distance=Q(1, 10),
                feedback_basis_seed=None):
    if (type(precision) is not int or precision < 128
            or type(variance_bins) is not int or not 1 <= variance_bins <= 64
            or not 0 < Q(distance) < Q(1, 2)
            or feedback_basis_seed is not None and type(feedback_basis_seed) is not int):
        raise ValueError('precision >=128, valid variance partition and exact distance required')
    source = Path(source).resolve()
    raw = source.read_bytes()
    old = json.loads(raw)['scope']
    if (old.get('K') != 1 << 20 or old.get('N') != 1 << 21
            or old.get('block_width') != 8 or old.get('minimum_groups') != 33
            or old.get('maximum_groups') != 2048
            or old.get('comparison') != 'direct-expected-shell-majorant'):
        raise ValueError('matching K20 canonical GL32 comparison-mixture proposal required')
    print('S16 freshly authenticating outer count premises', flush=True)
    caps, premises = authenticated_bch_cdf()
    comparison = transport_shells(premises['canonical_cdf'], full_block(8))
    if (premises != old['outer_premises']
            or cell_search.dense.fingerprint(caps) != old['expected_cdf_sha256']
            or cell_search.dense.fingerprint(comparison) != old['comparison_caps_sha256']):
        raise ValueError('fresh outer count premises differ from proposed mixture scope')
    mixture, verification = cell_search.dense.exact_mixture(comparison, old['mixture'])
    print('S16 outer mixture checked; preparing actual inner', feedback, refresh, flush=True)
    data, maps = kernel_maps.prepare(feedback=feedback, seed=seed,
                                     refresh=refresh, updates=updates)
    if feedback_basis_seed is not None:
        import basis_feedback
        basis, attempts = basis_feedback.seeded_basis(feedback_basis_seed)
        columns, inverse = basis_feedback.transform(data['columns'], basis)
        data = kernel_maps.prepare_maps(data['map_images'], columns, bits=16,
            updates=updates, distribution=data['distribution'])
        maps = dict(maps, feedback_columns=list(columns),
            feedback_definition='B times the original declared feedback map; expansion unchanged',
            feedback_basis=dict(seed=feedback_basis_seed, rows_hex=list(map(hex, basis)),
                inverse_rows_hex=list(map(hex, inverse)), sampling_attempts=attempts,
                original_feedback=feedback, kernel_and_row_code_preserved=True),
            map_sha256=data['map_sha256'])
    maps = json.loads(json.dumps(maps))
    ctx.prec = precision
    threshold = int(Q(distance)*(1 << 21))
    model = Model(mixture, data, threshold, variance_bins=variance_bins)
    sampling = dict(kind='uniform_gl', bits=16) if refresh == 'uniform' else dict(
        kind='transvections', bits=16, updates=updates)
    scope = dict(schema='packed-gl32-s16-dense-context-2',
        K=1 << 20, N=1 << 21, t=128, s=16, block_width=8,
        outer='fixed BCH[256,128]', block_rows=4, groups=2048,
        mixing='independent uniform GL32 on each canonical four-row/eight-column block',
        routing='shared uniform column shuffle per four-row group; independent uniform packet shuffle per region',
        inner_recurrence='y=x+Aq; next_q=Mq+Cx; zero initial state; no flush',
        sampling=sampling, maps=maps, maps_sha256=fingerprint(maps),
        birth_density=data.get('birth_density', 'classes'),
        feedback=feedback, feedback_seed=seed,
        feedback_basis_seed=feedback_basis_seed,
        distance=str(Q(distance)), threshold=threshold,
        minimum_groups=33, maximum_groups=2048,
        comparison='direct-expected-shell-majorant',
        base_tilt=str(model.tilt), variance_bins=variance_bins,
        root=list(map(str, model.root)), outer_premises=premises,
        expected_cdf_sha256=cell_search.dense.fingerprint(caps),
        comparison_caps_sha256=cell_search.dense.fingerprint(comparison),
        mixture=[dict(mass=str(c), activity=str(p)) for c, p in mixture],
        mixture_verification=verification, outer_authenticated=True)
    if source.read_bytes() != raw:
        raise ValueError('mixture proposal source changed during preparation')
    return model, scope, dict(path=str(source), sha256=hashlib.sha256(raw).hexdigest())
