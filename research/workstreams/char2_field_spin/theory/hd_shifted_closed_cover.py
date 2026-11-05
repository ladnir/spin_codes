"""Freeze the explicit research cover after observing its fresh dual replays.

This inventory builder is not a proof verifier. Review the emitted manifest
and publish its hash separately before using the trusted-cache checker.
The --fresh mode of that checker recomputes the arithmetic from the pinned
proposals and frozen sparse schedule, without using cached bound values.
"""
from hashlib import sha256
import json
from pathlib import Path

from hd_shifted_checkpoint_export import PLANS as TAIL_PLANS
from hd_shifted_dense_checkpoint_export import PLANS as DENSE_PLANS
from hd_shifted_outer_density import DEFAULT_OUTER_MEASURE
from hd_shifted_swap_certificate import SWAP_CONSTRUCTION
from hd_shifted_verified_cover import TRUST_BASIS, validate_manifest


# Deliberately explicit: overlapping earlier checkpoints are not discovered
# and added by a glob. Each listed receipt covers all profiles at its weights.
POTENTIALS = (
    'hd_shifted_potential_band_272_257_279',
    'hd_shifted_potential_band_320_280_349',
    'hd_shifted_potential_band_400_350_433',
    'hd_shifted_potential_band_448_434_457',
    'hd_shifted_potential_band_512_458_566',
    'hd_shifted_potential_band_608_567_633',
    'hd_shifted_potential_band_672_634_699',
    'hd_shifted_potential_band_768_700_809',
    'hd_shifted_potential_band_864_810_896',
    'hd_shifted_potential_band_944_897_974',
    'hd_shifted_potential_band_992_975_1007',
    'hd_shifted_potential_band_1024',
    'hd_shifted_potential_band_1064_dense_1041_1088',
    'hd_shifted_potential_band_1108_dense_1089_1129',
    'hd_shifted_potential_band_1168_1130_1199',
    'hd_shifted_potential_band_1232_1200_1255',
    'hd_shifted_potential_band_1280_toeplitz_1256_1325',
    'hd_shifted_midpoint_1352_1326_1383',
    'hd_shifted_midpoint_1412_1384_1441',
    'hd_shifted_midpoint_1472_1442_1499',
    'hd_shifted_potential_band_1536_1500_1571',
    'hd_shifted_midpoint_1608_1572_1645',
    'hd_shifted_midpoint_1728_1646_1824',
    'hd_shifted_midpoint_1856_1825_1881',
    'hd_shifted_potential_band_2048_1882_2153',
    'hd_shifted_potential_band_2208_2154_2274',
    'hd_shifted_potential_band_2496_2275_2683',
    'hd_shifted_potential_band_2752_2684_2827',
    'hd_shifted_potential_band_3072_2828_3285_toeplitz',
    'hd_shifted_potential_band_3304_3286_3325',
    'hd_shifted_potential_band_3584_3326_3768',
    'hd_shifted_potential_band_3936_3769_4050',
    'hd_shifted_potential_band_4080_4051_4067',
)
PROPOSAL_NAMES = {
    'hd_shifted_potential_band_1024': 'hd_shifted_potential_band_1024_frontier.json',
    'hd_shifted_potential_band_2496_2275_2683':
        'hd_shifted_potential_band_2496_choice0_2154_2827_selected.json',
}
LEGACY_WITHOUT_PROPOSAL_HASH = frozenset((
    'hd_shifted_potential_band_1280_toeplitz_1256_1325',
    'hd_shifted_potential_band_1536_1500_1571',
))


def freeze_inventory(directory):
    directory = Path(directory)
    def pin(name):
        return dict(path=name, sha256=sha256((directory/name).read_bytes()).hexdigest())
    def receipts(stem):
        return {str(bits): pin(f'{stem}_directed{bits}.json') for bits in (256, 384)}
    entries = []
    for stem in POTENTIALS:
        receipt = json.loads((directory/f'{stem}_directed256.json').read_bytes())
        entry = dict(id=stem, kind='potential-band', weight_range=receipt['weight_range'],
            observed_fresh_dual_replay=True, receipts=receipts(stem),
            proposal=pin(PROPOSAL_NAMES.get(stem, stem+'.json')),
            outer_measure=receipt.get('outer_measure', DEFAULT_OUTER_MEASURE))
        if stem in LEGACY_WITHOUT_PROPOSAL_HASH:
            entry['legacy_missing_proposal_pin'] = True
        entries.append(entry)
    for kind, interval in (('sparse', [1, 256]), ('dense', [4096, 5119]), ('tail', [5120, 8192])):
        entry = dict(id=kind+'-checkpoint', kind='checkpoint', checkpoint_kind=kind,
            weight_range=interval, observed_fresh_dual_replay=True,
            receipts=receipts('hd_shifted_'+kind+'_checkpoint'))
        if kind != 'sparse':
            entry['replay_source_pin'] = pin('hd_shifted_cover_ledger.py')
            specs = [('swap-plan', *row) for row in DENSE_PLANS] if kind == 'dense' else [
                ('unsigned-plan', name, low, high, expected) for _, name, low, high, expected in TAIL_PLANS]
            entry['plan_bands'] = []
            for plan_kind, name, low, high, expected in specs:
                reference = pin(name)
                if reference['sha256'] != expected:
                    raise ValueError(f'Established proposal changed: {name}')
                entry['plan_bands'].append(dict(label=name, kind=plan_kind,
                    weight_range=[low, high], proposal=reference))
        entries.append(entry)
    entries.append(dict(id='final-boundary-4068-4095', kind='swap-plan',
        weight_range=[4068, 4095], observed_fresh_dual_replay=True,
        proposal=pin('hd_shifted_swap_4068_4095_frozen_plan3.json'),
        receipts=receipts('hd_shifted_swap_4068_4095')))
    entries.sort(key=lambda row: row['weight_range'])
    manifest = dict(schema_version=1, construction=SWAP_CONSTRUCTION, trust_basis=TRUST_BASIS,
        target=dict(source_dimension=2**20, output_length=2**21, rate='1/2', field_order=256,
            forbidden_weight=734003, minimum_margin_bits=40, group=4, inner_rank=4,
            setup='ideal independent group-column mixers, route permutations, swaps, and inner diagonals',
            scalar_extension='same GF256 scalar map on all sixteen byte components; no factor16 failure union',
            native_setup_qualification='shared-seed MT is a benchmark implementation, not a certificate for seed113'),
        entries=entries,
        implementation_pins=[pin(name) for name in (
            'hd_shifted_cover_ledger.py', 'hd_shifted_swap_ledger.py',
            'hd_shifted_potential_band.py', 'hd_shifted_run_potential.py',
            'hd_shifted_outer_density.py', 'hd_shifted_swap_sparse.py',
            'hd_shifted_verified_cover.py', 'hd_shifted_checkpoint_export.py',
            'hd_shifted_dense_checkpoint_export.py', 'hd_shifted_closed_cover.py')])
    validate_manifest(manifest)
    intervals = [row['weight_range'] for row in entries]
    if intervals[0][0] != 1 or intervals[-1][1] != 8192 or any(
            left[1]+1 != right[0] for left, right in zip(intervals, intervals[1:])):
        raise ValueError('The frozen inventory must partition every weight1..8192')
    return manifest


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        parser.error('Frozen manifests are never overwritten')
    result = freeze_inventory(Path(__file__).parent)
    with args.output.open('x', encoding='utf-8') as target:
        json.dump(result, target, allow_nan=False, indent=2)
        target.write('\n')
    print(json.dumps(dict(path=str(args.output.resolve()), entries=len(result['entries']),
                          manifest_sha256=sha256(args.output.read_bytes()).hexdigest())), flush=True)
