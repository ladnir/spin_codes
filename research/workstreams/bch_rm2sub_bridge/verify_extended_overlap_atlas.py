"""Exact extended-domain ledger; optional sequential512-bit tile replay."""
import argparse
from fractions import Fraction as F
from pathlib import Path
import bridge as base
import verify_overlap_atlas as central
import certify_low_overlap_tile as low
import certify_mid_overlap_tile as mid


def run(stage, replay=False):
    assert stage in ('low', 'all')
    central.run(replay=False)
    central_path = base.HERE/'generated/overlap_atlas_replay.json'
    checked = base.read(central_path)
    assert checked['all_tiles_replayed_at512_bits']
    for name, digest in checked['local_sha256'].items():
        assert base.sha(base.HERE/name) == digest
    sources = [Path(__file__), central_path, Path(low.__file__)]
    ranges = [list(range(740, 819, 2)), list(range(820, 901, 2))]
    total = 793881
    for family, module, count in [('low', low, 6)]+([('mid', mid, 3)] if stage == 'all' else []):
        sources.append(Path(module.__file__))
        for tile in range(count):
            path = base.HERE/'generated'/f'{family}_overlap_tile{tile:02d}.json'
            saved = base.read(path)
            for name, digest in saved['local_sha256'].items():
                assert base.sha(base.HERE/name) == digest
            i, j = [(0, 0), (0, 1), (1, 1)][tile//2 if family == 'low' else tile]
            ks = list(range((tile % 2)*20, (tile % 2+1)*20)) if family == 'low' else list(range(161, 201))
            assert saved['first_weights'] == ranges[i] and saved['second_weights'] == ranges[j]
            assert saved['overlaps'] == ks and [r['overlap'] for r in saved['rows']] == ks
            assert saved['configuration'] == 't128_s15'
            assert base.decode(saved['quadratic_coefficient']) == F(3, 25000)
            intercept = max(base.decode(r['intercept_upper']) for r in saved['rows'])
            assert intercept == base.decode(saved['intercept_upper']) and intercept <= 0
            assert saved['common_quadratic_envelope_passed']
            count_types = len(ranges[i])*len(ranges[j])*len(ks)
            assert count_types == saved['type_count']
            total += count_types*(2 if i != j else 1)
            sources.append(path)
            if replay:
                module.run(tile, verify=True)
    maximum = 160 if stage == 'low' else 200
    assert total == 81*81*(maximum+1)
    result = dict(status='EXTENDED_ACTUAL_PAIR_ATLAS_REPLAY' if replay else 'EXTENDED_ACTUAL_PAIR_ATLAS_DOMAIN_LEDGER',
        stage=stage, actual_marginal_weights=list(range(740, 901, 2)), actual_overlaps=list(range(maximum+1)),
        type_count=total, quadratic_coefficient=base.encode(F(3, 25000)),
        central_replay_checked=True, extension_tiles_replayed_at512_bits=replay,
        full_second_moment_certified=False,
        local_sha256={str(p.relative_to(base.HERE)):base.sha(p) for p in sources})
    output = base.HERE/'generated'/f'extended_overlap_{stage}_{"replay" if replay else "coverage"}.json'
    if output.exists():
        assert result == base.read(output)
    else:
        base.write_new(output, result)
    print('Extended overlap atlas', total, 'types, k0..'+str(maximum), 'extension replay:', replay, flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--stage', choices=['low', 'all'], required=True)
    parser.add_argument('--replay', action='store_true')
    args = parser.parse_args()
    run(args.stage, args.replay)
