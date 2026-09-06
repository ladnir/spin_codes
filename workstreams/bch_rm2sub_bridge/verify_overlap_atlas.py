"""Exact domain ledger and optional sequential512-bit replay for the atlas."""
import argparse
from pathlib import Path
from fractions import Fraction as F
import bridge as base
import certify_overlap_atlas_tile as atlas


def run(replay=False):
    pair_path = base.HERE/'generated/t128_s15_pair_spectrum.json'
    pairs = base.read(pair_path)
    for name, digest in pairs['local_sha256'].items():
        assert base.sha(base.HERE/name) == digest
    counts = {(row['first_weight'], row['second_weight'], row['sum_weight']):row['ordered_pairs']
              for row in pairs['rows']}
    assert all(counts[a, b, c] == counts[b, a, c] for a, b, c in counts)
    ranges = [list(range(740, 819, 2)), list(range(820, 901, 2))]
    blocks = [list(range(40, 80)), list(range(80, 120)), list(range(120, 161))]
    sources, total = [Path(__file__), Path(atlas.__file__), pair_path], 0
    maxima = []
    for tile in range(9):
        i, j = [(0, 0), (0, 1), (1, 1)][tile//3]
        path = base.HERE/'generated'/f'overlap_atlas_tile{tile:02d}.json'
        saved = base.read(path)
        for name, digest in saved['local_sha256'].items():
            assert base.sha(base.HERE/name) == digest
        assert saved['status'] == 'OUTWARD_ACTUAL_PAIR_OVERLAP_ENVELOPE_TILE'
        assert saved['configuration'] == 't128_s15' and saved['tile'] == tile
        assert saved['first_weights'] == ranges[i] and saved['second_weights'] == ranges[j]
        assert saved['overlaps'] == blocks[tile % 3]
        assert base.decode(saved['quadratic_coefficient']) == F(3, 25000)
        maximum = max(base.decode(row['intercept_upper']) for row in saved['rows'])
        assert maximum == base.decode(saved['intercept_upper']) and maximum < 0
        maxima.append(maximum)
        count = len(ranges[i])*len(ranges[j])*len(blocks[tile % 3])
        assert saved['type_count'] == count
        total += count*(2 if i != j else 1)
        sources.append(path)
        if replay:
            atlas.run(tile, verify=True)
    assert total == 81*81*121 == 793881
    result = dict(status='FULL_793881_TYPE_ATLAS_REPLAY' if replay else 'EXACT_793881_TYPE_ATLAS_DOMAIN_LEDGER',
        configuration='t128_s15', first_weights=list(range(740, 901, 2)), second_weights=list(range(740, 901, 2)),
        overlaps=list(range(40, 161)), type_count=total, intercept_upper=base.encode(max(maxima)),
        common_relaxed_intercept=base.encode(F(0)), quadratic_coefficient=base.encode(F(3, 25000)),
        pair_swap_symmetry_checked=True, full_second_moment_certified=False,
        all_tiles_replayed_at512_bits=replay,
        local_sha256={str(path.relative_to(base.HERE)):base.sha(path) for path in sources})
    path = base.HERE/'generated'/('overlap_atlas_replay.json' if replay else 'overlap_atlas_coverage.json')
    if path.exists():
        assert result == base.read(path)
    else:
        base.write_new(path, result)
    print('Atlas domain', total, 'types: log pair ratio <=(3/25000)(k-xy/8192)^2;',
          'all512-bit replays:', replay, flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--replay', action='store_true')
    run(parser.parse_args().replay)
