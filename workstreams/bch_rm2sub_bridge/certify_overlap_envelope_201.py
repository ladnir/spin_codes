"""Join the atlas with overlap-dependent intercepts, including failed gates."""
import argparse
from fractions import Fraction as F
from pathlib import Path
import bridge as base
import verify_extended_overlap_atlas as low_ledger
import certify_low_overlap_tile as low
import certify_mid_overlap_tile as mid


def run(replay=False):
    low_ledger.run('low', replay=False)
    low_path = base.HERE/'generated/extended_overlap_low_coverage.json'
    sources = [Path(__file__), low_path, Path(low.__file__), Path(mid.__file__)]
    intercepts = [F(0)]*161+[None]*40
    ranges = [list(range(740, 819, 2)), list(range(820, 901, 2))]
    total = 1056321
    if replay:
        for tile in range(6):
            low.run(tile, verify=True)
    for tile, (i, j) in enumerate([(0, 0), (0, 1), (1, 1)]):
        path = base.HERE/'generated'/f'mid_overlap_tile{tile:02d}.json'
        saved = base.read(path)
        for name, digest in saved['local_sha256'].items():
            assert base.sha(base.HERE/name) == digest
        assert saved['configuration'] == 't128_s15'
        assert saved['first_weights'] == ranges[i] and saved['second_weights'] == ranges[j]
        assert saved['overlaps'] == [r['overlap'] for r in saved['rows']] == list(range(161, 201))
        assert base.decode(saved['quadratic_coefficient']) == F(3, 25000)
        assert base.decode(saved['intercept_upper']) == max(base.decode(r['intercept_upper']) for r in saved['rows'])
        for row in saved['rows']:
            k, value = row['overlap'], base.decode(row['intercept_upper'])
            intercepts[k] = max(F(0), value, intercepts[k] if intercepts[k] is not None else F(0))
        count = len(ranges[i])*len(ranges[j])*40
        assert saved['type_count'] == count
        total += count*(2 if i != j else 1)
        sources.append(path)
        if replay:
            mid.run(tile, verify=True)
    assert total == 81*81*201 == 1318761 and all(v is not None for v in intercepts)
    result = dict(status='OUTWARD_201_OVERLAP_ENVELOPE_REPLAY' if replay else 'OUTWARD_201_OVERLAP_ENVELOPE_DOMAIN_LEDGER',
        actual_marginal_weights=list(range(740, 901, 2)), actual_overlaps=list(range(201)),
        type_count=total, quadratic_coefficient=base.encode(F(3, 25000)),
        intercept_by_overlap=[base.encode(v) for v in intercepts],
        form='log R(x,y,k) <= intercept[k] + (3/25000)*(k-x*y/8192)^2',
        all_extension_tiles_replayed_at512_bits=replay, full_second_moment_certified=False,
        local_sha256={str(p.relative_to(base.HERE)):base.sha(p) for p in sources})
    output = base.HERE/'generated'/f'overlap_envelope_201_{"replay" if replay else "coverage"}.json'
    if output.exists():
        assert result == base.read(output)
    else:
        base.write_new(output, result)
    print('201-overlap envelope:', total, 'types; max intercept', float(max(intercepts)),
          'all extension replays', replay, flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--replay', action='store_true')
    run(parser.parse_args().replay)
