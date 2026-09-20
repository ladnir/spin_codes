"""Merge discovery witnesses with checked provenance; never certifies a bound."""
import argparse
from pathlib import Path
import bridge as base
import screen_exponential_modes as caps_loader


def combine(names, tag):
    assert tag.isidentifier() and names
    caps, sources = caps_loader.latest_caps()
    best = {}
    hashes = {}
    configuration = bands = None
    for name in names:
        assert Path(name).name == name
        path = base.HERE/'generated'/name
        saved = base.read(path)
        if configuration is None:
            configuration, bands = saved['configuration'], saved['bands']
        assert saved['configuration'] == configuration and saved['bands'] == bands
        assert saved['used_caps'] == {str(w):caps[w] for w in base.WEIGHTS}
        for source, digest in saved['local_sha256'].items():
            assert base.sha(base.HERE/source) == digest
            assert source not in hashes or hashes[source] == digest
            hashes[source] = digest
        hashes[str(path.relative_to(base.HERE))] = base.sha(path)
        for row in saved['rows']:
            q = row['occupation']
            if q not in best or row['margin_bits_diagnostic'] > best[q]['margin_bits_diagnostic']:
                best[q] = {**row, 'discovery_source':name}
    hashes.update({str(p.relative_to(base.HERE)):base.sha(p) for p in [Path(__file__),Path(caps_loader.__file__)]+sources})
    base.write_new(base.HERE/'generated'/f'larger_witnesses_{tag}_screen.json', dict(
        status='COMBINED_LARGER_STATE_SCREEN_ONLY', configuration=configuration, bands=bands,
        rows=[best[q] for q in sorted(best)], used_caps={str(w):caps[w] for w in base.WEIGHTS},
        cap_receipts=[str(p.relative_to(base.HERE)) for p in sources], local_sha256=hashes,
        all_occupations_certified=False))
    print([(q,round(best[q]['margin_bits_diagnostic'],2)) for q in sorted(best)], flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--screens', nargs='+', required=True)
    parser.add_argument('--tag', required=True)
    args = parser.parse_args()
    combine(args.screens, args.tag)
