"""Check the frozen K16 source identity, without running a benchmark or proof.

Uses Python 3.11+ only. Source hashes normalize CRLF to LF for fresh Windows
checkouts. The optional receipt check authenticates the original local archive;
it is not a substitute for the independent mathematical audit or fresh replay.
"""
import argparse
import hashlib
import json
from pathlib import Path
import tomllib


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]


def digest_source(path):
    return hashlib.sha256(path.read_bytes().replace(b'\r\n', b'\n')).hexdigest()


def verify(root=ROOT, receipt=None):
    manifest = tomllib.loads((HERE / 'source_manifest.toml').read_text(encoding='utf-8'))
    if manifest['checkpoint'] != 'k16-rs16-paired15-v1':
        raise ValueError('unexpected checkpoint')
    pins = manifest['sources']
    if len(pins) != manifest['source_count']:
        raise ValueError('incomplete source manifest')
    root = root.resolve()
    for name, expected in pins.items():
        path = (root / name).resolve()
        if not path.is_relative_to(root):
            raise ValueError('manifest path leaves repository')
        if digest_source(path) != expected:
            raise ValueError('source mismatch: ' + name)
    if receipt is not None:
        data = receipt.read_bytes()
        if hashlib.sha256(data).hexdigest() != manifest['receipt_sha256']:
            raise ValueError('archived receipt identity mismatch')
        record = json.loads(data)
        expected_geometry = dict(K=65536, N=131072, group_dimension=128,
                                 groups=512, packet_bits=4, regions=64,
                                 physical_t=64, state_bits=15,
                                 target_minimum_distance=13108,
                                 whole_code_certificate=True,
                                 all_occupancies_covered=True, target_met=True)
        for key, value in expected_geometry.items():
            if record[key] != value:
                raise ValueError('receipt geometry mismatch: ' + key)
        if set(record['occupancy_uppers']) != set(map(str, range(1, 513))):
            raise ValueError('receipt has incomplete occupancy coverage')
    return manifest


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root', type=Path, default=ROOT)
    parser.add_argument('--receipt', type=Path)
    args = parser.parse_args()
    manifest = verify(args.root, args.receipt)
    print(manifest['checkpoint'] + ': ' + str(manifest['source_count']) + ' source pins match')
    print('Archived receipt authenticated.' if args.receipt else
          'Source-only check; no numerical receipt or full proof replay performed.')
