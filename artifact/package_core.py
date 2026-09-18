"""Package only the selected encoder sources and their direct dependencies."""
import argparse
import hashlib
from pathlib import Path, PurePosixPath
import re
import subprocess
import zipfile

ROOT = Path(__file__).resolve().parents[1]
BARE = 'workstreams/bare_bch_rm2sub'
HALF = 'workstreams/inner_design/asymmetric/bch256/weight5/implementation'
QUARTER = 'workstreams/inner_design/routing_opt/deployment'
VENDOR = 'constructions/riffle_parityfanout31x33_bchperm_transpose_bitshuffle_splitstate_preaddmul_rm2sub_t128_s19/implementation/vendor'


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--hypercat-root', type=Path)
    parser.add_argument('--output', required=True, type=Path)
    parser.add_argument('--reviewer', action='store_true', help='Omit source provenance and scan for identifying text')
    parser.add_argument('--without-forward', action='store_true', help='Package only the two transposed encoders')
    parser.add_argument('--anonymous-owner-notice', action='store_true',
                        help='Use only with explicit rights-holder authorization for an anonymous review notice')
    args = parser.parse_args()
    if not args.without_forward and args.hypercat_root is None:
        parser.error('--hypercat-root is required unless --without-forward is selected')
    if args.anonymous_owner_notice and (not args.reviewer or args.without_forward):
        parser.error('--anonymous-owner-notice applies only to a reviewer package that includes forward code')
    if args.reviewer and not args.without_forward and not args.anonymous_owner_notice:
        parser.error('Forward code has an identifying copyright notice. Obtain rights-holder authorization '
                     'for --anonymous-owner-notice, or select --without-forward; do not silently strip the notice.')
    external = args.hypercat_root.resolve() if args.hypercat_root else None
    members, sources = {}, []

    def add(root, source, destination, label='SPIN'):
        path = (root / source).resolve(strict=True)
        if not path.is_relative_to(root.resolve()):
            raise ValueError(f'Source escapes root: {source}')
        relative = PurePosixPath(destination)
        if relative.is_absolute() or '..' in relative.parts or destination in members:
            raise ValueError(f'Invalid or duplicate archive path: {destination}')
        members[destination] = path.read_bytes()
        sources.append(f'{destination}\t{label}/{source}')

    for name in ('Spin.h', 'Inner.h', 'WorkspaceRouting.h'):
        add(ROOT, f'{BARE}/{name}', f'src/common/{name}')
    for name in ('BchCircuit.cpp', 'BchCircuit.h', 'SelectedMaps.h'):
        add(ROOT, f'{BARE}/generated/{name}', f'src/common/generated/{name}')
    for name in ('QuarterCircuit.cpp', 'QuarterCircuit.h'):
        add(ROOT, f'workstreams/rate_quarter_bch/implementation/generated/{name}', f'src/common/generated/{name}')
    for name in ('Weight5Spin.cpp', 'Weight5Inner.h', 'AsymmetricMap.h'):
        add(ROOT, f'{HALF}/{name}', f'src/half/{name}')
    for name in ('AsymmetricSpin.cpp', 'AsymmetricMap.h'):
        add(ROOT, f'{QUARTER}/{name}', f'src/routing/deployment/{name}')
    add(ROOT, 'workstreams/inner_design/asymmetric/AsymmetricInner.h', 'src/asymmetric/AsymmetricInner.h')
    for name in ('Defines.h', 'config.h', 'block.h', 'Bit.h'):
        add(ROOT, f'{VENDOR}/cryptoTools/Common/{name}', f'vendor/cryptoTools/Common/{name}')
    for source, destination in ((f'{HALF}/correctness.cpp', 'tests/half.cpp'),
                                (f'{HALF}/inner_identity.cpp', 'tests/inner.cpp'),
                                ('workstreams/rate_quarter_bch/implementation/correctness.cpp', 'tests/quarter.cpp'),
                                (f'{QUARTER}/page_test.cpp', 'tests/page.cpp'),
                                ('artifact/core_package/CMakeLists.txt', 'CMakeLists.txt')):
        add(ROOT, source, destination)
    if not args.without_forward:
        add(ROOT, 'artifact/core_package/bidirectional.cpp', 'tests/bidirectional.cpp')
        for name in ('Spin.cpp', 'Spin.h', 'Block.h', 'Inner.h', 'WorkspaceRouting.h',
                     'generated/BchCircuit.cpp', 'generated/BchCircuit.h',
                     'generated/BchForward.cpp', 'generated/BchForward.h', 'generated/SelectedMaps.h'):
            add(external, f'hypercat/native/spin/{name}', f'src/bidirectional/{name}', 'bidirectional-upstream')
        add(external, 'LICENSE', 'licenses/bidirectional-MIT.txt', 'bidirectional-upstream')
        if args.anonymous_owner_notice:
            license_name = 'licenses/bidirectional-MIT.txt'
            original = b'Copyright (c) 2026 Hypercat contributors'
            if members[license_name].count(original) != 1:
                raise ValueError('Unexpected license notice; explicit review is required')
            members[license_name] = members[license_name].replace(
                original, b'Copyright (c) 2026 the authors (identity withheld for anonymous review)')
    readme = 'README.reviewer.md' if args.reviewer or args.without_forward else 'README.md'
    add(ROOT, f'artifact/core_package/{readme}', 'README.md')
    if readme == 'README.reviewer.md':
        text = members['README.md'].decode('utf-8')
        if args.without_forward:
            text = re.sub(r'<!-- BEGIN FORWARD -->.*?<!-- END FORWARD -->\n?', '', text, flags=re.S)
        else:
            text = re.sub(r'<!-- (?:BEGIN|END) FORWARD -->\n?', '', text)
        members['README.md'] = text.encode('utf-8')
    if not args.reviewer:
        members['SOURCE_MAP.txt'] = ('Archive path\tOriginal source path\n' + '\n'.join(sorted(sources)) + '\n').encode()
    else:
        identifiers = re.compile(r'hypercat|ladnir|rindal|peceny|rachuri|raghuraman|'
                                 r'github\.com/|C:[/\\]Users|/home/|/Users/|\.codex|'
                                 r'peach48|devcore4|[\w.+-]+@[\w.-]+\.[A-Za-z]{2,}', re.I)
        for name, data in members.items():
            if identifiers.search(name) or identifiers.search(data.decode('utf-8')):
                raise ValueError(f'Identifying text in reviewer member: {name}')
    members['SHA256SUMS'] = ''.join(f'{hashlib.sha256(data).hexdigest()}  {name}\n'
                                  for name, data in sorted(members.items())).encode()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with zipfile.ZipFile(args.output, 'x', compression=zipfile.ZIP_DEFLATED, compresslevel=9) as archive:
        for name, data in sorted(members.items()):
            entry = zipfile.ZipInfo(f'spin-core/{name}', date_time=(1980,1,1,0,0,0))
            entry.create_system = 3
            entry.external_attr = 0o100644 << 16
            entry.compress_type = zipfile.ZIP_DEFLATED
            archive.writestr(entry, data)
    with zipfile.ZipFile(args.output) as archive:
        assert archive.testzip() is None
        assert set(archive.namelist()) == {f'spin-core/{name}' for name in members}
        for name, data in members.items():
            assert archive.read(f'spin-core/{name}') == data
    print(f'{len(members)} files; {args.output.stat().st_size:,} compressed bytes')
    print(f'SHA-256 {hashlib.sha256(args.output.read_bytes()).hexdigest()}')
    roots = [('SPIN', ROOT)]
    if not args.without_forward:
        roots.append(('bidirectional-upstream', external))
    for label, root in roots:
        revision = subprocess.check_output(['git', '-C', str(root), 'rev-parse', 'HEAD'], text=True).strip()
        print(f'{label} source revision: {revision}')
    print(args.output.resolve())


if __name__ == '__main__':
    main()
