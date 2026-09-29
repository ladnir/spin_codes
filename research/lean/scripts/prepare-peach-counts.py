"""Package the checked dependency closure for an isolated Peach assembly run.

Project .olean files reuse the local kernel-checked replay. The final assembly
is compiled on Peach; this is not a fresh replay of every imported certificate.
"""
import hashlib
import json
from pathlib import Path
import re
import subprocess
import tarfile

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / 'scripts/map_data'
REMOTE = '/tmp/spin-lean-peach'


def main():
    pending = ['SpinCodes.Structured.' + name for name in
               ['MapSpectrumBridge', 'FiberFrozen', 'ConcreteShells'] +
               [f'FiberNumericsData.Weight{j}' for j in range(129)]]
    seen = set()
    while pending:
        module = pending.pop()
        if module in seen:
            continue
        seen.add(module)
        source = ROOT / (module.replace('.', '/') + '.lean')
        pending.extend(re.findall(r'^import (SpinCodes\.\S+)', source.read_text(encoding='utf-8'), re.M))
    records = []
    archive = DATA / 'peach-count-dependencies.tar.gz'
    with tarfile.open(archive, 'w:gz') as tar:
        for module in sorted(seen):
            relative = Path(module.replace('.', '/'))
            source = ROOT / relative.with_suffix('.lean')
            obj = ROOT / '.lake/build/lib/lean' / relative.with_suffix('.olean')
            assert obj.exists(), module
            files = [source] + sorted(obj.parent.glob(obj.name + '*'))
            for path in files:
                rel = path.relative_to(ROOT).as_posix()
                tar.add(path, arcname=rel)
                records.append(dict(path=rel, sha256=hashlib.sha256(path.read_bytes()).hexdigest()))
    report = dict(scope='Reuse locally checked dependencies for remote final assembly; not a fresh numerical replay.',
                  remote=REMOTE, modules=len(seen), files=records,
                  archive_sha256=hashlib.sha256(archive.read_bytes()).hexdigest())
    (DATA / 'peach_dependencies.json').write_text(json.dumps(report, indent=2) + '\n')
    # The installed WinSCP rejects the skill's SHA256 fingerprint syntax. Plink
    # accepts that exact pinned fingerprint; stream bytes without shell conversion.
    command = [r'C:\Users\peter\OneDrive\tools\plink.exe', '-ssh', '-batch',
               '-P', '9022', '-l', 'prindal', '-i',
               r'C:\Users\peter\OneDrive\keys\key_cat-9-3-2025.ppk', '-hostkey',
               'ssh-ed25519 255 SHA256:6ZRPOoZl4NGXghNQToErE0anQ0PkFggvYu7eBjXbK9o',
               'peach48.devcore4.com', f'cat > {REMOTE}/{archive.name}']
    with archive.open('rb') as stream:
        subprocess.run(command, stdin=stream, check=True, timeout=300)
    print(f'Uploaded {len(seen)} checked dependency modules; archive SHA256 {report["archive_sha256"]}')


if __name__ == '__main__':
    main()
