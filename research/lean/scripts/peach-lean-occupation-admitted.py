"""Compile one project module on the prepared Peach workspace and fetch its objects.

Dependencies reuse existing local checked .olean files, transferred only when
their hashes differ from the recorded remote snapshot. All SSH calls pin the
configured host key. This command does not constitute a full dependency replay.
"""
import argparse
import hashlib
import json
from pathlib import Path
import re
import shlex
import subprocess
import sys
import tarfile
import time

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT/'scripts/map_data'
REMOTE = '/tmp/spin-lean-peach/project'
SSH = [r'C:\Users\peter\OneDrive\tools\plink.exe', '-ssh', '-batch', '-P', '9022',
       '-l', 'prindal', '-i', r'C:\Users\peter\OneDrive\keys\key_cat-9-3-2025.ppk',
       '-hostkey', 'ssh-ed25519 255 SHA256:6ZRPOoZl4NGXghNQToErE0anQ0PkFggvYu7eBjXbK9o',
       'peach48.devcore4.com']


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def admitted_lean(relative, output):
    # Only occupation semantic imports enter the two-slot semaphore.
    lean = f'lake env lean -j1 -M 20000 -o {shlex.quote(output)} {shlex.quote(relative)}'
    if not re.fullmatch(r'SpinCodes/Structured/DenseOccupationFixed(?:W[0-9]{3}|B[0-9]{3}(?:Rate)?|First)\.lean', relative):
        return lean
    gate = """import fcntl, os, subprocess, sys, time
os.makedirs('/tmp/spin-occupation-semantic-admission', exist_ok=True)
while True:
 for slot in range(2):
  f=open('/tmp/spin-occupation-semantic-admission/slot'+str(slot),'a')
  try: fcntl.flock(f,fcntl.LOCK_EX|fcntl.LOCK_NB)
  except BlockingIOError: f.close();continue
  sys.exit(subprocess.call(sys.argv[1:]))
 time.sleep(.05)
"""
    return 'python3 -c '+shlex.quote(gate)+' '+lean


def main():
    sys.stdout.reconfigure(encoding='utf-8')
    parser = argparse.ArgumentParser()
    parser.add_argument('source')
    args = parser.parse_args()
    source = (ROOT / args.source).resolve()
    relative = source.relative_to(ROOT).as_posix()
    assert relative.startswith('SpinCodes/') and source.suffix == '.lean'
    cachefile = DATA/'peach_extra_files.json'
    known = {r['path']: r['sha256'] for r in json.loads((DATA/'peach_dependencies.json').read_text())['files']}
    if cachefile.exists():
        known.update(json.loads(cachefile.read_text()))
    pending = re.findall(r'^import (SpinCodes\.\S+)', source.read_text(encoding='utf-8'), re.M)
    seen, files = set(), [source]
    while pending:
        module = pending.pop()
        if module in seen:
            continue
        seen.add(module)
        src = ROOT/(module.replace('.', '/')+'.lean')
        obj = ROOT/'.lake/build/lib/lean'/src.relative_to(ROOT).with_suffix('.olean')
        assert obj.exists(), f'Compile dependency first: {src}'
        pending.extend(re.findall(r'^import (SpinCodes\.\S+)', src.read_text(encoding='utf-8'), re.M))
        files += [src] + sorted(obj.parent.glob(obj.name+'*'))
    updates = {}
    archive = DATA/f'peach-module-{source.stem}.tar.gz'
    remote_archive = f'/tmp/spin-lean-peach/{archive.name}'
    with tarfile.open(archive, 'w:gz') as tar:
        for path in files:
            rel, digest = path.relative_to(ROOT).as_posix(), sha(path)
            if known.get(rel) != digest or path == source:
                tar.add(path, arcname=rel)
                updates[rel] = digest
    with archive.open('rb') as stream:
        subprocess.run(SSH+[f'cat > {shlex.quote(remote_archive)}'], stdin=stream, check=True, timeout=300)
    output = '.lake/build/lib/lean/'+relative.removesuffix('.lean')+'.olean'
    command = (f'cd {REMOTE} && tar -xzf {shlex.quote(remote_archive)} && '
               f'mkdir -p {shlex.quote(str(Path(output).parent).replace(chr(92), "/"))} && '
               'export PATH=/tmp/spin-lean-peach/toolchain/lean-4.34.0-linux/bin:$PATH && '
               + admitted_lean(relative, output))
    result = subprocess.run(SSH+[command], stdout=subprocess.PIPE, stderr=subprocess.STDOUT, timeout=900)
    text = result.stdout.decode('utf-8', errors='replace')
    (DATA/f'peach_{source.stem}.log').write_text(text, encoding='utf-8')
    if cachefile.exists():
        known.update(json.loads(cachefile.read_text()))
    known.update(updates)
    cachetmp = cachefile.with_suffix(f'.{source.stem}.tmp')
    cachetmp.write_text(json.dumps(known, indent=2)+'\n')
    # Windows readers can briefly hold the destination against replacement.
    # The cache is only a transfer optimization; concurrent records may cause
    # redundant uploads, but a transient sharing violation must not discard a
    # successful kernel check before its objects are fetched.
    for attempt in range(40):
        try:
            cachetmp.replace(cachefile)
            break
        except PermissionError:
            if attempt == 39:
                raise
            time.sleep(0.05)
    print(text, end='')
    if result.returncode:
        raise SystemExit(result.returncode)
    assert sha(source) == updates[relative], 'Source changed during compilation.'
    archive = DATA/f'peach-module-result-{source.stem}.tar.gz'
    with archive.open('wb') as stream:
        subprocess.run(SSH+[f'cd {REMOTE} && tar -czf - {shlex.quote(output)}*'], stdout=stream, check=True)
    with tarfile.open(archive) as tar:
        assert all(m.name == output or m.name.startswith(output+'.') for m in tar.getmembers())
        tar.extractall(ROOT, filter='data')
    receipt={'status':'PASS','name':source.stem,'source_sha256':sha(source),'olean_sha256':sha(ROOT/output),'admission':'two semantic slots; five matrix workers'}
    receipt_path=DATA/f'occupation_receipt_{source.stem}.json'
    receipt_temp=receipt_path.with_suffix('.tmp');receipt_temp.write_text(json.dumps(receipt)+'\n');receipt_temp.replace(receipt_path)
    print(f'PASS {relative}; fetched kernel-checked objects from Peach.')


if __name__ == '__main__':
    main()

