"""Remote half of the isolated package-metadata evidence check."""
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys
import time

PROJECT = Path('/tmp/spin-lean-peach/project')
RUN = Path(__file__).resolve().parent
EXPECTED = json.loads((RUN / 'expected.json').read_text())


def sha(path):
    with path.open('rb') as stream:
        return hashlib.file_digest(stream, 'sha256').hexdigest()


def free_gib():
    rows = dict(line.split(':', 1) for line in Path('/proc/meminfo').read_text().splitlines())
    return int(rows['MemAvailable'].split()[0]) / 1024**2


def main():
    module = sys.argv[1]
    objects = RUN / 'objects'
    if module == '--prepare':
        for name, expected in EXPECTED['dependencies'].items():
            original = PROJECT / '.lake/build/lib/lean' / (name.replace('.', '/') + '.olean')
            assert sha(original) == expected, f'Imported object hash mismatch: {name}'
            for artifact in original.parent.glob(original.name + '*'):
                dest = objects / artifact.relative_to(PROJECT / '.lake/build/lib/lean')
                dest.parent.mkdir(parents=True, exist_ok=True)
                os.link(artifact, dest)
        print(json.dumps({'status': 'PREPARED', 'free_gib': free_gib()}), flush=True)
        return
    source = RUN / 'sources' / (module.replace('.', '/') + '.lean')
    audit = module == 'EncoderPackageEvidencePin'
    entry = EXPECTED['modules'][module] if not audit else {}
    if not audit:
        assert sha(source) == entry['source_sha256']
    environment = os.environ.copy()
    environment['PATH'] = '/tmp/spin-lean-peach/toolchain/lean-4.34.0-linux/bin:' + environment['PATH']
    result = subprocess.run(['lake', 'env', 'python3', '-c',
        "import os,json; print(json.dumps({k:os.environ.get(k,'') for k in ['LEAN_PATH','PATH']}))"],
        cwd=PROJECT, env=environment, capture_output=True, text=True, check=True)
    environment.update(json.loads(result.stdout))
    environment['LEAN_PATH'] = str(objects) + ':' + environment['LEAN_PATH']
    output = objects / (module.replace('.', '/') + '.olean')
    output.parent.mkdir(parents=True, exist_ok=True)
    log = RUN / (module.rsplit('.', 1)[-1] + '.log')
    setup = RUN / (module.rsplit('.', 1)[-1] + '.setup.json')
    setup_data = {'name': module, 'isModule': False,
                  'importArts': {}, 'plugins': [], 'dynlibs': [], 'options': {}}
    if not EXPECTED.get('no_package', False):
        setup_data['package'] = 'spincodes'
    setup.write_text(json.dumps(setup_data))
    command = ['lean', '-j1', '-M', '8000', '--root', str(RUN / 'sources'),
               '--setup', str(setup), '-o', str(output), str(source)]
    use_stdin = EXPECTED.get('original_filenames', False) and not audit
    if use_stdin:
        command = ['lean', '-j1', '-M', '8000', '--setup', str(setup), '-o', str(output),
                   '--stdin', entry['original_source_filename']]
    minimum = free_gib()
    assert minimum >= 20, 'Remote headroom below 20 GiB before start'
    start = time.monotonic()
    with log.open('wb') as stream, source.open('rb') as source_stream:
        proc = subprocess.Popen(command, cwd=PROJECT, env=environment, stdout=stream, stderr=subprocess.STDOUT,
                                stdin=source_stream if use_stdin else None)
        while proc.poll() is None:
            minimum = min(minimum, free_gib())
            if minimum < 20:
                proc.terminate()
                proc.wait()
                raise RuntimeError('Remote headroom below 20 GiB; own diagnostic stopped')
            time.sleep(1)
    record = {'module': module, 'exit_code': proc.returncode, 'command': command,
              'minimum_free_gib': minimum, 'seconds': round(time.monotonic() - start, 3),
              'setup_sha256': sha(setup), 'log_sha256': sha(log), 'remote_log': str(log),
              'source_sha256': sha(source), 'remote_output': str(output)}
    record['lean_version'] = subprocess.check_output(['lean', '--version'], cwd=PROJECT, env=environment, text=True).strip()
    assert '293d5d0c0c3f3dded4688b3ccd6a33939ac5102b' in record['lean_version']
    if proc.returncode == 0:
        record.update(status='PASS', fresh_object_sha256=sha(output))
        if not audit:
            assert sha(source) == entry['source_sha256']
            record['original_object_sha256'] = entry['original_object_sha256']
            record['fresh_object_bytes_match_original'] = sha(output) == entry['original_object_sha256']
            if record['fresh_object_bytes_match_original']:
                record['object_sha256'] = sha(output)
    else:
        record['status'] = 'FAIL'
    print(json.dumps(record), flush=True)


if __name__ == '__main__':
    main()
