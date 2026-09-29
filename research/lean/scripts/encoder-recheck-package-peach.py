"""One-at-a-time isolated Peach source checks, preserving Lake package metadata."""
from datetime import datetime, timezone
import argparse
import importlib.util
import json
from pathlib import Path
import re
import shlex
import subprocess
import sys
import tarfile
import time

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / 'scripts/map_data'


def load(name, file):
    spec = importlib.util.spec_from_file_location(name, ROOT / 'scripts' / file)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def main():
    sys.stdout.reconfigure(encoding='utf-8')
    helper = load('peach', 'peach-lean.py')
    foundations = load('foundations', 'encoder-recheck-foundations-isolated.py')
    replay = load('replay', 'replay-native-closure.py')
    parser = argparse.ArgumentParser()
    parser.add_argument('--original-filenames', action='store_true')
    parser.add_argument('--modules', nargs='+')
    parser.add_argument('--plan')
    parser.add_argument('--optional', action='store_true')
    parser.add_argument('--report-name')
    parser.add_argument('--no-package', action='store_true')
    parser.add_argument('--relative-filenames', action='store_true')
    args = parser.parse_args()
    if args.plan:
        plan = json.loads((ROOT / args.plan).read_text(encoding='utf-8'))
        foundations.MODULES = list(plan)
        foundations.THEOREMS = list(plan.values())
    prefix = 'encoder_package_peach_filename_' if args.original_filenames else 'encoder_package_peach_'
    run_name = prefix + time.strftime('%Y%m%d_%H%M%S')
    run_dir = DATA / run_name
    run_dir.mkdir()
    remote = '/tmp/spin-lean-peach/' + run_name
    selected = set(args.modules or foundations.MODULES)
    audit_theorems = [t for m, t in zip(foundations.MODULES, foundations.THEOREMS) if m in selected]
    assert len(audit_theorems) == len(selected)
    imported, order = set(), []
    for module in foundations.MODULES:
        if module not in selected:
            continue
        _, chain = replay.graph(module)
        imported.update(chain)
        order += [m for m in chain if m in selected and m not in order]
    if not args.original_filenames:
        assert order[0] == 'SpinCodes.Framework'
    expected = {'modules': {}, 'dependencies': {}, 'original_filenames': args.original_filenames,
                'no_package': args.no_package}
    for module in order:
        src = ROOT / (module.replace('.', '/') + '.lean')
        obj = ROOT / '.lake/build/lib/lean' / (module.replace('.', '/') + '.olean')
        expected['modules'][module] = {'source_sha256': foundations.sha(src), 'original_object_sha256': foundations.sha(obj),
                                       'original_source_filename': src.relative_to(ROOT).as_posix() if args.relative_filenames else str(src)}
    for module in imported - selected:
        expected['dependencies'][module] = foundations.sha(ROOT / '.lake/build/lib/lean' / (module.replace('.', '/') + '.olean'))
    (run_dir / 'expected.json').write_text(json.dumps(expected, indent=2))
    pin = run_dir / 'EncoderPackageEvidencePin.lean'
    pin.write_text('\n'.join('import ' + m for m in order) + '\n' +
                   '\n'.join('#print axioms ' + t for t in audit_theorems) + '\n', encoding='utf-8')
    archive = run_dir / 'input.tar.gz'
    with tarfile.open(archive, 'w:gz') as tar:
        tar.add(run_dir / 'expected.json', arcname='expected.json')
        tar.add(ROOT / 'scripts/encoder-peach-package-worker.py', arcname='worker.py')
        tar.add(pin, arcname='sources/EncoderPackageEvidencePin.lean')
        for module in order:
            rel = module.replace('.', '/') + '.lean'
            tar.add(ROOT / rel, arcname='sources/' + rel)
    report_path = DATA / ('encoder_semantic_package_peach_verification.json' if args.plan else
                         'encoder_foundations_filename_peach_verification.json' if args.original_filenames else
                         'encoder_foundations_package_peach_verification.json')
    if args.report_name:
        assert re.fullmatch(r'encoder_[A-Za-z_0-9]+\.json', args.report_name)
        report_path = DATA / args.report_name
    report = {'status': 'RUNNING', 'started_utc': datetime.now(timezone.utc).isoformat(),
              'scope': 'Isolated current-source rechecks with package metadata; live objects unchanged',
              'run_directory': str(run_dir), 'remote_directory': remote, 'compiler_jobs': 1,
              'compiler_cap_mb': 8000, 'minimum_remote_free_gib': 20, 'modules': []}

    def save():
        text = json.dumps(report, indent=2) + '\n'
        report_path.write_text(text)
        (run_dir / 'report.json').write_text(text)

    def ssh(command):
        result = subprocess.run(helper.SSH + [command], capture_output=True, timeout=1800)
        assert result.returncode == 0, result.stdout.decode(errors='replace') + result.stderr.decode(errors='replace')
        return result.stdout

    save()
    try:
        with archive.open('rb') as stream:
            subprocess.run(helper.SSH + [f'mkdir -p {shlex.quote(remote)} && tar -xzf - -C {shlex.quote(remote)}'], stdin=stream, check=True)
        report['preparation'] = json.loads(ssh(f'python3 {shlex.quote(remote + "/worker.py")} --prepare'))
        for module in order + ['EncoderPackageEvidencePin']:
            if args.optional:
                base = json.loads((DATA / 'native_theorem_verification.json').read_text())
                if base['status'] in {'RUNNING', 'PASS'}:
                    report.update(status='PARTIAL', stopped_for_final_assembly=True,
                                  note='Optional checks stopped scheduling new modules when base assembly began.')
                    print('PARTIAL: final assembly started; no further optional compilers scheduled.', flush=True)
                    return
            record = json.loads(ssh(f'python3 {shlex.quote(remote + "/worker.py")} {shlex.quote(module)}'))
            log = run_dir / (module.rsplit('.', 1)[-1] + '.log')
            log.write_bytes(ssh('cat ' + shlex.quote(record['remote_log'])))
            record['log'] = log.relative_to(ROOT).as_posix()
            if module in selected:
                report['modules'].append(record)
            else:
                record['audit_source_sha256'] = record.pop('source_sha256')
                record['audit_object_sha256'] = record.pop('fresh_object_sha256', None)
                report['audit_pin'] = record
            save()
            assert record['exit_code'] == 0, log.read_text(encoding='utf-8')[-4000:]
            output = run_dir / (module.rsplit('.', 1)[-1] + '.olean')
            output.write_bytes(ssh('cat ' + shlex.quote(record['remote_output'])))
            assert foundations.sha(output) == record.get('fresh_object_sha256', record.get('audit_object_sha256'))
            record['fresh_object'] = output.relative_to(ROOT).as_posix()
            if module in selected:
                assert foundations.sha(ROOT / (module.replace('.', '/') + '.lean')) == record['source_sha256']
                assert foundations.sha(ROOT / '.lake/build/lib/lean' / (module.replace('.', '/') + '.olean')) == record['original_object_sha256']
                print(f'PASS {module}; matches live bytes={record["fresh_object_bytes_match_original"]}; min remote free={record["minimum_free_gib"]:.2f} GiB', flush=True)
                if module == 'SpinCodes.Framework':
                    assert record['fresh_object_bytes_match_original'], 'Framework gate: no exact byte match; stop before remaining 11'
            else:
                audits = re.findall(r"'([^']+)' depends on axioms: \[([^]]*)\]", log.read_text(encoding='utf-8'))
                audits += [(name, '') for name in re.findall(r"'([^']+)' does not depend on any axioms", log.read_text(encoding='utf-8'))]
                assert len(audits) == len(audit_theorems)
                assert all({a.strip() for a in names.split(',') if a.strip()} <= foundations.ALLOWED for _, names in audits)
                report.update(axiom_audit_count=len(audits), axiom_theorems=[name for name, _ in audits])
            save()
        for module, before in expected['modules'].items():
            assert foundations.sha(ROOT / (module.replace('.', '/') + '.lean')) == before['source_sha256']
            assert foundations.sha(ROOT / '.lake/build/lib/lean' / (module.replace('.', '/') + '.olean')) == before['original_object_sha256']
        report.update(status='PASS', finished_utc=datetime.now(timezone.utc).isoformat(), all_live_objects_unchanged=True,
                      all_fresh_objects_match_live=all(m.get('fresh_object_bytes_match_original', True) for m in report['modules']))
        print('PASS package-aware isolated source checks and axiom audits.', flush=True)
    except BaseException as error:
        report.update(status='FAIL', error=str(error))
        raise
    finally:
        save()


if __name__ == '__main__':
    main()
