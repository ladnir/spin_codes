"""Check paper-milestone provenance without running searches or benchmarks.

This authenticates retained audit evidence; it is not an interval replay.
Missing bulk evidence is reported and can optionally be made an error.
"""
import argparse
from fractions import Fraction
import hashlib
import json
from pathlib import Path
import statistics
import subprocess

ROOT = Path(__file__).resolve().parents[2]
HERE = Path(__file__).resolve().parent
AUDIT = HERE / 'generated/t128_s19_ladder_completion_audit_v1.json'
AUDIT_SHA = 'dc06c1e9770a33ca37c8ea8824ae554d3b2ce84011477e55b210d39a38446aa9'


def require(condition, message):
    if not condition:
        raise ValueError(message)


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def read(path):
    return json.loads(path.read_text(encoding='utf-8'))


def rational(value):
    return Fraction(int(value['numerator']), int(value['denominator']))


def run(require_local=False, check_index=False):
    require(sha(AUDIT) == AUDIT_SHA, 'Completion audit changed')
    audit = read(AUDIT)
    require(audit['status'] == 'FULL_M20_M22_M24_LADDER_COMPLETION_AUDIT_PASSED',
            'Completion audit did not pass')
    require([r['message_exponent'] for r in audit['instances']] == [20, 22, 24],
            'Incomplete ladder')
    require(sum(r['dense_rectangles'] for r in audit['instances']) == 1497,
            'Unexpected dense coverage count')
    pins = dict(audit['source_sha256'])
    pins[AUDIT.relative_to(ROOT).as_posix()] = AUDIT_SHA
    for row in audit['instances']:
        m = row['message_exponent']
        path = HERE / f'generated/t128_s19_m{m}_ladder_full_v1.json'
        require(sha(path) == row['ledger_sha256'], f'Changed K=2^{m} ledger')
        ledger = read(path)
        total = rational(ledger['union_upper'])
        require(total == rational(row['union_upper']) and 0 < total < Fraction(1, 2**40),
                f'Wrong final union at K=2^{m}')
        require(total == sum((rational(ledger[k]) for k in
                             ('q1_upper', 'sparse_rest_upper', 'dense_upper')), Fraction()),
                f'Component sum mismatch at K=2^{m}')
        require(ledger['covered_intervals'] == [[1, 2**(m-7)]] and
                ledger['complete_coverage'] and row['all_occupancies_covered'],
                f'Incomplete recorded coverage at K=2^{m}')

    perf = read(ROOT / 'workstreams/bare_bch_rm2sub/PERFORMANCE.json')
    for name, digest in perf['source_sha256'].items():
        require(name not in pins or pins[name] == digest, 'Conflicting source pin: ' + name)
        pins[name] = digest
    raw = perf['raw_final_samples']
    require(len(raw) == 36, 'Missing performance samples')
    for row in perf['summary']:
        values = [r['median_ms'] for r in raw if
                  (r['configuration'], r['m']) == (row['configuration'], row['m'])]
        require(len(values) == 3 and statistics.median(values) == row['median_ms'],
                'Performance summary mismatch')

    tracked = set(subprocess.check_output(['git', 'ls-files', '-z'], cwd=ROOT)
                  .decode().split('\0'))
    missing, external, index_errors = [], [], []
    for name, expected in sorted(pins.items()):
        path = (ROOT / name).resolve()
        require(path.is_relative_to(ROOT), 'Dependency escapes repository: ' + name)
        if path.is_file():
            require(sha(path) == expected, 'Pinned dependency changed: ' + name)
        else:
            missing.append(name)
        if name not in tracked:
            external.append(name)
        elif check_index:
            # Existing calibration files use explicit eol=crlf attributes.
            # Authenticate the bytes a checkout produces, not normalized blobs.
            blob = subprocess.check_output(['git', 'cat-file', '--filters', ':' + name], cwd=ROOT)
            if hashlib.sha256(blob).hexdigest() != expected:
                index_errors.append(name)
    require(not index_errors, 'Staged checkout bytes differ from pins: ' + ', '.join(index_errors))
    require(not require_local or not missing, 'Missing local evidence: ' + ', '.join(missing))
    print(json.dumps(dict(status='PAPER_MILESTONE_PROVENANCE_PASSED',
                          numerical_replay_performed=False,
                          authenticated_pins=len(pins)-len(missing),
                          missing_local_dependencies=len(missing),
                          dependencies_not_in_git=len(external),
                          staged_checkout_bytes_checked=check_index), indent=2))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--require-local-evidence', action='store_true')
    parser.add_argument('--check-index', action='store_true')
    args = parser.parse_args()
    run(args.require_local_evidence, args.check_index)
