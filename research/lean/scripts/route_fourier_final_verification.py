"""Verify completed Fourier proof records; never rebuild dependencies or boxes."""
from pathlib import Path
from datetime import datetime, timezone
import argparse
import hashlib
import json
import re
import subprocess
import sys
import time

root = Path(__file__).resolve().parents[1]
data = root / 'scripts/map_data'
parser = argparse.ArgumentParser()
parser.add_argument('--wait', action='store_true')
args = parser.parse_args()
names = ['dense_fourier_box_replay.json', 'dense_fourier_box_second_replay.json',
         'dense_fourier_box_prefix_tail_replay.json', 'dense_fourier_box_suffix_tail_replay.json']
last_count = None
while True:
    try:
        states = [json.loads((data / name).read_text()) for name in names]
    except (FileNotFoundError, json.JSONDecodeError):
        if not args.wait:
            raise
        time.sleep(10)
        continue
    count = sum(state['checked'] for state in states)
    if count == 307 and all(state['status'] == 'PASS' for state in states):
        break
    assert count <= 307
    if not args.wait:
        raise RuntimeError(f'Fourier boxes still pending: {count}/307')
    if count != last_count:
        print(f'WAIT {count}/307 indexed Fourier boxes', flush=True)
        last_count = count
    time.sleep(30)

def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()

lanes = json.loads((data / 'dense_fourier_box_lanes.json').read_text())
tail = lanes['suffix_tail']['indices']
ownership = [list(range(115)), sorted(set(range(154, 307)) - set(tail)),
             lanes['prefix_tail']['indices'], tail]
all_boxes = []
for name, state, owned in zip(names, states, ownership):
    assert state['total'] == len(owned), name
    assert sorted(record['index'] for record in state['boxes']) == sorted(owned), name
    assert state['checked'] == len(state['boxes']), name
    assert state['common_prefactor'] == 24000000000000, name
    all_boxes.extend(state['boxes'])
assert sorted(record['index'] for record in all_boxes) == list(range(307))

witness_state = json.loads((data / 'dense_fourier_witness_replay.json').read_text())
assert witness_state['status'] == 'PASS'
assert witness_state['checked'] == witness_state['total'] == 105
assert sorted(record['index'] for record in witness_state['witnesses']) == list(range(105))

allowed = {'propext', 'Classical.choice', 'Quot.sound'}
verified_modules = []
audit_records = []
for family, records, expected_audits in [('W', witness_state['witnesses'], 3), ('B', all_boxes, 4)]:
    for record in sorted(records, key=lambda record: record['index']):
        tag = f"{family}{record['index']:03d}"
        assert record['audits'] == expected_audits
        assert len(record['modules']) == 2
        for module, suffix in zip(record['modules'], ['Data', '']):
            stem = f'DenseFourierExact{tag}{suffix}'
            relative = f'SpinCodes/Structured/{stem}.lean'
            assert module['source'].replace('\\', '/') == relative
            source = root / relative
            obj = root / '.lake/build/lib/lean' / Path(relative).with_suffix('.olean')
            assert sha(source) == module['source_sha256'], source
            assert sha(obj) == module['object_sha256'], obj
            source_text = source.read_text(encoding='utf-8-sig')
            assert not re.search(r'\b(sorry|admit|native_decide|axiom)\b', source_text), source
            compile_log = data / f'dense_fourier_batch_{tag}{suffix}.log'
            assert f'PASS {relative}; fetched kernel-checked objects from Peach.' in compile_log.read_text(encoding='utf-8')
            verified_modules.append({**module, 'source': relative,
                                     'compile_log': compile_log.name, 'compile_log_sha256': sha(compile_log)})
        audit_log = data / f'peach_DenseFourierExact{tag}.log'
        audits = re.findall(r"'([^']+)' depends on axioms: \[([^]]*)\]", audit_log.read_text(encoding='utf-8'))
        assert len(audits) == expected_audits, audit_log
        for theorem, axioms in audits:
            assert theorem.startswith(f'Spin.Structured.DenseFourierExact.{tag}.')
            actual = [name.strip() for name in axioms.split(',') if name.strip()]
            assert set(actual) <= allowed, (theorem, actual)
            audit_records.append({'theorem': theorem, 'axioms': actual})

check = subprocess.run([sys.executable, str(root / 'scripts/route_fourier_index_audit.py')],
                       cwd=root, capture_output=True, text=True, encoding='utf-8', check=True)
structure = json.loads(check.stdout)
assert structure['checked'] == structure['emitted'] == structure['total'] == 307
assert not structure['missing_are_pending']
assert len(verified_modules) == 824 and len(audit_records) == 1543
report = {
    'status': 'PASS', 'verified_at_utc': datetime.now(timezone.utc).isoformat(),
    'scope': 'All 105 exact Fourier matrix witnesses and 307 original-indexed Fourier box certificates. Individual modules were kernel-checked on Peach with previously checked dependencies reused. This validates current source/object/log hashes and all stored axiom audits; it is not a fresh full dependency-closure replay. Mixed-family aggregation and the final native theorem are separate checks.',
    'witnesses': 105, 'boxes': 307, 'modules_count': len(verified_modules),
    'audit_count': len(audit_records), 'allowed_axioms': sorted(allowed),
    'common_prefactor': 24000000000000, 'eta': '4/10000000',
    'structural_audit': structure,
    'input_records': [{'file': name, 'sha256': sha(data / name)} for name in names + ['dense_fourier_witness_replay.json']],
    'modules': verified_modules, 'audits': audit_records,
}
destination = data / 'dense_fourier_complete_verification.json'
destination.write_text(json.dumps(report, indent=2) + '\n')
print(f'PASS {destination.name}: 105 witnesses, 307 boxes, 824 current source/object pairs, 1543 standard-only audits.', flush=True)
