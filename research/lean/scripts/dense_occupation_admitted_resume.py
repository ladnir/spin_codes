"""Hash-gated occupation replay: five matrices, three local, three indexed workers.

Each mode has one manifest writer. All semantic compilers share the dedicated
helper's two remote permits. Start only after predecessor controllers drain.
"""
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor
from datetime import datetime, timezone
import hashlib, json, re, subprocess, sys, threading, time

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / 'scripts/map_data'
MODE = sys.argv[1]
assert MODE in ('matrices', 'boxes', 'rates')
WORKERS = 5 if MODE == 'matrices' else 3
FILE = {'matrices': 'dense_occupation_fixed_remaining_verification.json',
        'boxes': 'dense_occupation_fixed_all_boxes_verification.json',
        'rates': 'dense_occupation_fixed_rates_verification.json'}[MODE]
KEY = 'witnesses' if MODE == 'matrices' else 'boxes'
COUNT = {'matrices': 'numerical_witnesses_checked',
         'boxes': 'numerical_local_boxes_checked',
         'rates': 'indexed_rate_boxes_checked'}[MODE]
OUT = DATA / FILE
LOCK = threading.Lock()

def read(path):
    for _ in range(100):
        try:
            return json.loads(path.read_text())
        except json.JSONDecodeError:
            time.sleep(.05)
    raise RuntimeError(f'Unstable JSON: {path}')

def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()

def valid_module(rec):
    name = rec.get('name')
    if not name:
        return False
    source = ROOT / f'SpinCodes/Structured/{name}.lean'
    obj = ROOT / f'.lake/build/lib/lean/SpinCodes/Structured/{name}.olean'
    return (source.exists() and obj.exists() and
            sha(source) == rec.get('source_sha256') and
            sha(obj) == rec.get('olean_sha256'))

def records(rec, tag, mode):
    return ([dict(rec, name='DenseOccupationFixed' + tag + 'Rate')]
            if mode == 'rates' else rec.get('modules', []))

def valid_pass(rec, tag, mode):
    modules = records(rec, tag, mode)
    return (rec.get('status') == 'PASS' and
            len(modules) == (1 if mode == 'rates' else 2) and
            all(valid_module(m) for m in modules))

REPORT = read(OUT)
for tag, rec in REPORT[KEY].items():
    if rec.get('status') == 'PASS':
        assert valid_pass(rec, tag, MODE), f'Invalid existing PASS: {tag}'

if MODE == 'matrices':
    ITEMS = [dict(box=f'W{i:03d}') for i in range(9, 133)]
else:
    candidate = 'dense_occupation_fixed_' + ('all_boxes' if MODE == 'boxes' else 'rates') + '_candidates.json'
    ITEMS = read(DATA / candidate)['boxes' if MODE == 'boxes' else 'entries']

def save():
    REPORT[COUNT] = sum(rec.get('status') == 'PASS' for rec in REPORT[KEY].values())
    tmp = OUT.with_suffix('.admitted.tmp')
    tmp.write_text(json.dumps(REPORT, indent=2) + '\n')
    for attempt in range(100):
        try:
            tmp.replace(OUT)
            return
        except PermissionError:
            if attempt == 99:
                raise
            time.sleep(.02)

def wait_parent(item):
    if MODE == 'matrices':
        return
    tag = item['witness'] if MODE == 'boxes' else item['box']
    while True:
        if MODE == 'rates':
            parent = read(DATA / 'dense_occupation_fixed_all_boxes_verification.json')
            rec = parent['boxes'].get(tag, {})
            ok = valid_pass(rec, tag, 'boxes')
        elif tag == 'First':
            parent = read(DATA / 'dense_occupation_fixed_verification.json')
            mods = [m for m in parent['modules'] if m['name'] in
                    ('DenseOccupationFixedFirst', 'DenseOccupationFixedFirstData')]
            rec = parent
            ok = parent.get('status') == 'PASS' and len(mods) == 2 and all(valid_module(m) for m in mods)
        else:
            file = ('dense_occupation_fixed_batch_verification.json' if int(tag[1:]) <= 8
                    else 'dense_occupation_fixed_remaining_verification.json')
            parent = read(DATA / file)
            rec = parent['witnesses'].get(tag, {})
            ok = valid_pass(rec, tag, 'matrices')
        if ok:
            return
        if rec.get('status') == 'PASS':
            raise RuntimeError(f'Parent PASS hash mismatch: {tag}')
        if rec.get('status') in ('FAIL', 'FAILED') or parent.get('status') in ('FAIL', 'FAILED'):
            raise RuntimeError(f'Parent replay failed: {tag}')
        time.sleep(2)

def audit(text, expected):
    found = re.findall(r"'([^']+)' depends on axioms: \[([^\]]*)\]", text, re.S)
    if len(found) != expected or any(set(x.strip() for x in ax.split(',')) -
            {'propext', 'Classical.choice', 'Quot.sound'} for _, ax in found):
        raise RuntimeError('Unexpected axiom audit')
    return dict(found)

def run(item):
    tag = item['box']
    start = time.monotonic()
    old = REPORT[KEY].get(tag, {})
    # Only reuse individually recorded, matching hashes. Existence plus an old
    # log is deliberately insufficient; unrecorded drained work is rechecked.
    previous = {m['name']: m for m in records(old, tag, MODE) if valid_module(m)}
    rec = {**old, **item, 'status': 'WAITING', 'modules': []}
    with LOCK:
        REPORT[KEY][tag] = rec
        save()
    try:
        wait_parent(item)
        with LOCK:
            rec.update(status='RUNNING', started=datetime.now(timezone.utc).isoformat())
            save()
        print(tag, MODE, 'RUNNING', flush=True)
        for suffix in (['Rate'] if MODE == 'rates' else ['Data', '']):
            name = 'DenseOccupationFixed' + tag + suffix
            receipt = DATA / f'occupation_receipt_{name}.json'
            cached = previous.get(name)
            if cached is None and receipt.exists():
                candidate = read(receipt)
                if candidate.get('status') == 'PASS' and valid_module(candidate):
                    cached = candidate
            if cached is not None:
                text = (DATA / f'peach_{name}.log').read_text(encoding='utf-8')
            else:
                result = subprocess.run([sys.executable, str(ROOT / 'scripts/peach-lean-occupation-admitted.py'),
                    f'SpinCodes/Structured/{name}.lean'], cwd=ROOT, text=True, encoding='utf-8',
                    errors='replace', stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
                text = result.stdout
                (DATA / f'dense_occupation_admitted_{name}.log').write_text(text, encoding='utf-8')
                if result.returncode:
                    raise RuntimeError(text[-3500:])
                cached = read(receipt)
                assert cached.get('status') == 'PASS' and valid_module(cached)
            with LOCK:
                if MODE == 'rates':
                    rec.update(source_sha256=cached['source_sha256'], olean_sha256=cached['olean_sha256'])
                else:
                    rec['modules'].append(cached)
                if suffix != 'Data':
                    rec['axioms'] = audit(text, 3 if MODE == 'matrices' else 2)
                save()
        with LOCK:
            rec.update(status='PASS', elapsed_seconds=time.monotonic()-start)
            save()
        print(tag, MODE, 'PASS', flush=True)
        return True
    except Exception as exc:
        with LOCK:
            rec.update(status='FAIL', error=str(exc))
            REPORT['status'] = 'FAIL'
            save()
        print(tag, MODE, 'FAIL', str(exc), flush=True)
        return False

REPORT.update(status='RUNNING', concurrency=WORKERS, semantic_permits=2,
              admitted_helper='peach-lean-occupation-admitted.py', resumed=True)
save()
pending = [item for item in ITEMS if REPORT[KEY].get(item['box'], {}).get('status') != 'PASS']
if MODE != 'matrices':
    # Follow matrix completion order, instead of blocking all workers on an
    # early rectangle while later rectangles already have checked witnesses.
    pending.sort(key=lambda item: (0 if item['witness'] == 'First' else int(item['witness'][1:]), item['box']))
with ThreadPoolExecutor(max_workers=WORKERS) as pool:
    results = list(pool.map(run, pending))
REPORT['status'] = 'PASS' if all(results) and REPORT[COUNT] == len(ITEMS) else 'FAIL'
save()
print(REPORT['status'], MODE, REPORT[COUNT], flush=True)
