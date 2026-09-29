"""Read-only consolidation of existing Lean verification evidence; no compilation.

Writes INCOMPLETE until every required final gate and current hash pairing passes.
This is an evidence index for human inspection, not an independent proof checker.
"""
from datetime import datetime, timezone
from pathlib import Path
import hashlib
import json
import re
import sys
import traceback

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / 'scripts/map_data'
ALLOWED = {'propext', 'Classical.choice', 'Quot.sound'}
PIN_HASH = '447d959bec1591dd678a9696262ceb7a3df818b59cdda38086c2feb259f7ceca'
gates, inputs, pairs, issues, pair_kinds, extra_artifacts = {}, {}, {}, [], {}, {}
DESTINATION = DATA / 'final_closure_evidence.json'


def publish(record):
    temporary = DESTINATION.with_suffix('.route.tmp')
    temporary.write_text(json.dumps(record, indent=2) + '\n', encoding='utf-8')
    temporary.replace(DESTINATION)


def unhandled_error(kind, error, tb):
    # Invalidate an older PASS even if a malformed input or I/O error interrupts
    # this run before its normal diagnostic report can be assembled.
    try:
        publish({'status': 'INCOMPLETE', 'phase': 'EXCEPTION',
                 'checked_at_utc': datetime.now(timezone.utc).isoformat(),
                 'error_type': kind.__name__, 'error': str(error),
                 'diagnostic': ''.join(traceback.format_exception(kind, error, tb)),
                 'scope': 'Read-only evidence consolidation failed; no final-closure claim.'})
    finally:
        sys.__excepthook__(kind, error, tb)


sys.excepthook = unhandled_error
publish({'status': 'INCOMPLETE', 'phase': 'SCANNING',
         'started_at_utc': datetime.now(timezone.utc).isoformat(),
         'scope': 'Read-only scan in progress; any previous PASS is invalidated.'})


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def issue(gate, message):
    gates[gate]['problems'].append(message)
    issues.append({'gate': gate, 'message': message})


def need(gate, condition, message):
    if not condition:
        issue(gate, message)
    return condition


def report(gate, filename):
    gates.setdefault(gate, {'problems': [], 'reports': []})
    gates[gate]['reports'].append(filename)
    path = DATA / filename
    try:
        raw = path.read_bytes()
        record = json.loads(raw)
    except (OSError, json.JSONDecodeError) as exc:
        issue(gate, f'Unavailable/incomplete JSON {filename}: {exc}')
        return {}
    inputs[filename] = hashlib.sha256(raw).hexdigest()
    need(gate, record.get('status') == 'PASS', f'{filename} status={record.get("status")}')
    return record


def pair(gate, relative, record):
    relative = relative.replace('\\', '/')
    if not relative.endswith('.lean'):
        relative = f'SpinCodes/Structured/{relative}.lean'
    try:
        src = ROOT / relative
        obj = (ROOT / '.lake/build/lib/lean' / relative).with_suffix('.olean')
        expected = (record['source_sha256'], record.get('object_sha256', record.get('olean_sha256')))
        actual = (digest(src), digest(obj))
        need(gate, actual == expected, f'Current source/object hash mismatch: {relative}')
        if relative in pairs:
            need(gate, pairs[relative] == expected, f'Conflicting recorded pair: {relative}')
        pairs[relative] = expected
        kind = 'individual_check_record_pair'
        if gate == 'actual_native_theorem_and_full_statement_pin' and Path(relative).stem not in {
                'ConcreteNativeTheorem', 'ConcreteNativeTheoremPin'}:
            kind = 'final_assembly_dependency_snapshot'
        pair_kinds.setdefault(relative, set()).add(kind)
    except (OSError, KeyError, TypeError) as exc:
        issue(gate, f'Missing source/object evidence for {relative}: {exc}')


def module_list(gate, records):
    if isinstance(records, dict):
        for relative, entry in records.items():
            pair(gate, relative, entry)
    else:
        for entry in records:
            name = entry.get('source') or entry.get('name') or entry.get('module')
            if name:
                pair(gate, name, entry)
            else:
                issue(gate, 'Module record has no source/name/module')


def exact_modules(gate, records, expected):
    actual = [entry.get('source') or entry.get('name') or entry.get('module') for entry in records]
    actual = [name.replace('\\', '/') if isinstance(name, str) else name for name in actual]
    need(gate, len(actual) == len(expected) and set(actual) == set(expected),
         'Exact module identity/count mismatch')


def audit_map(gate, records, count=None):
    if count is not None:
        need(gate, len(records) == count, f'Expected {count} axiom audits, found {len(records)}')
    for theorem, names in records.items():
        actual = set(names if isinstance(names, list) else [s.strip() for s in names.split(',') if s.strip()])
        need(gate, actual <= ALLOWED, f'Nonstandard axiom in {theorem}: {sorted(actual - ALLOWED)}')


def audit_log(gate, filename, count, required=()):
    try:
        raw = (DATA / filename).read_text(encoding='utf-8')
        entries = re.findall(r"'([^']+)' depends on axioms: \[([^]]*)\]", raw)
        need(gate, len(entries) == count, f'{filename}: expected {count} audits, found {len(entries)}')
        audit_map(gate, dict(entries))
        need(gate, set(required) <= {name for name, _ in entries}, f'{filename}: missing required theorem audits')
        need(gate, 'sorryAx' not in raw and 'Lean.ofReduceBool' not in raw,
             f'{filename}: disallowed trust marker')
    except OSError as exc:
        issue(gate, f'Missing audit log {filename}: {exc}')


def indexed(gate, filename, count, family, suffix='', local=False):
    data = report(gate, filename)
    boxes = data.get('boxes', {})
    expected = {f'B{i:03d}' for i in range(count)}
    need(gate, set(boxes) == expected, f'Expected exact B000..B{count-1:03d} box set')
    passed = 0
    for tag, entry in boxes.items():
        if not need(gate, entry.get('status') == 'PASS', f'{tag}: status={entry.get("status")}'):
            continue
        passed += 1
        stem = f'{family}{tag}{suffix}'
        if local:
            module_list(gate, entry.get('modules', []))
            exact_modules(gate, entry.get('modules', []), [stem + 'Data', stem])
        else:
            pair(gate, stem, entry)
        count_audits = 1 if family == 'DenseScalarGeometry' else 2
        audit_log(gate, f'peach_{stem}.log', count_audits)
        if 'axioms' in entry:
            audit_map(gate, entry['axioms'], count_audits)
    gates[gate]['checked_count'] = passed
    gates[gate]['expected_count'] = count
    return data


# Original statement and paper preservation are independent of all active workers.
gate = 'original_pin_and_papers'
baseline = report(gate, 'parallel_finite_verification.json')
need(gate, digest(ROOT / 'SpinCodes/Pin.lean') == PIN_HASH, 'Original SpinCodes/Pin.lean changed')
expected_papers = {path.replace('\\', '/'): value for path, value in baseline.get('paper_sha256', {}).items()}
current_papers = {path.relative_to(ROOT.parent).as_posix(): digest(path)
                  for path in sorted((ROOT.parent / 'paper').glob('*.tex'))}
need(gate, bool(expected_papers) and current_papers == expected_papers, 'Original paper file set/hashes changed')

# Three independently checked numerical families, then their aggregate objects.
indexed('scalar_local_boxes', 'dense_scalar_exact_boxes_verification.json', 283, 'DenseScalarExact', local=True)
indexed('scalar_indexed_boxes', 'dense_scalar_geometry_verification.json', 283, 'DenseScalarGeometry')

gate = 'fourier_complete'
fourier = report(gate, 'dense_fourier_complete_verification.json')
for field, expected in [('witnesses', 105), ('boxes', 307), ('modules_count', 824), ('audit_count', 1543)]:
    need(gate, fourier.get(field) == expected, f'Expected {field}={expected}')
module_list(gate, fourier.get('modules', []))
exact_modules(gate, fourier.get('modules', []),
              [f'SpinCodes/Structured/DenseFourierExact{family}{i:03d}{suffix}.lean'
               for family, count in [('W', 105), ('B', 307)]
               for i in range(count) for suffix in ['Data', '']])
audit_map(gate, {entry['theorem']: entry['axioms'] for entry in fourier.get('audits', [])}, 1543)
expected_fourier_audits = {
    f'Spin.Structured.DenseFourierExact.{family}{i:03d}.{theorem}'
    for family, count, theorems in [
        ('W', 105, ['collatz', 'witness_floor', 'witness_prefactor']),
        ('B', 307, ['exponent_bound', 'collatz', 'certified', 'certified_uniform'])]
    for i in range(count) for theorem in theorems}
need(gate, len(fourier.get('audits', [])) == 1543 and
     {entry['theorem'] for entry in fourier.get('audits', [])} == expected_fourier_audits,
     'Exact Fourier theorem-audit identity/count mismatch')
for entry in fourier.get('input_records', []):
    try:
        need(gate, digest(DATA / entry['file']) == entry['sha256'], f'Fourier input record changed: {entry["file"]}')
    except OSError as exc:
        issue(gate, str(exc))
for entry in fourier.get('modules', []):
    try:
        need(gate, digest(DATA / entry['compile_log']) == entry['compile_log_sha256'],
             f'Fourier compile log changed: {entry["compile_log"]}')
    except (KeyError, OSError) as exc:
        issue(gate, f'Fourier compile log evidence missing: {exc}')
structure = fourier.get('structural_audit', {})
need(gate, structure.get('checked') == structure.get('emitted') == 307 and
     structure.get('duplicate_indices') == 0 and not structure.get('missing_are_pending'),
     'Fourier complete index/ownership audit missing')

gate = 'occupation_witnesses'
first = report(gate, 'dense_occupation_fixed_verification.json')
module_list(gate, first.get('modules', []))
exact_modules(gate, first.get('modules', []), [f'DenseOccupationFixed{suffix}' for suffix in [
    'Defs', 'Basic', 'Moments', 'ColumnDefs', 'Column', 'Mixture', 'Check', 'Pin', 'FirstData', 'First']])
audit_map(gate, first.get('axioms', {}), 9)
need(gate, first.get('numerical_witnesses_checked') == 1, 'First occupation witness missing')
for filename, expected in [('dense_occupation_fixed_batch_verification.json', range(1, 9)),
                           ('dense_occupation_fixed_remaining_verification.json', range(9, 133))]:
    data = report(gate, filename)
    entries = data.get('witnesses', {})
    need(gate, set(entries) == {f'W{i:03d}' for i in expected}, f'{filename}: witness index set mismatch')
    for tag, entry in entries.items():
        if need(gate, entry.get('status') == 'PASS', f'{tag}: status={entry.get("status")}'):
            module_list(gate, entry.get('modules', []))
            exact_modules(gate, entry.get('modules', []),
                          [f'DenseOccupationFixed{tag}Data', f'DenseOccupationFixed{tag}'])
            audit_map(gate, entry.get('axioms', {}), 3)

indexed('occupation_local_boxes', 'dense_occupation_fixed_all_boxes_verification.json',
        433, 'DenseOccupationFixed', local=True)
occupation = indexed('occupation_indexed_boxes', 'dense_occupation_fixed_rates_verification.json',
                     433, 'DenseOccupationFixed', suffix='Rate')

gate = 'global_geometry_and_indices'
geometry = report(gate, 'dense_occupation_geometry_verification.json')
need(gate, geometry.get('geometry_boxes_covered') == 1023, 'Expected 1023 covered rectangles')
module_list(gate, geometry.get('modules', []))
audit_map(gate, geometry.get('axioms', {}))
global_indices = report(gate, 'dense_occupation_global_indices_verification.json')
pair(gate, 'DenseOccupationGlobalIndices', global_indices)
audit_map(gate, global_indices.get('axioms', {}), 2)
try:
    text = (ROOT / 'SpinCodes/Structured/DenseOccupationGlobalIndices.lean').read_text(encoding='utf-8-sig')
    lists = {family: json.loads(re.search(rf'def {family}Indices : List ℕ := (\[[^\n]*\])', text).group(1))
             for family in ['scalar', 'occupation', 'fourier']}
    need(gate, [len(lists[name]) for name in ['scalar', 'occupation', 'fourier']] == [283, 433, 307],
         'Family index counts differ')
    need(gate, sorted(sum(lists.values(), [])) == list(range(1023)), 'Global families do not partition indices 0..1022')
    for tag, entry in occupation.get('boxes', {}).items():
        need(gate, entry.get('global_index') == lists['occupation'][int(tag[1:])], f'{tag}: occupation global index mismatch')
except (OSError, AttributeError, IndexError, ValueError) as exc:
    issue(gate, f'Global index parsing failed: {exc}')

gate = 'occupation_aggregate'
aggregate = report(gate, 'dense_occupation_aggregate_verification.json')
if aggregate.get('status') == 'PASS':
    need(gate, aggregate.get('indexed_rate_boxes_checked') == 433, 'Expected 433 indexed occupation boxes')
    pair(gate, 'DenseOccupationAllCertified', aggregate)
    audit_map(gate, aggregate.get('axioms', {}), 1)

gate = 'occupation_complete'
complete = report(gate, 'dense_occupation_complete_verification.json')
for field, expected in [('witnesses', 133), ('local_boxes', 433), ('indexed_boxes', 433),
                        ('module_count', 1566), ('audit_count', 2132)]:
    need(gate, complete.get(field) == expected, f'Expected {field}={expected}')
expected_occupation_modules = (
    {'DenseOccupationFixed' + tag + suffix
     for tag in ['First'] + [f'W{i:03d}' for i in range(1, 133)] for suffix in ['Data', '']} |
    {f'DenseOccupationFixedB{i:03d}{suffix}' for i in range(433) for suffix in ['Data', '', 'Rate']} |
    {'DenseOccupationAllCertified'})
need(gate, set(complete.get('modules', {})) == expected_occupation_modules,
     'Exact complete occupation module identity/count mismatch')
module_list(gate, complete.get('modules', {}))
expected_occupation_audits = (
    {f'Spin.Structured.DenseOccupationFixed.{tag}.{theorem}'
     for tag in ['First'] + [f'W{i:03d}' for i in range(1, 133)]
     for theorem in ['collatz', 'iterate', 'routed_probability']} |
    {f'Spin.Structured.DenseOccupationFixed.B{i:03d}.{theorem}' for i in range(433)
     for theorem in ['exponent_bound', 'collatz', 'certified', 'certified_uniform']} |
    {'Spin.Structured.DenseGeometry.occupation_certified'})
need(gate, set(complete.get('axioms', {})) == expected_occupation_audits,
     'Exact complete occupation theorem-audit identity/count mismatch')
audit_map(gate, complete.get('axioms', {}), 2132)
expected_producers = {
    'dense_occupation_fixed_verification.json', 'dense_occupation_fixed_batch_verification.json',
    'dense_occupation_fixed_remaining_verification.json', 'dense_occupation_fixed_all_boxes_verification.json',
    'dense_occupation_fixed_rates_verification.json', 'dense_occupation_aggregate_verification.json'}
producer_hashes = complete.get('producer_manifests', {})
need(gate, set(producer_hashes) == expected_producers, 'Complete occupation producer set mismatch')
for filename, expected in producer_hashes.items():
    need(gate, inputs.get(filename) == expected, f'Occupation producer snapshot mismatch: {filename}')
    need(gate, digest(DATA / filename) == expected, f'Occupation producer changed: {filename}')
if complete.get('status') == 'PASS':
    checker = ROOT / 'scripts/dense_occupation_complete_wait.py'
    need(gate, complete.get('script_sha256') == digest(checker), 'Occupation complete checker hash mismatch')
    need(gate, bool(complete.get('verified_at_utc')), 'Occupation complete verification timestamp missing')
    extra_artifacts[checker] = complete.get('script_sha256')

gate = 'mixed_aggregate'
mixed = report(gate, 'dense_mixed_aggregate_verification.json')
for stem, count in [('DenseOccupationScalarCertified', 283), ('DenseOccupationFourierCertified', 307),
                    ('DenseOccupationMixedCertified', 1023)]:
    entry = mixed.get('groups', {}).get(stem, {})
    if need(gate, entry.get('status') == 'PASS', f'{stem}: aggregate PASS missing'):
        need(gate, entry.get('indexed_boxes') == count, f'{stem}: expected {count} indices')
        pair(gate, stem, entry)
        audit_map(gate, entry.get('axioms', {}), 1)
need(gate, mixed.get('eta') == '4/10000000' and mixed.get('common_constant') == 24000000000000,
     'Mixed rate constants mismatch')

gate = 'native_fixed_all'
fixed = report(gate, 'native_fixed_all_verification.json')
module_list(gate, fixed.get('modules', {}))
need(gate, 'SpinCodes/Structured/ConcreteNativeFixedAll.lean' in fixed.get('modules', {}), 'FixedAll object missing')
need(gate, fixed.get('axiom_audits') == 2, 'FixedAll audit count mismatch')
audit_log(gate, 'peach_ConcreteNativeFixedAll.log', 2)

gate = 'dependency_timestamp_provenance'
provenance_filename = 'native_dependency_timestamp_verification.json'
provenance = report(gate, provenance_filename)
expected_provenance_modules = (
    {'SpinCodes/Structured/LowCancellationData/Weight1Block0.lean'} |
    {f'SpinCodes/Structured/LowCancellationData/Weight2Block{i}.lean' for i in range(62)} |
    {'SpinCodes/Structured/MapSpectrumData/Basis.lean', 'SpinCodes/Structured/SparseProgramDefs.lean',
     'SpinCodes/Structured/SparseMaxima.lean', 'SpinCodes/Structured/SparseContributionSound.lean',
     'SpinCodes/Majorant/RefinedData.lean'})
need(gate, provenance.get('checked_modules') == 68 and
     set(provenance.get('modules', {})) == expected_provenance_modules,
     'Exact 68 dependency-provenance module identities missing')
module_list(gate, provenance.get('modules', {}))
expected_provenance_inputs = {'scripts/map_data/low_cancellation_verification.json',
                              'scripts/map_data/encoder_mtime_gap_peach_verification.json'}
provenance_inputs = provenance.get('input_report_sha256', {})
need(gate, set(provenance_inputs) == expected_provenance_inputs,
     'Exact historical/fresh provenance input-report set missing')
for relative, expected in provenance_inputs.items():
    producer = report(gate, Path(relative).name)
    need(gate, inputs.get(Path(relative).name) == expected, f'Provenance producer snapshot mismatch: {relative}')
    entries = producer.get('modules', {})
    if isinstance(entries, list):
        entries = {entry['module'].replace('.', '/') + '.lean': entry for entry in entries}
    for module, entry in provenance.get('modules', {}).items():
        if entry.get('verification_report') != relative:
            continue
        compiled = entries.get(module, {})
        need(gate, compiled.get('exit_code') == 0 and
             (compiled.get('source_sha256'), compiled.get('object_sha256')) ==
             (entry.get('source_sha256'), entry.get('object_sha256')),
             f'Provenance lacks matching successful compiler record: {module}')
        if relative.endswith('encoder_mtime_gap_peach_verification.json'):
            need(gate, compiled.get('status') == 'PASS' and compiled.get('fresh_object_bytes_match_original') is True and
                 compiled.get('fresh_object_sha256') == compiled.get('original_object_sha256') == entry.get('object_sha256'),
                 f'Fresh isolated object failed byte-match evidence: {module}')
            if compiled.get('fresh_object'):
                artifact = ROOT / compiled['fresh_object']
                need(gate, digest(artifact) == entry.get('object_sha256'), f'Fresh isolated object changed: {module}')
                extra_artifacts[artifact] = entry.get('object_sha256')
for module, entry in provenance.get('modules', {}).items():
    expected_producer = ('scripts/map_data/low_cancellation_verification.json'
                         if '/LowCancellationData/' in module else 'scripts/map_data/encoder_mtime_gap_peach_verification.json')
    need(gate, entry.get('verification_report') == expected_producer, f'Wrong provenance producer: {module}')
if provenance.get('status') == 'PASS':
    checker = ROOT / 'scripts/check-native-dependency-provenance.py'
    need(gate, provenance.get('script_sha256') == digest(checker), 'Dependency provenance checker hash mismatch')
    extra_artifacts[checker] = provenance.get('script_sha256')

gate = 'actual_native_theorem_and_full_statement_pin'
base = report(gate, 'native_theorem_verification.json')
if base.get('status') == 'PASS':
    need(gate, base.get('dependency_provenance_report_sha256') == inputs.get(provenance_filename),
         'Base theorem consumes a different dependency-provenance report')
    module_list(gate, base.get('modules', {}))
    closure, pending = set(), ['SpinCodes.Structured.ConcreteNativeTheoremPin']
    while pending:
        module = pending.pop()
        relative = module.replace('.', '/') + '.lean'
        if relative in closure:
            continue
        closure.add(relative)
        try:
            source_text = (ROOT / relative).read_text(encoding='utf-8-sig')
            pending.extend(re.findall(r'^import (SpinCodes\S*)', source_text, re.M))
        except OSError as exc:
            issue(gate, f'Base import graph source missing: {exc}')
    need(gate, set(base.get('modules', {})) == closure, 'Base hash record is not the complete current project import graph')
    need(gate, base.get('project_module_count') == len(closure), 'Base project-module count mismatch')
    checks = {entry['module']: entry for entry in base.get('checks', [])}
    for stem in ['ConcreteNativeTheorem', 'ConcreteNativeTheoremPin']:
        need(gate, checks.get(stem, {}).get('exit_code') == 0, f'{stem}: successful compile record missing')
        need(gate, f'SpinCodes/Structured/{stem}.lean' in base.get('modules', {}), f'{stem}: hash pair missing')
    required = ['Spin.Structured.ConcreteNativeTheoremPin.actual_minimum_distance',
                'Spin.Structured.ConcreteNativeTheoremPin.actual_rate']
    audit_log(gate, 'native_theorem_check_ConcreteNativeTheorem.log', 2)
    audit_log(gate, 'native_theorem_check_ConcreteNativeTheoremPin.log', 2, required)
    need(gate, base.get('axiom_audits') == 4, 'Base theorem/pin audit count mismatch')
    need(gate, base.get('original_pin_sha256') == PIN_HASH, 'Base preserved-pin record mismatch')
    need(gate, {p.replace('\\', '/'): h for p, h in base.get('paper_sha256', {}).items()} == expected_papers,
         'Base preserved-paper record mismatch')
    input_hashes = base.get('input_report_sha256', {})
    need(gate, set(input_hashes) == {'native_fixed_all_verification.json', 'dense_mixed_aggregate_verification.json'},
         'Base input-report snapshot is incomplete')
    for filename, expected in input_hashes.items():
        try:
            need(gate, digest(DATA / filename) == expected, f'Base input record changed: {filename}')
        except OSError as exc:
            issue(gate, f'Base input record missing: {exc}')

gate = 'unconditional_corollaries'
consequences = report(gate, 'encoder_native_distance_consequences_verification.json')
if consequences.get('status') == 'PASS':
    module_list(gate, consequences.get('modules', []))
    stems = ['ConcreteNativeDistanceConsequences', 'ConcreteNativeDistanceConsequencesPin',
             'ConcreteNativeDistanceConsequencesFinal', 'ConcreteNativeDistanceConsequencesFinalPin']
    need(gate, {entry['module'] for entry in consequences.get('modules', [])} == set(stems),
         'Expected all four consequence module hash pairs')
    checks = {entry['module']: entry for entry in consequences.get('checks', [])}
    for stem in stems:
        need(gate, checks.get(stem, {}).get('exit_code') == 0, f'{stem}: fresh owned-module compile missing')
    need(gate, consequences.get('axiom_audit_count') == 7, 'Consequence audit count mismatch')
    audit_log(gate, 'encoder_consequence_check_ConcreteNativeDistanceConsequencesPin.log', 5)
    audit_log(gate, 'encoder_consequence_check_ConcreteNativeDistanceConsequencesFinalPin.log', 2,
              ['Spin.Structured.ConcreteNativeFamily.relative_distance_success_tendsto',
               'Spin.Structured.ConcreteNativeFamily.eventually_exists_rate_half_distance_gt_eleven_percent'])
    need(gate, consequences.get('base_report_sha256') == digest(DATA / 'native_theorem_verification.json'),
         'Corollaries consume a different base theorem report')

gate = 'fresh_final_default_invariants'
invariants = report(gate, 'native_final_invariants_verification.json')
if invariants.get('status') == 'PASS':
    command = invariants.get('command')
    need(gate, isinstance(command, list) and len(command) == 2 and
         Path(command[0]).name.lower() == 'bash.exe' and command[1] == 'scripts/check.sh' and
         invariants.get('exit_code') == 0,
         'Successful final invariant command missing')
    need(gate, invariants.get('checker_sha256') == digest(ROOT / 'scripts/check.sh'),
         'Final invariant checker source hash mismatch')
    extra_artifacts[ROOT / 'scripts/check.sh'] = invariants.get('checker_sha256')
    need(gate, invariants.get('original_pin_sha256') == PIN_HASH,
         'Final invariant preserved-pin hash mismatch')
    need(gate, {path.replace('\\', '/'): value for path, value in invariants.get('paper_sha256', {}).items()}
         == expected_papers, 'Final invariant preserved-paper hashes mismatch')
    for field, filename in [('base_report_sha256', 'native_theorem_verification.json'),
                            ('corollaries_report_sha256', 'encoder_native_distance_consequences_verification.json')]:
        need(gate, invariants.get(field) == digest(DATA / filename), f'Final invariants reference changed {filename}')
    try:
        log = Path(invariants['log'])
        if not log.is_absolute():
            log = ROOT / log if str(log).replace('\\', '/').startswith('scripts/') else DATA / log
        need(gate, digest(log) == invariants.get('log_sha256'), 'Final invariant log hash mismatch')
        need(gate, 'ALL CHECKS PASS' in log.read_text(encoding='utf-8'), 'Final invariant log lacks success marker')
        extra_artifacts[log] = invariants.get('log_sha256')
    except (OSError, KeyError, TypeError) as exc:
        issue(gate, f'Final invariant log missing: {exc}')

# Refuse a PASS if any producer record or checked artifact changed during this scan.
if not issues:
    gates['stable_final_snapshot'] = {'problems': [], 'reports': []}
    need('stable_final_snapshot', digest(ROOT / 'SpinCodes/Pin.lean') == PIN_HASH,
         'Original pin changed during scan')
    need('stable_final_snapshot',
         {path.relative_to(ROOT.parent).as_posix(): digest(path)
          for path in sorted((ROOT.parent / 'paper').glob('*.tex'))} == expected_papers,
         'Original paper set/hashes changed during scan')
    for filename, expected in inputs.items():
        need('stable_final_snapshot', digest(DATA / filename) == expected, f'Input changed during scan: {filename}')
    for relative, expected in pairs.items():
        obj = (ROOT / '.lake/build/lib/lean' / relative).with_suffix('.olean')
        need('stable_final_snapshot', (digest(ROOT / relative), digest(obj)) == expected,
             f'Source/object changed during scan: {relative}')
    for path, expected in extra_artifacts.items():
        need('stable_final_snapshot', digest(path) == expected, f'Final invariant artifact changed during scan: {path.name}')

for entry in gates.values():
    entry['status'] = 'PASS' if not entry['problems'] else 'INCOMPLETE'
individual_pairs = sorted(path for path, kinds in pair_kinds.items() if 'individual_check_record_pair' in kinds)
snapshot_only_pairs = sorted(path for path, kinds in pair_kinds.items() if kinds == {'final_assembly_dependency_snapshot'})
coverage = {'scope': 'Classification within the mandatory reports consumed by this run. A final assembly dependency snapshot is not evidence of a fresh individual compile. Other individual records may exist outside these mandatory gates.',
            'individual_check_record_pairs': len(individual_pairs),
            'final_assembly_snapshot_only_pairs': len(snapshot_only_pairs),
            'final_assembly_snapshot_only_modules': snapshot_only_pairs}
coverage_file = DATA / 'encoder_native_evidence_coverage.json'
if coverage_file.exists():
    try:
        coverage_raw = coverage_file.read_bytes()
        external = json.loads(coverage_raw)
        coverage['external_evidence_review'] = {
            'record': coverage_file.name, 'record_sha256': hashlib.sha256(coverage_raw).hexdigest(),
            'status': external.get('status'), 'finished_utc': external.get('finished_utc'),
            'counts': external.get('counts'),
            'qualification': 'Historical read-only coverage snapshot, not a mandatory PASS gate and not revalidated by this run. Its pass_record_pair and completed_fourier_batch_pair classes are distinct from initial/cache/final-assembly snapshots; source-only records remain source-only.'}
    except (OSError, json.JSONDecodeError) as exc:
        coverage['external_evidence_review'] = {'status': 'UNAVAILABLE', 'diagnostic': str(exc)}
result = {
    'status': 'PASS' if not issues else 'INCOMPLETE',
    'checked_at_utc': datetime.now(timezone.utc).isoformat(),
    'scope': 'Read-only consolidation of existing individually kernel-checked verification reports. Cached project and mathlib dependencies were reused; this is not a fresh full source/dependency replay, a new Lean check, or independent proof authority.',
    'script_sha256': digest(Path(__file__)),
    'expected_counts': {'scalar_indexed_boxes': 283, 'occupation_indexed_boxes': 433,
                        'fourier_indexed_boxes': 307, 'global_boxes': 1023,
                        'fourier_witnesses': 105, 'occupation_witnesses': 133},
    'gates': gates, 'missing_or_invalid_gates': [name for name, entry in gates.items() if entry['problems']],
    'issues': issues, 'current_recorded_pairs_checked': len(pairs),
    'evidence_classification': coverage,
    'input_reports_sha256': inputs, 'original_pin_sha256': PIN_HASH, 'paper_sha256': current_papers,
}
publish(result)
print(json.dumps({'status': result['status'], 'gates': len(gates),
                  'current_recorded_pairs_checked': len(pairs),
                  'missing_or_invalid_gates': result['missing_or_invalid_gates'],
                  'issues': len(issues), 'record': DESTINATION.name}))
