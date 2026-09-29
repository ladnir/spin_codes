"""Read-only structural audit of the incremental Fourier box replay."""
from pathlib import Path
import hashlib
import json
import re

root = Path(__file__).resolve().parents[1]
source = json.loads((root.parent / 'workstreams/inner_design/imt_asymptotic/d11/DENSE_REPLAY.json').read_text())
selected = [(i, b) for i, b in enumerate(source['leaves']) if b['rational_witness']['family'] == 'fourier']
assert len(selected) == 307
unique = {}
for _, box in selected:
    key = json.dumps(box['rational_witness'], sort_keys=True)
    if key not in unique:
        unique[key] = len(unique)
names = ['dense_fourier_box_replay.json', 'dense_fourier_box_second_replay.json',
         'dense_fourier_box_prefix_tail_replay.json', 'dense_fourier_box_suffix_tail_replay.json']
records = []
for name in names:
    state = json.loads((root / 'scripts/map_data' / name).read_text())
    assert state['checked'] == len(state['boxes']), name
    records.extend(state['boxes'])
indices = [record['index'] for record in records]
assert len(indices) == len(set(indices)), 'duplicate ownership'
for record in records:
    index = record['index']
    global_index, box = selected[index]
    witness = unique[json.dumps(box['rational_witness'], sort_keys=True)]
    assert (record['global_index'], record['witness_index']) == (global_index, witness)
    assert record['audits'] == 4
    for module in record['modules']:
        path = root / module['source'].replace('\\', '/')
        assert hashlib.sha256(path.read_bytes()).hexdigest() == module['source_sha256'], path

emitted = []
for index, (global_index, box) in enumerate(selected):
    tag = f'B{index:03d}'
    path = root / f'SpinCodes/Structured/DenseFourierExact{tag}.lean'
    if not path.exists():
        continue
    emitted.append(index)
    text = path.read_text(encoding='utf-8-sig')
    witness = unique[json.dumps(box['rational_witness'], sort_keys=True)]
    assert f'namespace Spin.Structured.DenseFourierExact.{tag}' in text, path
    assert f'import SpinCodes.Structured.DenseFourierExact{tag}Data' in text, path
    assert f'import SpinCodes.Structured.DenseFourierExactW{witness:03d}' in text, path
    assert f'geometry_match : boxes.getD {global_index} zeroRect' in text, path
    assert f'certified_uniform : CertifiedBox {global_index} (4/10000000) 24000000000000' in text, path

aggregate = (root / 'SpinCodes/Structured/DenseOccupationFourierCertified.lean').read_text(encoding='utf-8-sig')
imports = [int(i) for i in re.findall(r'^import SpinCodes.Structured.DenseFourierExactB(\d{3})$', aggregate, re.M)]
terms = [int(i) for i in re.findall(r'DenseFourierExact.B(\d{3}).certified_uniform', aggregate)]
assert imports == list(range(307)), 'aggregate imports'
assert terms == list(range(307)), 'aggregate proof terms'
global_text = (root / 'SpinCodes/Structured/DenseOccupationGlobalIndices.lean').read_text(encoding='utf-8-sig')
global_list = json.loads(re.search(r'def fourierIndices : List ℕ := (\[[^\n]*\])', global_text).group(1))
assert global_list == [i for i, _ in selected], 'global Fourier index list'
print(json.dumps({'status': 'PASS', 'checked': len(indices), 'emitted': len(emitted),
                  'total': 307, 'witnesses': len(unique), 'duplicate_indices': 0,
                  'missing_are_pending': sorted(set(range(307)) - set(emitted)),
                  'scope': 'Structural source/index/import/hash audit; does not compile pending boxes.'}))
