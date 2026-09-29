"""Generate untrusted Fourier rows and integer cap certificates for Lean replay."""
import hashlib
import json
from pathlib import Path

from map_fiber_preview import choose

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT/'scripts/map_data'
OUT = ROOT/'SpinCodes/Structured/FiberNumericsData'


def emit(path, lines):
    path.write_text('\n'.join(lines)+'\n',encoding='utf-8')
    return dict(path=path.relative_to(ROOT).as_posix(),
                sha256=hashlib.sha256(path.read_bytes()).hexdigest())


def main():
    OUT.mkdir(exist_ok=True)
    spectra = json.loads((DATA/'spectra.json').read_text())
    model = json.loads((ROOT/'scripts/sparse_data/model.json').read_text())
    candidates = json.loads((DATA/'fiber_candidates.json').read_text())['records']
    table = ['import SpinCodes.Structured.FiberNumericsDefs',
             'namespace Spin.Structured.FiberNumerics.Data']
    for name,values in [('spectrum',spectra['ct']),('kernels',model['kernel']),('caps',model['caps'])]:
        table.append(f'def {name} : List Nat := ['+', '.join(map(str,values))+']')
    table.append('end Spin.Structured.FiberNumerics.Data')
    table_record = emit(OUT/'Tables.lean',table)
    (DATA/'fiber_tables_manifest.json').write_text(json.dumps(table_record,indent=2)+'\n')
    manifest = []
    for j,record in enumerate(candidates):
        assert record['weight'] == j and j < 129
        vals = [sum((-1)**h * choose(w,h) * choose(128-w,j-h)
                    for h in range(j+1)) for w in range(129)]
        method = record['methods'][0]
        lines = ['import SpinCodes.Structured.FiberNumericsData.Tables',
                 'namespace Spin.Structured.FiberNumerics.Data',
                 'set_option maxRecDepth 1000000', 'set_option maxHeartbeats 0',
                 f'def values{j} : List Int := ['+', '.join(map(str,vals))+']',
                 f'theorem values{j}_checked : (List.range 129).map (kraw {j}) = values{j} := by decide +kernel',
                 f'theorem signed{j}_checked : signedSum spectrum values{j} = 524288 * (kernels.getD {j} 0 : Int) := by decide +kernel',
                 f'theorem square{j}_checked : squareSum spectrum values{j} = 524288 * ({record["pairs"]} : Int) := by decide +kernel',
                 f'theorem absolute{j}_checked : absSum spectrum values{j} = {record["absolute"]} := by decide +kernel']
        if method == 'fourier':
            condition = f'({record["absolute"]} : Int) < 524288 * ((caps.getD {j} 0 : Int) + 1)'
        elif method == 'complement':
            condition = f'polyChoose 128 {j} - kernels.getD {j} 0 ≤ caps.getD {j} 0'
        elif method == 'packing':
            condition = f'polyChoose 128 ({j} - 1) < (caps.getD {j} 0 + 1) * polyChoose {j} ({j} - 1)'
        elif method == 'packing_compl':
            condition = f'polyChoose 128 ((128 - {j}) - 1) < (caps.getD {j} 0 + 1) * polyChoose (128 - {j}) ((128 - {j}) - 1)'
        else:
            raise AssertionError(method)
        lines += [f'theorem cap{j}_checked : {condition} := by decide +kernel',
                  'end Spin.Structured.FiberNumerics.Data']
        manifest.append(dict(weight=j,method=method,**emit(OUT/f'Weight{j}.lean',lines)))
    (DATA/'fiber_manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
    print('Generated 129 Fourier rows and cap checks; these require kernel replay.')


if __name__ == '__main__':
    main()
