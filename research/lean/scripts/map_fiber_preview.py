"""Prepare untrusted integer witnesses for the concrete fiber-cap connection.

This uses proposed spectrum totals. Its output is planning data, not a Lean
certificate; both the spectrum and these arithmetic identities need replay.
"""
from collections import Counter
import json
from math import comb
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT/'scripts/map_data'


def choose(n, k):
    return comb(n, k) if 0 <= k <= n else 0


def main():
    spectra = json.loads((DATA/'spectra.json').read_text())
    model = json.loads((ROOT/'scripts/sparse_data/model.json').read_text())
    order, nonzero = 2**19, 2**19-1
    records = []
    for j in range(129):
        values = [sum((-1)**h * choose(w,h) * choose(128-w,j-h)
                      for h in range(j+1)) for w in range(129)]
        fourier = sum(b*k for b,k in zip(spectra['ct'],values))
        parseval = sum(b*k*k for b,k in zip(spectra['ct'],values))
        absolute = sum(b*abs(k) for b,k in zip(spectra['ct'],values))
        assert fourier % order == parseval % order == 0
        kernel, pairs = fourier//order, parseval//order
        assert kernel == model['kernel'][j]
        cap, a = model['caps'][j], choose(128,j)-kernel
        squares = pairs-kernel*kernel
        checks = dict(complement=a <= cap,
            fourier=absolute < order*(cap+1),
            packing=choose(128,max(0,j-1)) < (cap+1)*choose(j,max(0,j-1)),
            packing_compl=choose(128,max(0,127-j)) < (cap+1)*choose(128-j,max(0,127-j)),
            variance=a <= nonzero*(cap+1) and
                2*a*(cap+1)+(nonzero-1)*squares < nonzero*(cap+1)**2+a*a)
        assert any(checks.values()), (j,cap,checks)
        records.append(dict(weight=j,kernel=kernel,pairs=pairs,absolute=absolute,
            cap=cap,methods=[method for method,holds in checks.items() if holds]))
    report = dict(status='UNTRUSTED_PREVIEW',
        scope='All 129 proposed kernel counts match the frozen model and each proposed cap has an integer witness. Requires Lean replay and spectrum identification.',
        records=records)
    (DATA/'fiber_candidates.json').write_text(json.dumps(report,indent=2)+'\n')
    print('129/129 proposed kernel counts match; 129/129 proposed caps have candidate witnesses.')
    print(dict(Counter(r['methods'][0] for r in records)))
    print('UNTRUSTED_PREVIEW: no Lean proof claimed.')


if __name__ == '__main__':
    main()
