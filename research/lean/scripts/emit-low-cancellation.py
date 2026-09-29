"""Emit untrusted low-layer witnesses, checked by Lean's kernel.

The maps are read from the authoritative Lean hexadecimal tables. Every
syndrome, weight, pattern, coverage and shell histogram is subsequently checked.
"""
import argparse
import ast
from collections import defaultdict
import itertools
from pathlib import Path
import re

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / 'SpinCodes/Structured/LowCancellationData'


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--chunk', type=int, default=128)
    args = parser.parse_args()
    OUT.mkdir(exist_ok=True)
    maps = (ROOT / 'SpinCodes/Structured/ConcreteMapData.lean').read_text(encoding='utf-8')
    def rows(name):
        body = re.search(r'def ' + name + r' : List Nat :=\s*\[([^]]*)\]', maps).group(1)
        return [int(x.strip(), 0) for x in body.split(',') if x.strip()]
    ar, ct = rows('aRows'), rows('cTransposeRows')
    cr = [sum(((r >> i) & 1) << k for k, r in enumerate(ct)) for i in range(128)]
    def expand(q):
        ans = 0
        for i, a in enumerate(ar):
            if q >> i & 1:
                ans ^= a
        return ans
    sparse = (ROOT / 'SpinCodes/Structured/SparseModelData.lean').read_text(encoding='utf-8')
    for j in (1, 2):
        data = sparse.split(f'def weight{j} :')[1].split('def weight')[0]
        patterns = ast.literal_eval(re.search(r'some (.*),\n', data).group(1))
        shells = ast.literal_eval(data.splitlines()[2].strip().removesuffix(','))
        groups = defaultdict(list)
        for inds in itertools.combinations(range(128), j):
            x, q = 0, 0
            for i in inds:
                x |= 1 << i
                q ^= cr[i]
            a = expand(q)
            groups[q].append((x, (x ^ a).bit_count()))
        grouped = [(q, expand(q).bit_count(), sorted(es, key=lambda r: (r[1], r[0])))
                   for q, es in sorted(groups.items())]
        assert all(q for q, _, _ in grouped)
        expanded = [[e for e, c in p for _ in range(c)] for p in patterns]
        assert all([e for _, e in es] in expanded for _, _, es in grouped)
        for i, w in enumerate((48, 56, 64, 72, 80)):
            assert sorted(e for _, s, es in grouped if s == w for _, e in es) == [e for e,c in shells[i] for _ in range(c)]
        chunks = []
        for k, start in enumerate(range(0, len(grouped), args.chunk)):
            name = f'Weight{j}Block{k}'
            chunks.append(name)
            lines = ['import SpinCodes.Structured.LowCancellationDefs',
                     'import SpinCodes.Structured.SparseModelData',
                     'namespace Spin.Structured.LowCancellation.Data',
                     'set_option maxRecDepth 100000', 'set_option maxHeartbeats 0',
                     f'def {name} : List Group := [']
            batch = grouped[start:start+args.chunk]
            for b, (q,w,es) in enumerate(batch):
                entries = ', '.join(f'(0x{x:x}, {e})' for x,e in es)
                lines.append(f'  ⟨{q}, {w}, [{entries}]⟩' + (',' if b+1 < len(batch) else ']'))
            lines += [f'theorem {name}_valid : ∀ g ∈ {name}, valid {j} g := by decide +kernel',
                      f'theorem {name}_patterns : ∀ g ∈ {name}, g.entries.map Prod.snd ∈',
                      f'    ((SparsePolynomial.Data.weight{j}).lowPatterns.getD []).map expand := by decide +kernel',
                      'end Spin.Structured.LowCancellation.Data', '']
            (OUT / f'{name}.lean').write_text('\n'.join(lines), encoding='utf-8')
        lines = [f'import SpinCodes.Structured.LowCancellationData.{name}' for name in chunks]
        lines += ['import SpinCodes.Structured.LowCancellationBridge',
                  'namespace Spin.Structured.LowCancellation.Data',
                  'set_option maxRecDepth 100000', 'set_option maxHeartbeats 0',
                  f'def groups{j} : List Group := ' + ' ++ '.join(chunks),
                  f'theorem groups{j}_valid : ∀ g ∈ groups{j}, valid {j} g := by',
                  f'  simp only [groups{j}, List.mem_append, or_imp, forall_and, and_assoc]',
                  '  exact ⟨' + ', '.join(name+'_valid' for name in chunks) + '⟩' if len(chunks)>1 else f'  exact {chunks[0]}_valid',
                  f'theorem groups{j}_patterns : ∀ g ∈ groups{j}, g.entries.map Prod.snd ∈',
                  f'    ((SparsePolynomial.Data.weight{j}).lowPatterns.getD []).map expand := by',
                  f'  simp only [groups{j}, List.mem_append, or_imp, forall_and, and_assoc]',
                  '  exact ⟨' + ', '.join(name+'_patterns' for name in chunks) + '⟩' if len(chunks)>1 else f'  exact {chunks[0]}_patterns',
                  f'theorem groups{j}_length : (inputs groups{j}).length = Nat.choose 128 {j} := by decide +kernel',
                  f'theorem groups{j}_inputs : (inputs groups{j}).Nodup :=',
                  f'  sorted_nodup _ (by decide +kernel)',
                  f'theorem groups{j}_keys : (groups{j}.map Group.syndrome).Nodup := by',
                  '  apply List.Pairwise.imp (fun h => Nat.ne_of_lt h)',
                  '  apply List.isChain_iff_pairwise.mp',
                  '  decide +kernel']
        for i,w in enumerate((48,56,64,72,80)):
            lines += [f'theorem groups{j}_shell{i} : certSort (shellExponents groups{j} {w}) =',
                      f'    expand ((SparsePolynomial.Data.weight{j}).lowShells.getD {i} []) := by decide +kernel']
        lines += ['end Spin.Structured.LowCancellation.Data', '']
        (OUT/f'Weight{j}.lean').write_text('\n'.join(lines), encoding='utf-8')
        print(f'Weight {j}: {sum(len(es) for _,_,es in grouped)} inputs, {len(grouped)} syndromes, {len(chunks)} blocks.')


if __name__ == '__main__':
    main()
