"""Transcribe the paper's map table and generate exact inverse witnesses.

This script is untrusted. Every emitted finite claim is checked with Lean's
kernel, and PackedMap.eval_inverse proves its consequence for all inputs.
"""
from pathlib import Path
import re

ROOT = Path(__file__).resolve().parents[1]
PAPER = ROOT.parent / 'paper/structured_imt_appendix.tex'


def encode(rows, q):
    result = 0
    for row in rows:
        if q & 1:
            result ^= row
        q >>= 1
    return result


def transpose(rows, width):
    return [sum(((row >> i) & 1) << j for j, row in enumerate(rows))
            for i in range(width)]


def elimination(rows):
    basis = {}
    for j, row in enumerate(rows):
        v, coefficient = row, 1 << j
        while v:
            pivot = v.bit_length()-1
            if pivot not in basis:
                basis[pivot] = (v, coefficient)
                break
            r, c = basis[pivot]
            v ^= r
            coefficient ^= c
    def solve(v):
        coefficient = 0
        for pivot in sorted(basis, reverse=True):
            if (v >> pivot) & 1:
                r, c = basis[pivot]
                v ^= r
                coefficient ^= c
        return coefficient, v
    return solve


def main():
    entries = re.findall(r'(?m)^(\d+) & \\texttt\{([0-9a-f]+)\} & \\texttt\{([0-9a-f]+)\}',
                         PAPER.read_text(encoding='utf-8'))
    assert [int(j) for j, _, _ in entries] == list(range(19))
    a = [int(x, 16) for _, x, _ in entries]
    ct = [int(x, 16) for _, _, x in entries]
    c = transpose(ct, 128)
    adecode = [elimination(a)(1 << i)[0] for i in range(128)]
    cright = [elimination(c)(1 << i)[0] for i in range(19)]
    ctdecode = transpose(cright, 128)
    assert [encode(adecode, r) for r in a] == [1 << i for i in range(19)]
    assert [encode(c, r) for r in cright] == [1 << i for i in range(19)]
    assert [encode(ctdecode, r) for r in ct] == [1 << i for i in range(19)]
    def literal(name, values):
        return f'def {name} : List Nat :=\n  [' + ',\n   '.join(hex(x) for x in values) + ']\n'
    source = ['import SpinCodes.Structured.PackedMapDefs',
              '', 'namespace Spin.Structured.ConcreteMaps',
              'open PackedMap', 'set_option maxRecDepth 100000',
              'set_option maxHeartbeats 0',
              '/-- Basis images of A from `tab:imt-maps`, least significant bit first. -/',
              literal('aRows', a),
              '/-- Basis images of C transpose from the same table. -/',
              literal('cTransposeRows', ct),
              'def cRows : List Nat := transpose 128 cTransposeRows',
              literal('aLeftInverse', adecode),
              literal('cRightInverse', cright),
              literal('cTransposeLeftInverse', ctdecode)]
    checks = {
        'a_length': 'aRows.length = 19',
        'cTranspose_length': 'cTransposeRows.length = 19',
        'c_length': 'cRows.length = 128',
        'a_rows_bounded': 'aRows.all (fun r => r < 2 ^ 128) = true',
        'cTranspose_rows_bounded': 'cTransposeRows.all (fun r => r < 2 ^ 128) = true',
        'c_rows_bounded': 'cRows.all (fun r => r < 2 ^ 19) = true',
        'cRightInverse_bounded': 'cRightInverse.all (fun r => r < 2 ^ 128) = true',
        'a_inverse_checked': 'aRows.map (eval aLeftInverse) = identityRows 19',
        'c_inverse_checked': 'cRightInverse.map (eval cRows) = identityRows 19',
        'cTranspose_inverse_checked': 'cTransposeRows.map (eval cTransposeLeftInverse) = identityRows 19',
        'a_coordinates_nonzero_checked': '(transpose 128 aRows).all (fun r => r != 0) = true',
        'c_columns_weight_checked': 'cRows.all (fun r => weight 19 r == 5) = true',
        'c_columns_distinct_checked': 'cRows.Nodup',
    }
    source += [f'theorem {name} : {claim} := by decide' for name, claim in checks.items()]
    source += ['end Spin.Structured.ConcreteMaps', '']
    target = ROOT/'SpinCodes/Structured/ConcreteMapData.lean'
    target.write_text('\n'.join(source), encoding='utf-8')
    print(f'Wrote {target.name}: exact paper rows and {len(checks)} kernel checks.')


if __name__ == '__main__':
    main()
