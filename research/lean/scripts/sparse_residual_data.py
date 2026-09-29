"""Emit exact identities summing the 903 contributions into the seven rows."""
import hashlib
import json
import math
from fractions import Fraction
from pathlib import Path
import re
from flint import fmpq, fmpq_poly
from sparse_bridge_data import ratpoly

ROOT = Path(__file__).resolve().parents[1]
CHUNK = 8


def polynomials(path):
    text = path.read_text(encoding='utf-8')
    return {name: fmpq_poly([fmpq(int(c),int(den)) for c in nums.split(',') if c.strip()])
            for name,nums,den in re.findall(r'def (\w+) : RatPoly := ⟨\[(.*?)\], (\d+)⟩',text)}


def denominator(p):
    return math.lcm(*(Fraction(str(c)).denominator for c in p.coeffs())) if p else 1


def main():
    directory = ROOT/'SpinCodes/Structured/SparseBridge'
    data = [polynomials(directory/f'Weight{j}.lean') for j in range(129)]
    manifest = []
    sums = []
    for row in range(7):
        pieces = []
        for part,start in enumerate(range(0,129,CHUNK)):
            weights = list(range(start,min(129,start+CHUNK)))
            poly = sum((data[j]['prob']*data[j][f'a{row}'] for j in weights),fmpq_poly())
            pieces.append(poly)
            common = math.lcm(denominator(poly), *(denominator(data[j]['prob'])*denominator(data[j][f'a{row}']) for j in weights))
            lines = ['import SpinCodes.Structured.PolyPackedDefs']
            lines += [f'import SpinCodes.Structured.SparseBridge.Weight{j}' for j in weights]
            lines += [f'namespace Spin.Structured.SparsePolynomial.Row{row}Part{part}',
                      'set_option maxRecDepth 1000000', 'set_option maxHeartbeats 0',
                      'def terms : List PackedTerm := [' +
                      ', '.join(f'⟨Weight{j}.prob, Weight{j}.a{row}, {common//(denominator(data[j]["prob"])*denominator(data[j][f"a{row}"]))}⟩' for j in weights) + ']',
                      f'def expected : RatPoly := {ratpoly(poly)}',
                      f'theorem checked : checkSumProducts {common} terms expected {common//denominator(poly)} = true := by decide',
                      f'end Spin.Structured.SparsePolynomial.Row{row}Part{part}']
            output = directory/f'Row{row}Part{part}.lean'
            output.write_text('\n'.join(lines)+'\n',encoding='utf-8')
            manifest.append(dict(weight=len(manifest), row=row, part=part, input_weights=weights,
                                 path=str(output.relative_to(ROOT)),
                                 sha256=hashlib.sha256(output.read_bytes()).hexdigest()))
        sums.append(sum(pieces,fmpq_poly()))
    model = json.loads((ROOT/'scripts/sparse_data/model.json').read_text())
    x = fmpq_poly([0,1])
    vector = [fmpq_poly([fmpq(c) for c in cs]) for cs in model['vector']]
    root = [f'import SpinCodes.Structured.SparseBridge.Row{r}Part{p}'
            for r in range(7) for p in range((129+CHUNK-1)//CHUNK)]
    root += [f'import SpinCodes.Structured.SparseData.Row{r}' for r in range(7)]
    root += ['namespace Spin.Structured.SparsePolynomial.Residuals',
             'set_option maxRecDepth 1000000', 'set_option maxHeartbeats 0',
             'def negativeContraction : RatPoly := ⟨[-10000, 96], 10000⟩']
    for row in range(7):
        residual = sums[row] - (1-96*x/10000)*vector[row]
        source = (ROOT/f'SpinCodes/Structured/SparseData/Row{row}.lean').read_text(encoding='utf-8')
        nums = re.search(r'def row\d+\s*:\s*List (?:ℤ|Int)\s*:=\s*\[(.*?)\]',source,re.S)[1]
        original_den = int(model['normalizations'][row]['denominator'])
        coeffs = [0]+[int(s.strip().replace('(', '').replace(')', '')) for s in nums.split(',') if s.strip()]
        assert fmpq_poly([fmpq(c,original_den) for c in coeffs]) == residual
        parts = [sum((data[j]['prob']*data[j][f'a{row}'] for j in range(start,min(129,start+CHUNK))),fmpq_poly())
                 for start in range(0,129,CHUNK)]
        common = math.lcm(original_den, 10000*denominator(vector[row]), *(denominator(p) for p in parts))
        terms = [f'⟨Row{row}Part{p}.expected, RatPoly.constant 1, {common//denominator(poly)}⟩' for p,poly in enumerate(parts)]
        terms.append(f'⟨negativeContraction, vector{row}, {common//(10000*denominator(vector[row]))}⟩')
        root += [f'def vector{row} : RatPoly := {ratpoly(vector[row])}',
                 f'theorem vector{row}_checked : RatPoly.checkEq (witness {row}) vector{row} = true := by decide',
                 f'def terms{row} : List PackedTerm := [' + ', '.join(terms) + ']',
                 f'def original{row} : RatPoly := ⟨0 :: SparseData.row{row}, {original_den}⟩',
                 f'theorem row{row}_checked : checkSumProducts {common} terms{row} original{row} {common//original_den} = true := by decide']
    root += ['end Spin.Structured.SparsePolynomial.Residuals']
    (ROOT/'SpinCodes/Structured/SparseResiduals.lean').write_text('\n'.join(root)+'\n',encoding='utf-8')
    (ROOT/'scripts/sparse_data/contribution_manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
    print(f'Emitted {len(manifest)} contribution sums and seven exact residual identities.')


if __name__ == '__main__':
    main()
