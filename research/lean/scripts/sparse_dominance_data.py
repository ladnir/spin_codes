"""Emit kernel checks for every maximum replacement in the sparse program."""
import hashlib
import json
import math
from pathlib import Path
from flint import fmpq, fmpq_poly
import sparse_bridge_data as bridge

ROOT = bridge.LEAN


def main():
    model = bridge.verifier.screen.Model()
    e = model.engine
    z = fmpq_poly([1, fmpq(-1,6250)])
    powers = [z**i for i in range(129)]
    directory = ROOT/'SpinCodes/Structured/SparseBridge'
    manifest = []
    comparisons = 0
    for j in range(129):
        total = math.comb(128,j)
        moments = [sum((math.comb(w,h)*math.comb(128-w,j-h)*powers[w+j-2*h]
                        for h in range(max(0,j-128+w),min(w,j)+1)), fmpq_poly()) / total
                   for w in e.levels]
        selected = max(range(5),key=lambda i: moments[i](fmpq(1,2)))
        lines = ['import SpinCodes.Structured.SparseModelData',
                 'import SpinCodes.Structured.PolySignIdentityDefs',
                 f'namespace Spin.Structured.SparsePolynomial.Dominance{j}',
                 'set_option maxRecDepth 1000000', 'set_option maxHeartbeats 0']
        for i,(w,m) in enumerate(zip(e.levels,moments)):
            lines += [f'def m{i} : RatPoly := {bridge.ratpoly(m)}',
                      f'theorem m{i}_checked : RatPoly.checkEq (moment {w} {j}) m{i} = true := by decide']
        lines += [f'theorem selected_checked : Data.weight{j}.choices.momentMax = {selected} := by decide']
        for i in range(5):
            bridge.verifier.sign_certificate(moments[i]-moments[selected])
            lines += [f'theorem m{i}_le_checked : RatPoly.checkLE m{i} m{selected} = true := by decide']
            comparisons += 1
        if j in e.low:
            patterns = e.low[j]['patterns']
            ps = [sum((ct*powers[w] for w,ct in p),fmpq_poly())/total for p in patterns]
            if ps:
                selected_low = max(range(len(ps)),key=lambda i: ps[i](fmpq(1,2)))
                lines += [f'theorem selected_low_checked : Data.weight{j}.choices.lowMax = {selected_low} := by decide']
                for i,(p,polynomial) in enumerate(zip(patterns,ps)):
                    lines += [f'def low{i} : RatPoly := {bridge.ratpoly(polynomial)}',
                              f'theorem low{i}_checked : RatPoly.checkEq (pattern {j} {bridge.lean_pattern(p)}) low{i} = true := by decide']
                for i,p in enumerate(ps):
                    bridge.verifier.sign_certificate(p-ps[selected_low])
                    lines += [f'theorem low{i}_le_checked : RatPoly.checkLE low{i} low{selected_low} = true := by decide']
                    comparisons += 1
        lines.append(f'end Spin.Structured.SparsePolynomial.Dominance{j}')
        output = directory/f'Dominance{j}.lean'
        output.write_text('\n'.join(lines)+'\n',encoding='utf-8')
        manifest.append(dict(weight=j,path=str(output.relative_to(ROOT)),
                             sha256=hashlib.sha256(output.read_bytes()).hexdigest()))
    (ROOT/'scripts/sparse_data/dominance_manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
    validation = ['import SpinCodes.Structured.SparseModelData',
                  'import SpinCodes.Structured.SparseValidationDefs',
                  'namespace Spin.Structured.SparsePolynomial.Data',
                  'set_option maxRecDepth 100000', 'set_option maxHeartbeats 0']
    for j in range(129):
        validation.append(f'theorem weight{j}_valid : checkData {j} weight{j} = true := by decide')
    validation.append('end Spin.Structured.SparsePolynomial.Data')
    (ROOT/'SpinCodes/Structured/SparseDataChecks.lean').write_text('\n'.join(validation)+'\n',encoding='utf-8')
    print(f'Emitted {comparisons} sign comparisons and their polynomial identities.')


if __name__ == '__main__':
    main()
