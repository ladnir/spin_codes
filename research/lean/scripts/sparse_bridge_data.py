"""Recover exact sparse model data and positive row normalizations.

This is an untrusted generator. It checks its transcription against the
existing artifacts; Lean must separately verify every polynomial identity.
"""
from fractions import Fraction
import hashlib
import json
import math
from pathlib import Path
import re
import sys
from flint import fmpq, fmpq_poly

LEAN = Path(__file__).resolve().parents[1]
SOURCE = LEAN.parent / 'workstreams/inner_design/imt_asymptotic'
sys.path.insert(0, str(SOURCE))
import certify_imt_sparse as verifier


def ratpoly(p):
    coeffs = [Fraction(str(c)) for c in p.coeffs()]
    denominator = math.lcm(*(c.denominator for c in coeffs)) if coeffs else 1
    return '⟨[' + ', '.join(str(int(c * denominator)) for c in coeffs) + f'], {denominator}⟩'


def lean_list(xs):
    return '[' + ', '.join(map(str, xs)) + ']'


def lean_pattern(p):
    return '[' + ', '.join(f'({w}, {c})' for w, c in p) + ']'


def emit_powers():
    x = fmpq_poly([0, 1])
    lines = ['import SpinCodes.Structured.PolyIdentityDefs',
             'namespace Spin.Structured.SparsePowers',
             'set_option maxRecDepth 1000000', 'set_option maxHeartbeats 0']
    for name, base in [('z', 1-x/6250), ('beta', x/12500), ('complement', 1-x/12500)]:
        for n in range(129):
            lines.append(f'def {name}{n} : RatPoly := {ratpoly(base**n)}')
        lines.append(f'def {name}Powers : List RatPoly := ' + lean_list(f'{name}{i}' for i in range(129)))
        lines.append(f'def {name}Pow (n : Nat) : RatPoly := {name}Powers.getD n ({name}1.pow n)')
        lines.append(f'theorem {name}_zero_checked : RatPoly.checkEq {name}0 (RatPoly.constant 1) = true := by decide')
        for n in range(128):
            lines.append(f'theorem {name}_step{n}_checked : RatPoly.checkMul {name}{n} {name}1 {name}{n+1} = true := by decide')
    lines.append('end Spin.Structured.SparsePowers')
    (LEAN / 'SpinCodes/Structured/SparsePowers.lean').write_text('\n'.join(lines)+'\n',encoding='utf-8')
    proof = ['import SpinCodes.Structured.SparsePowers',
             'import SpinCodes.Structured.PolyIdentity',
             'namespace Spin.Structured.SparsePowers',
             'set_option maxRecDepth 100000', 'set_option maxHeartbeats 0']
    for name, expression in [('z', '(1 - x / 6250)'), ('beta', '(x / 12500)'),
                             ('complement', '(1 - x / 12500)')]:
        proof += [f'theorem {name}0_eval (x : ℝ) : {name}0.eval x = {expression} ^ 0 := by',
                  f'  norm_num [{name}0, RatPoly.eval, listEval]',
                  f'theorem {name}1_eval (x : ℝ) : {name}1.eval x = {expression} := by',
                  f'  simp only [{name}1, RatPoly.eval, listEval]', '  push_cast', '  ring']
        for n in range(2,129):
            proof += [f'theorem {name}{n}_eval (x : ℝ) : {name}{n}.eval x = {expression} ^ {n} := by',
                      f'  rw [RatPoly.checkMul_sound _ _ _ {name}_step{n-1}_checked,',
                      f'    {name}{n-1}_eval, {name}1_eval, ← pow_succ]']
        # At n=2 the preceding theorem has exponent one simplified.
        idx = proof.index(f'    {name}1_eval, {name}1_eval, ← pow_succ]')
        proof[idx] = f'    {name}1_eval, pow_two]'
        proof += [f'theorem {name}Pow_eval (n : ℕ) (x : ℝ) : ({name}Pow n).eval x = {expression} ^ n := by',
                  '  by_cases h : n < 129', '  · interval_cases n']
        for n in range(129):
            proof += [f'    · '+(f'rw [pow_one]; exact {name}1_eval x' if n==1 else f'exact {name}{n}_eval x')]
        proof += [f'  · rw [{name}Pow, List.getD_eq_default _ _ (by simpa only [{name}Powers, List.length_cons, List.length_nil] using Nat.le_of_not_gt h)]',
                  f'    rw [RatPoly.eval_pow, {name}1_eval]',
                  f'theorem {name}Pow_den_pos (n : ℕ) : 0 < ({name}Pow n).den := by',
                  '  by_cases h : n < 129', '  · interval_cases n <;> decide',
                  f'  · rw [{name}Pow, List.getD_eq_default _ _ (by simpa only [{name}Powers, List.length_cons, List.length_nil] using Nat.le_of_not_gt h)]',
                  f'    exact RatPoly.pow_den_pos _ (by decide) n']
    proof += ['end Spin.Structured.SparsePowers']
    (LEAN / 'SpinCodes/Structured/SparsePowersSound.lean').write_text('\n'.join(proof)+'\n',encoding='utf-8')


def emit_model(model, vector, original_rows, out):
    """Record the finite branch choices and replay the source formulas."""
    e = model.engine
    x = fmpq_poly([0, 1]); z = 1 - x / 6250; beta = x / 12500
    zp = [z**i for i in range(129)]
    total_rows = [fmpq_poly() for _ in range(7)]
    header = ['import SpinCodes.Structured.SparsePolynomialDefs',
              'namespace Spin.Structured.SparsePolynomial.Data', '']
    checks = []
    bridge = LEAN / 'SpinCodes/Structured/SparseBridge'
    bridge.mkdir(exist_ok=True)
    manifest = []
    for j in range(129):
        t = math.comb(128, j); live = fmpq_poly([fmpq(t - e.kernel[j], t)])
        moments = [sum((math.comb(w,h) * math.comb(128-w,j-h) * zp[w+j-2*h]
                        for h in range(max(0,j-128+w),min(w,j)+1)),fmpq_poly()) / t
                   for w in e.levels]
        def pick(ps, largest=False):
            return (max if largest else min)(range(len(ps)),key=lambda i: ps[i](fmpq(1,2)))
        am = pick(moments, True); a = moments[am]
        bounds = [a,live,int(e.caps[j]['cap']) * zp[min(abs(w-j) for w in e.levels)] / t]
        lowmax = 0
        if j in e.low:
            patterns = e.low[j]['patterns']
            ps = [sum((c*zp[w] for w,c in p),fmpq_poly()) / t for p in patterns]
            if ps: lowmax = pick(ps, True)
            bounds.append(ps[lowmax] if ps else fmpq_poly())
            low_text = 'some ' + '[' + ', '.join(lean_pattern(p) for p in patterns) + ']'
            shells = [sorted(e.low[j]['by_weight'].get(w,{}).items()) for w in e.levels]
        else:
            low_text = 'none'; shells = []
        cd = pick(bounds); fd = pick([a,live])
        actions = [fmpq(e.kernel[j],t)*zp[j] + live*zp[j]*vector[1],
                   bounds[cd]/2 + [a,live][fd]/(2*e.m) + a*(vector[1]+fmpq(1,1024))/2]
        cs=[]; fs=[]
        for i,w in enumerate(e.levels):
            m = moments[i]; c = e.spectrum[w]
            sb = [m, fmpq(min(t-e.kernel[j], c*int(e.caps[j]['cap'])),c*t)*zp[abs(w-j)]]
            if j in e.low: sb.append(sum((ct*zp[ow] for ow,ct in shells[i]), fmpq_poly())/(c*t))
            cs.append(pick(sb)); fs.append(pick([m,live]))
            lazy = vector[i+2] if t == e.kernel[j] else vector[1]
            actions.append(sb[cs[-1]]/2 + [m,live][fs[-1]]/(2*e.m) + m*(fmpq(1,1024)+lazy)/2)
        header.extend([f'def weight{j} : WeightData :=',
                       f'  ⟨{e.kernel[j]}, {int(e.caps[j]["cap"])}, {low_text},',
                       '   [' + ', '.join(lean_pattern(p) for p in shells) + '],',
                       f'   ⟨{am}, {lowmax}, {cd}, {fd}, {lean_list(cs)}, {lean_list(fs)}⟩⟩', ''])
        p = t * beta**j * (1-beta)**(128-j)
        module = ['import SpinCodes.Structured.SparseModelData',
                  f'namespace Spin.Structured.SparsePolynomial.Weight{j}',
                  'set_option maxRecDepth 1000000', 'set_option maxHeartbeats 0',
                  f'def prob : RatPoly := {ratpoly(p)}',
                  f'theorem prob_checked : RatPoly.checkEq (probability {j}) prob = true := by decide']
        for i,action in enumerate(actions):
            module += [f'def a{i} : RatPoly := {ratpoly(action)}',
                       f'theorem a{i}_checked : RatPoly.checkEq (action {j} Data.weight{j} {i}) a{i} = true := by decide']
            contribution = p*action
            total_rows[i] += contribution
            # Small first case exercises the exact checker before expanding the cover.
            if j == 0 and i == 0:
                checks = ['import SpinCodes.Structured.SparseModelData',
                          'namespace Spin.Structured.SparsePolynomial.Data',
                          'set_option maxRecDepth 1000000',
                          'set_option maxHeartbeats 0',
                          f'def contribution0_0 : RatPoly := {ratpoly(contribution)}',
                          'theorem contribution0_0_checked :',
                          '    RatPoly.checkEq (contribution 0 weight0 0) contribution0_0 = true := by decide',
                          'end Spin.Structured.SparsePolynomial.Data']
        module.append(f'end Spin.Structured.SparsePolynomial.Weight{j}')
        file = bridge / f'Weight{j}.lean'
        file.write_text('\n'.join(module)+'\n',encoding='utf-8')
        manifest.append(dict(weight=j, path=str(file.relative_to(LEAN)),
                             sha256=hashlib.sha256(file.read_bytes()).hexdigest()))
    assert total_rows == original_rows
    header.append('end Spin.Structured.SparsePolynomial.Data')
    (LEAN / 'SpinCodes/Structured/SparseModelData.lean').write_text('\n'.join(header)+'\n',encoding='utf-8')
    (LEAN / 'SpinCodes/Structured/SparseData/Contribution0.lean').write_text('\n'.join(checks)+'\n',encoding='utf-8')
    (out / 'witness_polynomials.txt').write_text('\n'.join(ratpoly(p) for p in vector)+'\n',encoding='utf-8')
    (out / 'manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
    print('All 903 reconstructed contribution polynomials sum to the seven original rows.')


def main():
    model = verifier.screen.Model()
    vector, rows, residuals, dominance = verifier.build(model)
    saved = json.loads((SOURCE / 'SPARSE_EXACT.json').read_text())
    out = LEAN / 'scripts/sparse_data'
    out.mkdir(exist_ok=True)
    normalizations = []
    for i, poly in enumerate(residuals):
        digest = hashlib.sha256(str(poly).encode()).hexdigest()
        assert digest == saved['row_checks'][i]['coefficients_sha256']
        assert poly[0] == 0
        coefficients = [Fraction(str(c)) for c in poly.coeffs()[1:]]
        denominator = math.lcm(*(c.denominator for c in coefficients))
        numerator = [int(c * denominator) for c in coefficients]
        source = (LEAN / f'SpinCodes/Structured/SparseData/Row{i}.lean').read_text()
        literal = re.search(r'def row\d+\s*:\s*List (?:ℤ|Int)\s*:=\s*\[(.*?)\]', source, re.S)
        assert literal is not None
        existing = [int(s.strip().replace('(', '').replace(')', ''))
                    for s in literal[1].split(',') if s.strip()]
        assert numerator == existing, i
        normalizations.append(dict(row=i, denominator=str(denominator),
                                   numerator_sha256=hashlib.sha256(str(numerator).encode()).hexdigest(),
                                   residual_sha256=digest))
    e = model.engine
    data = dict(kernel=list(map(int, e.kernel)),
                caps=[int(c['cap']) for c in e.caps],
                levels=list(map(int, e.levels)), counts=list(map(int, model.counts)),
                low=e.low, vector=[[str(c) for c in p.coeffs()] for p in vector],
                normalizations=normalizations,
                dominance_methods={method: sum(c['method'] == method for c in dominance)
                                   for method in sorted({c['method'] for c in dominance})},
                status='EXACT_TRANSCRIPTION_ONLY_IDENTITIES_REQUIRE_LEAN')
    (out / 'model.json').write_text(json.dumps(data, indent=2) + '\n')
    print('Recovered seven positive denominators; all seven integer lists match.')
    print('Dominance methods:', data['dominance_methods'])
    emit_model(model, vector, rows, out)
    emit_powers()


if __name__ == '__main__':
    main()


