from pathlib import Path
import re

ROOT = Path('.')
T = ROOT / 'SpinCodes/Turbo'
GEN = T / 'Generated/K18Tail'
NL = chr(10)

# 1. generic assembly lemmas (Mathlib-level file)
p = T / 'PairedFracPacked22.lean'
s = p.read_text(encoding='utf-8')
if 'rpowCheck_of' not in s:
    lemmas = '''
theorem rpowCheck_of (den sy : Int) (R : List (List (List Int)))
    (Ys : List (List (List (List Int)))) (h0 : 0 < den)
    (h : ∀ a, a < 17 → (sparseMatShapeL22 (R.getD a []) && (List.range 17).all fun b =>
      sparseRpow45CheckL22 (den * den) sy
        (sparseMatMulL22 (R.getD a []) (R.getD b [])) ((Ys.getD a []).getD b [])) = true) :
    pairedRpowCheckL22 den sy R Ys = true := by
  unfold pairedRpowCheckL22
  rw [decide_eq_true h0, Bool.true_and, List.all_eq_true]
  intro a ha
  exact h a (List.mem_range.mp ha)

theorem baseGo_append (Ys : List (List (List (List Int)))) :
    ∀ (l1 l2 : List (List (List Int))) (d : Nat),
      pairedBaseGoL22 Ys d (l1 ++ l2) =
        (pairedBaseGoL22 Ys d l1 && pairedBaseGoL22 Ys (d + l1.length) l2)
  | [], l2, d => by simp [pairedBaseGoL22]
  | m :: ms, l2, d => by
    rw [List.cons_append, pairedBaseGoL22, pairedBaseGoL22, baseGo_append Ys ms l2 (d + 1),
      List.length_cons, show d + (ms.length + 1) = d + 1 + ms.length by omega, Bool.and_assoc]

'''
    s = s.replace('end PairedFracChain22\n', lemmas + 'end PairedFracChain22\n', 1)
    p.write_text(s, encoding='utf-8', newline=NL)

ns = 'Spin.Turbo.Hi.ProfileFamilyK18Tail722.Fr'
pre = 'SpinCodes.Turbo.Generated.K18Tail.'
chain = GEN / 'Tail7FrChain.lean'
c = chain.read_text(encoding='utf-8')

# 2. Ys module
ys_line = re.search(r'^def Ys : .*$', c, re.M).group(0)
(GEN / 'Tail7FrYs.lean').write_text(NL.join(
    [f'import {pre}Tail7FrY{i}' for i in range(17)] +
    ['', f'namespace {ns}', ys_line, f'end {ns}', '']), encoding='utf-8', newline=NL)

DEN = '(16 ^ 16 * 2 ^ 384)'
# 3. rpow per a
for a in range(17):
    body = (f'(sparseMatShapeL22 (R.getD {a} []) && (List.range 17).all fun b =>' + NL +
            f'      sparseRpow45CheckL22 ({DEN} * {DEN}) (2 ^ 704)' + NL +
            f'        (sparseMatMulL22 (R.getD {a} []) (R.getD b [])) ((Ys.getD {a} []).getD b [])) = true')
    (GEN / f'Tail7FrRpow{a}.lean').write_text(NL.join([
        f'import {pre}Tail7FrR', f'import {pre}Tail7FrYs', 'import SpinCodes.Turbo.PairedFracChain22', '',
        f'namespace {ns}', 'open Spin.Turbo', '',
        'set_option maxHeartbeats 0 in',
        f'theorem rpow_{a} :', f'    {body} := by', '  decide +kernel', '',
        f'end {ns}', '']), encoding='utf-8', newline=NL)

# 4. base chunks
CH = 12
nparts = 144
for ci in range(nparts // CH):
    parts = ', '.join(f'base_part{t}' for t in range(ci * CH, (ci + 1) * CH))
    (GEN / f'Tail7FrBaseC{ci}.lean').write_text(NL.join([
        f'import {pre}Tail7FrBase', f'import {pre}Tail7FrYs', 'import SpinCodes.Turbo.PairedFracChain22', '',
        f'namespace {ns}', 'open Spin.Turbo', '',
        f'def baseC{ci} : List (List (List Int)) := [{parts}]', '',
        'set_option maxHeartbeats 0 in',
        f'theorem baseC{ci}_ok : pairedBaseGoL22 Ys {ci * CH} baseC{ci} = true := by', '  decide +kernel', '',
        f'end {ns}', '']), encoding='utf-8', newline=NL)

# 5. rewrite FrChain
c = c.replace(ys_line + NL, '')
old_rpow = re.search(r'set_option maxHeartbeats 0 in\ntheorem rpowChecked : .*\n', c).group(0)
old_base = re.search(r'set_option maxHeartbeats 0 in\ntheorem baseChecked : .*\n', c).group(0)
arms = NL.join(f'    | {a}, _ => rpow_{a}' for a in range(17))
new_rpow = (f'theorem rpowChecked : pairedRpowCheckL22 {DEN} (2 ^ 704) R Ys = true :=' + NL +
            f'  PairedFracChain22.rpowCheck_of _ _ _ _ (by decide) fun a ha =>' + NL +
            '    match a, ha with' + NL + arms + NL + '    | _ + 17, h => absurd h (by omega)' + NL)
nch = nparts // CH
appended = ' ++ '.join(f'baseC{ci}' for ci in range(nch))
new_base = (f'theorem baseChecked : pairedBaseCheckDL22 143 Ys base = true := by' + NL +
            f'  have hb : base = {appended} := rfl' + NL +
            '  simp only [pairedBaseCheckDL22, Bool.and_eq_true, decide_eq_true_eq]' + NL +
            '  refine ⟨rfl, ?_⟩' + NL +
            '  rw [hb]' + NL +
            '  simp only [List.append_assoc, PairedFracChain22.baseGo_append, Bool.and_eq_true]' + NL +
            '  exact ⟨' + ', '.join(f'baseC{ci}_ok' for ci in range(nch)) + '⟩' + NL)
c = c.replace(old_rpow, new_rpow).replace(old_base, new_base)
imports = [f'import {pre}Tail7FrYs', 'import SpinCodes.Turbo.PairedFracPacked22'] + \
    [f'import {pre}Tail7FrRpow{a}' for a in range(17)] + [f'import {pre}Tail7FrBaseC{ci}' for ci in range(nch)]
c = NL.join(imports) + NL + c
chain.write_text(c, encoding='utf-8', newline=NL)
print('ok')
