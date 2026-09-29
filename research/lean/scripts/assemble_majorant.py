"""Pin the paper's exact 39 supports and assemble their active intervals."""
from fractions import Fraction as Q
import json
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'SpinCodes/Majorant'
data=json.loads((ROOT.parent/'workstreams/inner_design/imt_asymptotic/d11/OUTER_REFINED.json').read_text())
segments=[(Q(s['slope']),Q(s['intercept']),*map(Q,s['weight'])) for s in data['left_segments']]
central=Q(data['central_constant_upper'])
segments.append((Q(0),central,segments[-1][3],Q(1,2)))
left=[(s,i) for s,i,_,_ in segments]
all_lines=left+[(-s,s+i) for s,i in reversed(left[:-1])]


def rat(q):
    return f'({q.numerator} / {q.denominator})'


def pair(p):
    return '('+rat(p[0])+', '+rat(p[1])+')'


source='''import SpinCodes.Structured.Majorant
set_option maxRecDepth 100000
set_option maxHeartbeats 0
namespace Spin.Majorant
open Spin.Structured
'''
source+='def supports : Finset (ℚ × ℚ) := {'+', '.join(map(pair,all_lines))+'}\n'
source+='''def refined : ConcaveMajorant := ⟨supports, by exact ⟨_, Finset.mem_insert_self _ _⟩⟩
def reflect (p : ℚ × ℚ) : ℚ × ℚ := (-p.1, p.1+p.2)
'''
for j,line in enumerate(all_lines):
    proof='True.intro' if j==38 else '(Or.inl True.intro)'
    for _ in range(j): proof=f'(Or.inr {proof})'
    source+=f'''lemma member{j} : {pair(line)} ∈ supports := by
  simp only [supports, Finset.mem_insert, Finset.mem_singleton]
  exact {proof}
'''
source+='''theorem reflect_mem {p : ℚ × ℚ} (hp : p ∈ supports) : reflect p ∈ supports := by
  simp only [supports, Finset.mem_insert, Finset.mem_singleton] at hp
  rcases hp with '''+' | '.join(['rfl']*39)+'\n'
for j,line in enumerate(all_lines):
    k=all_lines.index((-line[0],line[0]+line[1]))
    source+=f'''  · have he : reflect {pair(line)} = {pair(all_lines[k])} := by norm_num [reflect]
    rw [he]
    exact member{k}
'''
source+='''
theorem refined_reflect (w : ℝ) : refined.toFun (1-w) = refined.toFun w := by
  apply le_antisymm
  · apply Finset.le_inf'
    intro p hp
    have h := refined.toFun_le_support (reflect_mem hp) (1-w)
    convert h using 1 <;> simp only [reflect, Rat.cast_neg, Rat.cast_add] <;> ring
  · apply Finset.le_inf'
    intro p hp
    have h := refined.toFun_le_support (reflect_mem hp) w
    convert h using 1 <;> simp only [reflect, Rat.cast_neg, Rat.cast_add] <;> ring

lemma affine_le_between {s i t j l u x : ℝ} (hx : x ∈ Set.Icc l u)
    (hl : s*l+i ≤ t*l+j) (hu : s*u+i ≤ t*u+j) : s*x+i ≤ t*x+j := by
  rcases le_total s t with h | h
  · nlinarith [mul_nonneg (sub_nonneg.mpr h) (sub_nonneg.mpr hx.1)]
  · nlinarith [mul_nonneg (sub_nonneg.mpr h) (sub_nonneg.mpr hx.2)]

lemma active_le_refined {s i l u : ℚ} {w : ℝ} (hw : w ∈ Set.Icc (l : ℝ) (u : ℝ))
    (h : ∀ p ∈ supports, s*l+i ≤ p.1*l+p.2 ∧ s*u+i ≤ p.1*u+p.2) :
    (s : ℝ)*w+(i : ℝ) ≤ refined.toFun w := by
  apply Finset.le_inf'
  intro p hp
  obtain ⟨hl, hu⟩ := h p hp
  apply affine_le_between hw
  · exact_mod_cast hl
  · exact_mod_cast hu
'''
for j,(s,i,l,u) in enumerate(segments):
    source+=f'''theorem active{j} : ∀ p ∈ supports,
    ({rat(s)}:ℚ)*{rat(l)}+{rat(i)} ≤ p.1*{rat(l)}+p.2 ∧
    ({rat(s)}:ℚ)*{rat(u)}+{rat(i)} ≤ p.1*{rat(u)}+p.2 := by
  simp only [supports, Finset.forall_mem_insert, Finset.mem_singleton, forall_eq]
  norm_num
'''
for j,(s,i,l,u) in enumerate(segments):
    source+=f'''theorem active_real{j} {{w : ℝ}} (hw : w ∈ Set.Icc ({rat(l)}:ℝ) {rat(u)}) :
    ({rat(s)}:ℝ)*w+{rat(i)} ≤ refined.toFun w := by
  have hwQ : w ∈ Set.Icc (({rat(l)}:ℚ):ℝ) (({rat(u)}:ℚ):ℝ) := by
    norm_num
    exact hw
  have h := active_le_refined hwQ active{j}
  norm_num at h
  linarith
'''
source+='end Spin.Majorant\n'
(OUT/'RefinedData.lean').write_text(source,encoding='utf-8')

source='\n'.join(f'import SpinCodes.Majorant.Segment{j}' for j in range(19))
source+='''
import SpinCodes.Majorant.RefinedData
import SpinCodes.Numeric.BACentral
namespace Spin.Majorant
open Spin.Numeric Spin.Numeric.Fix
'''
for j,(s,i,l,u) in enumerate(segments[:-1]):
    d=json.loads((ROOT/f'scripts/majorant_data/segment{j}.json').read_text())
    source+=f'''theorem segment{j}_real {{a b w : ℝ}}
    (ha : a ∈ Set.Icc (0:ℝ) 1) (hb : b ∈ Set.Icc (0:ℝ) 1)
    (hw : w ∈ Set.Icc ({rat(l)}:ℝ) {rat(u)}) :
    MajorantClaim {rat(s)} {rat(i)} a b w := by
  apply S{j}.cover (s := {rat(s)}) (i := {rat(i)})
  · norm_num [Mem, scale]
  · norm_num [Mem, scale]
  · simpa [Mem, scale] using ha
  · simpa [Mem, scale] using hb
  · constructor <;> norm_num [scale] <;> linarith [hw.1, hw.2]
'''
source+='''
theorem majorant_left {a b w : ℝ}
    (ha : a ∈ Set.Icc (0:ℝ) 1) (hb : b ∈ Set.Icc (0:ℝ) 1)
    (hw : w ∈ Set.Icc ((13:ℝ)/125) (1/2))
    (hf : a/2 ≤ b ∧ b ≤ 1-a/2 ∧ b/2 ≤ w ∧ w ≤ 1-b/2 ∧ 0 ≤ a ∧ a ≤ 1) :
    gBA a + piBA a b + piBA b w ≤ refined.toFun w := by
'''
for j,(s,i,l,u) in enumerate(segments[:-1]):
    source+=f'''  by_cases h{j} : w ≤ ({rat(u)}:ℝ)
  · have hwj : w ∈ Set.Icc ({rat(l)}:ℝ) {rat(u)} := by constructor <;> linarith [hw.1]
    exact le_trans (segment{j}_real ha hb hwj hf) (active_real{j} hwj)
'''
s,i,l,u=segments[-1]
source+=f'''  have hwj : w ∈ Set.Icc ({rat(l)}:ℝ) {rat(u)} := by constructor <;> linarith [hw.2]
  have hc := le_trans (ba_central ha hf.1 hf.2.1 hf.2.2.1 hf.2.2.2.1) central_rounding
  have hactive := active_real19 hwj
  norm_num at hactive hc
  exact le_trans hc hactive

theorem majorant_pointwise {{a b w : ℝ}}
    (ha : a ∈ Set.Icc (0:ℝ) 1) (hb : b ∈ Set.Icc (0:ℝ) 1)
    (hw : w ∈ Set.Icc ((13:ℝ)/125) (112/125))
    (hf : a/2 ≤ b ∧ b ≤ 1-a/2 ∧ b/2 ≤ w ∧ w ≤ 1-b/2) :
    gBA a + piBA a b + piBA b w ≤ refined.toFun w := by
  by_cases h : w ≤ (1:ℝ)/2
  · exact majorant_left ha hb ⟨hw.1, h⟩ ⟨hf.1, hf.2.1, hf.2.2.1, hf.2.2.2, ha.1, ha.2⟩
  · have hwr : 1-w ∈ Set.Icc ((13:ℝ)/125) (1/2) := by constructor <;> linarith [hw.2]
    have hfr : a/2 ≤ b ∧ b ≤ 1-a/2 ∧ b/2 ≤ 1-w ∧ 1-w ≤ 1-b/2 ∧ 0 ≤ a ∧ a ≤ 1 := by
      exact ⟨hf.1, hf.2.1, by linarith [hf.2.2.2], by linarith [hf.2.2.1], ha.1, ha.2⟩
    have hr := majorant_left ha hb hwr hfr
    rwa [piBA_reflect, refined_reflect] at hr

/-- The paper's variational exponent. Restricting to feasible pairs implements
the paper's value -∞ for infeasible paths. -/
noncomputable def baExponent (w : ℝ) : ℝ := sSup
  ((fun p : ℝ × ℝ => gBA p.1 + piBA p.1 p.2 + piBA p.2 w) ''
    {{p | p.1 ∈ Set.Icc (0:ℝ) 1 ∧ p.2 ∈ Set.Icc (0:ℝ) 1 ∧
      p.1/2 ≤ p.2 ∧ p.2 ≤ 1-p.1/2 ∧ p.2/2 ≤ w ∧ w ≤ 1-p.2/2}})

/-- Equation ba-spectrum-majorant, with every affine support checked by the
kernel, including all feasible boundary points. -/
theorem baExponent_le_refined {{w : ℝ}} (hw : w ∈ Set.Icc ((13:ℝ)/125) (112/125)) :
    baExponent w ≤ refined.toFun w := by
  apply csSup_le
  · refine ⟨_, ⟨(0,0), ?_, rfl⟩⟩
    norm_num
    constructor <;> linarith [hw.1, hw.2]
  · rintro y ⟨⟨a,b⟩, ⟨ha,hb,hf⟩, rfl⟩
    exact majorant_pointwise ha hb hw hf

end Spin.Majorant
'''
(OUT/'Refined.lean').write_text(source,encoding='utf-8')
