/-
Fixed-point arithmetic: *soundness against `ℝ`*.

`FixedDefs` gives the kernel-reducible operations; this file proves, once and
for all, that each of them returns a genuine enclosure.  Downstream a box check
is then a `decide` over `Int` comparisons, and the bridge to the real statement
is a single application of the lemmas here.

The shape of every soundness lemma is the same:

    Fix.Mem a x → Fix.Mem b y → Fix.Mem (op a b) (opR x y)

so enclosures compose by `exact`, with no side conditions except positivity of
an explicit divisor.
-/
import Mathlib
import SpinCodes.Numeric.FixedDefs

set_option linter.unusedSectionVars false

namespace Spin.Numeric

/-! ## The denominator -/

lemma scale_pos : (0 : Int) < scale := by unfold scale; norm_num

lemma scaleR_pos : (0 : ℝ) < ((scale : Int) : ℝ) := by exact_mod_cast scale_pos

lemma scaleR_ne : ((scale : Int) : ℝ) ≠ 0 := ne_of_gt scaleR_pos

/-! ## Directed division

`Int`'s division is Euclidean, so for a positive divisor it is floor division. -/

lemma ediv_le_div (a : Int) {d : Int} (hd : 0 < d) :
    ((a / d : Int) : ℝ) ≤ (a : ℝ) / (d : ℝ) := by
  have hid := Int.mul_ediv_add_emod a d
  have hm := Int.emod_nonneg a (ne_of_gt hd)
  have key : (a / d) * d ≤ a := by linarith
  have hdR : (0 : ℝ) < (d : ℝ) := by exact_mod_cast hd
  rw [le_div_iff₀ hdR]
  exact_mod_cast key

lemma fdiv_le (n : Int) {d : Int} (hd : 0 < d) :
    ((fdiv n d : Int) : ℝ) ≤ (n : ℝ) / (d : ℝ) := ediv_le_div n hd

lemma le_cdiv (n : Int) {d : Int} (hd : 0 < d) :
    (n : ℝ) / (d : ℝ) ≤ ((cdiv n d : Int) : ℝ) := by
  have h := ediv_le_div (-n) hd
  unfold cdiv
  push_cast at h ⊢
  rw [neg_div] at h
  linarith

/-! ## Order helpers -/

lemma imin_le_left (a b : Int) : imin a b ≤ a := by unfold imin; split_ifs <;> omega
lemma imin_le_right (a b : Int) : imin a b ≤ b := by unfold imin; split_ifs <;> omega
lemma left_le_imax (a b : Int) : a ≤ imax a b := by unfold imax; split_ifs <;> omega
lemma right_le_imax (a b : Int) : b ≤ imax a b := by unfold imax; split_ifs <;> omega

/-! ## Enclosure -/

namespace Fix

/-- `a` encloses the real `x`. -/
def Mem (a : Fix) (x : ℝ) : Prop :=
  (a.lo : ℝ) / ((scale : Int) : ℝ) ≤ x ∧ x ≤ (a.hi : ℝ) / ((scale : Int) : ℝ)

lemma sc_mem (n : Int) : Mem (sc n) ((n : ℝ) / ((scale : Int) : ℝ)) :=
  ⟨le_refl _, le_refl _⟩

lemma ofInt_mem (n : Int) : Mem (ofInt n) (n : ℝ) := by
  unfold ofInt Mem
  push_cast
  rw [mul_div_assoc, div_self scaleR_ne, mul_one]
  exact ⟨le_refl _, le_refl _⟩

lemma div_le_div_right_of_le {a b c : ℝ} (h : a ≤ b) (hc : 0 < c) : a / c ≤ b / c := by
  gcongr

lemma ofFrac_mem {p q : Int} (hq : 0 < q) :
    Mem (ofFrac p q) ((p : ℝ) / (q : ℝ)) := by
  have hqR : (0 : ℝ) < (q : ℝ) := by exact_mod_cast hq
  have hS := scaleR_pos
  have hlo := fdiv_le (p * scale) hq
  have hhi := le_cdiv (p * scale) hq
  have hrw : (((p * scale : Int) : ℝ) / (q : ℝ)) / ((scale : Int) : ℝ)
      = (p : ℝ) / (q : ℝ) := by
    push_cast
    field_simp
  refine ⟨?_, ?_⟩
  · show ((fdiv (p * scale) q : Int) : ℝ) / ((scale : Int) : ℝ) ≤ (p : ℝ) / (q : ℝ)
    calc ((fdiv (p * scale) q : Int) : ℝ) / ((scale : Int) : ℝ)
        ≤ (((p * scale : Int) : ℝ) / (q : ℝ)) / ((scale : Int) : ℝ) :=
          div_le_div_right_of_le hlo hS
      _ = (p : ℝ) / (q : ℝ) := hrw
  · show (p : ℝ) / (q : ℝ) ≤ ((cdiv (p * scale) q : Int) : ℝ) / ((scale : Int) : ℝ)
    calc (p : ℝ) / (q : ℝ)
        = (((p * scale : Int) : ℝ) / (q : ℝ)) / ((scale : Int) : ℝ) := hrw.symm
      _ ≤ ((cdiv (p * scale) q : Int) : ℝ) / ((scale : Int) : ℝ) :=
          div_le_div_right_of_le hhi hS

/-! ## Arithmetic -/

lemma add_mem {a b : Fix} {x y : ℝ} (hx : Mem a x) (hy : Mem b y) :
    Mem (add a b) (x + y) := by
  obtain ⟨hx1, hx2⟩ := hx
  obtain ⟨hy1, hy2⟩ := hy
  unfold add Mem
  push_cast
  constructor <;> rw [add_div] <;> linarith

lemma neg_mem {a : Fix} {x : ℝ} (hx : Mem a x) : Mem (neg a) (-x) := by
  obtain ⟨hx1, hx2⟩ := hx
  unfold neg Mem
  push_cast
  constructor <;> rw [neg_div] <;> linarith

lemma sub_mem {a b : Fix} {x y : ℝ} (hx : Mem a x) (hy : Mem b y) :
    Mem (sub a b) (x - y) := by
  obtain ⟨hx1, hx2⟩ := hx
  obtain ⟨hy1, hy2⟩ := hy
  unfold sub Mem
  push_cast
  constructor <;> rw [sub_div] <;> linarith

/-! ### Multiplication

The product of two intervals is the hull of the four corner products.  Proving
that needs a case split on the sign of the multiplier; once the sign is fixed,
`nlinarith` closes each branch. -/

lemma mul_le_of_corners {al ah bl bh x y M : ℝ}
    (hx : al ≤ x) (hx' : x ≤ ah) (hy : bl ≤ y) (hy' : y ≤ bh)
    (h1 : al * bl ≤ M) (h2 : al * bh ≤ M) (h3 : ah * bl ≤ M) (h4 : ah * bh ≤ M) :
    x * y ≤ M := by
  by_cases hy0 : (0 : ℝ) ≤ y
  · by_cases hah : (0 : ℝ) ≤ ah
    · nlinarith
    · nlinarith
  · by_cases hal : (0 : ℝ) ≤ al
    · nlinarith
    · nlinarith

lemma le_mul_of_corners {al ah bl bh x y m : ℝ}
    (hx : al ≤ x) (hx' : x ≤ ah) (hy : bl ≤ y) (hy' : y ≤ bh)
    (h1 : m ≤ al * bl) (h2 : m ≤ al * bh) (h3 : m ≤ ah * bl) (h4 : m ≤ ah * bh) :
    m ≤ x * y := by
  have h := mul_le_of_corners (al := al) (ah := ah) (bl := -bh) (bh := -bl)
    (x := x) (y := -y) (M := -m) hx hx' (by linarith) (by linarith)
    (by nlinarith) (by nlinarith) (by nlinarith) (by nlinarith)
  nlinarith [h]

private lemma corner_eq (u v : Int) :
    ((u : ℝ) / ((scale : Int) : ℝ)) * ((v : ℝ) / ((scale : Int) : ℝ))
      = ((u * v : Int) : ℝ) / (((scale : Int) : ℝ) * ((scale : Int) : ℝ)) := by
  push_cast
  ring

lemma mul_mem {a b : Fix} {x y : ℝ} (hx : Mem a x) (hy : Mem b y) :
    Mem (mul a b) (x * y) := by
  obtain ⟨hx1, hx2⟩ := hx
  obtain ⟨hy1, hy2⟩ := hy
  have hSpos : (0 : ℝ) < ((scale : Int) : ℝ) := scaleR_pos
  have hSS : (0 : ℝ) < ((scale : Int) : ℝ) * ((scale : Int) : ℝ) := by positivity
  set p1 := a.lo * b.lo with hp1
  set p2 := a.lo * b.hi with hp2
  set p3 := a.hi * b.lo with hp3
  set p4 := a.hi * b.hi with hp4
  set m : Int := imin (imin p1 p2) (imin p3 p4) with hmdef
  set M : Int := imax (imax p1 p2) (imax p3 p4) with hMdef
  have hm1 : m ≤ p1 := le_trans (imin_le_left _ _) (imin_le_left _ _)
  have hm2 : m ≤ p2 := le_trans (imin_le_left _ _) (imin_le_right _ _)
  have hm3 : m ≤ p3 := le_trans (imin_le_right _ _) (imin_le_left _ _)
  have hm4 : m ≤ p4 := le_trans (imin_le_right _ _) (imin_le_right _ _)
  have hM1 : p1 ≤ M := le_trans (left_le_imax _ _) (left_le_imax _ _)
  have hM2 : p2 ≤ M := le_trans (right_le_imax _ _) (left_le_imax _ _)
  have hM3 : p3 ≤ M := le_trans (left_le_imax _ _) (right_le_imax _ _)
  have hM4 : p4 ≤ M := le_trans (right_le_imax _ _) (right_le_imax _ _)
  have cast_le : ∀ {u v : Int}, u ≤ v →
      ((u : ℝ)) / (((scale : Int) : ℝ) * ((scale : Int) : ℝ))
        ≤ ((v : ℝ)) / (((scale : Int) : ℝ) * ((scale : Int) : ℝ)) := by
    intro u v h
    exact div_le_div_right_of_le (by exact_mod_cast h) hSS
  have hlow : ((m : Int) : ℝ) / (((scale : Int) : ℝ) * ((scale : Int) : ℝ)) ≤ x * y := by
    refine le_mul_of_corners hx1 hx2 hy1 hy2 ?_ ?_ ?_ ?_ <;> rw [corner_eq]
    · exact cast_le hm1
    · exact cast_le hm2
    · exact cast_le hm3
    · exact cast_le hm4
  have hhigh : x * y ≤ ((M : Int) : ℝ) / (((scale : Int) : ℝ) * ((scale : Int) : ℝ)) := by
    refine mul_le_of_corners hx1 hx2 hy1 hy2 ?_ ?_ ?_ ?_ <;> rw [corner_eq]
    · exact cast_le hM1
    · exact cast_le hM2
    · exact cast_le hM3
    · exact cast_le hM4
  refine ⟨?_, ?_⟩
  · show ((fdiv m scale : Int) : ℝ) / ((scale : Int) : ℝ) ≤ x * y
    have h1 : ((fdiv m scale : Int) : ℝ) ≤ ((m : Int) : ℝ) / ((scale : Int) : ℝ) :=
      fdiv_le m scale_pos
    calc ((fdiv m scale : Int) : ℝ) / ((scale : Int) : ℝ)
        ≤ (((m : Int) : ℝ) / ((scale : Int) : ℝ)) / ((scale : Int) : ℝ) :=
          div_le_div_right_of_le h1 hSpos
      _ = ((m : Int) : ℝ) / (((scale : Int) : ℝ) * ((scale : Int) : ℝ)) := div_div _ _ _
      _ ≤ x * y := hlow
  · show x * y ≤ ((cdiv M scale : Int) : ℝ) / ((scale : Int) : ℝ)
    have h1 : ((M : Int) : ℝ) / ((scale : Int) : ℝ) ≤ ((cdiv M scale : Int) : ℝ) :=
      le_cdiv M scale_pos
    calc x * y ≤ ((M : Int) : ℝ) / (((scale : Int) : ℝ) * ((scale : Int) : ℝ)) := hhigh
      _ = (((M : Int) : ℝ) / ((scale : Int) : ℝ)) / ((scale : Int) : ℝ) :=
          (div_div _ _ _).symm
      _ ≤ ((cdiv M scale : Int) : ℝ) / ((scale : Int) : ℝ) :=
          div_le_div_right_of_le h1 hSpos

/-! ### Division by an integer -/

lemma divInt_mem {a : Fix} {x : ℝ} {k : Int} (hx : Mem a x) (hk : 0 < k) :
    Mem (divInt a k) (x / (k : ℝ)) := by
  obtain ⟨hx1, hx2⟩ := hx
  have hkR : (0 : ℝ) < (k : ℝ) := by exact_mod_cast hk
  have hS := scaleR_pos
  refine ⟨?_, ?_⟩
  · show ((fdiv a.lo k : Int) : ℝ) / ((scale : Int) : ℝ) ≤ x / (k : ℝ)
    have h1 : ((fdiv a.lo k : Int) : ℝ) ≤ (a.lo : ℝ) / (k : ℝ) := fdiv_le a.lo hk
    calc ((fdiv a.lo k : Int) : ℝ) / ((scale : Int) : ℝ)
        ≤ ((a.lo : ℝ) / (k : ℝ)) / ((scale : Int) : ℝ) := div_le_div_right_of_le h1 hS
      _ = ((a.lo : ℝ) / ((scale : Int) : ℝ)) / (k : ℝ) := by
          field_simp
      _ ≤ x / (k : ℝ) := div_le_div_right_of_le hx1 hkR
  · show x / (k : ℝ) ≤ ((cdiv a.hi k : Int) : ℝ) / ((scale : Int) : ℝ)
    have h1 : (a.hi : ℝ) / (k : ℝ) ≤ ((cdiv a.hi k : Int) : ℝ) := le_cdiv a.hi hk
    calc x / (k : ℝ) ≤ ((a.hi : ℝ) / ((scale : Int) : ℝ)) / (k : ℝ) :=
          div_le_div_right_of_le hx2 hkR
      _ = ((a.hi : ℝ) / (k : ℝ)) / ((scale : Int) : ℝ) := by
          field_simp
      _ ≤ ((cdiv a.hi k : Int) : ℝ) / ((scale : Int) : ℝ) := div_le_div_right_of_le h1 hS

/-! ### Reciprocal and division -/

lemma pos_hi {a : Fix} {y : ℝ} (hy : Mem a y) (ha : 0 < a.lo) : 0 < a.hi := by
  have haR : (0 : ℝ) < (a.lo : ℝ) := by exact_mod_cast ha
  have hpos : (0 : ℝ) < (a.hi : ℝ) / ((scale : Int) : ℝ) :=
    lt_of_lt_of_le (div_pos haR scaleR_pos) (le_trans hy.1 hy.2)
  have hR : (0 : ℝ) < (a.hi : ℝ) := by
    rcases div_pos_iff.mp hpos with ⟨h1, _⟩ | ⟨_, h2⟩
    · exact h1
    · exact absurd h2 (not_lt.mpr scaleR_pos.le)
  exact_mod_cast hR

lemma pos_of_mem {a : Fix} {y : ℝ} (hy : Mem a y) (ha : 0 < a.lo) : 0 < y := by
  have haR : (0 : ℝ) < (a.lo : ℝ) := by exact_mod_cast ha
  exact lt_of_lt_of_le (div_pos haR scaleR_pos) hy.1

lemma inv_mem {a : Fix} {y : ℝ} (hy : Mem a y) (ha : 0 < a.lo) :
    Mem (inv a) (1 / y) := by
  have hS := scaleR_pos
  have haR : (0 : ℝ) < (a.lo : ℝ) := by exact_mod_cast ha
  have hhi := pos_hi hy ha
  have hhiR : (0 : ℝ) < (a.hi : ℝ) := by exact_mod_cast hhi
  have hypos : 0 < y := pos_of_mem hy ha
  have hSS : ((scale * scale : Int) : ℝ) = ((scale : Int) : ℝ) * ((scale : Int) : ℝ) := by
    push_cast; ring
  refine ⟨?_, ?_⟩
  · show ((fdiv (scale * scale) a.hi : Int) : ℝ) / ((scale : Int) : ℝ) ≤ 1 / y
    have h1 : ((fdiv (scale * scale) a.hi : Int) : ℝ)
        ≤ ((scale * scale : Int) : ℝ) / (a.hi : ℝ) := fdiv_le _ hhi
    have h2 : ((scale * scale : Int) : ℝ) / (a.hi : ℝ) / ((scale : Int) : ℝ)
        = ((scale : Int) : ℝ) / (a.hi : ℝ) := by rw [hSS]; field_simp
    have h3 : ((scale : Int) : ℝ) / (a.hi : ℝ) ≤ 1 / y := by
      rw [div_le_div_iff₀ hhiR hypos, one_mul]
      have := hy.2
      rw [le_div_iff₀ hS] at this
      linarith
    calc ((fdiv (scale * scale) a.hi : Int) : ℝ) / ((scale : Int) : ℝ)
        ≤ ((scale * scale : Int) : ℝ) / (a.hi : ℝ) / ((scale : Int) : ℝ) :=
          div_le_div_right_of_le h1 hS
      _ = ((scale : Int) : ℝ) / (a.hi : ℝ) := h2
      _ ≤ 1 / y := h3
  · show 1 / y ≤ ((cdiv (scale * scale) a.lo : Int) : ℝ) / ((scale : Int) : ℝ)
    have h1 : ((scale * scale : Int) : ℝ) / (a.lo : ℝ)
        ≤ ((cdiv (scale * scale) a.lo : Int) : ℝ) := le_cdiv _ ha
    have h2 : ((scale * scale : Int) : ℝ) / (a.lo : ℝ) / ((scale : Int) : ℝ)
        = ((scale : Int) : ℝ) / (a.lo : ℝ) := by rw [hSS]; field_simp
    have h3 : 1 / y ≤ ((scale : Int) : ℝ) / (a.lo : ℝ) := by
      rw [div_le_div_iff₀ hypos haR, one_mul]
      have := hy.1
      rw [div_le_iff₀ hS] at this
      linarith
    calc 1 / y ≤ ((scale : Int) : ℝ) / (a.lo : ℝ) := h3
      _ = ((scale * scale : Int) : ℝ) / (a.lo : ℝ) / ((scale : Int) : ℝ) := h2.symm
      _ ≤ ((cdiv (scale * scale) a.lo : Int) : ℝ) / ((scale : Int) : ℝ) :=
          div_le_div_right_of_le h1 hS

lemma div_mem {a b : Fix} {x y : ℝ} (hx : Mem a x) (hy : Mem b y) (hb : 0 < b.lo) :
    Mem (div a b) (x / y) := by
  have h := mul_mem hx (inv_mem hy hb)
  rw [show x * (1 / y) = x / y by ring] at h
  exact h

/-! ### Powers -/

lemma pow_mem {a : Fix} {x : ℝ} (hx : Mem a x) : ∀ n : Nat, Mem (pow a n) (x ^ n)
  | 0 => by simpa [pow] using ofInt_mem 1
  | n + 1 => by
      have h := mul_mem (pow_mem hx n) hx
      simpa [pow, pow_succ] using h

/-! ## Reading a verdict off an enclosure -/

lemma neg_of_isNeg {a : Fix} {x : ℝ} (hx : Mem a x) (h : isNeg a = true) : x < 0 := by
  have hhi : a.hi < 0 := by simpa [isNeg] using h
  have hhiR : ((a.hi : Int) : ℝ) < 0 := by exact_mod_cast hhi
  have hlt : (a.hi : ℝ) / ((scale : Int) : ℝ) < 0 :=
    div_neg_of_neg_of_pos hhiR scaleR_pos
  linarith [hx.2]

lemma nonneg_of_isNonneg {a : Fix} {x : ℝ} (hx : Mem a x) (h : isNonneg a = true) :
    0 ≤ x := by
  have hlo : 0 ≤ a.lo := by simpa [isNonneg] using h
  have hloR : (0 : ℝ) ≤ ((a.lo : Int) : ℝ) := by exact_mod_cast hlo
  have hge : (0 : ℝ) ≤ (a.lo : ℝ) / ((scale : Int) : ℝ) :=
    div_nonneg hloR scaleR_pos.le
  linarith [hx.1]

end Fix

end Spin.Numeric
