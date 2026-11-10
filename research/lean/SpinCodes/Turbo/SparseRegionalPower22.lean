import SpinCodes.Turbo.SparseRegionalPowerDefs22
import SpinCodes.Turbo.SparseFixedSound22
import SpinCodes.Turbo.SparseActiveKernel22

/-! Soundness of the regional binary-powering certificate.

`regionalCountDP22 B n d` is the unnormalized ordered occupancy coefficient
inside `sparseRegionalDP22`. Each stored stage is checked once against the
exact convolution `orderedMatrixDP_add`; an induction over the trace (not a
kernel reduction of a recursive power) gives the R-step bound, and
`normalize_sound` divides by choose(16R,d). -/
noncomputable section
namespace Spin.Turbo
open Finset Matrix SparseFixedSound22

def regionalCountDP22 (B : Fin 17 → Matrix (Fin 11) (Fin 11) ℝ) (n d : ℕ) :
    Matrix (Fin 11) (Fin 11) ℝ :=
  orderedMatrixDP Fin.val (fun j : Fin 17 => (Nat.choose 16 j.val : ℝ) • B j) n d

namespace SparseRegionalPower22

theorem count_nonneg (B : Fin 17 → Matrix (Fin 11) (Fin 11) ℝ)
    (hB : ∀ a i j, 0 ≤ B a i j) (n d : ℕ) (i j : Fin 11) :
    0 ≤ regionalCountDP22 B n d i j := by
  unfold regionalCountDP22
  rw [← orderedMatrixCoefficient_eq_dp]
  unfold orderedMatrixCoefficient
  rw [Matrix.sum_apply]
  apply Finset.sum_nonneg
  intro xs _
  split_ifs
  · exact listNonneg _ (fun M hM => by
      obtain ⟨step, rfl⟩ := List.mem_ofFn.mp hM
      intro a b
      exact mul_nonneg (Nat.cast_nonneg _) (hB _ _ _)) i j
  · exact le_rfl

theorem count_one (B : Fin 17 → Matrix (Fin 11) (Fin 11) ℝ) (d : ℕ) (hd : d < 17) :
    regionalCountDP22 B 1 d = (Nat.choose 16 d : ℝ) • B ⟨d, hd⟩ := by
  unfold regionalCountDP22
  simp only [orderedMatrixDP]
  rw [Finset.sum_eq_single (⟨d, hd⟩ : Fin 17)]
  · simp
  · intro j _ hj
    have hjd : j.val ≠ d := fun h => hj (Fin.ext h)
    split_ifs with h1 h2
    · exact absurd (by omega : j.val = d) hjd
    · simp
    · rfl
  · intro h
    exact absurd (Finset.mem_univ _) h

theorem count_add_bound (B : Fin 17 → Matrix (Fin 11) (Fin 11) ℝ)
    (hB : ∀ a i j, 0 ≤ B a i j) (scale : ℤ) (hs : 0 < scale)
    (a b u : Array ℤ) (m n : ℕ)
    (hcheck : sparseIntPolyMulUpperCheck22 scale a b u = true)
    (ha : ∀ d, d ≤ 2 → ∀ i j,
      regionalCountDP22 B m d i j ≤ (sparseIntPolyMatrix22 a d i j : ℝ) / scale)
    (hb : ∀ d, d ≤ 2 → ∀ i j,
      regionalCountDP22 B n d i j ≤ (sparseIntPolyMatrix22 b d i j : ℝ) / scale) :
    ∀ d, d ≤ 2 → ∀ i j,
      regionalCountDP22 B (m + n) d i j ≤ (sparseIntPolyMatrix22 u d i j : ℝ) / scale := by
  intro d hd i j
  unfold regionalCountDP22
  rw [orderedMatrixDP_add]
  exact polyMul_sound scale hs a b u hcheck (regionalCountDP22 B m) (regionalCountDP22 B n)
    (count_nonneg B hB m) (count_nonneg B hB n) ha hb d hd i j

theorem chain_sound (B : Fin 17 → Matrix (Fin 11) (Fin 11) ℝ)
    (hB : ∀ a i j, 0 ≤ B a i j) (scale : ℤ) (hs : 0 < scale) (base : Array ℤ)
    (hbase : ∀ d, d ≤ 2 → ∀ i j,
      regionalCountDP22 B 1 d i j ≤ (sparseIntPolyMatrix22 base d i j : ℝ) / scale) :
    ∀ (trace : List (Bool × Array ℤ × Array ℤ)) (prev : Array ℤ) (e : ℕ),
      (∀ d, d ≤ 2 → ∀ i j,
        regionalCountDP22 B e d i j ≤ (sparseIntPolyMatrix22 prev d i j : ℝ) / scale) →
      sparseRegionalPowerCheck22 scale base prev trace = true →
      ∀ d, d ≤ 2 → ∀ i j,
        regionalCountDP22 B (trace.foldl sparseRegionalStep22 e) d i j ≤
          (sparseIntPolyMatrix22 (sparseRegionalLast22 prev trace) d i j : ℝ) / scale
  | [], _, _, h, _ => h
  | (bit, sq, next) :: tail, prev, e, h, hc => by
    have hc' : ((sparseIntPolyMulUpperCheck22 scale prev prev sq &&
        (if bit then sparseIntPolyMulUpperCheck22 scale sq base next else true)) &&
        sparseRegionalPowerCheck22 scale base (if bit then next else sq) tail) = true := hc
    rw [Bool.and_eq_true, Bool.and_eq_true] at hc'
    obtain ⟨⟨hsq, hbit⟩, hrest⟩ := hc'
    have hsqB := count_add_bound B hB scale hs prev prev sq e e hsq h h
    cases bit
    · have hrest' : sparseRegionalPowerCheck22 scale base sq tail = true := hrest
      have he : sparseRegionalStep22 e (false, sq, next) = e + e := by
        unfold sparseRegionalStep22; simp <;> omega
      have hsq' : ∀ d, d ≤ 2 → ∀ i j,
          regionalCountDP22 B (sparseRegionalStep22 e (false, sq, next)) d i j ≤
            (sparseIntPolyMatrix22 sq d i j : ℝ) / scale := by
        rw [he]; exact hsqB
      exact chain_sound B hB scale hs base hbase tail sq _ hsq' hrest'
    · have hbit' : sparseIntPolyMulUpperCheck22 scale sq base next = true := hbit
      have hrest' : sparseRegionalPowerCheck22 scale base next tail = true := hrest
      have hnext := count_add_bound B hB scale hs sq base next (e + e) 1 hbit' hsqB hbase
      have he : sparseRegionalStep22 e (true, sq, next) = e + e + 1 := by
        unfold sparseRegionalStep22; simp <;> omega
      have hnext' : ∀ d, d ≤ 2 → ∀ i j,
          regionalCountDP22 B (sparseRegionalStep22 e (true, sq, next)) d i j ≤
            (sparseIntPolyMatrix22 next d i j : ℝ) / scale := by
        rw [he]; exact hnext
      exact chain_sound B hB scale hs base hbase tail next _ hnext' hrest'

/-- Base-stage bound from bounds on the three retained physical rows. -/
theorem base_of_rows (B : Fin 17 → Matrix (Fin 11) (Fin 11) ℝ) (scale : ℤ)
    (base : Array ℤ)
    (h : ∀ d (hd : d ≤ 2), ∀ i j, (Nat.choose 16 d : ℝ) * B ⟨d, by omega⟩ i j ≤
      (sparseIntPolyMatrix22 base d i j : ℝ) / scale) :
    ∀ d, d ≤ 2 → ∀ i j,
      regionalCountDP22 B 1 d i j ≤ (sparseIntPolyMatrix22 base d i j : ℝ) / scale := by
  intro d hd i j
  rw [count_one B d (by omega), Matrix.smul_apply, smul_eq_mul]
  exact h d hd i j

/-- A passing regional certificate bounds every retained normalized
occupancy kernel entry. -/
theorem certificate_sound (B : Fin 17 → Matrix (Fin 11) (Fin 11) ℝ)
    (hB : ∀ a i j, 0 ≤ B a i j) (R : ℕ) (scale : ℤ) (base result : Array ℤ)
    (trace : List (Bool × Array ℤ × Array ℤ))
    (hcheck : sparseRegionalCertificateCheck22 R scale base trace result = true)
    (hbase : ∀ d, d ≤ 2 → ∀ i j,
      regionalCountDP22 B 1 d i j ≤ (sparseIntPolyMatrix22 base d i j : ℝ) / scale) :
    ∀ d, d ≤ 2 → ∀ i j,
      sparseRegionalDP22 B R d i j ≤ (sparseIntPolyMatrix22 result d i j : ℝ) / scale := by
  simp only [sparseRegionalCertificateCheck22, Bool.and_eq_true, decide_eq_true_eq] at hcheck
  obtain ⟨⟨⟨hs, hR, hexp⟩, hchain⟩, hnorm⟩ := hcheck
  have hc := chain_sound B hB scale hs base hbase trace base 1 hbase hchain
  rw [hexp] at hc
  intro d hd i j
  exact normalize_sound R hR scale hs _ result hnorm (regionalCountDP22 B R) hc d hd i j

end SparseRegionalPower22

open SparseRegionalPower22 in
#print axioms chain_sound
open SparseRegionalPower22 in
#print axioms certificate_sound
end Spin.Turbo
