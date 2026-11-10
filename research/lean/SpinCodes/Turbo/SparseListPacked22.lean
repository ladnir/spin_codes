import SpinCodes.Turbo.SparseListPolyD22
import SpinCodes.Turbo.SparseListPackedDefs22

/-!
# Packed (Kronecker-substitution) degree-D polynomial-matrix checks

`SparseListPolyD22` checks a truncated product of degree-`D` polynomial
matrices by `(D+1)(D+2)/2` interpreted 11×11 matrix products. Kernel
evaluation of that is far too slow for `D` in the hundreds.

Here each entry `(i,k)` of a polynomial matrix `P = [P_0, …, P_D]` is packed
into the natural number `Σ_t P_t[i][k] · 2^(w t)` (Horner, shifts). If all
entries are in `[0, 2^v)` and `11 (2D+1) 4^v < 2^w`, the base-`2^w` digits of
`Σ_k pack(A_ik) · pack(B_kj)` are exactly the convolution coefficients
`Σ_k Σ_{x+y=e} A_x[i][k] B_y[k][j]` (no carries). One check is then 121·11
GMP multiplications plus a linear digit walk against the claimed stage `u`.

Soundness yields the same `PolyBoundD` interface as `SparseListPolyD22`, so
the remaining tail pipeline is unchanged.
-/

/-! ## Soundness -/
noncomputable section
namespace Spin.Turbo
open Finset Matrix Polynomial SparseFixedSound22 SparseListSound22 SparseRegionalPower22

namespace SparseListPacked22

def polyOf (α : ℕ → ℕ) (n : ℕ) : ℕ[X] := ∑ t ∈ range n, C (α t) * X ^ t

theorem coeff_polyOf (α : ℕ → ℕ) (n e : ℕ) :
    (polyOf α n).coeff e = if e < n then α e else 0 := by
  simp only [polyOf, finset_sum_coeff, coeff_C_mul_X_pow, Finset.sum_ite_eq, Finset.mem_range]

theorem polyOf_succ (α : ℕ → ℕ) (n : ℕ) :
    polyOf α (n + 1) = C (α 0) + X * polyOf (fun t => α (t + 1)) n := by
  ext e
  rw [coeff_add, coeff_polyOf]
  cases e with
  | zero => simp [coeff_C]
  | succ e =>
    rw [coeff_X_mul, coeff_polyOf, coeff_C]
    simp only [Nat.succ_ne_zero, if_false, zero_add]
    by_cases h : e < n
    · rw [if_pos h, if_pos (by omega)]
    · rw [if_neg h, if_neg (by omega)]

theorem lpoly_succ (M : List (List Int)) (P : List (List (List Int))) (t i j : ℕ) :
    sparseLPoly22 (M :: P) (t + 1) i j = sparseLPoly22 P t i j := by
  simp [sparseLPoly22]

theorem lpoly_zero (M : List (List Int)) (P : List (List (List Int))) (i j : ℕ) :
    sparseLPoly22 (M :: P) 0 i j = sparseLMat22 M i j := by
  simp [sparseLPoly22]

theorem pkEntry_eq (w i j : ℕ) : ∀ P : List (List (List Int)),
    pkEntryL22 w i j P =
      (polyOf (fun t => (sparseLPoly22 P t i j).toNat) P.length).eval (2 ^ w)
  | [] => by simp [pkEntryL22, polyOf]
  | M :: P => by
    rw [pkEntryL22, pkEntry_eq w i j P, List.length_cons, polyOf_succ]
    simp only [lpoly_succ, lpoly_zero, eval_add, eval_C, eval_mul, eval_X, Nat.shiftLeft_eq,
      sparseLMat22]
    ring

theorem getD_mem_or {α : Type*} (l : List α) (n : ℕ) (d : α) : l.getD n d ∈ l ∨ l.getD n d = d := by
  by_cases h : n < l.length
  · left; rw [gd l d h]; exact List.getElem_mem h
  · right; exact SparseListPolyD22.getD_ge d l n (by omega)

theorem bound_sound (v : ℕ) (P : List (List (List Int))) (h : pkBoundL22 v P = true)
    (t i k : ℕ) : 0 ≤ sparseLPoly22 P t i k ∧ (sparseLPoly22 P t i k).toNat < 2 ^ v := by
  simp only [pkBoundL22, List.all_eq_true, Bool.and_eq_true, decide_eq_true_eq] at h
  unfold sparseLPoly22 sparseLMat22
  rcases getD_mem_or P t [] with hM | hM
  · rcases getD_mem_or (P.getD t []) i [] with hr | hr
    · rcases getD_mem_or ((P.getD t []).getD i []) k 0 with hx | hx
      · exact h _ hM _ hr _ hx
      · rw [hx]; exact ⟨le_rfl, by simp⟩
    · rw [hr]; exact ⟨le_rfl, by simp⟩
  · rw [hM]; exact ⟨le_rfl, by simp⟩

theorem pkMat_get (w : ℕ) (P : List (List (List Int))) (i k : ℕ) (hi : i < 11) (hk : k < 11) :
    ((pkMatL22 w P).getD i []).getD k 0 = pkEntryL22 w i k P := by
  simp [pkMatL22, List.getD_eq_getElem?_getD, List.getElem?_map, List.getElem?_range, hi, hk]

theorem pkDot_eq (A B : List (List Nat)) (i j : ℕ) : ∀ n,
    pkDotL22 A B i j n = ∑ k ∈ range n, (A.getD i []).getD k 0 * (B.getD k []).getD j 0
  | 0 => by simp [pkDotL22]
  | n + 1 => by rw [pkDotL22, pkDot_eq A B i j n, Finset.sum_range_succ]

/-- Base-`W` digit extraction from a polynomial with coefficients `< W`. -/
theorem digit_eval (W : ℕ) (hW : 0 < W) : ∀ (d : ℕ) (p : ℕ[X]), (∀ e, p.coeff e < W) →
    p.eval W / W ^ d % W = p.coeff d := by
  have hsplit : ∀ p : ℕ[X], p.eval W = p.coeff 0 + W * (divX p).eval W := by
    intro p
    conv_lhs => rw [← divX_mul_X_add p]
    simp only [eval_add, eval_mul, eval_X, eval_C]
    ring
  intro d
  induction d with
  | zero =>
    intro p hp
    rw [pow_zero, Nat.div_one, hsplit p, Nat.add_mul_mod_self_left, Nat.mod_eq_of_lt (hp 0)]
  | succ d ih =>
    intro p hp
    have hq : ∀ e, (divX p).coeff e < W := fun e => by rw [coeff_divX]; exact hp (e + 1)
    rw [pow_succ', ← Nat.div_div_eq_div_mul, hsplit p, Nat.add_mul_div_left _ _ hW,
      Nat.div_eq_of_lt (hp 0), zero_add, ih (divX p) hq, coeff_divX]

theorem digits_sound (w : ℕ) (s : ℤ) (i j : ℕ) : ∀ (u : List (List (List Int))) (N : ℕ),
    pkDigitsLeL22 w s i j N u = true → ∀ e, e < u.length →
      ((N >>> (w * e) % 2 ^ w : ℕ) : ℤ) ≤ s * sparseLMat22 (u.getD e []) i j
  | [], _, _, e, he => by simp at he
  | M :: us, N, h, e, he => by
    simp only [pkDigitsLeL22, Bool.and_eq_true, decide_eq_true_eq] at h
    cases e with
    | zero => simpa [sparseLMat22] using h.1
    | succ e =>
      have ih := digits_sound w s i j us (N >>> w) h.2 e (by simp at he; omega)
      rw [List.getD_cons_succ]
      have hsh : N >>> (w * (e + 1)) = (N >>> w) >>> (w * e) := by
        rw [← Nat.shiftRight_add]; congr 1; ring
      rw [hsh]; exact ih

theorem entry_eval (D w : ℕ) (a b : List (List (List Int))) (hla : a.length = D + 1)
    (hlb : b.length = D + 1) (i j : ℕ) (hi : i < 11) (hj : j < 11) :
    pkDotL22 (pkMatL22 w a) (pkMatL22 w b) i j 11 =
      (∑ k ∈ range 11, polyOf (fun t => (sparseLPoly22 a t i k).toNat) (D + 1) *
        polyOf (fun t => (sparseLPoly22 b t k j).toNat) (D + 1)).eval (2 ^ w) := by
  rw [pkDot_eq, eval_finset_sum]
  refine Finset.sum_congr rfl fun k hk => ?_
  have hk' := Finset.mem_range.mp hk
  rw [pkMat_get w a i k hi hk', pkMat_get w b k j hk' hj, pkEntry_eq, pkEntry_eq, hla, hlb,
    eval_mul]

theorem coeff_prod (α β : ℕ → ℕ → ℕ) (n e : ℕ) :
    (∑ k ∈ range 11, polyOf (α k) n * polyOf (β k) n).coeff e =
      ∑ k ∈ range 11, ∑ x ∈ antidiagonal e, (polyOf (α k) n).coeff x.1 * (polyOf (β k) n).coeff x.2 := by
  rw [finset_sum_coeff]
  exact Finset.sum_congr rfl fun k _ => coeff_mul _ _ _

theorem coeff_prod_lt (D v w : ℕ) (α β : ℕ → ℕ → ℕ) (hα : ∀ k t, α k t < 2 ^ v)
    (hβ : ∀ k t, β k t < 2 ^ v) (hw : 11 * (2 * D + 1) * (2 ^ v * 2 ^ v) < 2 ^ w) (e : ℕ) :
    (∑ k ∈ range 11, polyOf (α k) (D + 1) * polyOf (β k) (D + 1)).coeff e < 2 ^ w := by
  rw [coeff_prod]
  by_cases he : e ≤ 2 * D
  · calc _ ≤ ∑ k ∈ range 11, ∑ x ∈ antidiagonal e, 2 ^ v * 2 ^ v := by
          apply Finset.sum_le_sum; intro k _
          apply Finset.sum_le_sum; intro x _
          rw [coeff_polyOf, coeff_polyOf]
          apply Nat.mul_le_mul
          · split_ifs <;> [exact (hα k _).le; exact Nat.zero_le _]
          · split_ifs <;> [exact (hβ k _).le; exact Nat.zero_le _]
      _ = 11 * ((e + 1) * (2 ^ v * 2 ^ v)) := by
          simp [Finset.sum_const, Finset.Nat.card_antidiagonal]
      _ ≤ 11 * (2 * D + 1) * (2 ^ v * 2 ^ v) := by
          rw [mul_assoc]; apply Nat.mul_le_mul_left; apply Nat.mul_le_mul_right; omega
      _ < 2 ^ w := hw
  · have hz : ∑ k ∈ range 11, ∑ x ∈ antidiagonal e,
        (polyOf (α k) (D + 1)).coeff x.1 * (polyOf (β k) (D + 1)).coeff x.2 = 0 := by
      apply Finset.sum_eq_zero; intro k _
      apply Finset.sum_eq_zero; intro x hx
      have hx : x.1 + x.2 = e := by simpa using hx
      rw [coeff_polyOf, coeff_polyOf]
      by_cases h1 : x.1 < D + 1
      · rw [if_neg (show ¬ x.2 < D + 1 by omega), mul_zero]
      · rw [if_neg h1, zero_mul]
    rw [hz]; positivity

theorem lmat_toNat (P : List (List (List Int))) (v t i k : ℕ) (h : pkBoundL22 v P = true) :
    (((sparseLPoly22 P t i k).toNat : ℕ) : ℤ) = sparseLMat22 (P.getD t []) i k :=
  Int.toNat_of_nonneg (bound_sound v P h t i k).1

theorem trunc_digit (x w n d : ℕ) (hd : d < n) :
    (x % 2 ^ (w * n)) >>> (w * d) % 2 ^ w = x >>> (w * d) % 2 ^ w := by
  rw [Nat.shiftRight_eq_div_pow, Nat.shiftRight_eq_div_pow,
    show w * n = w * d + w * (n - d) by rw [← mul_add]; congr 1; omega, pow_add,
    Nat.mod_mul_right_div_self]
  exact Nat.mod_mod_of_dvd _ (pow_dvd_pow 2 (by nlinarith [show 1 ≤ n - d by omega]))

theorem prod_sound (D v w : ℕ) (s : ℤ) (a b u : List (List (List Int)))
    (ha : ∀ t, t ≤ D → MatShape (a.getD t [])) (hb : ∀ t, t ≤ D → MatShape (b.getD t []))
    (hcheck : pkPolyMulCheckL22 D v w s a b u = true) :
    (∀ d, d ≤ D → MatShape (u.getD d [])) ∧ ∀ d, d ≤ D → ∀ i j, i < 11 → j < 11 →
      (∑ t ∈ range (d + 1),
        sparseLMat22 (sparseMatMulL22 (a.getD t []) (b.getD (d - t) [])) i j) ≤
        s * sparseLMat22 (u.getD d []) i j := by
  simp only [pkPolyMulCheckL22, pkHdrCheckL22, Bool.and_eq_true, decide_eq_true_eq] at hcheck
  obtain ⟨⟨⟨⟨hla, hlb⟩, hlu, hus⟩, ⟨hba, hbb⟩, hw⟩, hprod⟩ := hcheck
  rw [List.all_eq_true] at hus
  refine ⟨fun d hd => ?_, fun d hd i j hi hj => ?_⟩
  · have hdl : d < u.length := by omega
    rw [gd u [] hdl]
    exact matShape_of_check _ (hus _ (List.getElem_mem hdl))
  simp only [pkProdCheckL22, pkRowCheckL22, List.all_eq_true, List.mem_range] at hprod
  have hwalk := digits_sound w s i j u _ (hprod i hi j hj) d (by omega)
  rw [trunc_digit _ w (D + 1) d (by omega)] at hwalk
  have hα : ∀ k t, (fun k t => (sparseLPoly22 a t i k).toNat) k t < 2 ^ v :=
    fun k t => (bound_sound v a hba t i k).2
  have hβ : ∀ k t, (fun k t => (sparseLPoly22 b t k j).toNat) k t < 2 ^ v :=
    fun k t => (bound_sound v b hbb t k j).2
  have hN := entry_eval D w a b hla hlb i j hi hj
  have hdig := digit_eval (2 ^ w) (by positivity) d _ (coeff_prod_lt D v w _ _ hα hβ hw)
  beta_reduce at hdig
  rw [hN, Nat.shiftRight_eq_div_pow, pow_mul, hdig] at hwalk
  simp only [coeff_prod] at hwalk
  refine le_trans (le_of_eq ?_) hwalk
  push_cast
  simp_rw [Finset.Nat.sum_antidiagonal_eq_sum_range_succ_mk]
  rw [Finset.sum_comm]
  refine Finset.sum_congr rfl fun t ht => ?_
  have ht' := Finset.mem_range.mp ht
  rw [(matMul_spec _ _ (ha t (by omega)) (hb (d - t) (by omega))).2 i j hi hj,
    Fin.sum_univ_eq_sum_range
      (fun k => sparseLMat22 (a.getD t []) i k * sparseLMat22 (b.getD (d - t) []) k j) 11]
  refine Finset.sum_congr rfl fun k _ => ?_
  rw [coeff_polyOf, coeff_polyOf, if_pos (by omega), if_pos (by omega)]
  rw [lmat_toNat a v t i k hba, lmat_toNat b v (d - t) k j hbb]

open SparseListPolyD22 in
theorem polyMulPk_sound (D v w : ℕ) (scale : ℤ) (hs : 0 < scale) (a b u : List (List (List Int)))
    (ha : ∀ d, d ≤ D → MatShape (a.getD d [])) (hb : ∀ d, d ≤ D → MatShape (b.getD d []))
    (hcheck : pkPolyMulCheckL22 D v w scale a b u = true)
    (F G : ℕ → Matrix (Fin 11) (Fin 11) ℝ)
    (hF0 : ∀ d i j, 0 ≤ F d i j) (hG0 : ∀ d i j, 0 ≤ G d i j)
    (hF : ∀ d, d ≤ D → ∀ i j : Fin 11, F d i j ≤ (sparseLPoly22 a d i j : ℝ) / scale)
    (hG : ∀ d, d ≤ D → ∀ i j : Fin 11, G d i j ≤ (sparseLPoly22 b d i j : ℝ) / scale) :
    (∀ d, d ≤ D → MatShape (u.getD d [])) ∧
      ∀ d, d ≤ D → ∀ i j : Fin 11,
        (∑ p ∈ antidiagonal d, F p.1 * G p.2) i j ≤ (sparseLPoly22 u d i j : ℝ) / scale := by
  have hP := prod_sound D v w scale a b u ha hb hcheck
  refine ⟨hP.1, fun d hd i j => ?_⟩
  have hI := hP.2 d hd i j i.isLt j.isLt
  have hsR : (0 : ℝ) < scale := by exact_mod_cast hs
  have hpair : ∀ r t, r ≤ D → t ≤ D → (F r * G t) i j ≤
      (sparseLMat22 (sparseMatMulL22 (a.getD r []) (b.getD t [])) i j : ℝ) /
        ((scale : ℝ) * scale) := by
    intro r t hr ht
    rw [(matMul_spec _ _ (ha r hr) (hb t ht)).2 i j i.isLt j.isLt,
      Matrix.mul_apply, Int.cast_sum]
    simp only [Int.cast_mul]
    exact dot_le _ hsR _ _ _ _ (fun k => hF0 r i k) (fun k => hG0 t k j)
      (fun k => hF r hr i k) (fun k => hG t ht k j)
  rw [Finset.Nat.sum_antidiagonal_eq_sum_range_succ_mk, Matrix.sum_apply]
  refine scale_close scale _ _ hs _ ?_ hI
  rw [Int.cast_sum, Finset.sum_div]
  apply Finset.sum_le_sum
  intro t ht
  have ht' := Finset.mem_range.mp ht
  exact hpair t (d - t) (by omega) (by omega)

open SparseListPolyD22 in
theorem add_boundPk (D v w : ℕ) (B : Fin 17 → Matrix (Fin 11) (Fin 11) ℝ)
    (hB : ∀ a i j, 0 ≤ B a i j) (scale : ℤ) (hs : 0 < scale)
    (a b u : List (List (List Int))) (m n : ℕ)
    (hcheck : pkPolyMulCheckL22 D v w scale a b u = true)
    (ha : PolyBoundD B scale m a D) (hb : PolyBoundD B scale n b D) :
    PolyBoundD B scale (m + n) u D := by
  have h := polyMulPk_sound D v w scale hs a b u ha.1 hb.1 hcheck
    (regionalCountDP22 B m) (regionalCountDP22 B n)
    (count_nonneg B hB m) (count_nonneg B hB n) ha.2 hb.2
  refine ⟨h.1, fun d hd i j => ?_⟩
  unfold regionalCountDP22
  rw [orderedMatrixDP_add]
  exact h.2 d hd i j

open SparseListPolyD22 in
theorem chainPk_sound (D v w : ℕ) (B : Fin 17 → Matrix (Fin 11) (Fin 11) ℝ)
    (hB : ∀ a i j, 0 ≤ B a i j) (scale : ℤ) (hs : 0 < scale)
    (base : List (List (List Int))) (hbase : PolyBoundD B scale 1 base D) :
    ∀ (trace : List (Bool × List (List (List Int)) × List (List (List Int))))
      (prev : List (List (List Int))) (e : ℕ),
      PolyBoundD B scale e prev D →
      pkRegionalPowerCheckL22 D v w scale base prev trace = true →
      PolyBoundD B scale (trace.foldl sparseRegionalStepL22 e)
        (sparseRegionalLastL22 prev trace) D
  | [], _, _, h, _ => h
  | (bit, sq, next) :: tail, prev, e, h, hc => by
    have hc' : ((pkPolyMulCheckL22 D v w scale prev prev sq &&
        (if bit then pkPolyMulCheckL22 D v w scale sq base next else true)) &&
        pkRegionalPowerCheckL22 D v w scale base (if bit then next else sq) tail) = true := hc
    rw [Bool.and_eq_true, Bool.and_eq_true] at hc'
    obtain ⟨⟨hsq, hbit⟩, hrest⟩ := hc'
    have hsqB := add_boundPk D v w B hB scale hs prev prev sq e e hsq h h
    cases bit
    · have hrest' : pkRegionalPowerCheckL22 D v w scale base sq tail = true := hrest
      have he : sparseRegionalStepL22 e (false, sq, next) = e + e := by
        unfold sparseRegionalStepL22; simp <;> omega
      have hsq' : PolyBoundD B scale (sparseRegionalStepL22 e (false, sq, next)) sq D := by
        rw [he]; exact hsqB
      exact chainPk_sound D v w B hB scale hs base hbase tail sq _ hsq' hrest'
    · have hbit' : pkPolyMulCheckL22 D v w scale sq base next = true := hbit
      have hrest' : pkRegionalPowerCheckL22 D v w scale base next tail = true := hrest
      have hnext := add_boundPk D v w B hB scale hs sq base next (e + e) 1 hbit' hsqB hbase
      have he : sparseRegionalStepL22 e (true, sq, next) = e + e + 1 := by
        unfold sparseRegionalStepL22; simp <;> omega
      have hnext' : PolyBoundD B scale (sparseRegionalStepL22 e (true, sq, next)) next D := by
        rw [he]; exact hnext
      exact chainPk_sound D v w B hB scale hs base hbase tail next _ hnext' hrest'

/-- Kernel smoke: a 2-stage degree-1 product with `v = 3`, `w = 12`. -/
example : pkPolyMulCheckL22 1 3 12 1
    [List.replicate 11 (List.replicate 11 1), List.replicate 11 (List.replicate 11 1)]
    [List.replicate 11 (List.replicate 11 1), List.replicate 11 (List.replicate 11 1)]
    [List.replicate 11 (List.replicate 11 11), List.replicate 11 (List.replicate 11 22)] = true := by
  decide +kernel

#print axioms chainPk_sound

end SparseListPacked22
end Spin.Turbo
