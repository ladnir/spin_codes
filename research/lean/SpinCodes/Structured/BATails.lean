/-
The sparse half of `lem:structured-ba-tails`.

The paper's chain for the one-accumulator lower tail
(`eq:structured-ba-one-acc-tail`) is

    p_{2ℓ} ≤ C(b,ℓ)·C(D,ℓ) / C(b,2ℓ) = C(2ℓ,ℓ)·C(D,ℓ) / C(b-ℓ,ℓ) ≤ c_o^ℓ,

with `D = ⌊δ_o b⌋`.  This file proves the combinatorial part of that chain —
everything up to and including the middle equality — entirely in `ℕ`, with no
division, so nothing here depends on the numeric value of `δ_o`.

A trap worth recording: the sum must start at `h = 1`.  Truncated subtraction
makes `h - 1 = 0` at `h = 0`, so the bounding term `C(h-1, ℓ-1)` contributes a
spurious `C(0,0) = 1` there.  `accT` itself is correctly zero at `h = 0`
(its `c = 0` guard), so the fix is to split the summand off rather than to
weaken anything.
-/
import SpinCodes.Structured.Accumulator

set_option maxRecDepth 100000

set_option linter.unusedSectionVars false

namespace Spin.Structured

open Finset

/-- `T_b(2ℓ, h)` in closed form for an even input weight and `1 ≤ h ≤ b`. -/
lemma accT_even {b ℓ h : ℕ} (hℓ : 1 ≤ ℓ) (hh : 1 ≤ h) (hb : h ≤ b) :
    accT b (2 * ℓ) h = (h - 1).choose (ℓ - 1) * (b - h).choose ℓ := by
  obtain ⟨e, rfl⟩ := Nat.exists_eq_succ_of_ne_zero (by omega : h ≠ 0)
  obtain ⟨d, rfl⟩ := Nat.exists_eq_add_of_le hb
  have h1 : (2 * ℓ + 1) / 2 - 1 = ℓ - 1 := by omega
  have h2 : 2 * ℓ - (2 * ℓ + 1) / 2 = ℓ := by omega
  rw [accT_eval (by omega : (2 * ℓ) ≠ 0) e d, h1, h2]
  simp

/-- Hockey stick: `Σ_{j<D} C(j,k) = C(D, k+1)`. -/
lemma sum_range_choose_eq (D k : ℕ) : ∑ j ∈ range D, j.choose k = D.choose (k + 1) := by
  induction D with
  | zero => simp
  | succ n ih =>
      rw [Finset.sum_range_succ, ih, Nat.choose_succ_succ]
      simp only [Nat.succ_eq_add_one]
      omega

/-- **The lower-tail sum of the accumulator transition, at even input weight.**

`Σ_{h ≤ D} T_b(2ℓ, h) ≤ C(b,ℓ)·C(D,ℓ)`. -/
theorem sum_accT_even_le {b ℓ D : ℕ} (hℓ : 1 ≤ ℓ) :
    ∑ h ∈ range (D + 1), accT b (2 * ℓ) h ≤ b.choose ℓ * D.choose ℓ := by
  have hzero : accT b (2 * ℓ) 0 = 0 := accT_zero_right _ (by omega)
  rw [Finset.sum_range_succ', hzero, Nat.add_zero]
  have hterm : ∀ i ∈ range D,
      accT b (2 * ℓ) (i + 1) ≤ b.choose ℓ * i.choose (ℓ - 1) := by
    intro i _
    rcases Nat.lt_or_ge b (i + 1) with hb | hb
    · rw [accT_of_lt _ _ _ (by omega) hb]
      exact Nat.zero_le _
    · rw [accT_even hℓ (by omega) hb]
      have hi : i + 1 - 1 = i := by omega
      rw [hi, Nat.mul_comm]
      exact Nat.mul_le_mul_right _ (Nat.choose_le_choose ℓ (by omega))
  calc ∑ i ∈ range D, accT b (2 * ℓ) (i + 1)
      ≤ ∑ i ∈ range D, b.choose ℓ * i.choose (ℓ - 1) := Finset.sum_le_sum hterm
    _ = b.choose ℓ * ∑ i ∈ range D, i.choose (ℓ - 1) := by rw [Finset.mul_sum]
    _ = b.choose ℓ * D.choose ℓ := by
        rw [sum_range_choose_eq]
        congr 2
        omega

/-- The paper's reindexing `C(b,ℓ)/C(b,2ℓ) = C(2ℓ,ℓ)/C(b-ℓ,ℓ)`, cross-multiplied. -/
lemma choose_split (b ℓ : ℕ) :
    b.choose (2 * ℓ) * (2 * ℓ).choose ℓ = b.choose ℓ * (b - ℓ).choose ℓ := by
  have h := Nat.choose_mul (n := b) (k := 2 * ℓ) (s := ℓ) (by omega)
  rw [h]
  congr 2
  omega

/-- **`eq:structured-ba-one-acc-tail`, combinatorial form.**

`p_{2ℓ} ≤ C(2ℓ,ℓ)·C(D,ℓ)/C(b-ℓ,ℓ)`, stated without division: the lower-tail
mass times `C(b-ℓ,ℓ)` is at most `C(2ℓ,ℓ)·C(D,ℓ)` times the normalisation
`C(b,2ℓ)`. -/
theorem sum_accT_even_mul_le {b ℓ D : ℕ} (hℓ : 1 ≤ ℓ) :
    (∑ h ∈ range (D + 1), accT b (2 * ℓ) h) * (b - ℓ).choose ℓ
      ≤ (2 * ℓ).choose ℓ * D.choose ℓ * b.choose (2 * ℓ) := by
  calc (∑ h ∈ range (D + 1), accT b (2 * ℓ) h) * (b - ℓ).choose ℓ
      ≤ (b.choose ℓ * D.choose ℓ) * (b - ℓ).choose ℓ :=
        Nat.mul_le_mul_right _ (sum_accT_even_le hℓ)
    _ = (b.choose ℓ * (b - ℓ).choose ℓ) * D.choose ℓ := by ring
    _ = (b.choose (2 * ℓ) * (2 * ℓ).choose ℓ) * D.choose ℓ := by rw [choose_split]
    _ = (2 * ℓ).choose ℓ * D.choose ℓ * b.choose (2 * ℓ) := by ring

/-- The central binomial factor is at most `4^ℓ` — the source of the factor of
four in the paper's `c_o = 4δ_o/(1-δ_o)`. -/
lemma central_le (ℓ : ℕ) : (2 * ℓ).choose ℓ ≤ 4 ^ ℓ := by
  rw [← Nat.centralBinom_eq_two_mul_choose]
  exact Nat.centralBinom_le_four_pow ℓ

theorem sum_accT_even_mul_le_pow {b ℓ D : ℕ} (hℓ : 1 ≤ ℓ) :
    (∑ h ∈ range (D + 1), accT b (2 * ℓ) h) * (b - ℓ).choose ℓ
      ≤ 4 ^ ℓ * D.choose ℓ * b.choose (2 * ℓ) :=
  le_trans (sum_accT_even_mul_le hℓ)
    (Nat.mul_le_mul_right _ (Nat.mul_le_mul_right _ (central_le ℓ)))

/-- The paper's "if `ℓ > D`, the probability is zero": below the support
threshold there is no room for `ℓ` runs of ones in the first `D` positions. -/
theorem sum_accT_even_eq_zero {b ℓ D : ℕ} (hℓ : 1 ≤ ℓ) (hD : D < ℓ) :
    ∑ h ∈ range (D + 1), accT b (2 * ℓ) h = 0 := by
  refine Nat.le_zero.mp (le_trans (sum_accT_even_le hℓ) ?_)
  rw [Nat.choose_eq_zero_of_lt hD]
  simp



/-! ## The numeric step: `≤ c_o^ℓ`

`δ_o = 13/125` and `c_o = 4δ_o/(1-δ_o) = 13/28`.  Everything stays in `ℕ` by
cross-multiplying, so no rational arithmetic is needed. -/

/-- `(D-k)·m ≤ D·(m-k)` when `D ≤ m`: the step that makes the binomial ratio
`C(D,ℓ)/C(m,ℓ)` decrease factor by factor. -/
private lemma sub_mul_le_mul_sub {D m k : ℕ} (hDm : D ≤ m) :
    (D - k) * m ≤ D * (m - k) := by
  rcases Nat.lt_or_ge k D with hk | hk
  · have h1 : (D - k) * m = D * m - k * m := by rw [Nat.sub_mul]
    have h2 : D * (m - k) = D * m - D * k := by rw [Nat.mul_sub]
    have h3 : D * k ≤ k * m := by
      rw [Nat.mul_comm D k]
      exact Nat.mul_le_mul_left k hDm
    omega
  · have hz : D - k = 0 := by omega
    rw [hz]
    simp

/-- **`C(D,ℓ)/C(m,ℓ) ≤ (D/m)^ℓ` for `D ≤ m`**, cross-multiplied. -/
theorem choose_mul_pow_le {D m : ℕ} (hDm : D ≤ m) (ℓ : ℕ) :
    D.choose ℓ * m ^ ℓ ≤ D ^ ℓ * m.choose ℓ := by
  induction ℓ with
  | zero => simp
  | succ k ih =>
      have hstep : (D - k) * m ≤ D * (m - k) := sub_mul_le_mul_sub hDm
      have key : (D.choose (k + 1) * m ^ (k + 1)) * (k + 1)
          ≤ (D ^ (k + 1) * m.choose (k + 1)) * (k + 1) := by
        calc (D.choose (k + 1) * m ^ (k + 1)) * (k + 1)
            = (D.choose (k + 1) * (k + 1)) * m ^ (k + 1) := by ring
          _ = (D.choose k * (D - k)) * m ^ (k + 1) := by
              rw [Nat.choose_succ_right_eq]
          _ = (D.choose k * m ^ k) * ((D - k) * m) := by ring
          _ ≤ (D ^ k * m.choose k) * ((D - k) * m) := Nat.mul_le_mul_right _ ih
          _ ≤ (D ^ k * m.choose k) * (D * (m - k)) := Nat.mul_le_mul_left _ hstep
          _ = D ^ (k + 1) * (m.choose k * (m - k)) := by ring
          _ = D ^ (k + 1) * (m.choose (k + 1) * (k + 1)) := by
              rw [Nat.choose_succ_right_eq]
          _ = (D ^ (k + 1) * m.choose (k + 1)) * (k + 1) := by ring
      exact Nat.le_of_mul_le_mul_right key (by omega)

/-- The schedule constraint in the form the tail bound uses: from `ℓ ≤ D` and
`D ≤ δ_o b` one gets `D/(b-ℓ) ≤ δ_o/(1-δ_o) = 13/112`. -/
lemma tail_ratio {b ℓ D : ℕ} (hℓD : ℓ ≤ D) (hDb : 125 * D ≤ 13 * b) :
    112 * D ≤ 13 * (b - ℓ) := by omega

/-- **`eq:structured-ba-one-acc-tail`.**

For a uniformly interleaved accumulator on a weight-`2ℓ` input, the mass below
`D` satisfies `p_{2ℓ} ≤ (13/28)^ℓ`, stated as
`Σ_{h ≤ D} T_b(2ℓ,h) · 28^ℓ ≤ 13^ℓ · C(b,2ℓ)`.

The hypotheses are exactly the paper's: `ℓ ≤ D` (otherwise the tail is empty,
`sum_accT_even_eq_zero`) and `D ≤ δ_o b` with `δ_o = 13/125`. -/
theorem sum_accT_even_le_pow {b ℓ D : ℕ} (hℓ : 1 ≤ ℓ) (hℓD : ℓ ≤ D)
    (hDb : 125 * D ≤ 13 * b) :
    (∑ h ∈ range (D + 1), accT b (2 * ℓ) h) * 28 ^ ℓ ≤ 13 ^ ℓ * b.choose (2 * ℓ) := by
  have hℓb : 2 * ℓ ≤ b := by omega
  have hmpos : 0 < (b - ℓ).choose ℓ := Nat.choose_pos (by omega)
  have hbpos : 0 < (b - ℓ) ^ ℓ := pow_pos (by omega) ℓ
  have hbℓ : 112 * D ≤ 13 * (b - ℓ) := tail_ratio hℓD hDb
  have hratio : D.choose ℓ * (b - ℓ) ^ ℓ ≤ D ^ ℓ * (b - ℓ).choose ℓ :=
    choose_mul_pow_le (by omega) ℓ
  have hpow : 112 ^ ℓ * D ^ ℓ ≤ 13 ^ ℓ * (b - ℓ) ^ ℓ := by
    rw [← Nat.mul_pow, ← Nat.mul_pow]
    exact Nat.pow_le_pow_left hbℓ ℓ
  -- the ratio bound, with the powers of `(b-ℓ)` cancelled
  have hkey : 112 ^ ℓ * D.choose ℓ ≤ 13 ^ ℓ * (b - ℓ).choose ℓ := by
    refine Nat.le_of_mul_le_mul_right ?_ hbpos
    calc 112 ^ ℓ * D.choose ℓ * (b - ℓ) ^ ℓ
        = 112 ^ ℓ * (D.choose ℓ * (b - ℓ) ^ ℓ) := by ring
      _ ≤ 112 ^ ℓ * (D ^ ℓ * (b - ℓ).choose ℓ) := Nat.mul_le_mul_left _ hratio
      _ = (112 ^ ℓ * D ^ ℓ) * (b - ℓ).choose ℓ := by ring
      _ ≤ (13 ^ ℓ * (b - ℓ) ^ ℓ) * (b - ℓ).choose ℓ := Nat.mul_le_mul_right _ hpow
      _ = 13 ^ ℓ * (b - ℓ).choose ℓ * (b - ℓ) ^ ℓ := by ring
  set S := ∑ h ∈ range (D + 1), accT b (2 * ℓ) h with hS
  have hmain : S * (b - ℓ).choose ℓ ≤ 4 ^ ℓ * D.choose ℓ * b.choose (2 * ℓ) :=
    sum_accT_even_mul_le_pow hℓ
  refine Nat.le_of_mul_le_mul_right ?_ hmpos
  calc S * 28 ^ ℓ * (b - ℓ).choose ℓ
      = (S * (b - ℓ).choose ℓ) * 28 ^ ℓ := by ring
    _ ≤ (4 ^ ℓ * D.choose ℓ * b.choose (2 * ℓ)) * 28 ^ ℓ := Nat.mul_le_mul_right _ hmain
    _ = (112 ^ ℓ * D.choose ℓ) * b.choose (2 * ℓ) := by
        rw [show (112 : ℕ) ^ ℓ = 4 ^ ℓ * 28 ^ ℓ by rw [← Nat.mul_pow]]
        ring
    _ ≤ (13 ^ ℓ * (b - ℓ).choose ℓ) * b.choose (2 * ℓ) := Nat.mul_le_mul_right _ hkey
    _ = 13 ^ ℓ * b.choose (2 * ℓ) * (b - ℓ).choose ℓ := by ring


/-! ## The odd-weight case

For `v = 2ℓ+1` the run count is `r_v = ℓ+1`, so the transition is
`C(h-1, ℓ)·C(b-h, ℓ)` and the tail sum picks up `C(D, ℓ+1)` instead of
`C(D, ℓ)`.  The extra `D/(b-ℓ) ≤ 13/112 < 1` is the paper's "additional factor
at most `2δ_o/(1-δ_o) < 1`" — in fact a factor of two smaller than that, so
the same `c_o^{⌊v/2⌋}` bound holds. -/

/-- `T_b(2ℓ+1, h)` in closed form for `1 ≤ h ≤ b`. -/
lemma accT_odd {b ℓ h : ℕ} (hh : 1 ≤ h) (hb : h ≤ b) :
    accT b (2 * ℓ + 1) h = (h - 1).choose ℓ * (b - h).choose ℓ := by
  obtain ⟨e, rfl⟩ := Nat.exists_eq_succ_of_ne_zero (by omega : h ≠ 0)
  obtain ⟨d, rfl⟩ := Nat.exists_eq_add_of_le hb
  have h1 : (2 * ℓ + 1 + 1) / 2 - 1 = ℓ := by omega
  have h2 : 2 * ℓ + 1 - (2 * ℓ + 1 + 1) / 2 = ℓ := by omega
  rw [accT_eval (by omega : (2 * ℓ + 1) ≠ 0) e d, h1, h2]
  simp

/-- The odd analogue of `sum_accT_even_le`. -/
theorem sum_accT_odd_le {b ℓ D : ℕ} :
    ∑ h ∈ range (D + 1), accT b (2 * ℓ + 1) h ≤ b.choose ℓ * D.choose (ℓ + 1) := by
  have hzero : accT b (2 * ℓ + 1) 0 = 0 := accT_zero_right _ (by omega)
  rw [Finset.sum_range_succ', hzero, Nat.add_zero]
  have hterm : ∀ i ∈ range D,
      accT b (2 * ℓ + 1) (i + 1) ≤ b.choose ℓ * i.choose ℓ := by
    intro i _
    rcases Nat.lt_or_ge b (i + 1) with hb | hb
    · rw [accT_of_lt _ _ _ (by omega) hb]
      exact Nat.zero_le _
    · rw [accT_odd (by omega) hb]
      have hi : i + 1 - 1 = i := by omega
      rw [hi, Nat.mul_comm]
      exact Nat.mul_le_mul_right _ (Nat.choose_le_choose ℓ (by omega))
  calc ∑ i ∈ range D, accT b (2 * ℓ + 1) (i + 1)
      ≤ ∑ i ∈ range D, b.choose ℓ * i.choose ℓ := Finset.sum_le_sum hterm
    _ = b.choose ℓ * ∑ i ∈ range D, i.choose ℓ := by rw [Finset.mul_sum]
    _ = b.choose ℓ * D.choose (ℓ + 1) := by rw [sum_range_choose_eq]

/-- The odd reindexing: `C(b,ℓ)·C(b-ℓ,ℓ+1) = C(b,2ℓ+1)·C(2ℓ+1,ℓ)`. -/
lemma choose_split_odd (b ℓ : ℕ) :
    b.choose (2 * ℓ + 1) * (2 * ℓ + 1).choose ℓ = b.choose ℓ * (b - ℓ).choose (ℓ + 1) := by
  have h := Nat.choose_mul (n := b) (k := 2 * ℓ + 1) (s := ℓ) (by omega)
  rw [h]
  congr 2
  omega

/-- **`eq:structured-ba-one-acc-tail`, odd case**: `p_{2ℓ+1} ≤ (13/28)^ℓ`. -/
theorem sum_accT_odd_le_pow {b ℓ D : ℕ} (hℓ : 1 ≤ ℓ) (hℓD : ℓ ≤ D)
    (hDb : 125 * D ≤ 13 * b) :
    (∑ h ∈ range (D + 1), accT b (2 * ℓ + 1) h) * 28 ^ ℓ
      ≤ 13 ^ ℓ * b.choose (2 * ℓ + 1) := by
  have hb2 : 2 * ℓ + 1 ≤ b := by omega
  have hmpos : 0 < (b - ℓ).choose (ℓ + 1) := Nat.choose_pos (by omega)
  have hbpos : 0 < (b - ℓ) ^ (ℓ + 1) := pow_pos (by omega) _
  have hbℓ : 112 * D ≤ 13 * (b - ℓ) := tail_ratio hℓD hDb
  have hDbℓ : D ≤ b - ℓ := by omega
  have hratio : D.choose (ℓ + 1) * (b - ℓ) ^ (ℓ + 1)
      ≤ D ^ (ℓ + 1) * (b - ℓ).choose (ℓ + 1) := choose_mul_pow_le hDbℓ (ℓ + 1)
  have hpow : 112 ^ ℓ * D ^ ℓ ≤ 13 ^ ℓ * (b - ℓ) ^ ℓ := by
    rw [← Nat.mul_pow, ← Nat.mul_pow]
    exact Nat.pow_le_pow_left hbℓ ℓ
  -- the ratio bound, one power higher than in the even case
  have hkey : 112 ^ ℓ * D.choose (ℓ + 1) ≤ 13 ^ ℓ * (b - ℓ).choose (ℓ + 1) := by
    refine Nat.le_of_mul_le_mul_right ?_ hbpos
    calc 112 ^ ℓ * D.choose (ℓ + 1) * (b - ℓ) ^ (ℓ + 1)
        = 112 ^ ℓ * (D.choose (ℓ + 1) * (b - ℓ) ^ (ℓ + 1)) := by ring
      _ ≤ 112 ^ ℓ * (D ^ (ℓ + 1) * (b - ℓ).choose (ℓ + 1)) := Nat.mul_le_mul_left _ hratio
      _ = (112 ^ ℓ * D ^ ℓ) * (D * (b - ℓ).choose (ℓ + 1)) := by ring
      _ ≤ (13 ^ ℓ * (b - ℓ) ^ ℓ) * (D * (b - ℓ).choose (ℓ + 1)) :=
          Nat.mul_le_mul_right _ hpow
      _ ≤ (13 ^ ℓ * (b - ℓ) ^ ℓ) * ((b - ℓ) * (b - ℓ).choose (ℓ + 1)) := by
          exact Nat.mul_le_mul_left _ (Nat.mul_le_mul_right _ hDbℓ)
      _ = 13 ^ ℓ * (b - ℓ).choose (ℓ + 1) * (b - ℓ) ^ (ℓ + 1) := by ring
  set S := ∑ h ∈ range (D + 1), accT b (2 * ℓ + 1) h with hS
  have hmain : S * (b - ℓ).choose (ℓ + 1)
      ≤ 4 ^ ℓ * D.choose (ℓ + 1) * b.choose (2 * ℓ + 1) := by
    calc S * (b - ℓ).choose (ℓ + 1)
        ≤ (b.choose ℓ * D.choose (ℓ + 1)) * (b - ℓ).choose (ℓ + 1) :=
          Nat.mul_le_mul_right _ sum_accT_odd_le
      _ = (b.choose ℓ * (b - ℓ).choose (ℓ + 1)) * D.choose (ℓ + 1) := by ring
      _ = (b.choose (2 * ℓ + 1) * (2 * ℓ + 1).choose ℓ) * D.choose (ℓ + 1) := by
          rw [choose_split_odd]
      _ ≤ (b.choose (2 * ℓ + 1) * 4 ^ ℓ) * D.choose (ℓ + 1) := by
          exact Nat.mul_le_mul_right _
            (Nat.mul_le_mul_left _ (Nat.choose_middle_le_pow ℓ))
      _ = 4 ^ ℓ * D.choose (ℓ + 1) * b.choose (2 * ℓ + 1) := by ring
  refine Nat.le_of_mul_le_mul_right ?_ hmpos
  calc S * 28 ^ ℓ * (b - ℓ).choose (ℓ + 1)
      = (S * (b - ℓ).choose (ℓ + 1)) * 28 ^ ℓ := by ring
    _ ≤ (4 ^ ℓ * D.choose (ℓ + 1) * b.choose (2 * ℓ + 1)) * 28 ^ ℓ :=
        Nat.mul_le_mul_right _ hmain
    _ = (112 ^ ℓ * D.choose (ℓ + 1)) * b.choose (2 * ℓ + 1) := by
        rw [show (112 : ℕ) ^ ℓ = 4 ^ ℓ * 28 ^ ℓ by rw [← Nat.mul_pow]]
        ring
    _ ≤ (13 ^ ℓ * (b - ℓ).choose (ℓ + 1)) * b.choose (2 * ℓ + 1) :=
        Nat.mul_le_mul_right _ hkey
    _ = 13 ^ ℓ * b.choose (2 * ℓ + 1) * (b - ℓ).choose (ℓ + 1) := by ring



/-- `eq:structured-ba-one-acc-tail` for odd weights **including `ℓ = 0`**.

At `ℓ = 0` the general bound does not apply (it needs `ℓ ≥ 1`), but the claim
is trivial there: `T_b(1,h) = 1` for `1 ≤ h ≤ b`, so the tail mass is at most
`D ≤ b = C(b,1)`.  Having the uniform statement removes the `v = 2` edge case
from the even upper-tail assembly. -/
theorem sum_accT_odd_le_pow' {b ℓ D : ℕ} (hℓD : ℓ ≤ D) (hDb : 125 * D ≤ 13 * b) :
    (∑ h ∈ range (D + 1), accT b (2 * ℓ + 1) h) * 28 ^ ℓ
      ≤ 13 ^ ℓ * b.choose (2 * ℓ + 1) := by
  rcases Nat.eq_zero_or_pos ℓ with rfl | hℓ
  · have h := sum_accT_odd_le (b := b) (ℓ := 0) (D := D)
    simp only [Nat.mul_zero, Nat.zero_add, Nat.choose_zero_right, Nat.one_mul,
      Nat.zero_add, Nat.choose_one_right] at h ⊢
    simp only [pow_zero, Nat.mul_one, Nat.one_mul, Nat.choose_one_right]
    omega
  · exact sum_accT_odd_le_pow hℓ hℓD hDb

/-! ## The sparse moment ratio

The last step of `eq:structured-ba-sparse-moment` is

    (C(b,ℓ)/C(b,2ℓ)) · (ζ/(1-ζ))^ℓ ≤ (κℓ/(b-2ℓ+1))^ℓ,   κ = 4ζ/(1-ζ).

Writing `κ = 4·(ζ/(1-ζ))`, the factor `(ζ/(1-ζ))^ℓ` appears on both sides and
cancels, so the content is the `ζ`-free statement

    C(b,ℓ) · (b-2ℓ+1)^ℓ ≤ 4^ℓ · ℓ^ℓ · C(b,2ℓ),

which is proved here in `ℕ`.  Nothing about `ζ` is needed. -/

/-- `C(m,ℓ) ≥ (m-ℓ+1)^ℓ / ℓ^ℓ`, cross-multiplied.  Chains
`Nat.pow_sub_le_descFactorial`, `Nat.descFactorial_eq_factorial_mul_choose`
and `Nat.factorial_le_pow`. -/
theorem pow_le_pow_mul_choose (m ℓ : ℕ) :
    (m + 1 - ℓ) ^ ℓ ≤ ℓ ^ ℓ * m.choose ℓ := by
  calc (m + 1 - ℓ) ^ ℓ ≤ m.descFactorial ℓ := Nat.pow_sub_le_descFactorial m ℓ
    _ = Nat.factorial ℓ * m.choose ℓ := Nat.descFactorial_eq_factorial_mul_choose m ℓ
    _ ≤ ℓ ^ ℓ * m.choose ℓ := Nat.mul_le_mul_right _ (Nat.factorial_le_pow ℓ)

/-- **The sparse moment ratio, `ζ`-free.**

`C(b,ℓ)·(b-2ℓ+1)^ℓ ≤ 4^ℓ·ℓ^ℓ·C(b,2ℓ)` — equivalently
`C(b,ℓ)/C(b,2ℓ) ≤ (4ℓ/(b-2ℓ+1))^ℓ / 4^ℓ · 4^ℓ`, which is the last step of
`eq:structured-ba-sparse-moment` once the common `(ζ/(1-ζ))^ℓ` is cancelled. -/
theorem choose_ratio_le {b ℓ : ℕ} (hℓ : 1 ≤ ℓ) (hb : 2 * ℓ ≤ b) :
    b.choose ℓ * (b + 1 - 2 * ℓ) ^ ℓ ≤ 4 ^ ℓ * ℓ ^ ℓ * b.choose (2 * ℓ) := by
  have hmpos : 0 < (b - ℓ).choose ℓ := Nat.choose_pos (by omega)
  have hidx : b - ℓ + 1 - ℓ = b + 1 - 2 * ℓ := by omega
  have hlow : (b + 1 - 2 * ℓ) ^ ℓ ≤ ℓ ^ ℓ * (b - ℓ).choose ℓ := by
    have := pow_le_pow_mul_choose (b - ℓ) ℓ
    rwa [hidx] at this
  refine Nat.le_of_mul_le_mul_right ?_ hmpos
  calc b.choose ℓ * (b + 1 - 2 * ℓ) ^ ℓ * (b - ℓ).choose ℓ
      = (b.choose ℓ * (b - ℓ).choose ℓ) * (b + 1 - 2 * ℓ) ^ ℓ := by ring
    _ = (b.choose (2 * ℓ) * (2 * ℓ).choose ℓ) * (b + 1 - 2 * ℓ) ^ ℓ := by
        rw [choose_split]
    _ ≤ (b.choose (2 * ℓ) * 4 ^ ℓ) * (b + 1 - 2 * ℓ) ^ ℓ := by
        exact Nat.mul_le_mul_right _ (Nat.mul_le_mul_left _ (central_le ℓ))
    _ ≤ (b.choose (2 * ℓ) * 4 ^ ℓ) * (ℓ ^ ℓ * (b - ℓ).choose ℓ) :=
        Nat.mul_le_mul_left _ hlow
    _ = 4 ^ ℓ * ℓ ^ ℓ * b.choose (2 * ℓ) * (b - ℓ).choose ℓ := by ring

/-! ## Sanity checks

The closed form and the vanishing case, against hand-computed values. -/

/-- `b = 10`, input weight `4` (so `ℓ = 2`), `h = 3`:
`C(2,1)·C(7,2) = 2·21 = 42`. -/
example : accT 10 4 3 = 42 := by decide

/-- With `D = 2 < ℓ = 3` the lower tail is empty. -/
example : (∑ h ∈ range 3, accT 20 6 h) = 0 := by decide

/-- And the bound is not vacuous there: `C(20,3)·C(2,3) = 0`. -/
example : (20).choose 3 * (2).choose 3 = 0 := by decide


/-- A concrete instance of the tail bound, at `b = 125` (so `D = ⌊δ_o b⌋ = 13`)
and `ℓ = 2`.  The tail mass is `517595`, and the bound holds with a factor of
about four to spare — so it is not vacuous. -/
example : (∑ h ∈ range 14, accT 125 4 h) = 517595 := by decide

example : (∑ h ∈ range 14, accT 125 4 h) * 28 ^ 2 ≤ 13 ^ 2 * Nat.choose 125 4 := by
  decide


/-- The odd case at `b = 125`, `ℓ = 1` (input weight 3), `D = 13`. -/
example : (∑ h ∈ range 14, accT 125 3 h) * 28 ≤ 13 * Nat.choose 125 3 := by decide


/-- The moment ratio at `b = 100, ℓ = 1`, where it is tightest (ratio `1/2`). -/
example : Nat.choose 100 1 * (100 + 1 - 2) ^ 1 ≤ 4 ^ 1 * 1 ^ 1 * Nat.choose 100 2 := by
  decide

end Spin.Structured
