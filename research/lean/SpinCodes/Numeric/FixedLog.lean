/-
The logarithm in fixed point, sound against `Real.log`.

`RatLog` did this over `ℚ` and is correct, but `DECISIONS.md` D2 records that
the kernel cannot reduce `ℚ`, so it cannot drive a box cover.  This file
rebuilds the same enclosure on `Fix`, where the kernel *can* reduce, and
chains the soundness back to `LogBounds.log_mem_Icc`.

Two things differ from `RatLog` beyond the arithmetic:

* The series is evaluated by Horner, not term by term.  The naive form costs
  `O(n^2)` multiplications, and at the fifty-odd terms wanted here that is
  seconds per argument rather than a tenth of a second.
* The reduction window is `[3/4, 3/2]`, so `|1 - z| ≤ 1/2`.  A narrower window
  would converge faster, but a power-of-two shift can only guarantee a window
  of multiplicative width `2`, and `[3/4, 3/2]` is the symmetric one.  The
  remainder is therefore `≤ 2^{-n}`, which is why `n` is set in the forties
  rather than the twenties.
-/
import SpinCodes.Numeric.Fixed
import SpinCodes.Numeric.LogBounds

set_option linter.unusedSectionVars false

namespace Spin.Numeric

open Finset Real

namespace Fix

/-! ## Two small enclosures -/

lemma zero_mem : Mem zero 0 := by
  constructor <;> simp [zero, Mem]

/-- A symmetric enclosure absorbs anything small enough. -/
lemma pm_mem {r : Fix} {d : ℝ} (h : |d| ≤ (r.hi : ℝ) / ((scale : Int) : ℝ)) :
    Mem (pm r) d := by
  rw [abs_le] at h
  refine ⟨?_, h.2⟩
  show ((-r.hi : Int) : ℝ) / ((scale : Int) : ℝ) ≤ d
  push_cast
  rw [neg_div]
  linarith [h.1]

/-- Anything in an enclosure is at most its upper end. -/
lemma le_hi {a : Fix} {x : ℝ} (h : Mem a x) : x ≤ (a.hi : ℝ) / ((scale : Int) : ℝ) :=
  h.2

lemma half_mem : Mem half ((1 : ℝ) / 2) := by
  have h := ofFrac_mem (p := 1) (q := 2) (by norm_num)
  simpa [half] using h

/-! ## The Horner evaluation

The real-side mirror of `hornerAux`, and the identity that identifies it with
the partial sum `logSeries`. -/

/-- Real-side Horner form. -/
noncomputable def hornerR (t : ℝ) : Nat → Nat → ℝ
  | 0, _ => 0
  | m + 1, i => t * (1 / (i : ℝ) + hornerR t m (i + 1))

lemma hornerR_eq_sum (t : ℝ) :
    ∀ (m i : Nat), hornerR t m i = ∑ j ∈ range m, t ^ (j + 1) / ((i : ℝ) + j)
  | 0, i => by simp [hornerR]
  | m + 1, i => by
      have hterm : ∀ j ∈ range m,
          t * (t ^ (j + 1) / (((i + 1 : Nat) : ℝ) + (j : ℝ)))
            = t ^ ((j + 1) + 1) / ((i : ℝ) + ((j + 1 : Nat) : ℝ)) := by
        intro j _
        push_cast
        ring
      rw [hornerR, hornerR_eq_sum t m (i + 1), mul_add, Finset.mul_sum,
        Finset.sum_congr rfl hterm,
        Finset.sum_range_succ' (fun j => t ^ (j + 1) / ((i : ℝ) + (j : ℝ))) m]
      push_cast
      ring

lemma hornerAux_mem {y : Fix} {t : ℝ} (hy : Mem y t) :
    ∀ (m i : Nat), 0 < i → Mem (hornerAux y m i) (hornerR t m i)
  | 0, _, _ => by simpa [hornerAux, hornerR] using zero_mem
  | m + 1, i, hi => by
      have hifrac : Mem (ofFrac 1 (i : Int)) ((1 : ℝ) / (i : ℝ)) := by
        have h := ofFrac_mem (p := 1) (q := ((i : Int))) (by omega)
        simpa using h
      have hrec := hornerAux_mem hy m (i + 1) (Nat.succ_pos i)
      have := mul_mem hy (add_mem hifrac hrec)
      simpa [hornerAux, hornerR] using this

/-- **The series is enclosed.** -/
theorem series_mem {y : Fix} {t : ℝ} (hy : Mem y t) (n : Nat) :
    Mem (series n y) (logSeries n t) := by
  have h := hornerAux_mem hy n 1 Nat.one_pos
  rw [hornerR_eq_sum] at h
  have : (∑ j ∈ range n, t ^ (j + 1) / (((1 : Nat) : ℝ) + j)) = logSeries n t := by
    unfold logSeries
    apply Finset.sum_congr rfl
    intro j _
    push_cast
    ring
  rwa [this] at h

/-! ## The remainder -/

/-- On the reduction window the remainder is at most `2 * rho^(n+1)`. -/
theorem logRem_le_half {t : ℝ} (ht : |t| ≤ 1 / 2) (n : Nat) :
    logRem n t ≤ 2 * ((1 : ℝ) / 2) ^ (n + 1) := by
  unfold logRem
  have habs : (0 : ℝ) ≤ |t| := abs_nonneg t
  have hden : (1 : ℝ) / 2 ≤ 1 - |t| := by linarith
  have hdenpos : (0 : ℝ) < 1 - |t| := by linarith
  have hnum : |t| ^ (n + 1) ≤ ((1 : ℝ) / 2) ^ (n + 1) :=
    pow_le_pow_left₀ habs ht _
  have hpos : (0 : ℝ) < ((1 : ℝ) / 2) ^ (n + 1) := by positivity
  rw [div_le_iff₀ hdenpos]
  nlinarith

lemma remB_mem {rho : Fix} {r : ℝ} (hr : Mem rho r) (n : Nat) :
    Mem (remB n rho) (2 * r ^ (n + 1)) := by
  have h2 : Mem (ofInt 2) (2 : ℝ) := by simpa using ofInt_mem 2
  have := mul_mem h2 (pow_mem hr (n + 1))
  simpa [remB] using this

/-- The computed remainder bound dominates the true remainder. -/
lemma logRem_le_remB {t : ℝ} (ht : |t| ≤ 1 / 2) (n : Nat) :
    logRem n t ≤ ((remB n half).hi : ℝ) / ((scale : Int) : ℝ) :=
  le_trans (logRem_le_half ht n) (le_hi (remB_mem half_mem n))

/-! ## The odd series

`2·∑_{k<m} z^(2k+1)/(2k+1)` is exactly `logSeries (2m) z - logSeries (2m) (-z)`
— the even terms cancel — so the enclosure comes from applying
`log_one_sub_mem_Icc` at `z` and at `-z` and adding the two remainders.  No new
analysis, and the convergence improves from `|z|` per term to `|z|^2`. -/

/-- Real-side Horner in `z^2`. -/
noncomputable def oddGoR (z2 : ℝ) : Nat → Nat → ℝ
  | 0, _ => 0
  | m + 1, i => 1 / (i : ℝ) + z2 * oddGoR z2 m (i + 2)

/-- `∑_{k<m} z^(2k+1)/(2k+1)`. -/
noncomputable def oddSeriesR (m : Nat) (z : ℝ) : ℝ :=
  ∑ k ∈ range m, z ^ (2 * k + 1) / (2 * k + 1)

lemma oddGoR_eq_sum (z2 : ℝ) :
    ∀ (m i : Nat), oddGoR z2 m i = ∑ k ∈ range m, z2 ^ k / ((i : ℝ) + 2 * k)
  | 0, i => by simp [oddGoR]
  | m + 1, i => by
      have hterm : ∀ k ∈ range m,
          z2 * (z2 ^ k / (((i + 2 : Nat) : ℝ) + 2 * (k : ℝ)))
            = z2 ^ (k + 1) / ((i : ℝ) + 2 * ((k + 1 : Nat) : ℝ)) := by
        intro k _
        push_cast
        rw [pow_succ]
        ring
      rw [oddGoR, oddGoR_eq_sum z2 m (i + 2), Finset.mul_sum,
        Finset.sum_congr rfl hterm,
        Finset.sum_range_succ' (fun k => z2 ^ k / ((i : ℝ) + 2 * (k : ℝ))) m]
      push_cast
      ring

lemma oddSeriesR_eq (m : Nat) (z : ℝ) : oddSeriesR m z = z * oddGoR (z * z) m 1 := by
  rw [oddGoR_eq_sum, oddSeriesR, Finset.mul_sum]
  refine Finset.sum_congr rfl fun k _ => ?_
  rw [show z * z = z ^ 2 by ring, ← pow_mul]
  push_cast
  ring

lemma oddGo_mem {z2 : Fix} {t2 : ℝ} (hz : Mem z2 t2) :
    ∀ (m i : Nat), 0 < i → Mem (oddGo z2 m i) (oddGoR t2 m i)
  | 0, _, _ => by simpa [oddGo, oddGoR] using zero_mem
  | m + 1, i, hi => by
      have hifrac : Mem (ofFrac 1 (i : Int)) ((1 : ℝ) / (i : ℝ)) := by
        have h := ofFrac_mem (p := 1) (q := ((i : Int))) (by omega)
        simpa using h
      have hrec := oddGo_mem hz m (i + 2) (by omega)
      have := add_mem hifrac (mul_mem hz hrec)
      simpa [oddGo, oddGoR] using this

theorem oddSeries_mem {z : Fix} {t : ℝ} (hz : Mem z t) (m : Nat) :
    Mem (oddSeries m z) (oddSeriesR m t) := by
  rw [oddSeriesR_eq]
  exact mul_mem hz (oddGo_mem (mul_mem hz hz) m 1 Nat.one_pos)

/-- The even terms cancel. -/
lemma two_oddSeriesR (z : ℝ) :
    ∀ m : Nat, 2 * oddSeriesR m z = logSeries (2 * m) z - logSeries (2 * m) (-z)
  | 0 => by simp [oddSeriesR, logSeries]
  | m + 1 => by
      have hstep : ∀ x : ℝ, logSeries (2 * (m + 1)) x
          = logSeries (2 * m) x + x ^ (2 * m + 1) / (2 * m + 1)
            + x ^ (2 * m + 2) / (2 * m + 2) := by
        intro x
        have h : 2 * (m + 1) = (2 * m + 1) + 1 := by ring
        rw [h]
        simp only [logSeries]
        rw [Finset.sum_range_succ, Finset.sum_range_succ]
        push_cast
        ring
      have hodd : (-z) ^ (2 * m + 1) = -z ^ (2 * m + 1) :=
        Odd.neg_pow ⟨m, by ring⟩ z
      have heven : (-z) ^ (2 * m + 2) = z ^ (2 * m + 2) :=
        Even.neg_pow ⟨m + 1, by ring⟩ z
      have ih := two_oddSeriesR z m
      rw [hstep, hstep, oddSeriesR, Finset.sum_range_succ, ← oddSeriesR, hodd, heven]
      push_cast
      linear_combination ih

/-- **The enclosure.** -/
theorem log_ratio_bound {z : ℝ} (hz : |z| < 1) (m : Nat) :
    |Real.log ((1 + z) / (1 - z)) - 2 * oddSeriesR m z| ≤ 2 * logRem (2 * m) z := by
  rw [abs_lt] at hz
  have hp : (0 : ℝ) < 1 + z := by linarith [hz.1]
  have hn : (0 : ℝ) < 1 - z := by linarith [hz.2]
  have habsz : |z| < 1 := by rw [abs_lt]; exact hz
  have habsnz : |(-z)| < 1 := by rwa [abs_neg]
  have hA := log_one_sub_mem_Icc habsnz (2 * m)
  have hB := log_one_sub_mem_Icc habsz (2 * m)
  rw [show (1 : ℝ) - -z = 1 + z by ring] at hA
  have hrem : logRem (2 * m) (-z) = logRem (2 * m) z := by
    unfold logRem
    rw [abs_neg]
  rw [hrem] at hA
  have hsplit : Real.log ((1 + z) / (1 - z)) = Real.log (1 + z) - Real.log (1 - z) :=
    Real.log_div (ne_of_gt hp) (ne_of_gt hn)
  have hser := two_oddSeriesR z m
  rw [abs_le, hsplit, hser]
  constructor <;> linarith [hA.1, hA.2, hB.1, hB.2]

/-- The computed remainder constant dominates, for `|z| ≤ 1/c` with `3 ≤ c`. -/
lemma two_logRem_le {z : ℝ} {c : Int} (hc : 3 ≤ c) (hz : |z| ≤ 1 / (c : ℝ))
    (m : Nat) : 2 * logRem (2 * m) z ≤ 3 * ((1 : ℝ) / (c : ℝ)) ^ (2 * m + 1) := by
  have hcR : (3 : ℝ) ≤ (c : ℝ) := by exact_mod_cast hc
  have hcpos : (0 : ℝ) < (c : ℝ) := by linarith
  have hinv : (1 : ℝ) / (c : ℝ) ≤ 1 / 3 := by
    rw [div_le_div_iff₀ hcpos (by norm_num)]
    linarith
  have habs : (0 : ℝ) ≤ |z| := abs_nonneg z
  have hden : (2 : ℝ) / 3 ≤ 1 - |z| := by linarith
  have hdpos : (0 : ℝ) < 1 - |z| := by linarith
  have hnum : |z| ^ (2 * m + 1) ≤ ((1 : ℝ) / (c : ℝ)) ^ (2 * m + 1) :=
    pow_le_pow_left₀ habs hz _
  have hpos : (0 : ℝ) < ((1 : ℝ) / (c : ℝ)) ^ (2 * m + 1) := by positivity
  unfold logRem
  rw [mul_div_assoc', div_le_iff₀ hdpos]
  nlinarith [hnum, hpos, hden, mul_le_mul_of_nonneg_left hden hpos.le]

lemma remOdd_mem {c : Int} (hc : 0 < c) (m : Nat) :
    Mem (remOdd c m) (3 * ((1 : ℝ) / (c : ℝ)) ^ (2 * m + 1)) := by
  have h3 : Mem (ofInt 3) (3 : ℝ) := by simpa using ofInt_mem 3
  have hc1 : Mem (ofFrac 1 c) ((1 : ℝ) / (c : ℝ)) := by
    have h := ofFrac_mem (p := 1) (q := c) hc
    simpa using h
  have := mul_mem h3 (pow_mem hc1 (2 * m + 1))
  simpa [remOdd] using this

/-! ## `log 2` -/

theorem log2_mem : Mem log2 (Real.log 2) := by
  have hz : Mem (ofFrac 1 3) ((1 : ℝ) / 3) := by
    have h := ofFrac_mem (p := 1) (q := 3) (by norm_num)
    simpa using h
  have hser := oddSeries_mem hz 21
  have h2 : Mem (ofInt 2) (2 : ℝ) := by simpa using ofInt_mem 2
  have hmul := mul_mem h2 hser
  have habs : |(1 : ℝ) / 3| < 1 := by
    rw [abs_of_nonneg (by norm_num : (0:ℝ) ≤ 1/3)]; norm_num
  have hbd := log_ratio_bound habs 21
  rw [show ((1 : ℝ) + 1 / 3) / ((1 : ℝ) - 1 / 3) = 2 by norm_num] at hbd
  have habsle : |(1 : ℝ) / 3| ≤ 1 / (((3 : Int)) : ℝ) := by
    rw [abs_of_nonneg (by norm_num : (0:ℝ) ≤ 1/3)]; norm_num
  have hrem := two_logRem_le (c := 3) (by norm_num) habsle 21
  have hrm := remOdd_mem (c := 3) (by norm_num) 21
  have hd : |Real.log 2 - 2 * oddSeriesR 21 ((1 : ℝ) / 3)|
      ≤ ((remOdd 3 21).hi : ℝ) / ((scale : Int) : ℝ) :=
    le_trans hbd (le_trans hrem (le_hi hrm))
  have hfinal := add_mem hmul (pm_mem hd)
  rw [show 2 * oddSeriesR 21 ((1 : ℝ) / 3)
      + (Real.log 2 - 2 * oddSeriesR 21 ((1 : ℝ) / 3)) = Real.log 2 by ring] at hfinal
  exact hfinal

/-! ## The odd-series logarithm -/

/-- **`flogA` encloses `Real.log (P/Q)` on the reduction window `[2/3, 3/2]`.**

The window hypotheses are exactly what makes `|z| ≤ 1/5`: `2Q ≤ 3P` gives
`5(Q-P) ≤ P+Q` and `2P ≤ 3Q` gives `5(P-Q) ≤ P+Q`. -/
theorem flogA_mem {P Q : Int} {m : Nat} (hP : 0 < P) (hQ : 0 < Q)
    (h1 : 2 * Q ≤ 3 * P) (h2 : 2 * P ≤ 3 * Q) :
    Mem (flogA P Q m) (Real.log ((P : ℝ) / (Q : ℝ))) := by
  have hPR : (0 : ℝ) < (P : ℝ) := by exact_mod_cast hP
  have hQR : (0 : ℝ) < (Q : ℝ) := by exact_mod_cast hQ
  have h1R : 2 * (Q : ℝ) ≤ 3 * (P : ℝ) := by exact_mod_cast h1
  have h2R : 2 * (P : ℝ) ≤ 3 * (Q : ℝ) := by exact_mod_cast h2
  have hsum : (0 : ℝ) < (P : ℝ) + (Q : ℝ) := by linarith
  set z : ℝ := ((P : ℝ) - (Q : ℝ)) / ((P : ℝ) + (Q : ℝ)) with hzdef
  have hzle : |z| ≤ 1 / (((5 : Int)) : ℝ) := by
    have habs : |(P : ℝ) - (Q : ℝ)| ≤ ((P : ℝ) + (Q : ℝ)) / 5 := by
      rw [abs_le]
      constructor <;> linarith
    rw [hzdef, abs_div, abs_of_pos hsum,
      div_le_div_iff₀ hsum (by norm_num : (0:ℝ) < (((5 : Int)) : ℝ))]
    push_cast
    linarith [habs]
  have hzlt : |z| < 1 := by
    have : (1 : ℝ) / (((5 : Int)) : ℝ) < 1 := by norm_num
    linarith [hzle]
  -- the enclosure of `z` itself
  have hzmem : Mem (ofFrac (P - Q) (P + Q)) z := by
    have h := ofFrac_mem (p := P - Q) (q := P + Q) (by omega)
    have hc : (((P - Q : Int)) : ℝ) / (((P + Q : Int)) : ℝ) = z := by
      rw [hzdef]; push_cast; ring
    rwa [hc] at h
  have hser := oddSeries_mem hzmem m
  have h2i : Mem (ofInt 2) (2 : ℝ) := by simpa using ofInt_mem 2
  have hmul := mul_mem h2i hser
  have hbd := log_ratio_bound hzlt m
  have hsne : ((P : ℝ) + (Q : ℝ)) ≠ 0 := ne_of_gt hsum
  have hQne : (Q : ℝ) ≠ 0 := ne_of_gt hQR
  have h1z : (1 : ℝ) + z = 2 * (P : ℝ) / ((P : ℝ) + (Q : ℝ)) := by
    rw [hzdef]; field_simp; ring
  have h2z : (1 : ℝ) - z = 2 * (Q : ℝ) / ((P : ℝ) + (Q : ℝ)) := by
    rw [hzdef]; field_simp; ring
  have hratio : ((1 : ℝ) + z) / ((1 : ℝ) - z) = (P : ℝ) / (Q : ℝ) := by
    rw [h1z, h2z]
    field_simp
  rw [hratio] at hbd
  have hrem := two_logRem_le (c := 5) (by norm_num) hzle m
  have hrm := remOdd_mem (c := 5) (by norm_num) m
  have hd : |Real.log ((P : ℝ) / (Q : ℝ)) - 2 * oddSeriesR m z|
      ≤ ((remOdd 5 m).hi : ℝ) / ((scale : Int) : ℝ) :=
    le_trans hbd (le_trans hrem (le_hi hrm))
  have hfinal := add_mem hmul (pm_mem hd)
  rw [show 2 * oddSeriesR m z + (Real.log ((P : ℝ) / (Q : ℝ)) - 2 * oddSeriesR m z)
      = Real.log ((P : ℝ) / (Q : ℝ)) by ring] at hfinal
  exact hfinal

/-! ## The shift -/

lemma ipow2_cast : ∀ k : Nat, ((ipow2 k : Int) : ℝ) = 2 ^ k
  | 0 => by simp [ipow2]
  | k + 1 => by
      rw [ipow2]
      push_cast
      rw [ipow2_cast k, pow_succ]
      ring

lemma ipow2_pos : ∀ k : Nat, 0 < ipow2 k
  | 0 => by decide
  | k + 1 => by
      have := ipow2_pos k
      rw [ipow2]
      omega

/-! ## The enclosure -/

/-- **`flog` encloses `Real.log (p / q)`.**

The two window hypotheses say `3/4 ≤ (p/q) · 2^k ≤ 3/2`.  They are `Int`
inequalities, so a certificate discharges them by `decide` together with
whatever bound it is really checking. -/
theorem flog_mem {p q : Int} {k n : Nat} (hp : 0 < p) (hq : 0 < q)
    (hlo : 3 * q ≤ 4 * (p * ipow2 k)) (hhi : 2 * (p * ipow2 k) ≤ 3 * q) :
    Mem (flog p q k n) (Real.log ((p : ℝ) / (q : ℝ))) := by
  have hqR : (0 : ℝ) < (q : ℝ) := by exact_mod_cast hq
  have hpR : (0 : ℝ) < (p : ℝ) := by exact_mod_cast hp
  have hkR : (0 : ℝ) < (2 : ℝ) ^ k := by positivity
  set Q : ℝ := (p : ℝ) / (q : ℝ) with hQ
  set z : ℝ := Q * 2 ^ k with hz
  -- the window, transported to `ℝ`
  have hloR : 3 * (q : ℝ) ≤ 4 * ((p : ℝ) * 2 ^ k) := by
    have : ((3 * q : Int) : ℝ) ≤ ((4 * (p * ipow2 k) : Int) : ℝ) := by exact_mod_cast hlo
    push_cast at this
    rwa [ipow2_cast] at this
  have hhiR : 2 * ((p : ℝ) * 2 ^ k) ≤ 3 * (q : ℝ) := by
    have : ((2 * (p * ipow2 k) : Int) : ℝ) ≤ ((3 * q : Int) : ℝ) := by exact_mod_cast hhi
    push_cast at this
    rwa [ipow2_cast] at this
  have hzeq : z = (p : ℝ) * 2 ^ k / (q : ℝ) := by
    rw [hz, hQ]; ring
  have hz34 : (3 : ℝ) / 4 ≤ z := by
    rw [hzeq, le_div_iff₀ hqR]; linarith
  have hz32 : z ≤ (3 : ℝ) / 2 := by
    rw [hzeq, div_le_iff₀ hqR]; linarith
  have hz0 : (0 : ℝ) < z := by linarith
  have hz2 : z < 2 := by linarith
  -- the series argument
  set t : ℝ := 1 - z with ht
  have habs : |t| ≤ 1 / 2 := by
    rw [abs_le]; constructor <;> (rw [ht]; linarith)
  have hIcc := log_mem_Icc hz0 hz2 n
  rw [← ht] at hIcc
  -- the `Fix` side of the argument
  set y : Fix := ofFrac (q - p * ipow2 k) q with hy
  have hymem : Mem y t := by
    have h := ofFrac_mem (p := q - p * ipow2 k) (q := q) hq
    have hcast : (((q - p * ipow2 k : Int)) : ℝ) / (q : ℝ) = t := by
      push_cast
      rw [ipow2_cast, ht, hzeq]
      field_simp
    rwa [hcast] at h
  have hser : Mem (series n y) (logSeries n t) := series_mem hymem n
  have hnegser : Mem (neg (series n y)) (-logSeries n t) := neg_mem hser
  have hrem := logRem_le_remB habs n
  have hd : |Real.log z - -logSeries n t|
      ≤ ((remB n half).hi : ℝ) / ((scale : Int) : ℝ) := by
    rw [abs_le]
    constructor <;> linarith [hIcc.1, hIcc.2, hrem]
  have hmain : Mem (add (neg (series n y)) (pm (remB n half))) (Real.log z) := by
    have := add_mem hnegser (pm_mem hd)
    simpa using this
  -- the shift term
  have hklog : Mem (mul (ofInt (k : Int)) log2) ((k : ℝ) * Real.log 2) := by
    have hk : Mem (ofInt (k : Int)) ((k : ℝ)) := by
      have := ofInt_mem ((k : Int))
      simpa using this
    exact mul_mem hk log2_mem
  have hsplit : Real.log z - (k : ℝ) * Real.log 2 = Real.log Q := by
    rw [hz, Real.log_mul (by positivity) (by positivity), Real.log_pow]
    push_cast
    ring
  have := sub_mem hmain hklog
  rw [hsplit] at this
  simpa [flog, hy] using this

/-! ## The logarithm of an enclosure

`Real.log` is monotone, so an enclosure of `log` on `[lo, hi]` is pinned by its
endpoints — and each endpoint is a *rational*, `lo / scale`, which `flog`
already handles.  Two `flog` calls per interval logarithm, not a fresh series
over intervals. -/

theorem flogI_mem {a : Fix} {x : ℝ} {klo khi n : Nat} (hx : Mem a x) (ha : 0 < a.lo)
    (hlo1 : 3 * scale ≤ 4 * (a.lo * ipow2 klo))
    (hlo2 : 2 * (a.lo * ipow2 klo) ≤ 3 * scale)
    (hhi1 : 3 * scale ≤ 4 * (a.hi * ipow2 khi))
    (hhi2 : 2 * (a.hi * ipow2 khi) ≤ 3 * scale) :
    Mem (flogI a klo khi n) (Real.log x) := by
  have hS := scaleR_pos
  have haR : (0 : ℝ) < (a.lo : ℝ) := by exact_mod_cast ha
  have hhipos : 0 < a.hi := pos_hi hx ha
  have hxpos : 0 < x := pos_of_mem hx ha
  have hlopos : (0 : ℝ) < (a.lo : ℝ) / ((scale : Int) : ℝ) := div_pos haR hS
  have hL := flog_mem (p := a.lo) (q := scale) (k := klo) (n := n) ha scale_pos hlo1 hlo2
  have hH := flog_mem (p := a.hi) (q := scale) (k := khi) (n := n) hhipos scale_pos hhi1 hhi2
  refine ⟨?_, ?_⟩
  · show ((flog a.lo scale klo n).lo : ℝ) / ((scale : Int) : ℝ) ≤ Real.log x
    exact le_trans hL.1 (Real.log_le_log hlopos hx.1)
  · show Real.log x ≤ ((flog a.hi scale khi n).hi : ℝ) / ((scale : Int) : ℝ)
    exact le_trans (Real.log_le_log hxpos hx.2) hH.2

end Fix

end Spin.Numeric
