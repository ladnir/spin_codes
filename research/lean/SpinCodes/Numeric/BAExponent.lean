/-
The BA outer exponent, on `ℝ`.

This is the real-side half of T7b: the definitions `h`, `π` and `g` exactly as
`structured_appendix.tex` writes them, plus the two facts the interval
verifier silently relies on.

**The `π` identity.**  The paper defines

    π(a,c) = c·h(a/(2c)) + (1-c)·h(a/(2(1-c))) - h(a)

but the verifier evaluates a different expression — one with no division, and
with the `log a` terms already cancelled.  `piBA_eq_piEval` proves the two
agree.  The cancellation is not cosmetic: `a/(2c)` is unstable as `c → a/2`,
and the `log a` terms are what make the naive form lose all precision near
`a = 0`, which is exactly the corner of the box cover that matters.

**`g` is an infimum, so it needs a floor.**  `gBA a ≤ gObj u a` is only
available once the family is bounded below, and it is not bounded below for
`a > 1`.  For `0 ≤ a ≤ 1` there is a clean floor: `gObj u a ≥ 0` for every
`u > 0`, because `GolayG u` is at least both `1` and `u^24`.

The Golay constants `759, 2576, 759` are not typed in again here: `GolayG_eq_sum`
derives the generating function from `Golay.enumerator_eq`, which is checked
against the actual code in `Structured/Golay.lean`.
-/
import Mathlib
import SpinCodes.Structured.GolayEnum

set_option linter.unusedSectionVars false

namespace Spin.Numeric

open Real

/-! ## Entropy -/

/-- Binary entropy, in nats.  `Real.log 0 = 0` gives the endpoint values the
paper asks for by continuity. -/
noncomputable def hEnt (x : ℝ) : ℝ := -x * Real.log x - (1 - x) * Real.log (1 - x)

@[simp] lemma hEnt_zero : hEnt 0 = 0 := by simp [hEnt]

@[simp] lemma hEnt_one : hEnt 1 = 0 := by simp [hEnt]

/-! ## `x log x`

The box cover needs this as a unit, not as a product of a factor and a
logarithm.  On the feasibility boundary `c = a/2` the factor is zero and the
logarithm is `log 0`; the paper and the verifier both read the product as `0`
by continuity, which only makes sense if the whole term is the primitive. -/

/-- `t log t`.  `Real.log 0 = 0` gives the continuous value `0` at the origin. -/
noncomputable def xlogx (t : ℝ) : ℝ := t * Real.log t

@[simp] lemma xlogx_zero : xlogx 0 = 0 := by simp [xlogx]

/-- **The tangent bound.**  `x log x` lies above its tangent at any `m > 0`:

    x·log m + x - m ≤ x·log x

This is Gibbs' inequality `x log(x/m) ≥ x - m` rearranged, and it is exact at
`x = m`.  Taking `m` at the middle of a box gives a lower bound good enough
that no separate treatment of the minimum at `1/e` is needed. -/
theorem xlogx_tangent {x m : ℝ} (hx : 0 ≤ x) (hm : 0 < m) :
    x * Real.log m + x - m ≤ xlogx x := by
  rcases eq_or_lt_of_le hx with hzero | hxpos
  · rw [← hzero]
    simp [xlogx]
    linarith
  · have h1 : Real.log (m / x) ≤ m / x - 1 := Real.log_le_sub_one_of_pos (by positivity)
    have h2 : Real.log (m / x) = Real.log m - Real.log x :=
      Real.log_div (ne_of_gt hm) (ne_of_gt hxpos)
    rw [h2] at h1
    have h3 : x * (Real.log m - Real.log x) ≤ x * (m / x - 1) :=
      mul_le_mul_of_nonneg_left h1 hxpos.le
    have h4 : x * (m / x - 1) = m - x := by field_simp
    unfold xlogx
    nlinarith [h3, h4]

/-- **The maximum principle.**  `x log x` is convex, so on a box its largest
value sits at an endpoint. -/
theorem xlogx_le_max {lo hi x : ℝ} (hlo : 0 ≤ lo) (hx : x ∈ Set.Icc lo hi) :
    xlogx x ≤ max (xlogx lo) (xlogx hi) := by
  have hhi : 0 ≤ hi := le_trans hlo (le_trans hx.1 hx.2)
  exact Real.convexOn_mul_log.le_max_of_mem_Icc hlo hhi hx

/-- Entropy in terms of `x log x`. -/
lemma hEnt_eq_xlogx (t : ℝ) : hEnt t = -(xlogx t + xlogx (1 - t)) := by
  unfold hEnt xlogx
  ring

/-! ## The accumulator exponent -/

/-- `π(a,c)`, as the paper defines it. -/
noncomputable def piBA (a c : ℝ) : ℝ :=
  c * hEnt (a / (2 * c)) + (1 - c) * hEnt (a / (2 * (1 - c))) - hEnt a

/-- `π(a,c)`, as the interval verifier evaluates it: division-free, with the
`log a` terms cancelled. -/
noncomputable def piEval (a c : ℝ) : ℝ :=
  a * Real.log 2 + c * Real.log c + (1 - c) * Real.log (1 - c)
    - (c - a / 2) * Real.log (c - a / 2)
    - (1 - c - a / 2) * Real.log (1 - c - a / 2)
    + (1 - a) * Real.log (1 - a)

/-- One half of the expansion: `c · h(a/(2c))` with the logarithm split. -/
private lemma half_expand {a c : ℝ} (ha : 0 < a) (hc : 0 < c) (hac : a / 2 < c) :
    c * hEnt (a / (2 * c))
      = -(a / 2) * Real.log a + (a / 2) * Real.log 2 + c * Real.log c
        - (c - a / 2) * Real.log (c - a / 2) := by
  have hca : 0 < c - a / 2 := by linarith
  have h2c : (0 : ℝ) < 2 * c := by linarith
  have hsplit : Real.log (a / (2 * c)) = Real.log a - Real.log 2 - Real.log c := by
    rw [Real.log_div (ne_of_gt ha) (ne_of_gt h2c), Real.log_mul (by norm_num) (ne_of_gt hc)]
    ring
  have hone : (1 : ℝ) - a / (2 * c) = (c - a / 2) / c := by
    field_simp
  have hsplit2 : Real.log ((c - a / 2) / c) = Real.log (c - a / 2) - Real.log c :=
    Real.log_div (ne_of_gt hca) (ne_of_gt hc)
  have hmul : c * (a / (2 * c)) = a / 2 := by field_simp
  have hmul2 : c * ((c - a / 2) / c) = c - a / 2 := by field_simp
  unfold hEnt
  rw [hone, hsplit, hsplit2]
  have : c * (-(a / (2 * c)) * (Real.log a - Real.log 2 - Real.log c)
      - (c - a / 2) / c * (Real.log (c - a / 2) - Real.log c))
      = -(c * (a / (2 * c))) * (Real.log a - Real.log 2 - Real.log c)
        - (c * ((c - a / 2) / c)) * (Real.log (c - a / 2) - Real.log c) := by
    ring
  rw [this, hmul, hmul2]
  ring

/-- **The paper's `π` is the verifier's `p`.** -/
theorem piBA_eq_piEval {a c : ℝ} (ha : 0 ≤ a) (ha1 : a < 1)
    (hlo : a / 2 < c) (hhi : c < 1 - a / 2) :
    piBA a c = piEval a c := by
  have hc : 0 < c := lt_of_le_of_lt (by linarith) hlo
  have hc1 : c < 1 := by linarith
  rcases eq_or_lt_of_le ha with hzero | hpos
  · -- `a = 0`: both sides collapse to zero.
    subst_vars
    simp only [piBA, piEval, hEnt_zero]
    rw [show (0 : ℝ) / (2 * c) = 0 by ring, show (0 : ℝ) / (2 * (1 - c)) = 0 by ring]
    simp only [hEnt_zero]
    norm_num
  · -- the generic case
    have hd : 0 < 1 - c := by linarith
    have hdlo : a / 2 < 1 - c := by linarith
    have h1 := half_expand hpos hc hlo
    have h2 := half_expand hpos hd hdlo
    unfold piBA piEval
    rw [h1, h2]
    unfold hEnt
    rw [show (1 : ℝ) - c - a / 2 = (1 - c) - a / 2 by ring]
    ring

/-- The verifier's `π` in terms of `x log x`. -/
lemma piEval_eq_xlogx (a c : ℝ) :
    piEval a c = a * Real.log 2 + xlogx c + xlogx (1 - c)
      - xlogx (c - a / 2) - xlogx (1 - c - a / 2) + xlogx (1 - a) := by
  unfold piEval xlogx
  ring

/-! ## The partial derivatives of `π`

The interval verifier's decisive bound is a centered (mean-value) form for
`p`, and iteration 36 measured that this one bound is enough on its own: with
it and no other mean-value form the cover closes in 13,512 boxes, against
6,747 with both and no termination at all with neither.

So only these two partial derivatives are needed — the three-variable form for
the whole objective is not. -/

/-- `∂π/∂a = ln 2 + (ln(c - a/2) + ln(1 - c - a/2))/2 - ln(1 - a)`. -/
noncomputable def Dpa (a c : ℝ) : ℝ :=
  Real.log 2 + (Real.log (c - a / 2) + Real.log (1 - c - a / 2)) / 2 - Real.log (1 - a)

/-- `∂π/∂c = ln c - ln(1-c) - ln(c - a/2) + ln(1 - c - a/2)`. -/
noncomputable def Dpc (a c : ℝ) : ℝ :=
  Real.log c - Real.log (1 - c) - Real.log (c - a / 2) + Real.log (1 - c - a / 2)

theorem hasDerivAt_piEval_fst {a c : ℝ} (h1 : 0 < c - a / 2) (h2 : 0 < 1 - c - a / 2)
    (h3 : a < 1) : HasDerivAt (fun t => piEval t c) (Dpa a c) a := by
  have hu : HasDerivAt (fun t : ℝ => c - t / 2) (-(1 / 2)) a :=
    ((hasDerivAt_id a).div_const 2).const_sub c
  have hv : HasDerivAt (fun t : ℝ => 1 - c - t / 2) (-(1 / 2)) a :=
    ((hasDerivAt_id a).div_const 2).const_sub (1 - c)
  have hw : HasDerivAt (fun t : ℝ => 1 - t) (-1 : ℝ) a := (hasDerivAt_id a).const_sub 1
  have hU : HasDerivAt (fun t : ℝ => (c - t / 2) * Real.log (c - t / 2))
      ((Real.log (c - a / 2) + 1) * (-(1 / 2))) a := by
    have h := HasDerivAt.comp a (Real.hasDerivAt_mul_log (x := c - a / 2) (ne_of_gt h1)) hu
    rwa [Function.comp_def] at h
  have hV : HasDerivAt (fun t : ℝ => (1 - c - t / 2) * Real.log (1 - c - t / 2))
      ((Real.log (1 - c - a / 2) + 1) * (-(1 / 2))) a := by
    have h := HasDerivAt.comp a
      (Real.hasDerivAt_mul_log (x := 1 - c - a / 2) (ne_of_gt h2)) hv
    rwa [Function.comp_def] at h
  have hW : HasDerivAt (fun t : ℝ => (1 - t) * Real.log (1 - t))
      ((Real.log (1 - a) + 1) * (-1)) a := by
    have h := HasDerivAt.comp a
      (Real.hasDerivAt_mul_log (x := 1 - a) (by linarith : (1 : ℝ) - a ≠ 0)) hw
    rwa [Function.comp_def] at h
  have hlin : HasDerivAt (fun t : ℝ => t * Real.log 2) (Real.log 2) a := by
    simpa using (hasDerivAt_id a).mul_const (Real.log 2)
  have hconst : HasDerivAt
      (fun _ : ℝ => c * Real.log c + (1 - c) * Real.log (1 - c)) 0 a :=
    hasDerivAt_const _ _
  have hsum := (((hlin.add hconst).sub hU).sub hV).add hW
  rw [show Dpa a c = Real.log 2 + 0 - (Real.log (c - a / 2) + 1) * (-(1 / 2))
      - (Real.log (1 - c - a / 2) + 1) * (-(1 / 2))
      + (Real.log (1 - a) + 1) * (-1) from by unfold Dpa; ring]
  convert hsum using 1
  funext t
  simp only [Pi.add_apply, Pi.sub_apply]
  unfold piEval
  ring

theorem hasDerivAt_piEval_snd {a c : ℝ} (h1 : 0 < c - a / 2) (h2 : 0 < 1 - c - a / 2)
    (h4 : 0 < c) (h5 : c < 1) : HasDerivAt (fun t => piEval a t) (Dpc a c) c := by
  have hs : HasDerivAt (fun t : ℝ => 1 - t) (-1 : ℝ) c := (hasDerivAt_id c).const_sub 1
  have hu : HasDerivAt (fun t : ℝ => t - a / 2) (1 : ℝ) c := by
    simpa using (hasDerivAt_id c).sub_const (a / 2)
  have hv : HasDerivAt (fun t : ℝ => 1 - t - a / 2) (-1 : ℝ) c := by
    simpa using ((hasDerivAt_id c).const_sub 1).sub_const (a / 2)
  have hC : HasDerivAt (fun t : ℝ => t * Real.log t) (Real.log c + 1) c :=
    Real.hasDerivAt_mul_log (ne_of_gt h4)
  have hS : HasDerivAt (fun t : ℝ => (1 - t) * Real.log (1 - t))
      ((Real.log (1 - c) + 1) * (-1)) c := by
    have h := HasDerivAt.comp c
      (Real.hasDerivAt_mul_log (x := 1 - c) (by linarith : (1 : ℝ) - c ≠ 0)) hs
    rwa [Function.comp_def] at h
  have hU : HasDerivAt (fun t : ℝ => (t - a / 2) * Real.log (t - a / 2))
      ((Real.log (c - a / 2) + 1) * 1) c := by
    have h := HasDerivAt.comp c (Real.hasDerivAt_mul_log (x := c - a / 2) (ne_of_gt h1)) hu
    rwa [Function.comp_def] at h
  have hV : HasDerivAt (fun t : ℝ => (1 - t - a / 2) * Real.log (1 - t - a / 2))
      ((Real.log (1 - c - a / 2) + 1) * (-1)) c := by
    have h := HasDerivAt.comp c
      (Real.hasDerivAt_mul_log (x := 1 - c - a / 2) (ne_of_gt h2)) hv
    rwa [Function.comp_def] at h
  have hconst : HasDerivAt
      (fun _ : ℝ => a * Real.log 2 + (1 - a) * Real.log (1 - a)) 0 c :=
    hasDerivAt_const _ _
  have hsum := ((((hconst.add hC).add hS).sub hU).sub hV)
  rw [show Dpc a c = 0 + (Real.log c + 1) + (Real.log (1 - c) + 1) * (-1)
      - (Real.log (c - a / 2) + 1) * 1
      - (Real.log (1 - c - a / 2) + 1) * (-1) from by unfold Dpc; ring]
  convert hsum using 1
  funext t
  simp only [Pi.add_apply, Pi.sub_apply]
  unfold piEval
  ring

/-! ## The centered bound

`π` on a box, expanded about any point of that box, with the two remainders
carried by mean values of the partial derivatives *taken inside the box*.
This is what the interval verifier calls `p_taylor_upper`, and iteration 36
measured that it is the one bound the cover cannot do without.

The route is two applications of the one-dimensional mean value theorem, not a
multivariate one: hold `c` fixed and move `a`, then hold `am` fixed and move
`c`.  Both intermediate points stay in the box, which is all the enclosure
needs. -/

/-- One-dimensional mean value theorem in the form the expansion wants: a
difference equals a derivative *somewhere in the enclosing interval* times the
step, with no assumption about which endpoint is larger. -/
private lemma mvt_between {f f' : ℝ → ℝ} {lo hi u v : ℝ}
    (hu : u ∈ Set.Icc lo hi) (hv : v ∈ Set.Icc lo hi)
    (hd : ∀ x ∈ Set.Icc lo hi, HasDerivAt f (f' x) x) :
    ∃ ξ ∈ Set.Icc lo hi, f v - f u = f' ξ * (v - u) := by
  rcases lt_trichotomy u v with hlt | heq | hgt
  · have hsub : Set.Icc u v ⊆ Set.Icc lo hi := Set.Icc_subset_Icc hu.1 hv.2
    have hcont : ContinuousOn f (Set.Icc u v) := fun x hx =>
      ((hd x (hsub hx)).continuousAt).continuousWithinAt
    obtain ⟨ξ, hξ, hslope⟩ := exists_hasDerivAt_eq_slope f f' hlt hcont
      (fun x hx => hd x (hsub (Set.Ioo_subset_Icc_self hx)))
    refine ⟨ξ, hsub (Set.Ioo_subset_Icc_self hξ), ?_⟩
    have hne : v - u ≠ 0 := by intro h; apply absurd hlt; linarith [sub_eq_zero.mp h]
    rw [hslope]
    field_simp
  · exact ⟨u, hu, by rw [heq]; ring⟩
  · have hsub : Set.Icc v u ⊆ Set.Icc lo hi := Set.Icc_subset_Icc hv.1 hu.2
    have hcont : ContinuousOn f (Set.Icc v u) := fun x hx =>
      ((hd x (hsub hx)).continuousAt).continuousWithinAt
    obtain ⟨ξ, hξ, hslope⟩ := exists_hasDerivAt_eq_slope f f' hgt hcont
      (fun x hx => hd x (hsub (Set.Ioo_subset_Icc_self hx)))
    refine ⟨ξ, hsub (Set.Ioo_subset_Icc_self hξ), ?_⟩
    have hne : u - v ≠ 0 := by intro h; apply absurd hgt; linarith [sub_eq_zero.mp h]
    rw [hslope]
    field_simp
    ring

/-- **The centered expansion of `π` on a box.**

The five hypotheses are exactly the verifier's feasibility guard: the box must
sit strictly inside the region where all five logarithms have positive
arguments. -/
theorem piEval_centered {alo ahi clo chi : ℝ}
    (h1 : 0 < clo - ahi / 2) (h2 : 0 < 1 - chi - ahi / 2) (h3 : ahi < 1)
    (h4 : 0 < clo) (h5 : chi < 1)
    {a c am cm : ℝ} (ha : a ∈ Set.Icc alo ahi) (hc : c ∈ Set.Icc clo chi)
    (ham : am ∈ Set.Icc alo ahi) (hcm : cm ∈ Set.Icc clo chi) :
    ∃ ξ ∈ Set.Icc alo ahi, ∃ η ∈ Set.Icc clo chi,
      piEval a c = piEval am cm + Dpa ξ c * (a - am) + Dpc am η * (c - cm) := by
  have hda : ∀ t ∈ Set.Icc alo ahi, HasDerivAt (fun s => piEval s c) (Dpa t c) t := by
    intro t ht
    exact hasDerivAt_piEval_fst (a := t) (c := c)
      (by linarith [ht.2, hc.1]) (by linarith [ht.2, hc.2]) (by linarith [ht.2])
  obtain ⟨ξ, hξ, hA⟩ := mvt_between ham ha hda
  have hdc : ∀ s ∈ Set.Icc clo chi,
      HasDerivAt (fun t => piEval am t) (Dpc am s) s := by
    intro s hs
    exact hasDerivAt_piEval_snd (a := am) (c := s)
      (by linarith [ham.2, hs.1]) (by linarith [ham.2, hs.2])
      (by linarith [hs.1]) (by linarith [hs.2])
  obtain ⟨η, hη, hC⟩ := mvt_between hcm hc hdc
  exact ⟨ξ, hξ, η, hη, by linarith [hA, hC]⟩

/-! ## The Golay generating function -/

/-- `G(u) = ∑_w A_w u^w` for the extended binary Golay code. -/
noncomputable def GolayG (u : ℝ) : ℝ :=
  1 + 759 * u ^ 8 + 2576 * u ^ 12 + 759 * u ^ 16 + u ^ 24

/-- The constants come from the code, not from the paper: this is
`Golay.enumerator_eq` evaluated at `u`. -/
lemma GolayG_eq_sum (u : ℝ) :
    GolayG u = ∑ m : Fin 12 → Bool, u ^ (Spin.Structured.Golay.W m) := by
  have h := congrArg (Polynomial.eval₂ (Nat.castRingHom ℝ) u)
    Spin.Structured.Golay.enumerator_eq
  simp only [Spin.Structured.enumerator, Polynomial.eval₂_finsetSum, Polynomial.eval₂_pow,
    Polynomial.eval₂_X, Polynomial.eval₂_add, Polynomial.eval₂_one, Polynomial.eval₂_mul,
    Polynomial.eval₂_C] at h
  rw [GolayG, h]
  norm_num

lemma one_le_GolayG {u : ℝ} (hu : 0 < u) : 1 ≤ GolayG u := by
  unfold GolayG
  have h8 : (0 : ℝ) < u ^ 8 := by positivity
  have h12 : (0 : ℝ) < u ^ 12 := by positivity
  have h16 : (0 : ℝ) < u ^ 16 := by positivity
  have h24 : (0 : ℝ) < u ^ 24 := by positivity
  linarith

lemma GolayG_pos {u : ℝ} (hu : 0 < u) : 0 < GolayG u :=
  lt_of_lt_of_le zero_lt_one (one_le_GolayG hu)

lemma pow_le_GolayG {u : ℝ} (hu : 0 < u) : u ^ 24 ≤ GolayG u := by
  unfold GolayG
  have h8 : (0 : ℝ) < u ^ 8 := by positivity
  have h12 : (0 : ℝ) < u ^ 12 := by positivity
  have h16 : (0 : ℝ) < u ^ 16 := by positivity
  linarith

/-! ## The Golay exponent `g` -/

/-- The objective whose infimum is `g`. -/
noncomputable def gObj (u a : ℝ) : ℝ := Real.log (GolayG u) / 24 - a * Real.log u

/-- `g(a) := inf_{u > 0} { (1/24) ln G(u) - a ln u }`. -/
noncomputable def gBA (a : ℝ) : ℝ := sInf ((fun u => gObj u a) '' Set.Ioi 0)

/-- **The objective has a floor on `[0,1]`.**

`G(u) ≥ 1` handles `u ≤ 1` and `G(u) ≥ u^24` handles `u ≥ 1`.  Without this
there is no infimum bound to use: for `a > 1` the objective really does run to
`-∞`. -/
theorem gObj_nonneg {u a : ℝ} (hu : 0 < u) (ha0 : 0 ≤ a) (ha1 : a ≤ 1) :
    0 ≤ gObj u a := by
  unfold gObj
  rcases le_or_gt u 1 with h | h
  · have hlog : Real.log u ≤ 0 := Real.log_nonpos hu.le h
    have hG : 0 ≤ Real.log (GolayG u) := Real.log_nonneg (one_le_GolayG hu)
    nlinarith
  · have hlog : 0 < Real.log u := Real.log_pos h
    have hG : 24 * Real.log u ≤ Real.log (GolayG u) := by
      have h1 : Real.log (u ^ 24) ≤ Real.log (GolayG u) :=
        Real.log_le_log (by positivity) (pow_le_GolayG hu)
      rwa [Real.log_pow] at h1

    nlinarith

lemma gBA_bddBelow {a : ℝ} (ha0 : 0 ≤ a) (ha1 : a ≤ 1) :
    BddBelow ((fun u => gObj u a) '' Set.Ioi 0) := by
  refine ⟨0, ?_⟩
  rintro y ⟨u, hu, rfl⟩
  exact gObj_nonneg hu ha0 ha1

/-- **The witness bound.**  Any `u > 0` gives an upper bound on `g(a)`; this is
what makes the certificate's single witness per box enough. -/
theorem gBA_le {a u : ℝ} (ha0 : 0 ≤ a) (ha1 : a ≤ 1) (hu : 0 < u) :
    gBA a ≤ gObj u a :=
  csInf_le (gBA_bddBelow ha0 ha1) ⟨u, hu, rfl⟩

lemma gBA_nonneg {a : ℝ} (ha0 : 0 ≤ a) (ha1 : a ≤ 1) : 0 ≤ gBA a := by
  refine le_csInf ⟨gObj 1 a, ⟨1, by norm_num, rfl⟩⟩ ?_
  rintro y ⟨u, hu, rfl⟩
  exact gObj_nonneg hu ha0 ha1

end Spin.Numeric
