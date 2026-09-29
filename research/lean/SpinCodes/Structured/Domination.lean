/-
Assembling `lem:structured-route-domination`.

The paper's lemma is the composition of two mechanisms, neither of which
mentions the route:

*  conditioning `b` independent Bernoulli variables on their sum costs a
   factor of at most the reciprocal of the conditioning probability, which
   `types_lower_bound_at` bounds by `b+1` when the target is the mode;
*  a uniform shuffle inside each weight class turns the fibre-mass bound of
   `poissonBinom_le_binom` into a *pointwise* bound, and pointwise bounds
   multiply over independent coordinates.

This file proves those mechanisms and composes them.  What stays a hypothesis
is only the modelling: that the route law really is the conditioned product,
and really is uniform inside each weight class.  That is the "certified
interface" boundary — the combinatorics is proved, the identification of the
model is named.
-/
import SpinCodes.Structured.RouteDomination

namespace Spin

open Finset Filter

variable {Ω : Type*} [Fintype Ω]

/-! ## The domination calculus

`Dominates P Q c` is `P ≤ c · Q` pointwise.  It is exactly the right notion
because the lemma quantifies over *all* nonnegative `F`: a bound on
expectations for every nonnegative `F` is equivalent to a pointwise bound, so
nothing is lost by proving the stronger statement. -/

/-- `P ≤ c · Q` pointwise. -/
def Dominates (P Q : FinPMF Ω) (c : ℝ) : Prop := ∀ ω, P.p ω ≤ c * Q.p ω

/-- The conclusion of the paper's lemma: expectations of nonnegative functions
transfer. -/
theorem Dominates.expect_le {P Q : FinPMF Ω} {c : ℝ} (h : Dominates P Q c)
    {F : Ω → ℝ} (hF : ∀ ω, 0 ≤ F ω) : P.expect F ≤ c * Q.expect F := by
  unfold FinPMF.expect
  rw [Finset.mul_sum]
  refine Finset.sum_le_sum fun ω _ => ?_
  calc P.p ω * F ω ≤ c * Q.p ω * F ω := mul_le_mul_of_nonneg_right (h ω) (hF ω)
    _ = c * (Q.p ω * F ω) := by ring

/-- Costs compose. -/
theorem Dominates.trans {P R S : FinPMF Ω} {c₁ c₂ : ℝ} (h₁ : Dominates P R c₁)
    (h₂ : Dominates R S c₂) (hc₁ : 0 ≤ c₁) : Dominates P S (c₁ * c₂) := by
  intro ω
  calc P.p ω ≤ c₁ * R.p ω := h₁ ω
    _ ≤ c₁ * (c₂ * S.p ω) := mul_le_mul_of_nonneg_left (h₂ ω) hc₁
    _ = c₁ * c₂ * S.p ω := by ring

/-! ## The shuffle step

A uniform shuffle inside each weight class makes both laws constant on the
fibres of the weight statistic.  Once that holds, the fibre-mass bound of
`poissonBinom_le_binom` *is* a pointwise bound: both sides of the fibre
equation carry the same fibre cardinality, which cancels. -/

theorem dominates_of_fiber_uniform {κ : Type*} [DecidableEq κ]
    (P Q : FinPMF Ω) (S : Ω → κ) (c : ℝ)
    (hP : ∀ ω ω', S ω = S ω' → P.p ω = P.p ω')
    (hQ : ∀ ω ω', S ω = S ω' → Q.p ω = Q.p ω')
    (hmass : ∀ k, P.prob (fun ω => S ω = k) ≤ c * Q.prob (fun ω => S ω = k)) :
    Dominates P Q c := by
  intro ω
  have key : ∀ R : FinPMF Ω, (∀ a b, S a = S b → R.p a = R.p b) →
      R.prob (fun ω' => S ω' = S ω)
        = ((univ.filter (fun ω' => S ω' = S ω)).card : ℝ) * R.p ω := by
    intro R hR
    rw [FinPMF.prob, Finset.sum_congr rfl (fun ω' hω' =>
      hR ω' ω (Finset.mem_filter.mp hω').2), Finset.sum_const, nsmul_eq_mul]
  have hcard : 0 < (univ.filter (fun ω' => S ω' = S ω)).card :=
    Finset.card_pos.mpr ⟨ω, by simp⟩
  have hcR : (0 : ℝ) < ((univ.filter (fun ω' => S ω' = S ω)).card : ℝ) := by
    exact_mod_cast hcard
  have hk := hmass (S ω)
  rw [key P hP, key Q hQ] at hk
  refine le_of_mul_le_mul_left ?_ hcR
  calc ((univ.filter (fun ω' => S ω' = S ω)).card : ℝ) * P.p ω
      ≤ c * (((univ.filter (fun ω' => S ω' = S ω)).card : ℝ) * Q.p ω) := hk
    _ = ((univ.filter (fun ω' => S ω' = S ω)).card : ℝ) * (c * Q.p ω) := by ring

/-! ## Independent coordinates

The `b` region shuffles are independent, so their costs multiply. -/

/-- The independent product of finitely many laws. -/
def piPMF {ι : Type*} [Fintype ι] [DecidableEq ι] {α : ι → Type*}
    [∀ i, Fintype (α i)] (P : ∀ i, FinPMF (α i)) : FinPMF (∀ i, α i) where
  p := fun f => ∏ i, (P i).p (f i)
  nonneg := fun f => Finset.prod_nonneg fun i _ => (P i).nonneg (f i)
  total := by
    have h := Finset.prod_univ_sum (fun i => (univ : Finset (α i)))
      (fun i x => (P i).p x)
    rw [Fintype.piFinset_univ] at h
    rw [← h]
    simp [FinPMF.total]

@[simp] lemma piPMF_apply {ι : Type*} [Fintype ι] [DecidableEq ι] {α : ι → Type*}
    [∀ i, Fintype (α i)] (P : ∀ i, FinPMF (α i)) (f : ∀ i, α i) :
    (piPMF P).p f = ∏ i, (P i).p (f i) := rfl

/-- **Costs multiply over independent coordinates.** -/
theorem dominates_piPMF {ι : Type*} [Fintype ι] [DecidableEq ι] {α : ι → Type*}
    [∀ i, Fintype (α i)] (P Q : ∀ i, FinPMF (α i)) (c : ι → ℝ)
    (h : ∀ i, Dominates (P i) (Q i) (c i)) :
    Dominates (piPMF P) (piPMF Q) (∏ i, c i) := by
  intro f
  simp only [piPMF_apply, ← Finset.prod_mul_distrib]
  exact Finset.prod_le_prod₀ (fun i _ => (P i).nonneg (f i)) (fun i _ => h i (f i))

/-! ## The conditioning step

Conditioning costs the reciprocal of the conditioning probability, and when
the conditioning event is "the sum equals a mode" that reciprocal is at most
`b+1` — which is `inv_mode_le` applied to the binomial. -/

theorem dominates_condition (P : FinPMF Ω) (s : Ω → Prop) [DecidablePred s]
    (h : 0 < P.prob s) : Dominates (P.condition s h) P (1 / P.prob s) := by
  intro ω
  rw [FinPMF.condition_p_apply]
  by_cases hs : s ω
  · rw [if_pos hs, one_div, inv_mul_eq_div]
  · rw [if_neg hs]
    exact mul_nonneg (by positivity) (P.nonneg ω)

/-- The conditioning cost is at most `b+1` when the event has probability at
least `1/(b+1)`. -/
theorem dominates_condition_le (P : FinPMF Ω) (s : Ω → Prop) [DecidablePred s]
    {b : ℕ} (hb : 1 / ((b : ℝ) + 1) ≤ P.prob s) (h : 0 < P.prob s) :
    Dominates (P.condition s h) P ((b : ℝ) + 1) := by
  intro ω
  refine (dominates_condition P s h ω).trans ?_
  refine mul_le_mul_of_nonneg_right ?_ (P.nonneg ω)
  rw [div_le_iff₀ h]
  have hb1 : (0 : ℝ) < (b : ℝ) + 1 := by positivity
  have hmul := mul_le_mul_of_nonneg_left hb hb1.le
  rwa [mul_one_div, div_self hb1.ne'] at hmul
/-! ## The iid Bernoulli law is the binomial

The dominating side of the lemma is `Ber(q)^L` on one region.  Its mass
depends only on the weight — so it *is* fibre-uniform, provably — and its
fibre mass is the binomial.  The routed side is not fibre-uniform; that is
precisely what the region shuffle has to supply, and why it stays a
hypothesis. -/

lemma poissonBinom_const_apply {L : ℕ} {q : ℝ} (hq0 : ∀ _ : Fin L, (0 : ℝ) ≤ q)
    (hq1 : ∀ _ : Fin L, q ≤ 1) (T : Finset (Fin L)) :
    (poissonBinom (fun _ => q) hq0 hq1).p T
      = q ^ T.card * (1 - q) ^ (L - T.card) := by
  rw [poissonBinom_apply, Finset.prod_const, Finset.prod_const,
    Finset.card_univ_diff, Fintype.card_fin]

lemma poissonBinom_const_prob {L : ℕ} (k : ℕ) {q : ℝ}
    (hq0 : ∀ _ : Fin L, (0 : ℝ) ≤ q) (hq1 : ∀ _ : Fin L, q ≤ 1) :
    (poissonBinom (fun _ => q) hq0 hq1).prob (fun T => T.card = k)
      = (L.choose k : ℝ) * q ^ k * (1 - q) ^ (L - k) := by
  have hfilter : (univ.filter (fun T : Finset (Fin L) => T.card = k))
      = Finset.powersetCard k univ := by
    ext T; simp [Finset.mem_powersetCard]
  rw [FinPMF.prob, hfilter,
    Finset.sum_congr rfl (fun T hT => by
      rw [poissonBinom_const_apply, (Finset.mem_powersetCard.mp hT).2]),
    Finset.sum_const, Finset.card_powersetCard, Finset.card_univ, Fintype.card_fin,
    nsmul_eq_mul]
  ring

/-! ## The weight bound at every weight

`poissonBinom_le_binom` needs `0 < k < L`, because the optimal `z` degenerates
at the ends.  Both ends are AM-GM directly: at `k = 0` the mass is
`∏ (1-pᵢ) ≤ (1-q)^L`, and at `k = L` it is `∏ pᵢ ≤ q^L`. -/

theorem poissonBinom_le_binom_all {L k : ℕ} (hL : 0 < L) (hkL : k ≤ L)
    (p : Fin L → ℝ) (hp0 : ∀ i, 0 ≤ p i) (hp1 : ∀ i, p i ≤ 1)
    {q : ℝ} (hqval : (∑ i, p i) / (L : ℝ) = q) (hq0 : 0 < q) (hq1 : q < 1) :
    (poissonBinom p hp0 hp1).prob (fun T => T.card = k)
      ≤ ((L : ℝ) + 1) * ((L.choose k : ℝ) * q ^ k * (1 - q) ^ (L - k)) := by
  have hL1 : (1 : ℝ) ≤ (L : ℝ) + 1 := by
    have : (0 : ℝ) ≤ (L : ℝ) := Nat.cast_nonneg L
    linarith
  rcases Nat.eq_zero_or_pos k with rfl | hk0
  · have hmass : (poissonBinom p hp0 hp1).prob (fun T => T.card = 0)
        = ∏ i, (1 - p i) := by
      have hfil : (univ.filter (fun T : Finset (Fin L) => T.card = 0)) = {∅} := by
        ext T; simp [Finset.card_eq_zero]
      rw [FinPMF.prob, hfil, Finset.sum_singleton, poissonBinom_apply]
      simp
    have hamgm : ∏ i, (1 - p i) ≤ (1 - q) ^ L := by
      have h := prod_one_add_mul_le hL p (-1) hp0 hp1 (le_refl _)
      rw [hqval] at h
      calc ∏ i, (1 - p i) = ∏ i, (1 + p i * (-1)) :=
            Finset.prod_congr rfl fun i _ => by ring
        _ ≤ (1 + q * (-1)) ^ L := h
        _ = (1 - q) ^ L := by ring_nf
    have hnn : (0 : ℝ) ≤ (1 - q) ^ L := pow_nonneg (by linarith) L
    rw [hmass]
    simp only [Nat.choose_zero_right, Nat.cast_one, pow_zero, one_mul, Nat.sub_zero]
    nlinarith
  rcases eq_or_lt_of_le hkL with rfl | hklt
  · have hmass : (poissonBinom p hp0 hp1).prob (fun T => T.card = k)
        = ∏ i, p i := by
      have hfil : (univ.filter (fun T : Finset (Fin k) => T.card = k)) = {univ} := by
        ext T
        simp only [Finset.mem_filter, Finset.mem_univ, true_and, Finset.mem_singleton]
        exact ⟨fun h => Finset.eq_univ_of_card T (by simpa using h),
               fun h => by simp [h]⟩
      rw [FinPMF.prob, hfil, Finset.sum_singleton, poissonBinom_apply]
      simp
    have hamgm : ∏ i, p i ≤ q ^ k := by
      have h := prod_le_pow_avg hL p hp0
      rwa [hqval] at h
    have hnn : (0 : ℝ) ≤ q ^ k := pow_nonneg (le_of_lt hq0) k
    rw [hmass]
    simp only [Nat.choose_self, Nat.cast_one, one_mul, Nat.sub_self, pow_zero, mul_one]
    nlinarith
  · exact poissonBinom_le_binom hL hk0 hklt p hp0 hp1 hqval hq0 hq1

/-! ## The region step -/

/-- **One region shuffle costs `L+1`.**  The two hypotheses are what the
uniform region shuffle supplies: the shuffled law is constant inside each
weight class, and permuting within weight classes leaves the weight
distribution alone. -/
theorem dominates_region {L : ℕ} (hL : 0 < L) (R : FinPMF (Finset (Fin L)))
    (p : Fin L → ℝ) (hp0 : ∀ i, 0 ≤ p i) (hp1 : ∀ i, p i ≤ 1)
    {q : ℝ} (hqval : (∑ i, p i) / (L : ℝ) = q) (hq0 : 0 < q) (hq1 : q < 1)
    (hRuniform : ∀ T T' : Finset (Fin L), T.card = T'.card → R.p T = R.p T')
    (hRmass : ∀ k, R.prob (fun T => T.card = k)
      = (poissonBinom p hp0 hp1).prob (fun T => T.card = k)) :
    Dominates R (poissonBinom (fun _ => q) (fun _ => hq0.le) (fun _ => hq1.le))
      ((L : ℝ) + 1) := by
  refine dominates_of_fiber_uniform R _ (fun T => T.card) _ hRuniform ?_ ?_
  · intro T T' h
    rw [poissonBinom_const_apply, poissonBinom_const_apply, h]
  · intro k
    rw [hRmass k, poissonBinom_const_prob]
    rcases le_or_gt k L with hk | hk
    · exact poissonBinom_le_binom_all hL hk p hp0 hp1 hqval hq0 hq1
    · have hzero : (poissonBinom p hp0 hp1).prob (fun T => T.card = k) = 0 := by
        have hfil : (univ.filter (fun T : Finset (Fin L) => T.card = k)) = ∅ := by
          ext T
          simp only [Finset.mem_filter, Finset.mem_univ, true_and,
            Finset.notMem_empty, iff_false]
          intro h
          have hle := Finset.card_le_univ T
          rw [Fintype.card_fin] at hle
          omega
        rw [FinPMF.prob, hfil, Finset.sum_empty]
      rw [hzero, Nat.choose_eq_zero_of_lt hk]
      simp

/-! ## Assembling `eq:structured-route-domination` -/

/-- The `Q` row conditionings, costing `(b+1)^Q` together. -/
theorem dominates_rows {Q : ℕ} {α : Fin Q → Type*} [∀ j, Fintype (α j)]
    (P : ∀ j, FinPMF (α j)) (s : ∀ j, α j → Prop) [∀ j, DecidablePred (s j)]
    (hpos : ∀ j, 0 < (P j).prob (s j)) {b : ℕ}
    (hmode : ∀ j, 1 / ((b : ℝ) + 1) ≤ (P j).prob (s j)) :
    Dominates (piPMF (fun j => (P j).condition (s j) (hpos j))) (piPMF P)
      (((b : ℝ) + 1) ^ Q) := by
  have h := dominates_piPMF (fun j => (P j).condition (s j) (hpos j)) P
    (fun _ => (b : ℝ) + 1)
    (fun j => dominates_condition_le (P j) (s j) (hmode j) (hpos j))
  simpa using h

/-- The `b` region shuffles, costing `(L+1)^b` together. -/
theorem dominates_regions {b : ℕ} {α : Fin b → Type*} [∀ i, Fintype (α i)]
    (P Q : ∀ i, FinPMF (α i)) {L : ℕ}
    (h : ∀ i, Dominates (P i) (Q i) ((L : ℝ) + 1)) :
    Dominates (piPMF P) (piPMF Q) (((L : ℝ) + 1) ^ b) := by
  simpa using dominates_piPMF P Q (fun _ => (L : ℝ) + 1) h

/-- **`eq:structured-route-domination`.**  For every nonnegative `F`,

    E_route[F] ≤ (b+1)^Q (L+1)^b E_Ber[F].

`hrow` is `dominates_rows`, `hregion` is `dominates_regions` fed by
`dominates_region`; the composition is `Dominates.trans`. -/
theorem route_domination {Route Mid Ber : FinPMF Ω} {b L Q : ℕ}
    (hrow : Dominates Route Mid (((b : ℝ) + 1) ^ Q))
    (hregion : Dominates Mid Ber (((L : ℝ) + 1) ^ b))
    {F : Ω → ℝ} (hF : ∀ ω, 0 ≤ F ω) :
    Route.expect F ≤ ((b : ℝ) + 1) ^ Q * ((L : ℝ) + 1) ^ b * Ber.expect F :=
  (hrow.trans hregion (by positivity)).expect_le hF
/-! ## The factor is `exp(o(N))`

    ln((b+1)^Q (L+1)^b) = Q ln(b+1) + b ln(L+1) = o(Lb).

Only `Q ≤ L` is needed for the upper bound; the paper's extra hypothesis that
`Q/L` is bounded below is what makes `Lb` the right normaliser elsewhere, not
what makes this estimate true. -/

private lemma log_succ_div_tendsto :
    Tendsto (fun n : ℕ => Real.log ((n : ℝ) + 1) / (n : ℝ)) atTop (nhds 0) := by
  have h1 : Tendsto (fun x : ℝ => Real.log x / x) atTop (nhds 0) := by
    simpa using Real.isLittleO_log_id_atTop.tendsto_div_nhds_zero
  have h2 : Tendsto (fun n : ℕ => ((n : ℝ) + 1)) atTop atTop :=
    tendsto_atTop_add_const_right _ 1 tendsto_natCast_atTop_atTop
  have h4 : Tendsto (fun n : ℕ => 2 * (Real.log ((n : ℝ) + 1) / ((n : ℝ) + 1)))
      atTop (nhds 0) := by
    simpa using (h1.comp h2).const_mul 2
  refine squeeze_zero_norm' ?_ h4
  filter_upwards [eventually_ge_atTop 1] with n hn
  have hn1 : (1 : ℝ) ≤ (n : ℝ) := by exact_mod_cast hn
  have hnp : (0 : ℝ) < (n : ℝ) := by linarith
  have hn1p : (0 : ℝ) < (n : ℝ) + 1 := by linarith
  have hlog : 0 ≤ Real.log ((n : ℝ) + 1) := Real.log_nonneg (by linarith)
  rw [Real.norm_eq_abs, abs_of_nonneg (div_nonneg hlog (le_of_lt hnp)),
    mul_div_assoc', div_le_div_iff₀ hnp hn1p]
  nlinarith [mul_nonneg hlog (by linarith : (0 : ℝ) ≤ (n : ℝ) - 1)]

/-- **The multiplicative factor is `exp(o(N))`.**  For every `ε > 0` there is a
threshold beyond which `ln((b+1)^Q (L+1)^b) ≤ ε·Lb`. -/
theorem route_factor_isLittleO {ε : ℝ} (hε : 0 < ε) :
    ∃ M : ℕ, 0 < M ∧ ∀ L b Q : ℕ, M ≤ L → M ≤ b → Q ≤ L →
      Real.log (((b : ℝ) + 1) ^ Q * ((L : ℝ) + 1) ^ b) ≤ ε * ((L : ℝ) * b) := by
  obtain ⟨M₀, hM₀⟩ := Metric.tendsto_atTop.mp log_succ_div_tendsto (ε / 2) (by linarith)
  refine ⟨max M₀ 1, lt_of_lt_of_le Nat.one_pos (le_max_right M₀ 1), ?_⟩
  intro L b Q hL hb hQ
  have hLR : (1 : ℝ) ≤ (L : ℝ) := by
    exact_mod_cast le_trans (le_max_right M₀ 1) hL
  have hbR : (1 : ℝ) ≤ (b : ℝ) := by
    exact_mod_cast le_trans (le_max_right M₀ 1) hb
  have hLpos : (0 : ℝ) < (L : ℝ) := by linarith
  have hbpos : (0 : ℝ) < (b : ℝ) := by linarith
  have hlb := hM₀ b (le_trans (le_max_left M₀ 1) hb)
  have hlL := hM₀ L (le_trans (le_max_left M₀ 1) hL)
  rw [Real.dist_eq, sub_zero] at hlb hlL
  have hb2 : Real.log ((b : ℝ) + 1) / (b : ℝ) < ε / 2 := lt_of_abs_lt hlb
  have hL2 : Real.log ((L : ℝ) + 1) / (L : ℝ) < ε / 2 := lt_of_abs_lt hlL
  have hlogb : 0 ≤ Real.log ((b : ℝ) + 1) := Real.log_nonneg (by linarith)
  have hQR : (Q : ℝ) ≤ (L : ℝ) := by exact_mod_cast hQ
  have hexp : Real.log (((b : ℝ) + 1) ^ Q * ((L : ℝ) + 1) ^ b)
      = (Q : ℝ) * Real.log ((b : ℝ) + 1) + (b : ℝ) * Real.log ((L : ℝ) + 1) := by
    rw [Real.log_mul (by positivity) (by positivity), Real.log_pow, Real.log_pow]
  rw [hexp]
  have t1 : (Q : ℝ) * Real.log ((b : ℝ) + 1) ≤ (L : ℝ) * (b : ℝ) * (ε / 2) := by
    calc (Q : ℝ) * Real.log ((b : ℝ) + 1)
        ≤ (L : ℝ) * Real.log ((b : ℝ) + 1) := mul_le_mul_of_nonneg_right hQR hlogb
      _ ≤ (L : ℝ) * ((ε / 2) * (b : ℝ)) := by
          refine mul_le_mul_of_nonneg_left ?_ (by linarith)
          rw [← div_le_iff₀ hbpos]
          exact hb2.le
      _ = (L : ℝ) * (b : ℝ) * (ε / 2) := by ring
  have t2 : (b : ℝ) * Real.log ((L : ℝ) + 1) ≤ (L : ℝ) * (b : ℝ) * (ε / 2) := by
    calc (b : ℝ) * Real.log ((L : ℝ) + 1)
        ≤ (b : ℝ) * ((ε / 2) * (L : ℝ)) := by
          refine mul_le_mul_of_nonneg_left ?_ (by linarith)
          rw [← div_le_iff₀ hLpos]
          exact hL2.le
      _ = (L : ℝ) * (b : ℝ) * (ε / 2) := by ring
  nlinarith [t1, t2]

end Spin
