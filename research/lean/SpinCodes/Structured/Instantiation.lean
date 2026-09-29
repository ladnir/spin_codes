/-
A concrete Structured SPIN family, and the certificate constructor that
*proves* the first-moment linkage instead of assuming it.

The message and codeword spaces are pinned to `𝔽₂^K` and `𝔽₂^N`: the realized
code is binary linear of dimension `K_m` and length `N_m`.  The randomness
spaces stay abstract (carrying their own `Fintype`), because
`Setup.prob_bad_le_cond_sum` is uniform in them and committing to a
parametrization of the BA permutations, route shuffles and IMT transvections
before the outer and inner maps are formalized would be guessing.  That
parametrization arrives with the maps themselves.
-/
import SpinCodes.Framework
import SpinCodes.Distance
import SpinCodes.Structured.Certificate
import SpinCodes.Structured.Schedule

set_option linter.unusedSectionVars false

namespace Spin.Structured

open Finset Filter

attribute [local instance] Classical.propDecidable

/-- A sampled Structured SPIN family: one realized encoder per length index.

`S m` carries the outer law (the two BA permutations) and the route/inner law
(block shuffles, region shuffles, IMT transvections), which are independent by
construction of `Setup`. -/
structure Family where
  /-- Message dimension `K_m`. -/
  K : ℕ → ℕ
  /-- Direct length `N_m`. -/
  N : ℕ → ℕ
  /-- Outer length `L_m`. -/
  L : ℕ → ℕ
  /-- Outer randomness. -/
  Out : ℕ → Type
  /-- Route and inner randomness. -/
  Inr : ℕ → Type
  outFin : ∀ m, Fintype (Out m)
  inrFin : ∀ m, Fintype (Inr m)
  S : ∀ m, @Setup (Fin (K m) → ZMod 2) (Fin (N m) → ZMod 2) (Out m) (Inr m)
        inferInstance (outFin m) (inrFin m)
  /-- The selection event `𝒢_b`, determined by the outer setup alone. -/
  G : ∀ m, Out m → Prop
  Gpos : ∀ m, 0 < @FinPMF.prob _ (outFin m) (@Setup.Pout _ _ _ _ _ (outFin m) (inrFin m) (S m))
    (G m) (fun _ => Classical.propDecidable _)
  /-- Number of active outer blocks of the outer word of a message. -/
  occ : ∀ m, Out m → (Fin (K m) → ZMod 2) → ℕ
  occ_mem : ∀ m ω₁ x, x ∈ nonzeroMsgs (Fin (K m) → ZMod 2) →
    occ m ω₁ x ∈ Ico 1 (L m + 1)
  /-- The distance target `d_m = ⌊0.11 N_m⌋`. -/
  d : ℕ → ℕ

namespace Family

variable (F : Family)

/-- `Pr[∃ a nonzero message whose encoding has weight at most `d_m`]`.
By `Distance.one_le_Z_iff` this is exactly the bad event of the framework,
and by `Distance.exists_nonzero_light` it dominates
"`E_m` is not injective or `d_min(E_m) ≤ d_m`". -/
noncomputable def probBad (m : ℕ) : ℝ :=
  letI := F.outFin m
  letI := F.inrFin m
  (F.S m).joint.prob (fun ω => 1 ≤ (F.S m).Z (F.d m) ω)

/-- `Pr[¬𝒢_b]`, over the outer setup only. -/
noncomputable def probNotGood (m : ℕ) : ℝ :=
  letI := F.outFin m
  letI := F.inrFin m
  (F.S m).Pout.prob (fun ω₁ => ¬ F.G m ω₁)

/-- `𝔼[Z_{d,Q} ∣ 𝒢_b]`, conditioned on the outer-only selection event. -/
noncomputable def EZ (m Q : ℕ) : ℝ :=
  letI := F.outFin m
  letI := F.inrFin m
  (((F.S m).Pout.condition (F.G m) (F.Gpos m)).prod (F.S m).Pin).expect
    (fun ω => ((F.S m).ZQ (F.d m) (F.occ m) Q ω : ℝ))

lemma probBad_nonneg (m : ℕ) : 0 ≤ F.probBad m := by
  letI := F.outFin m
  letI := F.inrFin m
  exact FinPMF.prob_nonneg _ _

lemma EZ_nonneg (m Q : ℕ) : 0 ≤ F.EZ m Q := by
  letI := F.outFin m
  letI := F.inrFin m
  exact FinPMF.expect_nonneg _ fun _ => by positivity

/-- **The linkage, proved.**  `probBad_le` is no longer an assumption: it is
`Setup.prob_bad_le_cond_sum` applied to the family. -/
theorem probBad_le (m : ℕ) :
    F.probBad m ≤ F.probNotGood m + ∑ Q ∈ Ico 1 (F.L m + 1), F.EZ m Q := by
  letI := F.outFin m
  letI := F.inrFin m
  exact (F.S m).prob_bad_le_cond_sum (F.d m) (F.occ m) (F.G m) (F.Gpos m) (F.L m)
    (F.occ_mem m)

end Family

/-- Build a `ScalableCertificate` from a concrete family.  The three fields
`probBad_nonneg`, `EZ_nonneg` and `probBad_le` are supplied by proof; only the
schedule, the selection event's vanishing probability, and the three regime
envelopes remain as hypotheses. -/
noncomputable def ScalableCertificate.ofFamily (F : Family) (b : ℕ → ℕ)
    (N_eq : ∀ m, F.N m = F.L m * b m)
    (b_pos : ∀ m, 0 < b m)
    (b_tendsto : Tendsto b atTop atTop)
    (N_tendsto : Tendsto (fun m => (F.N m : ℝ)) atTop atTop)
    (L_le_N : ∀ m, (F.L m : ℝ) ≤ (F.N m : ℝ))
    (probNotGood_tendsto : Tendsto F.probNotGood atTop (nhds 0))
    (fixed : ∀ Q ∈ Ico 1 4096, Tendsto (fun m => F.EZ m Q) atTop (nhds 0))
    (cut : ℕ → ℕ)
    (cut_spec : ∀ᶠ m in atTop, 4095 ≤ cut m ∧ cut m ≤ F.L m)
    (sparse : ∀ m, ∀ Q ∈ Ico 4096 (cut m + 1),
      F.EZ m Q ≤ Real.exp (-(0.006 : ℝ) * Q * b m))
    (eta : ℝ) (eta_pos : 0 < eta)
    (denseErr : ℕ → ℝ)
    (denseErr_tendsto : Tendsto denseErr atTop (nhds 0))
    (dense : ∀ m, ∑ Q ∈ Ico (cut m + 1) (F.L m + 1), F.EZ m Q
      ≤ (F.L m : ℝ) * Real.exp (-eta * F.N m + denseErr m * F.N m)) :
    ScalableCertificate where
  L := F.L
  b := b
  N := F.N
  N_eq := N_eq
  b_pos := b_pos
  b_tendsto := b_tendsto
  N_tendsto := N_tendsto
  L_le_N := L_le_N
  EZ := F.EZ
  EZ_nonneg := F.EZ_nonneg
  probNotGood := F.probNotGood
  probNotGood_tendsto := probNotGood_tendsto
  fixed := fixed
  cut := cut
  cut_spec := cut_spec
  sparse := sparse
  eta := eta
  eta_pos := eta_pos
  denseErr := denseErr
  denseErr_tendsto := denseErr_tendsto
  dense := dense
  probBad := F.probBad
  probBad_nonneg := F.probBad_nonneg
  probBad_le := F.probBad_le

/-- **Distance claim for a concrete family.**  The conclusion is about the
family's own `probBad`, which is `Pr[∃ nonzero message of encoded weight
≤ d_m]` — not an abstract number supplied by a hypothesis. -/
theorem Family.distance_whp (F : Family) (b : ℕ → ℕ)
    (N_eq : ∀ m, F.N m = F.L m * b m) (b_pos : ∀ m, 0 < b m)
    (b_tendsto : Tendsto b atTop atTop)
    (N_tendsto : Tendsto (fun m => (F.N m : ℝ)) atTop atTop)
    (L_le_N : ∀ m, (F.L m : ℝ) ≤ (F.N m : ℝ))
    (probNotGood_tendsto : Tendsto F.probNotGood atTop (nhds 0))
    (fixed : ∀ Q ∈ Ico 1 4096, Tendsto (fun m => F.EZ m Q) atTop (nhds 0))
    (cut : ℕ → ℕ)
    (cut_spec : ∀ᶠ m in atTop, 4095 ≤ cut m ∧ cut m ≤ F.L m)
    (sparse : ∀ m, ∀ Q ∈ Ico 4096 (cut m + 1),
      F.EZ m Q ≤ Real.exp (-(0.006 : ℝ) * Q * b m))
    (eta : ℝ) (eta_pos : 0 < eta) (denseErr : ℕ → ℝ)
    (denseErr_tendsto : Tendsto denseErr atTop (nhds 0))
    (dense : ∀ m, ∑ Q ∈ Ico (cut m + 1) (F.L m + 1), F.EZ m Q
      ≤ (F.L m : ℝ) * Real.exp (-eta * F.N m + denseErr m * F.N m)) :
    Tendsto F.probBad atTop (nhds 0) :=
  (ScalableCertificate.ofFamily F b N_eq b_pos b_tendsto N_tendsto L_le_N
    probNotGood_tendsto fixed cut cut_spec sparse eta eta_pos denseErr
    denseErr_tendsto dense).distance_whp

/-- **A family on the native schedule.**

`eq:structured-native-schedule` supplies `N_eq`, `b_pos`, `b_tendsto`,
`N_tendsto` and `L_le_N`, and `Instantiation` supplies `probBad_le`,
`probBad_nonneg` and `EZ_nonneg`.  What is left is exactly the mathematics of
the construction: the selection event and the three regime envelopes. -/
noncomputable def ScalableCertificate.ofNativeFamily (F : Family)
    (hL : ∀ m, F.L m = Lsched m) (hN : ∀ m, F.N m = Nsched m)
    (probNotGood_tendsto : Tendsto F.probNotGood atTop (nhds 0))
    (fixed : ∀ Q ∈ Ico 1 4096, Tendsto (fun m => F.EZ m Q) atTop (nhds 0))
    (cut : ℕ → ℕ)
    (cut_spec : ∀ᶠ m in atTop, 4095 ≤ cut m ∧ cut m ≤ F.L m)
    (sparse : ∀ m, ∀ Q ∈ Ico 4096 (cut m + 1),
      F.EZ m Q ≤ Real.exp (-(0.006 : ℝ) * Q * bsched m))
    (eta : ℝ) (eta_pos : 0 < eta)
    (denseErr : ℕ → ℝ)
    (denseErr_tendsto : Tendsto denseErr atTop (nhds 0))
    (dense : ∀ m, ∑ Q ∈ Ico (cut m + 1) (F.L m + 1), F.EZ m Q
      ≤ (F.L m : ℝ) * Real.exp (-eta * F.N m + denseErr m * F.N m)) :
    ScalableCertificate :=
  ScalableCertificate.ofFamily F bsched
    (fun m => by rw [hN, hL, Nsched_eq])
    bsched_pos bsched_tendsto
    (by
      have : (fun m => (F.N m : ℝ)) = fun m => (Nsched m : ℝ) := by
        funext m; rw [hN]
      rw [this]; exact Nsched_tendsto)
    (fun m => by rw [hL, hN]; exact Lsched_le_Nsched m)
    probNotGood_tendsto fixed cut cut_spec sparse eta eta_pos denseErr
    denseErr_tendsto dense

/-- The distance claim for a native-schedule family. -/
theorem Family.distance_whp_native (F : Family)
    (hL : ∀ m, F.L m = Lsched m) (hN : ∀ m, F.N m = Nsched m)
    (probNotGood_tendsto : Tendsto F.probNotGood atTop (nhds 0))
    (fixed : ∀ Q ∈ Ico 1 4096, Tendsto (fun m => F.EZ m Q) atTop (nhds 0))
    (cut : ℕ → ℕ)
    (cut_spec : ∀ᶠ m in atTop, 4095 ≤ cut m ∧ cut m ≤ F.L m)
    (sparse : ∀ m, ∀ Q ∈ Ico 4096 (cut m + 1),
      F.EZ m Q ≤ Real.exp (-(0.006 : ℝ) * Q * bsched m))
    (eta : ℝ) (eta_pos : 0 < eta) (denseErr : ℕ → ℝ)
    (denseErr_tendsto : Tendsto denseErr atTop (nhds 0))
    (dense : ∀ m, ∑ Q ∈ Ico (cut m + 1) (F.L m + 1), F.EZ m Q
      ≤ (F.L m : ℝ) * Real.exp (-eta * F.N m + denseErr m * F.N m)) :
    Tendsto F.probBad atTop (nhds 0) :=
  (ScalableCertificate.ofNativeFamily F hL hN probNotGood_tendsto fixed cut
    cut_spec sparse eta eta_pos denseErr denseErr_tendsto dense).distance_whp

end Spin.Structured
