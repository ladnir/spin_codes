import SpinCodes.Turbo.RS8ShuffledPlacement
import SpinCodes.Turbo.PacketInner22
import SpinCodes.Framework
import SpinCodes.FiniteMoment

/-!
# The actual sampled Turbo setup with the 22-state inner

One outer setup samples one coefficient tuple for every RS8 group. Those
coefficients are fixed across every message in the realized code.
The route and independent state-map updates form a separate setup factor.
This is the concrete instance of the old `Framework.Setup` first-moment
interface for the 22-state Turbo family.
-/

noncomputable section

namespace Spin.Turbo

open Spin.Structured.ConcreteRoute Spin.UniformMixer

abbrev TurboMessage22 (R : ℕ) :=
  Fin (R * 16) → RS8Group Field16

abbrev TurboOuterSeed22 (R : ℕ) :=
  Fin (R * 16) → (Fin 8 → {c : Field32 // c ≠ 0})

abbrev TurboOuterWord22 (R : ℕ) :=
  Fin (R * 16) → PacketWord64

def turboOuterLaw22 (R : ℕ) : FinPMF (TurboOuterSeed22 R) :=
  piPMF (fun _ : Fin (R * 16) =>
    piPMF (fun _ : Fin 8 =>
      FinPMF.uniform {c : Field32 // c ≠ 0}))

def turboOuterWord22 (points : Fin 8 → Field16)
    (seed : TurboOuterSeed22 R) (msg : TurboMessage22 R) :
    TurboOuterWord22 R :=
  fun i => sampledRS8RowWord points (msg i) (seed i)

/-- The setup law and bad-distance event use the realized encoder, with
one shared outer seed and independent route/inner setup. -/
def turboSetup22 (R : ℕ) (points : Fin 8 → Field16) :
    Spin.Setup (TurboMessage22 R) (TurboOuterWord22 R)
      (TurboOuterSeed22 R) (RouteInnerSeed22 R) where
  Pout := turboOuterLaw22 R
  Pin := routeInnerLaw22 R
  outer := turboOuterWord22 points
  innerWt := fun seed rows => outputWeight22 rows seed

/-- Exact first-moment identity for the actual sampled Turbo encoder.
No independent setup is sampled for individual messages in this sum. -/
theorem turboSetup22_expected_bad_count
    (R : ℕ) (points : Fin 8 → Field16) (d : ℕ) :
    (turboSetup22 R points).joint.expect
      (fun ω => ((turboSetup22 R points).Z d ω : ℝ)) =
    (turboOuterLaw22 R).expect
      (fun seed => ∑ msg ∈ Spin.nonzeroMsgs (TurboMessage22 R),
        (turboSetup22 R points).qd d
          (turboOuterWord22 points seed msg)) :=
  (turboSetup22 R points).expect_Z_eq d

/-- Markov converts the expected count into the probability that the
realized code has some nonzero message of weight at most `d`. -/
theorem turboSetup22_bad_prob_le_mean
    (R : ℕ) (points : Fin 8 → Field16) (d : ℕ) :
    (turboSetup22 R points).joint.prob
      (fun ω => 1 ≤ (turboSetup22 R points).Z d ω) ≤
    (turboSetup22 R points).joint.expect
      (fun ω => ((turboSetup22 R points).Z d ω : ℝ)) :=
  (turboSetup22 R points).prob_bad_le_expect_Z d

/-- The lower-tail probability for a fixed outer packet array is bounded
by the exact route-and-22-state kernel moment. Numerical certificate
operators must upper-bound this right-hand side. -/
theorem turboSetup22_qd_le_kernel_moment
    (R : ℕ) (points : Fin 8 → Field16)
    (rows : TurboOuterWord22 R) (d : ℕ) {z : ℝ}
    (hz0 : 0 < z) (hz1 : z ≤ 1) :
    (turboSetup22 R points).qd d rows ≤
      (seedLaw (R * 16) 64).expect (fun route =>
        ∑ q : State22,
          inputListKernel z
            (List.ofFn (flatInput (routePackets rows route))) 0 q) /
        z ^ d := by
  change (routeInnerLaw22 R).prob
    (fun seed => outputWeight22 rows seed ≤ d) ≤ _
  rw [← routeInnerMoment22 rows z]
  exact FinPMF.prob_weight_le
    (routeInnerLaw22 R) (outputWeight22 rows) d hz0 hz1

#print axioms turboSetup22_expected_bad_count
#print axioms turboSetup22_bad_prob_le_mean
#print axioms turboSetup22_qd_le_kernel_moment

end Spin.Turbo
