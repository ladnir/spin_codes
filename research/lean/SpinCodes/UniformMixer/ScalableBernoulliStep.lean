import SpinCodes.UniformMixer.ScalableRouteDomination

/-!
# The iid Bernoulli reference step

The positive-occupation argument compares the routed word with iid
Bernoulli input bits. This module defines the exact averaged 64-bit
step kernel and identifies the full iid moment with its matrix power.
-/

noncomputable section

namespace Spin.UniformMixer

open Spin.Structured.Routing

def supportInput (S : Finset (Fin 64)) : Input64 :=
  fun j => if j ∈ S then 1 else 0

def bernoulliInputLaw (β : ℝ) (hβ0 : 0 ≤ β) (hβ1 : β ≤ 1) :
    FinPMF Input64 :=
  (poissonBinom (fun _ : Fin 64 => β) (fun _ => hβ0)
    (fun _ => hβ1)).map supportInput

def bernoulliStepKernel (β z : ℝ)
    (hβ0 : 0 ≤ β) (hβ1 : β ≤ 1) :
    Matrix State22 State22 ℝ :=
  matrixExpect22 (bernoulliInputLaw β hβ0 hβ1)
    (fun x => Matrix.of (stepKernel z x))

theorem bernoulliStepKernel_support (β z : ℝ)
    (hβ0 : 0 ≤ β) (hβ1 : β ≤ 1) (q r : State22) :
    bernoulliStepKernel β z hβ0 hβ1 q r =
      (poissonBinom (fun _ : Fin 64 => β) (fun _ => hβ0)
        (fun _ => hβ1)).expect
        (fun S => stepKernel z (supportInput S) q r) := by
  change (bernoulliInputLaw β hβ0 hβ1).expect
    (fun x => stepKernel z x q r) = _
  rw [bernoulliInputLaw, FinPMF.expect_map]

theorem inputListKernel_matrix (z : ℝ) (xs : List Input64) :
    Matrix.of (inputListKernel z xs) =
      (xs.map (fun x => Matrix.of (stepKernel z x))).prod := by
  induction xs with
  | nil =>
      ext q r
      simp [inputListKernel, identity, Matrix.one_apply]
  | cons x tail ih =>
      change Matrix.of (compose (stepKernel z x) (inputListKernel z tail)) =
        Matrix.of (stepKernel z x) *
          (tail.map (fun y => Matrix.of (stepKernel z y))).prod
      rw [← ih]
      exact compose_eq_matrixMul
        (Matrix.of (stepKernel z x)) (Matrix.of (inputListKernel z tail))

theorem bernoulli_input_kernel_power (R : ℕ) (β z : ℝ)
    (hβ0 : 0 ≤ β) (hβ1 : β ≤ 1) :
    matrixExpect22
      (piPMF (fun _ : Fin R => bernoulliInputLaw β hβ0 hβ1))
      (fun xs => Matrix.of (inputListKernel z (List.ofFn xs))) =
    (bernoulliStepKernel β z hβ0 hβ1) ^ R := by
  simp only [inputListKernel_matrix, List.map_ofFn]
  change matrixExpect22
      (piPMF (fun _ : Fin R => bernoulliInputLaw β hβ0 hβ1))
      (fun xs => (List.ofFn (fun i : Fin R =>
        Matrix.of (stepKernel z (xs i)))).prod) =
    (bernoulliStepKernel β z hβ0 hβ1) ^ R
  rw [matrixExpect22_pi_product
    (fun _ : Fin R => bernoulliInputLaw β hβ0 hβ1)
    (fun _ x => Matrix.of (stepKernel z x))]
  simp only [bernoulliStepKernel, List.ofFn_const, List.prod_replicate]

end Spin.UniformMixer
