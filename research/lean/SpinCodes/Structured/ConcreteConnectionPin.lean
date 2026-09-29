import SpinCodes.Structured.FiberFrozen
import SpinCodes.Structured.FiberNumericsSample
import SpinCodes.Structured.ConcreteShells
import SpinCodes.Structured.SparseMomentBridge

/-! Pins for numerical soundness, actual weight shells, and the supported
nonnegative transfer invariant. Spectrum and one-step premises stay visible. -/

noncomputable section
namespace Spin.Structured.ConcreteConnectionPin
open ConcreteMaps

example (j w : ℕ) : FiberNumerics.kraw j w = Spin.krawtchouk 128 j w :=
  FiberNumerics.kraw_eq j w

example (j : Fin 129) : (SparsePolynomial.weightN j).kernel = FiberNumerics.Data.kernels.getD j 0 :=
  FiberNumerics.kernel_frozen_eq j

example (j : Fin 129) : (SparsePolynomial.weightN j).cap = FiberNumerics.Data.caps.getD j 0 :=
  FiberNumerics.cap_frozen_eq j

example (hs : ∀ w : Fin 129, weightCounts CtransposeSet w = FiberNumerics.Data.spectrum.getD w 0)
    (q : Finset (Fin 19)) : (syndromeFiber 2 q).card ≤ 61 :=
  FiberNumerics.sample_cap_two hs q

example (hc : ∀ i : Fin 5, weightCounts Aset (shellWeight i) = shellCount i) :
    (Finset.univ : Finset (Fin 5)).biUnion weightShell = Spin.nonzeroStates 19 :=
  weightShell_covers hc

example (hc : ∀ i : Fin 5, weightCounts Aset (shellWeight i) = shellCount i) (i : Fin 5) :
    (((shellSystem hc).shell i).card : ℝ) = Spin.Imt.Occupation.Sparse.count i :=
  shellSystem_card hc i

example {sys : Spin.Imt.ShellSystem 19 5} {c : Spin.Imt.Coords 5}
    {μ : Finset (Fin 19) → ℝ} (h : Spin.Imt.LiveDominates sys c μ) : μ ∅ ≤ c.Z := h.2

example {sys : Spin.Imt.ShellSystem 19 5} {α : ℝ} (h0 : 0 < α) (h1 : α ≤ 1 / 10000)
    (step : (Finset (Fin 19) → ℝ) → (Finset (Fin 19) → ℝ))
    (hstep : ∀ (c : Spin.Imt.Coords 5) μ, c.Nonneg → Spin.Imt.LiveDominates sys c μ →
      Spin.Imt.LiveDominates sys
        ((Spin.Imt.Occupation.Sparse.numericalMatrix ((4 / 5) * α) (1 - (8 / 5) * α)).apply c)
        (step μ)) (R : ℕ) :
    ∑ q, (step^[R] (fun q => if q = ∅ then (1 : ℝ) else 0)) q ≤
      2048 * (1 - 96 * α) ^ R :=
  Spin.Imt.Occupation.Sparse.sparse_live_moment_of_step h0 h1 step hstep R

end Spin.Structured.ConcreteConnectionPin
