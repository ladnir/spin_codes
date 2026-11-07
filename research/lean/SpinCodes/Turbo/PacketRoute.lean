import SpinCodes.Structured.ConcreteRoutePermutation

/-!
# The sampled packet route

The old Structured SPIN route already samples independent row and regional
permutations. Here the same seeds move complete packet values. Taking the
nonzero-packet support commutes with routing, so the old support-law theorem
applies without discarding the labels used by Turbo SPIN's inner encoder.
-/

noncomputable section

namespace Spin.Turbo

open Spin.Structured Spin.Structured.ConcreteRoute
  Spin.Structured.Routing Finset

variable {L b : ℕ} {α : Type*} [Zero α] [DecidableEq α]

/-- The nonzero packet positions in each outer group. -/
def packetSupport (rows : Fin L → Fin b → α) :
    Fin L → Finset (Fin b) :=
  fun i => Finset.univ.filter (fun j => rows i j ≠ 0)

/-- Row permutation, transpose, then regional permutation, retaining values. -/
def routePackets (rows : Fin L → Fin b → α)
    (seed : Seeds L b) : Fin b → Fin L → α :=
  fun j i => rows (seed.2 j i) (seed.1 (seed.2 j i) j)

/-- The exact support of the value-carrying route is the old checked route. -/
theorem routePackets_support (rows : Fin L → Fin b → α)
    (seed : Seeds L b) (j : Fin b) :
    Finset.univ.filter (fun i => routePackets rows seed j i ≠ 0) =
      routeEval (packetSupport rows) seed j := by
  ext i
  simp only [Finset.mem_filter, Finset.mem_univ, true_and]
  simp only [routeEval, mem_shuffleSupport, mem_transpose]
  simp only [packetSupport]
  simp only [Finset.mem_filter, Finset.mem_univ, true_and]
  rfl

/-- The packet-support experiment has exactly the old permutation law. -/
theorem routePackets_support_law (rows : Fin L → Fin b → α) :
    (seedLaw L b).map
      (fun seed j => Finset.univ.filter
        (fun i => routePackets rows seed j i ≠ 0)) =
      law (packetSupport rows) := by
  have h : (fun seed : Seeds L b =>
      fun j => Finset.univ.filter
        (fun i => routePackets rows seed j i ≠ 0)) =
      routeEval (packetSupport rows) := by
    funext seed j
    exact routePackets_support rows seed j
  rw [h]
  exact permutation_law_eq (packetSupport rows)

/-- Routing complete packets preserves the number of active packets. -/
theorem routePackets_totalSupport (rows : Fin L → Fin b → α)
    (seed : Seeds L b) :
    (∑ j : Fin b,
      (Finset.univ.filter (fun i => routePackets rows seed j i ≠ 0)).card) =
      ∑ i : Fin L, (packetSupport rows i).card := by
  calc
    _ = totalWeight (routeEval (packetSupport rows) seed) := by
      simp only [totalWeight, routePackets_support]
    _ = totalWeight (packetSupport rows) :=
      routeEval_totalWeight (packetSupport rows) seed
    _ = _ := rfl

#print axioms routePackets_support
#print axioms routePackets_support_law
#print axioms routePackets_totalSupport

end Spin.Turbo
