import SpinCodes.Turbo.SampledRS8PacketLaw
import SpinCodes.Turbo.PacketShuffleValues

/-!
# From RS8 columns to one 64-packet row

The sampled RS8 group outputs eight columns of eight four-bit packets.
This coordinate bijection puts them in the row format consumed by the
packet route. It preserves each packet value and the packet count used
in the expected-spectrum theorem.
-/

noncomputable section

namespace Spin.Turbo

open Finset

def packetColumnsWordEquiv :
    (Fin 8 → PacketColumn) ≃ PacketWord64 where
  toFun cols k :=
    cols ((finProdFinEquiv : Fin 8 × Fin 8 ≃ Fin 64).symm k).1
      ((finProdFinEquiv : Fin 8 × Fin 8 ≃ Fin 64).symm k).2
  invFun word j i :=
    word ((finProdFinEquiv : Fin 8 × Fin 8 ≃ Fin 64) (j, i))
  left_inv cols := by
    funext j i
    simp
  right_inv word := by
    funext k
    simp

theorem packetColumnsWord_support_mem
    (cols : Fin 8 → PacketColumn) (j i : Fin 8) :
    ((finProdFinEquiv : Fin 8 × Fin 8 ≃ Fin 64) (j, i)) ∈
      packetWordSupport (packetColumnsWordEquiv cols) ↔
        i ∈ packetColumnSupport (cols j) := by
  simp only [packetWordSupport, packetColumnSupport,
    Finset.mem_filter, Finset.mem_univ, true_and]
  change (cols ((finProdFinEquiv : Fin 8 × Fin 8 ≃ Fin 64).symm
      ((finProdFinEquiv : Fin 8 × Fin 8 ≃ Fin 64) (j, i))).1
      ((finProdFinEquiv : Fin 8 × Fin 8 ≃ Fin 64).symm
        ((finProdFinEquiv : Fin 8 × Fin 8 ≃ Fin 64) (j, i))).2 ≠ 0) ↔
          cols j i ≠ 0
  simp

/-- Flattening columns neither drops packets nor changes their weight. -/
theorem packetColumnsWord_support_card
    (cols : Fin 8 → PacketColumn) :
    (packetWordSupport (packetColumnsWordEquiv cols)).card =
      ∑ j : Fin 8, (packetColumnSupport (cols j)).card := by
  classical
  have hc {n : ℕ} (S : Finset (Fin n)) :
      S.card = ∑ i : Fin n, if i ∈ S then 1 else 0 := by
    rw [← Finset.sum_filter]
    simp
  rw [hc]
  rw [← Fintype.sum_equiv
    (finProdFinEquiv : Fin 8 × Fin 8 ≃ Fin 64)
    (fun p : Fin 8 × Fin 8 =>
      if p.2 ∈ packetColumnSupport (cols p.1) then 1 else 0)
    (fun k : Fin 64 =>
      if k ∈ packetWordSupport (packetColumnsWordEquiv cols) then 1 else 0)
    (fun p => by
      rcases p with ⟨j, i⟩
      by_cases hi : i ∈ packetColumnSupport (cols j)
      · have hw := (packetColumnsWord_support_mem cols j i).mpr hi
        simp [hi, hw]
      · have hw : ((finProdFinEquiv : Fin 8 × Fin 8 ≃ Fin 64) (j, i)) ∉
            packetWordSupport (packetColumnsWordEquiv cols) :=
          fun h => hi ((packetColumnsWord_support_mem cols j i).mp h)
        simp [hi, hw])]
  rw [Fintype.sum_prod_type]
  exact Finset.sum_congr rfl fun j _ => (hc (packetColumnSupport (cols j))).symm

def sampledRS8RowWord (points : Fin 8 → Field16)
    (rows : RS8Group Field16)
    (coeff : Fin 8 → {c : Field32 // c ≠ 0}) : PacketWord64 :=
  packetColumnsWordEquiv (sampledRS8PacketColumns points rows coeff)

def sampledRS8RowLaw (points : Fin 8 → Field16)
    (rows : RS8Group Field16) : FinPMF PacketWord64 :=
  (sampledRS8PacketColumnsLaw points rows).map packetColumnsWordEquiv

theorem sampledRS8RowWord_weight
    (points : Fin 8 → Field16) (rows : RS8Group Field16)
    (coeff : Fin 8 → {c : Field32 // c ≠ 0}) :
    (packetWordSupport (sampledRS8RowWord points rows coeff)).card =
      rs8SampledPacketWeight points rows coeff := by
  exact packetColumnsWord_support_card _

theorem packetColumnsWord_support_profile_of_eq
    (y z : PacketWord64)
    (hs : packetWordSupport y = packetWordSupport z)
    (j : Fin 8) :
    packetColumnSupport ((packetColumnsWordEquiv.symm y) j) =
      packetColumnSupport ((packetColumnsWordEquiv.symm z) j) := by
  ext i
  have h := congrArg
    (fun S : Finset (Fin 64) =>
      ((finProdFinEquiv : Fin 8 × Fin 8 ≃ Fin 64) (j, i)) ∈ S) hs
  simp only [packetColumnSupport, Finset.mem_filter,
    Finset.mem_univ, true_and]
  simpa [packetWordSupport, packetColumnSupport,
    packetColumnsWordEquiv] using h

/-- The sampled RS8 row assigns equal probability to complete words with
the same packet support, even before the independent row shuffle. -/
theorem sampledRS8RowLaw_fiber_uniform
    (points : Fin 8 → Field16) (rows : RS8Group Field16)
    (y z : PacketWord64)
    (hs : packetWordSupport y = packetWordSupport z) :
    (sampledRS8RowLaw points rows).p y =
      (sampledRS8RowLaw points rows).p z := by
  rw [sampledRS8RowLaw, FinPMF.map_equiv_apply,
    FinPMF.map_equiv_apply]
  exact sampledRS8PacketColumnsLaw_fiber_uniform points rows
    (packetColumnsWordEquiv.symm y)
    (packetColumnsWordEquiv.symm z)
    (packetColumnsWord_support_profile_of_eq y z hs)

#print axioms packetColumnsWord_support_card
#print axioms sampledRS8RowWord_weight
#print axioms sampledRS8RowLaw_fiber_uniform

end Spin.Turbo
