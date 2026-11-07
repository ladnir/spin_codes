import SpinCodes.Turbo.ColumnCoordinates
import Mathlib.Data.ZMod.Basic

/-!
# Four-bit packet support inside a randomized column

One 32-bit column contains eight four-bit packets. For a specified
support `s`, each active packet has fifteen possible nonzero labels;
all other packets are zero. The exact count is therefore `15^|s|`.
-/

noncomputable section

namespace Spin.Turbo

abbrev Packet4 := Fin 4 → ZMod 2
abbrev PacketColumn := Fin 8 → Packet4

def packetColumnSupport (x : PacketColumn) : Finset (Fin 8) :=
  Finset.univ.filter (fun i => x i ≠ 0)

noncomputable instance (s : Finset (Fin 8)) :
    Fintype {x : PacketColumn // packetColumnSupport x = s} :=
  Fintype.ofFinite _

/-- A column with exact packet support `s` is the choice of one nonzero
four-bit label at each member of `s`. -/
noncomputable def packetColumnSupportEquiv (s : Finset (Fin 8)) :
    {x : PacketColumn // packetColumnSupport x = s} ≃
      (s → {p : Packet4 // p ≠ 0}) := by
  classical
  refine ⟨(fun x i => ⟨x.val i.val, ?_⟩),
    (fun f => ⟨(fun i => if h : i ∈ s then (f ⟨i,h⟩).val else 0), ?_⟩),
    ?_, ?_⟩
  · have hi : i.val ∈ packetColumnSupport x.val := by
      rw [x.property]
      exact i.property
    exact (Finset.mem_filter.mp hi).2
  · ext i
    simp only [packetColumnSupport, Finset.mem_filter, Finset.mem_univ,
      true_and]
    by_cases hi : i ∈ s
    · simp [hi, (f ⟨i,hi⟩).property]
    · simp [hi]
  · intro x
    apply Subtype.ext
    funext i
    by_cases hi : i ∈ s
    · simp [hi]
    · have hz : x.val i = 0 := by
        by_contra hn
        have hm : i ∈ packetColumnSupport x.val :=
          Finset.mem_filter.mpr ⟨Finset.mem_univ _, hn⟩
        rw [x.property] at hm
        exact hi hm
      simp [hi, hz]
  · intro f
    funext i
    apply Subtype.ext
    simp

theorem packet4_nonzero_card :
    Fintype.card {p : Packet4 // p ≠ 0} = 15 := by
  classical
  have hzero : Fintype.card {p : Packet4 // p = 0} = 1 :=
    Fintype.card_subtype_eq 0
  have hcomp := Fintype.card_subtype_compl (fun p : Packet4 => p = 0)
  change Fintype.card {p : Packet4 // p ≠ 0} =
    Fintype.card Packet4 - Fintype.card {p : Packet4 // p = 0} at hcomp
  rw [hcomp, hzero]
  norm_num [Packet4]

/-- Exactly `15^|s|` columns have packet support `s`. -/
theorem packetColumnSupport_card (s : Finset (Fin 8)) :
    Fintype.card {x : PacketColumn // packetColumnSupport x = s} =
      15 ^ s.card := by
  rw [Fintype.card_congr (packetColumnSupportEquiv s)]
  simp

/-- Among all 32-bit columns, exactly `binom(8,v) * 15^v` have `v`
nonzero four-bit packets. This includes the zero column at `v=0`. -/
theorem packetColumnSupport_weight_card (v : ℕ) :
    (Finset.univ.filter (fun x : PacketColumn =>
      (packetColumnSupport x).card = v)).card =
      (Nat.choose 8 v) * 15 ^ v := by
  classical
  let supports := (Finset.univ : Finset (Fin 8)).powersetCard v
  have hsum := Finset.sum_card_fiberwise_eq_card_filter
    (Finset.univ : Finset PacketColumn) supports packetColumnSupport
  have hterm (s : Finset (Fin 8)) :
      (Finset.univ.filter (fun x : PacketColumn =>
        packetColumnSupport x = s)).card = 15 ^ s.card := by
    rw [← Fintype.card_subtype]
    exact packetColumnSupport_card s
  calc
    _ = (Finset.univ.filter (fun x : PacketColumn =>
        packetColumnSupport x ∈ supports)).card := by
      apply congrArg Finset.card
      apply Finset.filter_congr
      intro x _
      simp [supports, Finset.mem_powersetCard]
    _ = ∑ s ∈ supports,
        (Finset.univ.filter (fun x : PacketColumn =>
          packetColumnSupport x = s)).card := hsum.symm
    _ = ∑ s ∈ supports, 15 ^ v := by
      apply Finset.sum_congr rfl
      intro s hs
      rw [hterm, (Finset.mem_powersetCard.mp hs).2]
    _ = (Nat.choose 8 v) * 15 ^ v := by
      simp [supports, Finset.card_powersetCard]

theorem packetColumn_card : Fintype.card PacketColumn = 2 ^ 32 := by
  simp [PacketColumn, Packet4]

private theorem nonzero_card_generic
    {B : Type*} [Zero B] [Fintype B] [DecidableEq B] :
    Fintype.card {x : B // x ≠ 0} = Fintype.card B - 1 := by
  classical
  have hzero : Fintype.card {x : B // x = 0} = 1 :=
    Fintype.card_subtype_eq 0
  have hcomp := Fintype.card_subtype_compl (fun x : B => x = 0)
  change Fintype.card {x : B // x ≠ 0} =
    Fintype.card B - Fintype.card {x : B // x = 0} at hcomp
  rw [hcomp, hzero]

theorem packetColumn_nonzero_card :
    Fintype.card {x : PacketColumn // x ≠ 0} = 2 ^ 32 - 1 := by
  rw [nonzero_card_generic, packetColumn_card]

instance : Nonempty {x : PacketColumn // x ≠ 0} := by
  refine ⟨⟨(fun _ _ => (1 : ZMod 2)), ?_⟩⟩
  intro h
  have h00 := congrFun (congrFun h 0) 0
  norm_num at h00

/-- The packet-weight distribution of a uniformly sampled nonzero
32-bit column. The zero-weight case is excluded by the premise. -/
theorem uniform_nonzero_packet_weight_prob (v : ℕ) (hv : 1 ≤ v) :
    (FinPMF.uniform {x : PacketColumn // x ≠ 0}).prob
      (fun x => (packetColumnSupport x.val).card = v) =
      ((Nat.choose 8 v * 15 ^ v : ℕ) : ℝ) / (2 ^ 32 - 1 : ℝ) := by
  classical
  let nz := {x : PacketColumn // x ≠ 0}
  let snz : Finset nz :=
    Finset.univ.filter (fun x => (packetColumnSupport x.val).card = v)
  let sall : Finset PacketColumn :=
    Finset.univ.filter (fun x => (packetColumnSupport x).card = v)
  have himage : snz.image Subtype.val = sall := by
    ext x
    constructor
    · intro hx
      obtain ⟨y, hy, rfl⟩ := Finset.mem_image.mp hx
      exact Finset.mem_filter.mpr ⟨Finset.mem_univ _,
        (Finset.mem_filter.mp hy).2⟩
    · intro hx
      have hwt := (Finset.mem_filter.mp hx).2
      have hnz : x ≠ 0 := by
        intro hzero
        subst x
        have hcard : (packetColumnSupport (0 : PacketColumn)).card = 0 := by
          simp [packetColumnSupport]
        omega
      exact Finset.mem_image.mpr
        ⟨⟨x, hnz⟩, Finset.mem_filter.mpr ⟨Finset.mem_univ _, hwt⟩, rfl⟩
  have hcount : snz.card = Nat.choose 8 v * 15 ^ v := by
    calc
      snz.card = (snz.image Subtype.val).card :=
        (Finset.card_image_of_injective snz Subtype.val_injective).symm
      _ = sall.card := by rw [himage]
      _ = _ := packetColumnSupport_weight_card v
  rw [FinPMF.uniform_prob]
  simp only [packetColumn_nonzero_card]
  rw [show (Finset.univ.filter
      (fun x : nz => (packetColumnSupport x.val).card = v)).card =
      Nat.choose 8 v * 15 ^ v from hcount]
  norm_num

#print axioms packetColumnSupportEquiv
#print axioms packet4_nonzero_card
#print axioms packetColumnSupport_card
#print axioms packetColumnSupport_weight_card
#print axioms packetColumn_card
#print axioms packetColumn_nonzero_card
#print axioms uniform_nonzero_packet_weight_prob

end Spin.Turbo
