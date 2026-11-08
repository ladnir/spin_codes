import SpinCodes.Turbo.ActualSetup22
import SpinCodes.UniformMixer.ScalableCodeProperties

/-!
# Algebraic properties of the actual 22-state Turbo encoder

This follows the old `ConcreteNativeLinearCode` proof architecture: first
connect the vector encoder to the weight experiment, then prove additivity
and injectivity for every fixed setup. The unchanged inner uses the checked
Scalable vector recurrence. Only the outer and packet reshaping change.
No finite distance certificate or concrete C++ basis equivalence is assumed.
-/

noncomputable section
namespace Spin.Turbo

open Spin.UniformMixer Spin.Structured.ConcreteRoute

theorem coordinateAdjoint_injective
    (c : {c : Field32 // c ≠ 0}) :
    Function.Injective (coordinateAdjoint field32Coordinates c.val) := by
  intro x y h
  have hz : coordinateAdjoint field32Coordinates c.val (x - y) = 0 := by
    simp only [coordinateAdjoint, Matrix.mulVec_sub] at *
    rw [h, sub_self]
  by_contra hxy
  exact adjoint_ne_zero field32Coordinates.toAddEquiv binaryPair
    binaryPair_nondegenerate binaryPair_zero_left
    (coordinateAdjoint field32Coordinates)
    (coordinateAdjoint_pairing field32Coordinates) c.property
    (sub_ne_zero.mpr hxy) hz

theorem sampledRS8RowWord_injective
    (points : Fin 8 → Field16) (hp : Function.Injective points)
    (coeff : Fin 8 → {c : Field32 // c ≠ 0}) :
    Function.Injective (fun rows => sampledRS8RowWord points rows coeff) := by
  intro x y h
  apply rs8Columns_injective points hp
  funext j
  apply rs8ColumnBits.injective
  apply coordinateAdjoint_injective (coeff j)
  apply bits32PacketEquiv.injective
  exact congrFun (packetColumnsWordEquiv.injective h) j

theorem sampledRS8RowWord_add
    (points : Fin 8 → Field16)
    (coeff : Fin 8 → {c : Field32 // c ≠ 0})
    (x y : RS8Group Field16) :
    sampledRS8RowWord points (x + y) coeff =
      sampledRS8RowWord points x coeff + sampledRS8RowWord points y coeff := by
  have hcol (j : Fin 8) : rs8Column points (x + y) j =
      rs8Column points x j + rs8Column points y j := by
    funext r
    exact Polynomial.eval_add
  have hbits (j : Fin 8) : rs8Bits points (x + y) j =
      rs8Bits points x j + rs8Bits points y j := by
    rw [rs8Bits, hcol, map_add]
    rfl
  funext k i
  simp only [sampledRS8RowWord, packetColumnsWordEquiv,
    Equiv.coe_fn_mk,
    sampledRS8PacketColumns, sampledPacketColumn, hbits,
    coordinateAdjoint, Matrix.mulVec_add, bits32PacketEquiv, Pi.add_apply]

theorem routePackets_injective {L b : ℕ} {α : Type*}
    [Zero α] [DecidableEq α] (seed : Seeds L b) :
    Function.Injective (fun rows : Fin L → Fin b → α => routePackets rows seed) := by
  intro x y h
  funext i k
  let j := (seed.1 i).symm k
  have hh := congrFun (congrFun h j) ((seed.2 j).symm i)
  simpa only [routePackets, Equiv.apply_symm_apply, j] using hh

theorem flatInput_injective (R : ℕ) :
    Function.Injective (flatInput (R := R)) := by
  intro x y h
  funext j i k
  obtain ⟨p, rfl⟩ :=
    (finProdFinEquiv : Fin R × Fin 16 ≃ Fin (R * 16)).surjective i
  have hh := congrFun (congrFun h
    ((finProdFinEquiv : Fin 64 × Fin R ≃ Fin (64 * R)) (j, p.1)))
    ((finProdFinEquiv : Fin 16 × Fin 4 ≃ Fin 64) (p.2, k))
  simpa only [flatInput, stepInput, Equiv.symm_apply_apply] using hh

def turboCodeword22 (R : ℕ) (points : Fin 8 → Field16)
    (seed : TurboOuterSeed22 R × RouteInnerSeed22 R)
    (msg : TurboMessage22 R) : Fin (64 * R) → Input64 :=
  sampledOutput (flatInput (routePackets (turboOuterWord22 points seed.1 msg)
    seed.2.1)) seed.2.2 0

theorem turboCodeword22_weight (R : ℕ) (points : Fin 8 → Field16)
    (seed : TurboOuterSeed22 R × RouteInnerSeed22 R)
    (msg : TurboMessage22 R) :
    sampledWordWeight (turboCodeword22 R points seed msg) =
      outputWeight22 (turboOuterWord22 points seed.1 msg) seed.2 := by
  exact sampledOutput_weight _ _ _

theorem turboCodeword22_injective (R : ℕ) (points : Fin 8 → Field16)
    (hp : Function.Injective points)
    (seed : TurboOuterSeed22 R × RouteInnerSeed22 R) :
    Function.Injective (turboCodeword22 R points seed) := by
  intro x y h
  have hinner := sampledOutput_injective seed.2.2 0 h
  have hroute := flatInput_injective R hinner
  have hout := routePackets_injective seed.2.1 hroute
  funext i
  exact sampledRS8RowWord_injective points hp (seed.1 i) (congrFun hout i)

theorem turboCodeword22_add (R : ℕ) (points : Fin 8 → Field16)
    (seed : TurboOuterSeed22 R × RouteInnerSeed22 R)
    (x y : TurboMessage22 R) :
    turboCodeword22 R points seed (x + y) =
      turboCodeword22 R points seed x + turboCodeword22 R points seed y := by
  have hout : turboOuterWord22 points seed.1 (x + y) =
      turboOuterWord22 points seed.1 x + turboOuterWord22 points seed.1 y := by
    funext i
    exact sampledRS8RowWord_add points (seed.1 i) (x i) (y i)
  unfold turboCodeword22
  rw [hout]
  exact (show sampledOutput _ _ (0 + 0) = _ from sampledOutput_add _ _ _ 0 0)

theorem turboCodeword22_zero (R : ℕ) (points : Fin 8 → Field16)
    (seed : TurboOuterSeed22 R × RouteInnerSeed22 R) :
    turboCodeword22 R points seed 0 = 0 := by
  have h := turboCodeword22_add R points seed 0 0
  simp only [add_zero] at h
  exact add_left_cancel (h.symm.trans (add_zero _).symm)

def turboEncoderLinearMap22 (R : ℕ) (points : Fin 8 → Field16)
    (seed : TurboOuterSeed22 R × RouteInnerSeed22 R) :
    TurboMessage22 R →ₗ[ZMod 2] (Fin (64 * R) → Input64) where
  toFun := turboCodeword22 R points seed
  map_add' := turboCodeword22_add R points seed
  map_smul' c x := by
    have hc : c = 0 ∨ c = 1 := (by decide : ∀ c : ZMod 2, c = 0 ∨ c = 1) c
    rcases hc with rfl | rfl <;> simp [turboCodeword22_zero]

#print axioms turboCodeword22_weight
#print axioms turboCodeword22_injective
#print axioms turboEncoderLinearMap22

end Spin.Turbo
