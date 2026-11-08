import SpinCodes.Turbo.ColumnAdjoint
import SpinCodes.Turbo.PacketLocalMoment
import Mathlib.LinearAlgebra.Matrix.Dual
import Mathlib.FieldTheory.Finite.GaloisField

/-!
# A binary adjoint for field multiplication

The field basis is explicit as a linear equivalence. The adjoint is the
transpose of the multiplication matrix in that basis. This identifies a
mathematical instance of the hypotheses of `ColumnAdjoint`; the concrete
C++ tower representation still requires a separate basis check.
-/

noncomputable section

namespace Spin.Turbo

open scoped Matrix

abbrev Bits32 := Fin 32 → ZMod 2
abbrev Field32 := GaloisField 2 32

variable {F : Type*} [Field F] [Algebra (ZMod 2) F]

/-- Multiplication by `c`, expressed in the chosen binary coordinates. -/
def coordinateMul (e : F ≃ₗ[ZMod 2] Bits32) (c : F) :
    Bits32 →ₗ[ZMod 2] Bits32 :=
  e.toLinearMap.comp ((LinearMap.mulLeft (ZMod 2) c).comp
    e.symm.toLinearMap)

theorem coordinateMul_apply (e : F ≃ₗ[ZMod 2] Bits32)
    (c : F) (y : Bits32) :
    coordinateMul e c y = e (c * e.symm y) := rfl

/-- The binary transpose of coordinate multiplication. -/
def coordinateAdjoint (e : F ≃ₗ[ZMod 2] Bits32)
    (c : F) (x : Bits32) : Bits32 :=
  (LinearMap.toMatrix' (coordinateMul e c)).transpose *ᵥ x

/-- The standard binary dot product as an additive character. -/
def binaryPair (x : Bits32) : Bits32 →+ ZMod 2 where
  toFun y := x ⬝ᵥ y
  map_zero' := dotProduct_zero x
  map_add' y z := dotProduct_add x y z

theorem binaryPair_nondegenerate (x : Bits32)
    (h : ∀ y : Bits32, binaryPair x y = 0) : x = 0 := by
  apply dotProduct_eq_zero
  exact h

theorem binaryPair_zero_left (y : Bits32) : binaryPair 0 y = 0 := by
  simp [binaryPair]

theorem coordinateAdjoint_pairing (e : F ≃ₗ[ZMod 2] Bits32)
    (c : F) (x y : Bits32) :
    binaryPair (coordinateAdjoint e c x) y =
      binaryPair x (e (c * e.symm y)) := by
  let M := LinearMap.toMatrix' (coordinateMul e c)
  change (M.transpose *ᵥ x) ⬝ᵥ y = x ⬝ᵥ e (c * e.symm y)
  calc
    _ = y ⬝ᵥ M.transpose *ᵥ x := dotProduct_comm _ _
    _ = x ⬝ᵥ M *ᵥ y := Matrix.dotProduct_transpose_mulVec M y x
    _ = _ := by rw [LinearMap.toMatrix'_mulVec, coordinateMul_apply]

noncomputable instance : Fintype Field32 := Fintype.ofFinite Field32
noncomputable instance : DecidableEq Field32 := Classical.decEq Field32

/-- Choose a binary basis for the mathematical 32-bit field. The uniform
action theorem below holds for every such choice. -/
def field32Coordinates : Field32 ≃ₗ[ZMod 2] Bits32 :=
  LinearEquiv.ofFinrankEq Field32 Bits32 (by
    rw [GaloisField.finrank 2 (by decide), Module.finrank_fin_fun])

instance : Nonempty {x : Bits32 // x ≠ 0} := by
  refine ⟨⟨fun i => if i = 0 then 1 else 0, ?_⟩⟩
  intro hz
  have h := congrFun hz (0 : Fin 32)
  norm_num at h

/-- The ideal uniform nonzero GF(2^32) coefficient gives a uniform
nonzero binary column under the transposed multiplication map. -/
theorem field32_uniform_adjoint (x : Bits32) (hx : x ≠ 0) :
    (FinPMF.uniform {c : Field32 // c ≠ 0}).map
      (nonzeroAdjointEquiv field32Coordinates.toAddEquiv binaryPair
        binaryPair_nondegenerate binaryPair_zero_left
        (coordinateAdjoint field32Coordinates)
        (coordinateAdjoint_pairing field32Coordinates) x hx) =
      FinPMF.uniform {y : Bits32 // y ≠ 0} :=
  uniform_nonzero_adjoint field32Coordinates.toAddEquiv binaryPair
    binaryPair_nondegenerate binaryPair_zero_left
    (coordinateAdjoint field32Coordinates)
    (coordinateAdjoint_pairing field32Coordinates) x hx

/-- Regroup the 32 binary coordinates into eight consecutive four-bit
packets. This is a coordinate relabeling and preserves zero. -/
private def packetIndex : Fin 8 × Fin 4 ≃ Fin 32 := finProdFinEquiv

def bits32PacketEquiv : Bits32 ≃ PacketColumn where
  toFun x i j := x (packetIndex (i, j))
  invFun p k := p (packetIndex.symm k).1
    (packetIndex.symm k).2
  left_inv x := by
    funext k
    change x (packetIndex (packetIndex.symm k)) = x k
    rw [Equiv.apply_symm_apply]
  right_inv p := by
    funext i j
    change p (packetIndex.symm (packetIndex (i, j))).1
      (packetIndex.symm (packetIndex (i, j))).2 = p i j
    rw [Equiv.symm_apply_apply]

theorem bits32PacketEquiv_zero : bits32PacketEquiv 0 = 0 := by
  funext i j
  change (0 : Bits32) (packetIndex (i, j)) = 0
  rfl

def nonzeroBits32PacketEquiv :
    {x : Bits32 // x ≠ 0} ≃ {y : PacketColumn // y ≠ 0} where
  toFun x := ⟨bits32PacketEquiv x.val, by
    intro h
    exact x.property (bits32PacketEquiv.injective
      (h.trans bits32PacketEquiv_zero.symm))⟩
  invFun y := ⟨bits32PacketEquiv.symm y.val, by
    intro h
    apply y.property
    calc
      y.val = bits32PacketEquiv (bits32PacketEquiv.symm y.val) :=
        (bits32PacketEquiv.apply_symm_apply _).symm
      _ = bits32PacketEquiv 0 := by rw [h]
      _ = 0 := bits32PacketEquiv_zero⟩
  left_inv x := Subtype.ext (bits32PacketEquiv.symm_apply_apply _)
  right_inv y := Subtype.ext (bits32PacketEquiv.apply_symm_apply _)

/-- For any fixed nonzero column, the sampled field multiplier produces a
uniform nonzero eight-packet column. This is the value-level distribution
needed before conditioning on a packet support. -/
theorem field32_uniform_packet_column (x : Bits32) (hx : x ≠ 0) :
    (FinPMF.uniform {c : Field32 // c ≠ 0}).map
      ((nonzeroAdjointEquiv field32Coordinates.toAddEquiv binaryPair
        binaryPair_nondegenerate binaryPair_zero_left
        (coordinateAdjoint field32Coordinates)
        (coordinateAdjoint_pairing field32Coordinates) x hx).trans
          nonzeroBits32PacketEquiv) =
      FinPMF.uniform {y : PacketColumn // y ≠ 0} := by
  let e := (nonzeroAdjointEquiv field32Coordinates.toAddEquiv binaryPair
    binaryPair_nondegenerate binaryPair_zero_left
    (coordinateAdjoint field32Coordinates)
    (coordinateAdjoint_pairing field32Coordinates) x hx).trans
      nonzeroBits32PacketEquiv
  change (FinPMF.uniform {c : Field32 // c ≠ 0}).map e = _
  apply FinPMF.ext
  funext y
  rw [FinPMF.map_equiv_apply]
  change 1 / (Fintype.card {c : Field32 // c ≠ 0} : ℝ) =
    1 / (Fintype.card {y : PacketColumn // y ≠ 0} : ℝ)
  rw [Fintype.card_congr e]

/-- The value-level law also factors across any fixed family of active
columns. The coefficient tuple is one shared setup for all messages; the
independence used here is only across its column coordinates. -/
theorem field32_uniform_active_packet_columns {n : ℕ}
    (x : Fin n → Bits32) (hx : ∀ i, x i ≠ 0) :
    (piPMF (fun _ : Fin n =>
      FinPMF.uniform {c : Field32 // c ≠ 0})).map
        (fun coeff i =>
          ((nonzeroAdjointEquiv field32Coordinates.toAddEquiv binaryPair
            binaryPair_nondegenerate binaryPair_zero_left
            (coordinateAdjoint field32Coordinates)
            (coordinateAdjoint_pairing field32Coordinates)
            (x i) (hx i)).trans nonzeroBits32PacketEquiv) (coeff i)) =
      piPMF (fun _ : Fin n =>
        FinPMF.uniform {y : PacketColumn // y ≠ 0}) := by
  rw [piPMF_map]
  simp_rw [field32_uniform_packet_column]

/-- The mathematical sampled adjoint column has the exact packet-weight
generating function used by the finite Turbo certificate. -/
theorem field32_adjoint_packet_gf (x : Bits32) (hx : x ≠ 0) (z : ℝ) :
    (FinPMF.uniform {c : Field32 // c ≠ 0}).expect
      (fun c => z ^ (packetColumnSupport
        (bits32PacketEquiv (coordinateAdjoint field32Coordinates c.val x))).card) =
      ((1 + 15 * z) ^ 8 - 1) / (2 ^ 32 - 1 : ℝ) := by
  let f := nonzeroAdjointEquiv field32Coordinates.toAddEquiv binaryPair
    binaryPair_nondegenerate binaryPair_zero_left
    (coordinateAdjoint field32Coordinates)
    (coordinateAdjoint_pairing field32Coordinates) x hx
  let g := nonzeroBits32PacketEquiv
  have hdist : (FinPMF.uniform {c : Field32 // c ≠ 0}).map (f.trans g) =
      FinPMF.uniform {y : PacketColumn // y ≠ 0} :=
    field32_uniform_packet_column x hx
  have hmap : (FinPMF.uniform {c : Field32 // c ≠ 0}).expect
      (fun c => z ^ (packetColumnSupport ((f.trans g c).val)).card) =
      ((FinPMF.uniform {c : Field32 // c ≠ 0}).map (f.trans g)).expect
        (fun y => z ^ (packetColumnSupport y.val).card) := by
    exact (FinPMF.expect_map
      (FinPMF.uniform {c : Field32 // c ≠ 0})
      (fun c => f.trans g c)
      (fun y : {y : PacketColumn // y ≠ 0} =>
        z ^ (packetColumnSupport y.val).card)).symm
  have hresult : (FinPMF.uniform {c : Field32 // c ≠ 0}).expect
      (fun c => z ^ (packetColumnSupport ((f.trans g c).val)).card) =
      ((1 + 15 * z) ^ 8 - 1) / (2 ^ 32 - 1 : ℝ) := by
    rw [hmap, hdist, uniform_nonzero_packet_gf]
  have hfg (c : {c : Field32 // c ≠ 0}) :
      ((f.trans g) c).val =
        bits32PacketEquiv (coordinateAdjoint field32Coordinates c.val x) := by
    rfl
  simpa only [hfg] using hresult

/-- For `h` fixed active columns, independent ideal field coefficients
give the `h`th power of the one-column packet generating function. -/
theorem field32_active_columns_gf (h : ℕ)
    (x : Fin h → Bits32) (hx : ∀ i, x i ≠ 0) (z : ℝ) :
    (piPMF (fun _ : Fin h =>
      FinPMF.uniform {c : Field32 // c ≠ 0})).expect
      (fun coeff => z ^ (∑ i : Fin h,
        (packetColumnSupport (bits32PacketEquiv
          (coordinateAdjoint field32Coordinates (coeff i).val (x i)))).card)) =
      (((1 + 15 * z) ^ 8 - 1) / (2 ^ 32 - 1 : ℝ)) ^ h := by
  have hpow (coeff : Fin h → {c : Field32 // c ≠ 0}) :
      z ^ (∑ i : Fin h,
        (packetColumnSupport (bits32PacketEquiv
          (coordinateAdjoint field32Coordinates (coeff i).val (x i)))).card) =
      ∏ i : Fin h, z ^ (packetColumnSupport (bits32PacketEquiv
        (coordinateAdjoint field32Coordinates (coeff i).val (x i)))).card :=
    (Finset.prod_pow_eq_pow_sum Finset.univ
      (fun i : Fin h => (packetColumnSupport (bits32PacketEquiv
        (coordinateAdjoint field32Coordinates (coeff i).val (x i)))).card) z).symm
  simp_rw [hpow]
  calc
    _ = ∏ i : Fin h,
        (FinPMF.uniform {c : Field32 // c ≠ 0}).expect
          (fun c => z ^ (packetColumnSupport (bits32PacketEquiv
            (coordinateAdjoint field32Coordinates c.val (x i)))).card) :=
      piPMF_expect_prod
        (fun _ : Fin h => FinPMF.uniform {c : Field32 // c ≠ 0})
        (fun i c => z ^ (packetColumnSupport (bits32PacketEquiv
          (coordinateAdjoint field32Coordinates c.val (x i)))).card)
    _ = _ := by simp [field32_adjoint_packet_gf, hx, div_pow]

#print axioms coordinateAdjoint_pairing
#print axioms field32_uniform_adjoint
#print axioms field32_uniform_packet_column
#print axioms field32_uniform_active_packet_columns
#print axioms field32_adjoint_packet_gf
#print axioms field32_active_columns_gf

end Spin.Turbo
