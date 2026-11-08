import SpinCodes.Structured.ConcreteOuterNativeLinearity

/-! Binary vector interfaces for the actual shared Golay--BAA outer. -/
noncomputable section
namespace Spin.UniformMixer
open Spin.Structured Spin.Structured.ConcreteOuter
  Spin.Structured.ConcreteBinaryEncoding

/-- Transport an established Boolean linear map to binary field vectors. -/
def liftBinaryLinearMap {a b : ℕ} (f : (Fin a → Bool) → (Fin b → Bool))
    (hz : f 0 = 0)
    (hadd : ∀ x y, f (fun i => xor (x i) (y i)) =
      fun i => xor (f x i) (f y i)) :
    (Fin a → ZMod 2) →ₗ[ZMod 2] (Fin b → ZMod 2) where
  toFun x := wordEquiv b (f ((wordEquiv a).symm x))
  map_add' x y := by
    rw [wordEquiv_symm_add, hadd]
    funext i
    exact bitEquiv_xor _ _
  map_smul' c x := by
    have hc : c = 0 ∨ c = 1 := (by decide : ∀ c : ZMod 2, c = 0 ∨ c = 1) c
    rcases hc with rfl | rfl <;> simp [hz]

/-- The fixed checked Golay generator, transported without changing its map. -/
def golayBlockLinearMap :
    (Fin 12 → ZMod 2) →ₗ[ZMod 2] (Fin 24 → ZMod 2) :=
  liftBinaryLinearMap Golay.enc golay_zero golay_enc_xor

/-- Running parity, with zero initial parity, on actual binary vectors. -/
def accumulatorLinearMap (n : ℕ) :
    (Fin n → ZMod 2) →ₗ[ZMod 2] (Fin n → ZMod 2) :=
  liftBinaryLinearMap accF (accF_zero n) accF_xor

/-- Semantic running-parity recurrence; its cost comes from the shared scan circuit. -/
def accumulatorRun : {n : ℕ} → ZMod 2 → (Fin n → ZMod 2) → (Fin n → ZMod 2)
  | 0, _, _ => Fin.elim0
  | n + 1, q, x => Fin.cons (q + x 0)
      (accumulatorRun (q + x 0) (Fin.tail x))

theorem accumulatorRun_bool {n : ℕ} (q : ZMod 2) (x : Fin n → ZMod 2) :
    List.ofFn ((wordEquiv n).symm (accumulatorRun q x)) =
      accLAux (bitEquiv.symm q) (List.ofFn ((wordEquiv n).symm x)) := by
  induction n generalizing q with
  | zero => simp [accumulatorRun, accLAux]
  | succ n ih =>
    rw [List.ofFn_succ, List.ofFn_succ]
    change bitEquiv.symm (q + x 0) ::
      List.ofFn ((wordEquiv n).symm
        (accumulatorRun (q + x 0) (Fin.tail x))) = _
    rw [ih]
    simp only [accLAux, bitEquiv_symm_add]
    rfl

theorem accumulatorLinearMap_eq_run (n : ℕ) (x : Fin n → ZMod 2) :
    accumulatorLinearMap n x = accumulatorRun 0 x := by
  apply (wordEquiv n).symm.injective
  apply List.ofFn_injective
  change List.ofFn ((wordEquiv n).symm
    (wordEquiv n (accF ((wordEquiv n).symm x)))) = _
  rw [Equiv.symm_apply_apply, ofFn_accF, accumulatorRun_bool]
  rfl

def blockGolayBool (k : ℕ) (x : Fin (k * 12) → Bool) :
    Fin (k * 24) → Bool :=
  outerWord Golay.enc (fun i j => x (finProdFinEquiv (i, j)))

def blockGolayLinearMap (k : ℕ) :
    (Fin (k * 12) → ZMod 2) →ₗ[ZMod 2] (Fin (k * 24) → ZMod 2) :=
  liftBinaryLinearMap (blockGolayBool k) (by
    funext j
    change Golay.enc (0 : Fin 12 → Bool) (finProdFinEquiv.symm j).2 = 0
    exact congrFun golay_zero _) (by
    intro x y
    exact outerWord_xor _ _)

def outerConstituentBool {k : ℕ} (seed : Seed k)
    (x : Fin (k * 12) → Bool) : Fin (k * 24) → Bool :=
  encode seed (fun i j => x (finProdFinEquiv (i, j)))

/-- The actual growing constituent, with the sampled pair of permutations fixed. -/
def outerConstituentLinearMap {k : ℕ} (seed : Seed k) :
    (Fin (k * 12) → ZMod 2) →ₗ[ZMod 2] (Fin (k * 24) → ZMod 2) :=
  liftBinaryLinearMap (outerConstituentBool seed) (encode_zero seed) (by
    intro x y
    exact encode_xor seed _ _)

def permuteVectorLinearMap {n : ℕ} (σ : Equiv.Perm (Fin n)) :
    (Fin n → ZMod 2) →ₗ[ZMod 2] (Fin n → ZMod 2) where
  toFun x := x ∘ σ
  map_add' _ _ := rfl
  map_smul' _ _ := rfl

theorem wordEquiv_comp {n : ℕ} (σ : Equiv.Perm (Fin n))
    (x : Fin n → Bool) :
    wordEquiv n (x ∘ σ) = wordEquiv n x ∘ σ := rfl

theorem wordEquiv_symm_comp {n : ℕ} (σ : Equiv.Perm (Fin n))
    (x : Fin n → ZMod 2) :
    (wordEquiv n).symm (x ∘ σ) = (wordEquiv n).symm x ∘ σ := rfl

/-- The semantic factorization used by the costed circuit construction. -/
theorem outerConstituentLinearMap_factor {k : ℕ} (seed : Seed k)
    (x : Fin (k * 12) → ZMod 2) :
    outerConstituentLinearMap seed x =
      accumulatorLinearMap (k * 24)
        (permuteVectorLinearMap seed.2
          (accumulatorLinearMap (k * 24)
            (permuteVectorLinearMap seed.1 (blockGolayLinearMap k x)))) := by
  change wordEquiv _ (accF ((accF
    (blockGolayBool k ((wordEquiv _).symm x) ∘ seed.1)) ∘ seed.2)) = _
  simp only [accumulatorLinearMap, permuteVectorLinearMap, blockGolayLinearMap,
    liftBinaryLinearMap, LinearMap.coe_mk, AddHom.coe_mk, wordEquiv_symm_comp,
    Equiv.symm_apply_apply, wordEquiv_comp]

/-- All native outer rows, expressed directly as binary matrix coordinates. -/
def nativeOuterVectorLinearMap (m : ℕ) (seed : NativeSeed m) :
    (Fin (Nsched m / 2) → ZMod 2) →ₗ[ZMod 2]
      ((Fin (Lsched m) × Fin (bsched m)) → ZMod 2) where
  toFun x p := outerConstituentLinearMap seed
    (fun t => x (finCongr (native_dimension m) (finProdFinEquiv (p.1, t))))
    ((finCongr (native_width m)).symm p.2)
  map_add' x y := by
    funext p
    exact congrFun ((outerConstituentLinearMap seed).map_add _ _)
      ((finCongr (native_width m)).symm p.2)
  map_smul' c x := by
    funext p
    exact congrFun ((outerConstituentLinearMap seed).map_smul c _)
      ((finCongr (native_width m)).symm p.2)

theorem nativeOuterVectorLinearMap_encode (m : ℕ) (seed : NativeSeed m)
    (x : Fin (Nsched m / 2) → ZMod 2) (i : Fin (Lsched m)) (j : Fin (bsched m)) :
    nativeOuterVectorLinearMap m seed x (i, j) =
      bitEquiv (encode seed
        (nativeMessageEquiv m ((wordEquiv (Nsched m / 2)).symm x) i)
        ((finCongr (native_width m)).symm j)) := rfl

/-- The costed outer's vector carrier is exactly the old support-based row carrier. -/
theorem nativeOuterVectorLinearMap_rows (m : ℕ) (seed : NativeSeed m)
    (x : Fin (Nsched m / 2) → ZMod 2) (i : Fin (Lsched m)) (j : Fin (bsched m)) :
    nativeOuterVectorLinearMap m seed x (i, j) =
      rowsToWord (nativeRows m seed ((wordEquiv (Nsched m / 2)).symm x))
        (finProdFinEquiv (i, j)) := by
  rw [nativeOuterVectorLinearMap_encode]
  simp [rowsToWord, nativeRows, rowSupports, support]

#print axioms nativeOuterVectorLinearMap_rows
#print axioms outerConstituentLinearMap_factor
end Spin.UniformMixer
