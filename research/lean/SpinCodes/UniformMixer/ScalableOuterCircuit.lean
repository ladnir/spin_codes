import SpinCodes.LinearCircuitScan
import SpinCodes.UniformMixer.ScalableOuterAlgebra

/-! Costed compositional circuits for the actual shared Golay--BAA outer. -/
noncomputable section
namespace Spin.UniformMixer
open Spin.Structured Spin.Structured.ConcreteOuter
  Spin.Structured.ConcreteBinaryEncoding Spin.LinearCircuit

/-- Fixed-size dense macro for the checked Golay generator. -/
def golayBlockCircuit : Circuit (Fin 12) (Fin 24) :=
  .dense (LinearMap.toMatrix' golayBlockLinearMap)

theorem golayBlockCircuit_eval (x : Fin 12 → ZMod 2) :
    golayBlockCircuit.eval x = golayBlockLinearMap x :=
  LinearMap.toMatrix'_mulVec _ _

@[simp] theorem golayBlockCircuit_size : golayBlockCircuit.size = 576 := rfl

/-- Tiling retains constant work per Golay block as the constituent grows. -/
def blockGolayCircuit (k : ℕ) : Circuit (Fin (k * 12)) (Fin (k * 24)) :=
  .comp (.reindex finProdFinEquiv)
    (.comp (.tiled (σ := Fin k) golayBlockCircuit)
      (.reindex finProdFinEquiv.symm))

theorem blockGolayCircuit_eval (k : ℕ) (x : Fin (k * 12) → ZMod 2) :
    (blockGolayCircuit k).eval x = blockGolayLinearMap k x := by
  funext j
  simp only [blockGolayCircuit, Circuit.eval, Equiv.symm_symm]
  rw [golayBlockCircuit_eval]
  rfl

@[simp] theorem blockGolayCircuit_size (k : ℕ) :
    (blockGolayCircuit k).size = 576 * k := by
  simp [blockGolayCircuit, Circuit.size, golayBlockCircuit_size, Nat.mul_comm]

def outerPermutationCircuit {n : ℕ} (σ : Equiv.Perm (Fin n)) :
    Circuit (Fin n) (Fin n) := .reindex σ.symm

theorem outerPermutationCircuit_eval {n : ℕ} (σ : Equiv.Perm (Fin n))
    (x : Fin n → ZMod 2) :
    (outerPermutationCircuit σ).eval x = permuteVectorLinearMap σ x := rfl

@[simp] theorem outerPermutationCircuit_size {n : ℕ} (σ : Equiv.Perm (Fin n)) :
    (outerPermutationCircuit σ).size = 0 := rfl

def accumulatorStepLinearMap :
    ((Unit ⊕ Unit) → ZMod 2) →ₗ[ZMod 2] ((Unit ⊕ Unit) → ZMod 2) where
  toFun x := fun _ => x (.inl ()) + x (.inr ())
  map_add' x y := by funext z; simp only [Pi.add_apply]; abel
  map_smul' c x := by funext z; simp only [Pi.smul_apply, smul_add, RingHom.id_apply]

def accumulatorStepCircuit : Circuit (Unit ⊕ Unit) (Unit ⊕ Unit) :=
  .dense (LinearMap.toMatrix' accumulatorStepLinearMap)

theorem accumulatorStepCircuit_eval (x : (Unit ⊕ Unit) → ZMod 2) :
    accumulatorStepCircuit.eval x = fun _ => x (.inl ()) + x (.inr ()) :=
  LinearMap.toMatrix'_mulVec _ _

@[simp] theorem accumulatorStepCircuit_size : accumulatorStepCircuit.size = 8 := by
  norm_num [accumulatorStepCircuit, Circuit.size]

def accumulatorCoordinates (n : ℕ) : Fin n ≃ Fin n × Unit where
  toFun i := (i, ())
  invFun p := p.1
  left_inv _ := rfl
  right_inv p := by rcases p with ⟨i, u⟩; cases u; rfl

theorem accumulator_runScan {n : ℕ} (q : ZMod 2) (x : Fin n → ZMod 2) :
    (runScan (fun _ : Fin n => accumulatorStepCircuit)
      (fun i (_ : Unit) => x i) (fun _ : Unit => q)).1 =
      fun i (_ : Unit) => accumulatorRun q x i := by
  induction n generalizing q with
  | zero => funext i; exact Fin.elim0 i
  | succ n ih =>
    funext i u
    refine Fin.cases ?_ (fun i => ?_) i
    · simp only [runScan, accumulatorStepCircuit_eval, Sum.elim_inl, Sum.elim_inr,
        Fin.cons_zero, accumulatorRun]
    · have h := congrFun (congrFun (ih (q + x 0) (Fin.tail x)) i) u
      simpa only [runScan, accumulatorStepCircuit_eval, Sum.elim_inl, Sum.elim_inr,
        accumulatorRun, Fin.cons_succ, Fin.tail_def] using h

/-- Prefix parity implemented by the shared constant-state scan construction. -/
def accumulatorCircuit (n : ℕ) : Circuit (Fin n) (Fin n) :=
  .comp (.reindex (accumulatorCoordinates n).symm)
    (.comp (zeroScan (fun _ : Fin n => accumulatorStepCircuit))
      (.reindex (accumulatorCoordinates n)))

theorem accumulatorCircuit_eval (n : ℕ) (x : Fin n → ZMod 2) :
    (accumulatorCircuit n).eval x = accumulatorLinearMap n x := by
  rw [accumulatorLinearMap_eq_run]
  funext i
  have h := congrFun
    (zeroScan_eval (fun _ : Fin n => accumulatorStepCircuit)
      (fun i (_ : Unit) => x i)) (i, ())
  have hh := h.trans (congrFun (congrFun (accumulator_runScan (n := n) 0 x) i) ())
  change (zeroScan (fun _ : Fin n => accumulatorStepCircuit)).eval
    (fun ix => x ix.1) (i, ()) = accumulatorRun 0 x i
  exact hh

@[simp] theorem accumulatorCircuit_size (n : ℕ) :
    (accumulatorCircuit n).size = 8 * n := by
  simp [accumulatorCircuit, Circuit.size, zeroScan_size, Nat.mul_comm]

/-- No dense macro grows with the constituent width. -/
def outerConstituentCircuit {k : ℕ} (seed : Seed k) :
    Circuit (Fin (k * 12)) (Fin (k * 24)) :=
  .comp (accumulatorCircuit (k * 24))
    (.comp (outerPermutationCircuit seed.2)
      (.comp (accumulatorCircuit (k * 24))
        (.comp (outerPermutationCircuit seed.1) (blockGolayCircuit k))))

theorem outerConstituentCircuit_eval {k : ℕ} (seed : Seed k)
    (x : Fin (k * 12) → ZMod 2) :
    (outerConstituentCircuit seed).eval x = outerConstituentLinearMap seed x := by
  simp only [outerConstituentCircuit, Circuit.eval, accumulatorCircuit_eval,
    outerPermutationCircuit_eval, blockGolayCircuit_eval]
  exact (outerConstituentLinearMap_factor seed x).symm

@[simp] theorem outerConstituentCircuit_size {k : ℕ} (seed : Seed k) :
    (outerConstituentCircuit seed).size = 960 * k := by
  simp only [outerConstituentCircuit, Circuit.size, accumulatorCircuit_size,
    outerPermutationCircuit_size, blockGolayCircuit_size]
  omega

def nativeOuterInputEquiv (m : ℕ) :
    Fin (Nsched m / 2) ≃ Fin (Lsched m) × Fin (nativeBlocks m * 12) :=
  (finCongr (native_dimension m).symm).trans finProdFinEquiv.symm

def nativeOuterOutputEquiv (m : ℕ) :
    (Fin (Lsched m) × Fin (nativeBlocks m * 24)) ≃
      (Fin (Lsched m) × Fin (bsched m)) :=
  Equiv.prodCongr (Equiv.refl _) (finCongr (native_width m))

/-- One fixed constituent circuit is applied to every outer message row. -/
def nativeOuterCircuit (m : ℕ) (seed : NativeSeed m) :
    Circuit (Fin (Nsched m / 2)) (Fin (Lsched m) × Fin (bsched m)) :=
  .comp (.reindex (nativeOuterOutputEquiv m))
    (.comp (.tiled (σ := Fin (Lsched m)) (outerConstituentCircuit seed))
      (.reindex (nativeOuterInputEquiv m)))

theorem nativeOuterCircuit_eval (m : ℕ) (seed : NativeSeed m)
    (x : Fin (Nsched m / 2) → ZMod 2) :
    (nativeOuterCircuit m seed).eval x = nativeOuterVectorLinearMap m seed x := by
  funext p
  simp only [nativeOuterCircuit, Circuit.eval, outerConstituentCircuit_eval]
  rfl

@[simp] theorem nativeOuterCircuit_size (m : ℕ) (seed : NativeSeed m) :
    (nativeOuterCircuit m seed).size = 40 * Nsched m := by
  simp only [nativeOuterCircuit, Circuit.size, outerConstituentCircuit_size,
    Fintype.card_fin, Nat.add_zero, Nat.zero_add]
  calc Lsched m * (960 * nativeBlocks m) =
      40 * (Lsched m * (nativeBlocks m * 24)) := by ring
    _ = 40 * (Lsched m * bsched m) := by rw [native_width]
    _ = 40 * Nsched m := by rw [Nsched_eq]

#print axioms nativeOuterCircuit_eval

#print axioms accumulatorCircuit_eval
#print axioms outerConstituentCircuit_eval

#print axioms blockGolayCircuit_eval
end Spin.UniformMixer
