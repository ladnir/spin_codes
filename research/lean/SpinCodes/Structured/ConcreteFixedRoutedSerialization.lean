import SpinCodes.Structured.ConcreteFixedRegionComposition
import SpinCodes.Structured.ConcreteMarkedReferenceMoment
import SpinCodes.Structured.ConcreteSerialization

noncomputable section
namespace Spin.Structured.Placement
open Finset ConcreteRoute ConcreteEncoder

def regionWiring {b R : ℕ} : (Fin b × Fin (128*R)) ≃ (Fin (b*R) × Fin 128) :=
  regionMajor (by ring)

theorem regionWiring_block {b R : ℕ} (i : Fin b) (j : Fin R) (c : Fin 128) :
    regionWiring (i,blockBit j c) = (finProdFinEquiv (i,j),c) := by
  apply finProdFinEquiv.injective
  change finProdFinEquiv (regionMajor _ (i,blockBit j c)) = _
  rw [show finProdFinEquiv (regionMajor (L:=128*R) (b:=b) (R:=b*R) (by ring) (i,blockBit j c)) =
    finCongr (by ring) (finProdFinEquiv (i,blockBit j c)) from finProdFinEquiv.apply_symm_apply _]
  apply Fin.ext
  change 128*j.val+c.val+(128*R)*i.val = c.val+128*(j.val+R*i.val)
  ring

theorem reshape_regionWiring {b R : ℕ} (regions : Fin b → Finset (Fin (128*R)))
    (i : Fin b) (j : Fin R) :
    reshape regionWiring regions (finProdFinEquiv (i,j)) = serializedInput (regions i) j := by
  ext c
  rw [mem_reshape]
  have he : regionWiring.symm (finProdFinEquiv (i,j),c) = (i,blockBit j c) := by
    rw [←regionWiring_block,Equiv.symm_apply_apply]
  rw [he]
  simp [serializedInput]

/-- The actual region-major wiring is exactly the concatenated region stream. -/
theorem ofFn_reshape_regionWiring {b R : ℕ} (regions : Fin b → Finset (Fin (128*R))) :
    List.ofFn (reshape regionWiring regions) = regionStream regions := by
  rw [List.ofFn_mul]
  unfold regionStream
  congr 1
  apply congrArg List.ofFn
  funext i
  apply congrArg List.ofFn
  funext j
  convert reshape_regionWiring regions i j using 1
  congr 1
  apply Fin.ext
  change i.val*R+j.val=j.val+R*i.val
  ring

/-- Equality of the weighted actual IMT experiments, avoiding any length-cast assumption. -/
theorem inputMoment_regionWiring {b R : ℕ} (z : ℝ)
    (regions : Fin b → Finset (Fin (128*R))) (q : State) :
    inputMoment z (reshape regionWiring regions) q = inputMoment z (regionStream regions).get q := by
  simp only [←endpointKernel_total,endpointKernel_eq_path,List.ofFn_get,
    ofFn_reshape_regionWiring]

theorem inputMoment_regionMajor {b R T : ℕ} (h : b*(128*R)=T*128) (z : ℝ)
    (regions : Fin b → Finset (Fin (128*R))) (q : State) :
    inputMoment z (reshape (regionMajor h) regions) q = inputMoment z (regionStream regions).get q := by
  have he : T=b*R := by nlinarith
  subst T
  exact inputMoment_regionWiring z regions q

#print axioms regionWiring_block
#print axioms ofFn_reshape_regionWiring
#print axioms inputMoment_regionWiring
#print axioms inputMoment_regionMajor
end Spin.Structured.Placement

