import SpinCodes.Structured.ConcreteFixedOrderStatisticVolume
import Mathlib.Order.Fin.Tuple

noncomputable section
namespace Spin.Structured.Placement
open MeasureTheory Set

def appendMeasurableEquiv (a b : ℕ) :
    ((Fin a→ℝ) × (Fin b→ℝ)) ≃ᵐ (Fin (a+b)→ℝ) :=
  (MeasurableEquiv.sumPiEquivProdPi (fun _ : Fin a ⊕ Fin b => ℝ)).symm.trans
    (MeasurableEquiv.piCongrLeft (fun _ : Fin (a+b) => ℝ) finSumFinEquiv)

theorem appendMeasurableEquiv_apply (a b : ℕ) (u : Fin a→ℝ) (v : Fin b→ℝ) :
    appendMeasurableEquiv a b (u,v)=Fin.append u v := by
  ext i
  refine Fin.addCases (fun j => ?_) (fun j => ?_) i <;>
    simp [appendMeasurableEquiv, MeasurableEquiv.trans_apply,
      MeasurableEquiv.piCongrLeft, Equiv.piCongrLeft, Equiv.piCongrLeft',
      MeasurableEquiv.sumPiEquivProdPi]

theorem appendMeasurableEquiv_preserving (a b : ℕ) :
    MeasurePreserving (appendMeasurableEquiv a b) volume volume :=
  (volume_measurePreserving_piCongrLeft (fun _ : Fin (a+b) => ℝ) finSumFinEquiv).comp
    (volume_measurePreserving_sumPiEquivProdPi_symm (fun _ : Fin a ⊕ Fin b => ℝ))

theorem strictMono_append_iff {a b : ℕ} (u : Fin a→ℝ) (v : Fin b→ℝ) :
    StrictMono (Fin.append u v) ↔ StrictMono u ∧ StrictMono v ∧ ∀ i j, u i<v j := by
  constructor
  · intro h
    refine ⟨?_,?_,?_⟩
    · intro i j hij
      simpa using h (show i.castAdd b<j.castAdd b from hij)
    · intro i j hij
      simpa using h (show i.natAdd a<j.natAdd a by simpa using hij)
    · intro i j
      simpa using h (show i.castAdd b<j.natAdd a by simp only [Fin.lt_def,Fin.val_castAdd,Fin.val_natAdd]; omega)
  · rintro ⟨hu,hv,huv⟩ i j hij
    cases i using Fin.addCases <;> cases j using Fin.addCases
    · simpa using hu (by simpa only [Fin.lt_def,Fin.val_castAdd] using hij)
    · simpa using huv _ _
    · simp only [Fin.lt_def,Fin.val_castAdd,Fin.val_natAdd] at hij; omega
    · simpa using hv (by simpa using hij)

def orderStatisticJoin {a b : ℕ} (y : ℝ) (u : Fin a→ℝ) (v : Fin b→ℝ) :
    Fin (a+b+1)→ℝ := (⟨a,by omega⟩ : Fin (a+b+1)).insertNth y (Fin.append u v)

theorem orderStatisticJoin_mem {a b : ℕ} (y : ℝ) (u : Fin a→ℝ) (v : Fin b→ℝ) :
    orderStatisticJoin y u v ∈ orderedSiteDomain (a+b+1) ↔
      y∈Ico (0:ℝ) 1 ∧ u∈orderedBetween a 0 y ∧
      v∈orderedBetween b y 1 ∧ ∀j,y≠v j := by
  let k : Fin (a+b+1) := ⟨a,by omega⟩
  change ((∀i,0≤(k.insertNth y (Fin.append u v) : Fin (a+b+1)→ℝ) i ∧ (k.insertNth y (Fin.append u v) : Fin (a+b+1)→ℝ) i<1) ∧
    StrictMono (k.insertNth y (Fin.append u v))) ↔ _
  rw [Fin.strictMono_insertNth_iff,strictMono_append_iff]
  constructor
  · rintro ⟨hx,⟨hu,hv,huv⟩,hl,hr⟩
    have hy : 0≤y ∧ y<1 := by simpa using hx k
    have hleft (i : Fin a) : u i<y := by
      simpa using hl (i.castAdd b) (by simp [k,Fin.lt_def])
    have hright (j : Fin b) : y<v j := by
      simpa using hr (j.natAdd a) (by simp [k,Fin.le_def])
    refine ⟨hy,⟨?_,hu⟩,⟨?_,hv⟩,fun j => (hright j).ne⟩
    · intro i
      have hh := hx (k.succAbove (i.castAdd b))
      simp only [Fin.insertNth_apply_succAbove,Fin.append_left] at hh
      exact ⟨hh.1,hleft i⟩
    · intro j
      have hh := hx (k.succAbove (j.natAdd a))
      simp only [Fin.insertNth_apply_succAbove,Fin.append_right] at hh
      exact ⟨(hright j).le,hh.2⟩
  · rintro ⟨hy,⟨hu,hum⟩,⟨hv,hvm⟩,hne⟩
    have hr (j : Fin b) : y<v j := lt_of_le_of_ne (hv j).1 (hne j)
    refine ⟨?_,⟨hum,hvm,fun i j => (hu i).2.trans (hr j)⟩,?_,?_⟩
    · intro i
      cases i using Fin.succAboveCases k
      · simpa using hy
      · rename_i i
        simp only [Fin.insertNth_apply_succAbove]
        cases i using Fin.addCases
        · simpa using And.intro (hu _).1 ((hu _).2.trans hy.2)
        · simpa using And.intro (hy.1.trans (hr _).le) (hv _).2
    · intro i hi
      cases i using Fin.addCases
      · simpa using (hu _).2
      · simp only [Fin.lt_def,Fin.val_castSucc,Fin.val_natAdd,k] at hi; omega
    · intro i hi
      cases i using Fin.addCases
      · simp only [Fin.le_def,Fin.val_castSucc,Fin.val_castAdd,k] at hi; omega
      · simpa using hr _

end Spin.Structured.Placement
