import SpinCodes.Structured.ConcreteBinaryEncoding
import SpinCodes.Structured.ConcreteRoutedFirstMoment
import SpinCodes.Structured.Instantiation

noncomputable section
attribute [local instance] Classical.propDecidable
namespace Spin.Structured.ConcreteNativeFamily
open ConcreteOuter ConcreteRoute ConcreteBinaryEncoding

def rounds (m : ℕ) : ℕ := bsched m * Lsched m / 128

abbrev InnerSeed (m : ℕ) :=
  ConcreteRoutedFirstMoment.InnerSeeds (Lsched m) (bsched m) (rounds m)

def rows (m : ℕ) (seed : NativeSeed m) (message : Fin (Nsched m / 2) → ZMod 2) :
    Fin (Lsched m) → Finset (Fin (bsched m)) :=
  nativeRows m seed ((wordEquiv (Nsched m / 2)).symm message)

def nativeSetup (m : ℕ) : Setup (Fin (Nsched m / 2) → ZMod 2)
    (Fin (Nsched m) → ZMod 2) (NativeSeed m) (InnerSeed m) where
  Pout := nativeSeedLaw m
  Pin := ConcreteRoutedEncoder.experimentLaw (Lsched m) (bsched m) (rounds m)
  outer seed message := rowsToWord (rows m seed message)
  innerWt seed word := ConcreteRoutedEncoder.weight (streamWiring (native_round_divisibility m))
    (wordToRows word) seed

def occ (m : ℕ) (seed : NativeSeed m) (message : Fin (Nsched m / 2) → ZMod 2) : ℕ :=
  activeRows (rows m seed message)

theorem occ_mem (m : ℕ) (seed : NativeSeed m) (message : Fin (Nsched m / 2) → ZMod 2)
    (hne : message ≠ 0) : occ m seed message ∈ Finset.Ico 1 (Lsched m + 1) := by
  apply native_occupation_mem
  intro hz
  have hh := congrArg (wordEquiv (Nsched m / 2)) hz
  simp only [Equiv.apply_symm_apply, wordEquiv_zero] at hh
  exact hne hh

theorem setup_inner_weight (m : ℕ) (out : NativeSeed m) (message : Fin (Nsched m / 2) → ZMod 2)
    (seed : InnerSeed m) :
    (nativeSetup m).innerWt seed ((nativeSetup m).outer out message) =
      ConcreteRoutedEncoder.weight (streamWiring (native_round_divisibility m))
        (rows m out message) seed := by
  change ConcreteRoutedEncoder.weight _ (wordToRows (rowsToWord _)) _ = _
  rw [wordToRows_rowsToWord]

/-- The exact rational threshold, before any asymptotic approximation. -/
def threshold (m : ℕ) : ℕ := 11 * Nsched m / 100

theorem threshold_floor (m : ℕ) : threshold m = ⌊(0.11 : ℝ) * Nsched m⌋₊ := by
  have he : (0.11 : ℝ) * Nsched m = ((11 * Nsched m : ℕ) : ℝ) / (100 : ℕ) := by
    push_cast
    ring
  rw [he, Nat.floor_div_eq_div]
  rfl

/-- The actual native construction; the positive outer selection event remains an explicit parameter. -/
def family (G : ∀ m, NativeSeed m → Prop)
    (hG : ∀ m, 0 < (nativeSeedLaw m).prob (G m)) : Family where
  K := fun m => Nsched m / 2
  N := Nsched
  L := Lsched
  Out := NativeSeed
  Inr := InnerSeed
  outFin := fun _ => inferInstance
  inrFin := fun _ => inferInstance
  S := nativeSetup
  G := G
  Gpos := hG
  occ := occ
  occ_mem := fun m seed message hm => occ_mem m seed message
    (Finset.mem_filter.mp hm).2
  d := threshold

@[simp] theorem family_N (G : ∀ m, NativeSeed m → Prop)
    (hG : ∀ m, 0 < (nativeSeedLaw m).prob (G m)) (m : ℕ) :
    (family G hG).N m = Nsched m := rfl

@[simp] theorem family_L (G : ∀ m, NativeSeed m → Prop)
    (hG : ∀ m, 0 < (nativeSeedLaw m).prob (G m)) (m : ℕ) :
    (family G hG).L m = Lsched m := rfl

@[simp] theorem family_threshold (G : ∀ m, NativeSeed m → Prop)
    (hG : ∀ m, 0 < (nativeSeedLaw m).prob (G m)) (m : ℕ) :
    (family G hG).d m = 11 * Nsched m / 100 := rfl

end Spin.Structured.ConcreteNativeFamily
