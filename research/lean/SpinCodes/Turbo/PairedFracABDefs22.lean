import SpinCodes.Turbo.SparseRpowABDefs22

/-! Core-only pair certificates for the paired fractional chain at α = a/b
(generalizing `pairedRpowCheckL22`, α = 4/5). Row `a` is a separate
function so kernel checks can be split per `a` across processes. -/
namespace Spin.Turbo

def pairedRpowABRowL22 (ea eb : Nat) (den sy : Int) (R : List (List (List Int)))
    (Ys : List (List (List (List Int)))) (a : Nat) : Bool :=
  sparseMatShapeL22 (R.getD a []) && (List.range 17).all fun b =>
    sparseRpowABCheckL22 ea eb (den * den) sy
      (sparseMatMulL22 (R.getD a []) (R.getD b [])) ((Ys.getD a []).getD b [])

def pairedRpowABCheckL22 (ea eb : Nat) (den sy : Int) (R : List (List (List Int)))
    (Ys : List (List (List (List Int)))) : Bool :=
  decide (0 < den) && (List.range 17).all fun a => pairedRpowABRowL22 ea eb den sy R Ys a

theorem pairedRpowABCheck_of (ea eb : Nat) (den sy : Int) (R : List (List (List Int)))
    (Ys : List (List (List (List Int)))) (h0 : 0 < den)
    (h : ∀ a, a < 17 → pairedRpowABRowL22 ea eb den sy R Ys a = true) :
    pairedRpowABCheckL22 ea eb den sy R Ys = true := by
  simp only [pairedRpowABCheckL22, decide_eq_true h0, Bool.true_and, List.all_eq_true,
    List.mem_range]
  exact h

end Spin.Turbo
