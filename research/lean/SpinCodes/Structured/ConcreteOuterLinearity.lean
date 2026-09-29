import SpinCodes.Structured.ConcreteOuterNative

namespace Spin.Structured.ConcreteOuter

lemma accLAux_xor (p q : Bool) (xs ys : List Bool) :
    accLAux (xor p q) (List.zipWith xor xs ys) =
      List.zipWith xor (accLAux p xs) (accLAux q ys) := by
  induction xs generalizing ys p q with
  | nil => simp [accLAux]
  | cons x xs ih =>
    cases ys with
    | nil => simp [accLAux]
    | cons y ys =>
      have he : xor (xor p q) (xor x y)=xor (xor p x) (xor q y) := by
        cases p <;> cases q <;> cases x <;> cases y <;> rfl
      simp only [List.zipWith_cons_cons,accLAux,he]
      rw [ih]

lemma ofFn_xor {n : ℕ} (x y : Fin n→Bool) :
    List.ofFn (fun i => xor (x i) (y i))=List.zipWith xor (List.ofFn x) (List.ofFn y) := by
  apply List.ext_getElem
  · simp
  · intro i hi hj
    simp

theorem accF_xor {n : ℕ} (x y : Fin n→Bool) :
    accF (fun i => xor (x i) (y i))=fun i => xor (accF x i) (accF y i) := by
  apply List.ofFn_injective
  rw [ofFn_accF,ofFn_xor,ofFn_xor,ofFn_accF,ofFn_accF]
  exact accLAux_xor false false _ _

lemma golay_fold_step_xor (p q : Bool) (a b r : List Bool)
    (har : a.length=r.length) (hbr : b.length=r.length) :
    (if xor p q then Golay.xorW (Golay.xorW a b) r else Golay.xorW a b)=
      Golay.xorW (if p then Golay.xorW a r else a) (if q then Golay.xorW b r else b) := by
  cases p <;> cases q <;> simp only [Bool.false_xor,Bool.true_xor,Bool.not_false,Bool.not_true,Bool.false_eq_true,ite_false,ite_true]
  all_goals apply List.ext_getElem
  all_goals try {simp [Golay.xorW,har,hbr]}
  all_goals intro i hi hj
  all_goals have hir : i < r.length := by simpa [Golay.xorW,har,hbr] using hi
  all_goals have hia : i < a.length := by omega
  all_goals have hib : i < b.length := by omega
  all_goals simp only [Golay.xorW,List.getElem_zipWith]
  all_goals first | rfl | (cases a[i] <;> cases b[i] <;> cases r[i] <;> rfl)

end Spin.Structured.ConcreteOuter

