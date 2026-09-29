import SpinCodes.Structured.ConcreteOuterLinearity

namespace Spin.Structured.ConcreteOuter

lemma golay_fold_xor (n : ℕ) (rs : List (List Bool)) (mx my a b : List Bool)
    (hx : mx.length=rs.length) (hy : my.length=rs.length)
    (ha : a.length=n) (hb : b.length=n) (hr : ∀ r∈rs,r.length=n) :
    (((List.zipWith xor mx my).zip rs).foldl
      (fun acc p => if p.1 then Golay.xorW acc p.2 else acc) (Golay.xorW a b)) =
    Golay.xorW ((mx.zip rs).foldl (fun acc p => if p.1 then Golay.xorW acc p.2 else acc) a)
      ((my.zip rs).foldl (fun acc p => if p.1 then Golay.xorW acc p.2 else acc) b) := by
  induction rs generalizing mx my a b with
  | nil => simp
  | cons r rs ih =>
    cases mx with
    | nil => simp at hx
    | cons x mx =>
      cases my with
      | nil => simp at hy
      | cons y my =>
        have hrl : r.length=n := hr r (by simp)
        simp only [List.length_cons,Nat.add_right_cancel_iff] at hx hy
        simp only [List.zipWith_cons_cons,List.zip_cons_cons,List.foldl_cons]
        rw [golay_fold_step_xor x y a b r (ha.trans hrl.symm) (hb.trans hrl.symm)]
        apply ih _ _ _ _ hx hy
        · cases x <;> simp [Golay.length_xorW,ha,hrl]
        · cases y <;> simp [Golay.length_xorW,hb,hrl]
        · intro s hs
          exact hr s (List.mem_cons_of_mem r hs)

lemma golay_encList_xor (x y : List Bool) (hx : x.length=12) (hy : y.length=12) :
    Golay.encList (List.zipWith xor x y)=Golay.xorW (Golay.encList x) (Golay.encList y) := by
  have hh := golay_fold_xor 24 Golay.rows x y (List.replicate 24 false) (List.replicate 24 false)
    (by simpa only [show Golay.rows.length=12 by decide] using hx)
    (by simpa only [show Golay.rows.length=12 by decide] using hy)
    (by simp) (by simp) Golay.length_rows_mem
  have hz : Golay.xorW (List.replicate 24 false) (List.replicate 24 false)=List.replicate 24 false := by decide
  rw [hz] at hh
  exact hh

theorem golay_enc_xor (x y : Fin 12→Bool) :
    Golay.enc (fun i => xor (x i) (y i))=fun i => xor (Golay.enc x i) (Golay.enc y i) := by
  apply List.ofFn_injective
  rw [Golay.ofFn_enc,ofFn_xor,ofFn_xor,Golay.ofFn_enc,Golay.ofFn_enc]
  exact golay_encList_xor _ _ (by simp) (by simp)

theorem outerWord_xor {k : ℕ} (x y : LocalMessage k) :
    outerWord Golay.enc (fun i j => xor (x i j) (y i j))=
      fun j => xor (outerWord Golay.enc x j) (outerWord Golay.enc y j) := by
  funext j
  exact congrFun (golay_enc_xor _ _) _

/-- The actual two-permutation, two-accumulator constituent is binary linear. -/
theorem encode_xor {k : ℕ} (seed : Seed k) (x y : LocalMessage k) :
    encode seed (fun i j => xor (x i j) (y i j))=
      fun j => xor (encode seed x j) (encode seed y j) := by
  unfold encode
  rw [outerWord_xor]
  change accF ((accF (fun i => xor ((outerWord Golay.enc x ∘ seed.1) i)
    ((outerWord Golay.enc y ∘ seed.1) i))) ∘ seed.2)=_
  rw [accF_xor]
  exact accF_xor _ _

end Spin.Structured.ConcreteOuter

