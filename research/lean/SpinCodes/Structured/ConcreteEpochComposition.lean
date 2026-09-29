import SpinCodes.Structured.ConcreteEpochKernel
import SpinCodes.Structured.ConcreteEpochVariance

/-! Exact concatenation of the actual encoder, including its weighted endpoint law. -/
noncomputable section
namespace Spin.Structured.ConcreteEncoder
open ConcreteMaps Spin.Imt
open scoped symmDiff

def runState : List (Input × Transvection) → State → State
  | [], q => q
  | p :: ps, q => runState ps (nextState q p.1 p.2)

def runWeight : List (Input × Transvection) → State → Nat
  | [], _ => 0
  | p :: ps, q => (p.1 ∆ Aset q).card + runWeight ps (nextState q p.1 p.2)

theorem runState_append (as bs : List (Input × Transvection)) (q : State) :
    runState (as ++ bs) q = runState bs (runState as q) := by
  induction as generalizing q with
  | nil => rfl
  | cons p ps ih => simpa only [List.cons_append, runState] using ih (nextState q p.1 p.2)

theorem runWeight_append (as bs : List (Input × Transvection)) (q : State) :
    runWeight (as ++ bs) q = runWeight as q + runWeight bs (runState as q) := by
  induction as generalizing q with
  | nil => simp only [List.nil_append, runWeight, runState, zero_add]
  | cons p ps ih => simp only [List.cons_append, runWeight, runState, ih, Nat.add_assoc]

theorem runState_ofFn {R : Nat} (xs : Fin R → Input) (ts : Fin R → Transvection) (q : State) :
    runState (List.ofFn (fun i => (xs i, ts i))) q = terminalState xs ts q := by
  induction R generalizing q with
  | zero => rfl
  | succ R ih =>
    rw [List.ofFn_succ]
    simpa only [runState, terminalState, Fin.tail] using
      ih (Fin.tail xs) (Fin.tail ts) (nextState q (xs 0) (ts 0))

theorem runWeight_ofFn {R : Nat} (xs : Fin R → Input) (ts : Fin R → Transvection) (q : State) :
    runWeight (List.ofFn (fun i => (xs i, ts i))) q = outputWeight xs ts q := by
  induction R generalizing q with
  | zero => rfl
  | succ R ih =>
    rw [List.ofFn_succ]
    simp only [runWeight, outputWeight]
    simpa only [Fin.tail] using congrArg (fun n => (xs 0 ∆ Aset q).card + n)
      (ih (Fin.tail xs) (Fin.tail ts) (nextState q (xs 0) (ts 0)))

theorem paired_append {m n : Nat} (xs : Fin m → Input) (ys : Fin n → Input)
    (ts : Fin m → Transvection) (us : Fin n → Transvection) :
    (fun i => (Fin.append xs ys i, Fin.append ts us i)) =
      Fin.append (fun i => (xs i, ts i)) (fun i => (ys i, us i)) := by
  funext i
  refine Fin.addCases (fun j => ?_) (fun j => ?_) i <;> simp only [Fin.append_left, Fin.append_right]

theorem terminalState_append {m n : Nat} (xs : Fin m → Input) (ys : Fin n → Input)
    (ts : Fin m → Transvection) (us : Fin n → Transvection) (q : State) :
    terminalState (Fin.append xs ys) (Fin.append ts us) q =
      terminalState ys us (terminalState xs ts q) := by
  rw [← runState_ofFn, paired_append, List.ofFn_fin_append, runState_append,
    runState_ofFn, runState_ofFn]

theorem outputWeight_append {m n : Nat} (xs : Fin m → Input) (ys : Fin n → Input)
    (ts : Fin m → Transvection) (us : Fin n → Transvection) (q : State) :
    outputWeight (Fin.append xs ys) (Fin.append ts us) q =
      outputWeight xs ts q + outputWeight ys us (terminalState xs ts q) := by
  rw [← runWeight_ofFn, paired_append, List.ofFn_fin_append, runWeight_append,
    runWeight_ofFn, runState_ofFn, runWeight_ofFn]

def endpointKernel {R : Nat} (z : ℝ) (xs : Fin R → Input) : Matrix State State ℝ := fun q r =>
  (Spin.piPMF (fun _ : Fin R => transvectionLaw)).expect (fun ts =>
    z ^ outputWeight xs ts q * (if terminalState xs ts q = r then 1 else 0))

def roundKernel (z : ℝ) (x : Input) : Matrix State State ℝ := fun q r =>
  emitted z q x * actLaw q (r ∆ Cset x)

def pathKernel (z : ℝ) : List Input → Matrix State State ℝ
  | [] => 1
  | x :: xs => roundKernel z x * pathKernel z xs

theorem pathKernel_append (z : ℝ) (xs ys : List Input) :
    pathKernel z (xs ++ ys) = pathKernel z xs * pathKernel z ys := by
  induction xs with
  | nil => simp only [List.nil_append, pathKernel, one_mul]
  | cons x xs ih => simp only [List.cons_append, pathKernel, ih, Matrix.mul_assoc]

theorem endpointKernel_zero (z : ℝ) (xs : Fin 0 → Input) : endpointKernel z xs = 1 := by
  ext q r
  change (Spin.piPMF (fun _ : Fin 0 => transvectionLaw)).expect
    (fun _ => z ^ 0 * (if q = r then 1 else 0)) = _
  rw [Spin.FinPMF.expect_const]
  simp only [pow_zero, one_mul, Matrix.one_apply]

theorem endpointKernel_succ {R : Nat} (z : ℝ) (xs : Fin (R + 1) → Input) (q r : State) :
    endpointKernel z xs q r = emitted z q (xs 0) *
      transvectionLaw.expect (fun t => endpointKernel z (Fin.tail xs) (nextState q (xs 0) t) r) := by
  unfold endpointKernel
  rw [expect_iid_succ]
  simp only [outputWeight, terminalState, Fin.cons_zero, Fin.tail_cons, pow_add, emitted]
  have he (t : Transvection) (ts : Fin R → Transvection) :
      z ^ (xs 0 ∆ Aset q).card * z ^ outputWeight (Fin.tail xs) ts (nextState q (xs 0) t) *
        (if terminalState (Fin.tail xs) ts (nextState q (xs 0) t) = r then 1 else 0) =
      z ^ (xs 0 ∆ Aset q).card *
        (z ^ outputWeight (Fin.tail xs) ts (nextState q (xs 0) t) *
          (if terminalState (Fin.tail xs) ts (nextState q (xs 0) t) = r then 1 else 0)) := by ring
  simp only [he, expect_const_mul]

theorem endpointKernel_eq_path {R : Nat} (z : ℝ) (xs : Fin R → Input) :
    endpointKernel z xs = pathKernel z (List.ofFn xs) := by
  induction R with
  | zero => rw [List.ofFn_zero, pathKernel, endpointKernel_zero]
  | succ R ih =>
    rw [List.ofFn_succ, pathKernel]
    ext q r
    rw [endpointKernel_succ]
    have he : (fun t : Transvection => endpointKernel z (Fin.tail xs) (nextState q (xs 0) t) r) =
        fun t => (fun w => pathKernel z (List.ofFn (Fin.tail xs)) w r) (nextState q (xs 0) t) := by
      funext t
      rw [ih]
    rw [he, nextState_expect q (xs 0) (fun w => pathKernel z (List.ofFn (Fin.tail xs)) w r)]
    simp only [Matrix.mul_apply, roundKernel, Finset.mul_sum, mul_assoc]
    rw [show List.ofFn (Fin.tail xs) = List.ofFn (fun i => xs i.succ) from rfl]

/-- Chapman–Kolmogorov for the actual independently sampled encoder rounds. -/
theorem endpointKernel_append {m n : Nat} (z : ℝ) (xs : Fin m → Input) (ys : Fin n → Input) :
    endpointKernel z (Fin.append xs ys) = endpointKernel z xs * endpointKernel z ys := by
  simp only [endpointKernel_eq_path, List.ofFn_fin_append, pathKernel_append]

theorem endpointKernel_empty (z : ℝ) (g : Nat) :
    endpointKernel z (fun _ : Fin g => ∅) = emptyKernel z g := rfl

theorem pathKernel_empty (z : ℝ) (g : Nat) :
    pathKernel z (List.replicate g ∅) = emptyKernel z g := by
  rw [← List.ofFn_const, ← endpointKernel_eq_path, endpointKernel_empty]

theorem endpointKernel_nonneg {R : Nat} {z : ℝ} (hz : 0 ≤ z)
    (xs : Fin R → Input) (q r : State) : 0 ≤ endpointKernel z xs q r := by
  apply Spin.FinPMF.expect_nonneg
  intro ts
  split_ifs <;> positivity

theorem endpointKernel_total {R : Nat} (z : ℝ) (xs : Fin R → Input) (q : State) :
    ∑ r, endpointKernel z xs q r = inputMoment z xs q := by
  unfold endpointKernel inputMoment
  rw [← Spin.FinPMF.expect_sum]
  apply congrArg (Spin.FinPMF.expect _)
  funext ts
  rw [← Finset.mul_sum]
  simp

end Spin.Structured.ConcreteEncoder
