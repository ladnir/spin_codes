/-
Zero extension to a requested length (`eq:structured-requested-length-wrapper`).

    Ẽ_n(x) := E_{m(n)}(x) ‖ 0^{n - N_{m(n)}}

"This zero extension preserves injectivity and minimum distance."  That is the
content proved here, and it is where the wrapper's correctness actually lives.

The companion rate claim — `N_{m(n)}/n = 1 - o(1)`, hence rate `1/2 - o(1)` —
is a statement about *consecutive* schedule lengths having relative gap `o(1)`,
not about the extension.  It is queued separately (T2b); nothing here assumes
it.
-/
import Mathlib

set_option linter.unusedSectionVars false

namespace Spin

open Finset

/-- Hamming weight of a binary word. -/
def wt {n : ℕ} (w : Fin n → ZMod 2) : ℕ := (univ.filter (fun i => w i ≠ 0)).card

/-- Append `n - N` zeros. -/
def zeroExtend {N n : ℕ} (h : N ≤ n) (w : Fin N → ZMod 2) : Fin n → ZMod 2 :=
  fun i => if hi : (i : ℕ) < N then w ⟨i, hi⟩ else 0

lemma zeroExtend_of_lt {N n : ℕ} (h : N ≤ n) (w : Fin N → ZMod 2) {i : Fin n}
    (hi : (i : ℕ) < N) : zeroExtend h w i = w ⟨i, hi⟩ := dif_pos hi

lemma zeroExtend_of_ge {N n : ℕ} (h : N ≤ n) (w : Fin N → ZMod 2) {i : Fin n}
    (hi : ¬ (i : ℕ) < N) : zeroExtend h w i = 0 := dif_neg hi

/-- Zero extension preserves Hamming weight: the support is carried across by
`Fin.castLE`. -/
theorem wt_zeroExtend {N n : ℕ} (h : N ≤ n) (w : Fin N → ZMod 2) :
    wt (zeroExtend h w) = wt w := by
  classical
  unfold wt
  have hset : (univ.filter (fun j : Fin n => zeroExtend h w j ≠ 0))
      = (univ.filter (fun i : Fin N => w i ≠ 0)).map (Fin.castLEEmb h) := by
    ext j
    simp only [Finset.mem_filter, Finset.mem_univ, true_and, Finset.mem_map,
      Fin.castLEEmb_apply]
    constructor
    · intro hj
      by_cases hlt : (j : ℕ) < N
      · exact ⟨⟨j, hlt⟩, by rwa [zeroExtend_of_lt h w hlt] at hj, by
          apply Fin.ext; simp [Fin.castLE]⟩
      · exact absurd (zeroExtend_of_ge h w hlt) hj
    · rintro ⟨i, hi, rfl⟩
      have hlt : ((Fin.castLE h i : Fin n) : ℕ) < N := by simpa [Fin.castLE] using i.isLt
      rw [zeroExtend_of_lt h w hlt]
      simpa [Fin.castLE] using hi
  rw [hset, Finset.card_map]

theorem zeroExtend_injective {N n : ℕ} (h : N ≤ n) :
    Function.Injective (zeroExtend (N := N) (n := n) h) := by
  intro u v huv
  funext i
  have := congrFun huv (Fin.castLE h i)
  have hlt : ((Fin.castLE h i : Fin n) : ℕ) < N := by simpa [Fin.castLE] using i.isLt
  rw [zeroExtend_of_lt h u hlt, zeroExtend_of_lt h v hlt] at this
  simpa [Fin.castLE] using this

/-- **The wrapper is sound.**  Zero extension preserves injectivity of the
encoder and the weight of every codeword — hence the minimum distance. -/
theorem zeroExtend_preserves {M : Type*} {N n : ℕ} (h : N ≤ n)
    (E : M → (Fin N → ZMod 2)) :
    (Function.Injective (fun x => zeroExtend h (E x)) ↔ Function.Injective E)
      ∧ ∀ x, wt (zeroExtend h (E x)) = wt (E x) := by
  refine ⟨⟨fun hinj x y hxy => hinj (by simp [hxy]), fun hinj x y hxy => ?_⟩,
    fun x => wt_zeroExtend h (E x)⟩
  exact hinj (zeroExtend_injective h hxy)

end Spin
