import SpinCodes.UniformMixer.QuadraticRankHistogramSum

/-!
# Walsh behavior within each quadratic-radical class

For a fixed quadratic part, Parseval and the radical character identity
fix both the number of nonzero Walsh coefficients and their magnitude.
The constant coefficient later pairs each weight with its complement.
-/

namespace Spin.UniformMixer

theorem activeWalsh_card_radical_16 (c : Fin 15 → ZMod 2)
    (hr : (quadraticRadical c).card = 16) :
    (activeWalsh c).card = 4 := by
  have h := quadraticPart_activeWalsh_count c
  rw [hr] at h
  norm_num only [Nat.cast_ofNat] at h
  omega

theorem walsh_value_radical_16 (c : Fin 15 → ZMod 2)
    (hr : (quadraticRadical c).card = 16)
    (a : State6) (ha : a ∈ activeWalsh c) :
    walsh6 (quadraticPart c) a = 32 ∨
      walsh6 (quadraticPart c) a = -32 := by
  have h := quadraticPart_walsh_square_nonzero c ha
  rw [hr] at h
  have hs : walsh6 (quadraticPart c) a *
      walsh6 (quadraticPart c) a = (32 : ℤ) * 32 := by
    norm_num at h ⊢
    exact h
  exact mul_self_eq_mul_self_iff.mp hs

theorem activeWalsh_card_radical_64 (c : Fin 15 → ZMod 2)
    (hr : (quadraticRadical c).card = 64) :
    (activeWalsh c).card = 1 := by
  have h := quadraticPart_activeWalsh_count c
  rw [hr] at h
  norm_num only [Nat.cast_ofNat] at h
  omega

theorem walsh_value_radical_64 (c : Fin 15 → ZMod 2)
    (hr : (quadraticRadical c).card = 64)
    (a : State6) (ha : a ∈ activeWalsh c) :
    walsh6 (quadraticPart c) a = 64 ∨
      walsh6 (quadraticPart c) a = -64 := by
  have h := quadraticPart_walsh_square_nonzero c ha
  rw [hr] at h
  have hs : walsh6 (quadraticPart c) a *
      walsh6 (quadraticPart c) a = (64 : ℤ) * 64 := by
    norm_num at h ⊢
    exact h
  exact mul_self_eq_mul_self_iff.mp hs

theorem activeWalsh_card_radical_4 (c : Fin 15 → ZMod 2)
    (hr : (quadraticRadical c).card = 4) :
    (activeWalsh c).card = 16 := by
  have h := quadraticPart_activeWalsh_count c
  rw [hr] at h
  norm_num only [Nat.cast_ofNat] at h
  omega

theorem walsh_value_radical_4 (c : Fin 15 → ZMod 2)
    (hr : (quadraticRadical c).card = 4)
    (a : State6) (ha : a ∈ activeWalsh c) :
    walsh6 (quadraticPart c) a = 16 ∨
      walsh6 (quadraticPart c) a = -16 := by
  have h := quadraticPart_walsh_square_nonzero c ha
  rw [hr] at h
  have hs : walsh6 (quadraticPart c) a *
      walsh6 (quadraticPart c) a = (16 : ℤ) * 16 := by
    norm_num at h ⊢
    exact h
  exact mul_self_eq_mul_self_iff.mp hs

theorem activeWalsh_card_radical_1 (c : Fin 15 → ZMod 2)
    (hr : (quadraticRadical c).card = 1) :
    (activeWalsh c).card = 64 := by
  have h := quadraticPart_activeWalsh_count c
  rw [hr] at h
  norm_num only [Nat.cast_ofNat] at h
  omega

theorem walsh_value_radical_1 (c : Fin 15 → ZMod 2)
    (hr : (quadraticRadical c).card = 1)
    (a : State6) (ha : a ∈ activeWalsh c) :
    walsh6 (quadraticPart c) a = 8 ∨
      walsh6 (quadraticPart c) a = -8 := by
  have h := quadraticPart_walsh_square_nonzero c ha
  rw [hr] at h
  have hs : walsh6 (quadraticPart c) a *
      walsh6 (quadraticPart c) a = (8 : ℤ) * 8 := by
    norm_num at h ⊢
    exact h
  exact mul_self_eq_mul_self_iff.mp hs

theorem weight_of_radical_16 (c : Fin 15 → ZMod 2)
    (hr : (quadraticRadical c).card = 16) (a : State6) :
    boolWeight6 (fun x => quadraticPart c x + dot6Linear x a) = 16 ∨
      boolWeight6 (fun x => quadraticPart c x + dot6Linear x a) = 32 ∨
      boolWeight6 (fun x => quadraticPart c x + dot6Linear x a) = 48 := by
  let w := boolWeight6 (fun x => quadraticPart c x + dot6Linear x a)
  have hw : walsh6 (quadraticPart c) a = 64 - 2 * (w : ℤ) :=
    walsh6_weight (quadraticPart c) a
  by_cases ha : a ∈ activeWalsh c
  · rcases walsh_value_radical_16 c hr a ha with hp | hn
    · left
      omega
    · right; right
      omega
  · have hz : walsh6 (quadraticPart c) a = 0 := by
      simpa [activeWalsh] using ha
    right; left
    omega

theorem weight_of_radical_64 (c : Fin 15 → ZMod 2)
    (hr : (quadraticRadical c).card = 64) (a : State6) :
    boolWeight6 (fun x => quadraticPart c x + dot6Linear x a) = 0 ∨
      boolWeight6 (fun x => quadraticPart c x + dot6Linear x a) = 32 ∨
      boolWeight6 (fun x => quadraticPart c x + dot6Linear x a) = 64 := by
  let w := boolWeight6 (fun x => quadraticPart c x + dot6Linear x a)
  have hw : walsh6 (quadraticPart c) a = 64 - 2 * (w : ℤ) :=
    walsh6_weight (quadraticPart c) a
  by_cases ha : a ∈ activeWalsh c
  · rcases walsh_value_radical_64 c hr a ha with hp | hn
    · left
      omega
    · right; right
      omega
  · have hz : walsh6 (quadraticPart c) a = 0 := by
      simpa [activeWalsh] using ha
    right; left
    omega

theorem weight_of_radical_4 (c : Fin 15 → ZMod 2)
    (hr : (quadraticRadical c).card = 4) (a : State6) :
    boolWeight6 (fun x => quadraticPart c x + dot6Linear x a) = 24 ∨
      boolWeight6 (fun x => quadraticPart c x + dot6Linear x a) = 32 ∨
      boolWeight6 (fun x => quadraticPart c x + dot6Linear x a) = 40 := by
  let w := boolWeight6 (fun x => quadraticPart c x + dot6Linear x a)
  have hw : walsh6 (quadraticPart c) a = 64 - 2 * (w : ℤ) :=
    walsh6_weight (quadraticPart c) a
  by_cases ha : a ∈ activeWalsh c
  · rcases walsh_value_radical_4 c hr a ha with hp | hn
    · left
      omega
    · right; right
      omega
  · have hz : walsh6 (quadraticPart c) a = 0 := by
      simpa [activeWalsh] using ha
    right; left
    omega

theorem weight_of_radical_1 (c : Fin 15 → ZMod 2)
    (hr : (quadraticRadical c).card = 1) (a : State6) :
    boolWeight6 (fun x => quadraticPart c x + dot6Linear x a) = 28 ∨
      boolWeight6 (fun x => quadraticPart c x + dot6Linear x a) = 32 ∨
      boolWeight6 (fun x => quadraticPart c x + dot6Linear x a) = 36 := by
  let w := boolWeight6 (fun x => quadraticPart c x + dot6Linear x a)
  have hw : walsh6 (quadraticPart c) a = 64 - 2 * (w : ℤ) :=
    walsh6_weight (quadraticPart c) a
  by_cases ha : a ∈ activeWalsh c
  · rcases walsh_value_radical_1 c hr a ha with hp | hn
    · left
      omega
    · right; right
      omega
  · have hz : walsh6 (quadraticPart c) a = 0 := by
      simpa [activeWalsh] using ha
    right; left
    omega

end Spin.UniformMixer
