import SpinCodes.UniformMixer.QuadraticRankCountBridge

/-!
# One-pass bounded histogram of six-bit quadratic radicals

The four counters are checked together so each quadratic form needs only
one rank computation. Projection lemmas connect them to the separately
defined and semantically justified range counts.
-/

namespace Spin.UniformMixer

abbrev RankHistogram := Nat × Nat × Nat × Nat

def rankHistogramStep (q : Nat) (h : RankHistogram) : RankHistogram :=
  let size := fastRadicalSize q
  (h.1 + if size = 64 then 1 else 0,
   h.2.1 + if size = 16 then 1 else 0,
   h.2.2.1 + if size = 4 then 1 else 0,
   h.2.2.2 + if size = 1 then 1 else 0)

def rankHistogramList (qs : List Nat) : RankHistogram :=
  qs.foldr rankHistogramStep (0, 0, 0, 0)

def fastRankHistogramBlock (start length : Nat) : RankHistogram :=
  rankHistogramList (List.range' start length)

theorem rankHistogramList_first (qs : List Nat) :
    (rankHistogramList qs).1 =
      (qs.filter (fun q => fastRadicalSize q = 64)).length := by
  induction qs with
  | nil => rfl
  | cons q qs ih =>
      simp only [rankHistogramList, List.foldr_cons, rankHistogramStep]
      change (rankHistogramList qs).1 +
        (if fastRadicalSize q = 64 then 1 else 0) = _
      by_cases h : fastRadicalSize q = 64 <;> simp [h, ih, add_comm]

theorem rankHistogramList_second (qs : List Nat) :
    (rankHistogramList qs).2.1 =
      (qs.filter (fun q => fastRadicalSize q = 16)).length := by
  induction qs with
  | nil => rfl
  | cons q qs ih =>
      simp only [rankHistogramList, List.foldr_cons, rankHistogramStep]
      change (rankHistogramList qs).2.1 +
        (if fastRadicalSize q = 16 then 1 else 0) = _
      by_cases h : fastRadicalSize q = 16 <;> simp [h, ih, add_comm]

theorem rankHistogramList_third (qs : List Nat) :
    (rankHistogramList qs).2.2.1 =
      (qs.filter (fun q => fastRadicalSize q = 4)).length := by
  induction qs with
  | nil => rfl
  | cons q qs ih =>
      simp only [rankHistogramList, List.foldr_cons, rankHistogramStep]
      change (rankHistogramList qs).2.2.1 +
        (if fastRadicalSize q = 4 then 1 else 0) = _
      by_cases h : fastRadicalSize q = 4 <;> simp [h, ih, add_comm]

theorem rankHistogramList_fourth (qs : List Nat) :
    (rankHistogramList qs).2.2.2 =
      (qs.filter (fun q => fastRadicalSize q = 1)).length := by
  induction qs with
  | nil => rfl
  | cons q qs ih =>
      simp only [rankHistogramList, List.foldr_cons, rankHistogramStep]
      change (rankHistogramList qs).2.2.2 +
        (if fastRadicalSize q = 1 then 1 else 0) = _
      by_cases h : fastRadicalSize q = 1 <;> simp [h, ih, add_comm]

theorem fastRankHistogramBlock_first (start length : Nat) :
    (fastRankHistogramBlock start length).1 =
      fastRadicalRangeCount start length 64 := by
  exact rankHistogramList_first _

theorem fastRankHistogramBlock_second (start length : Nat) :
    (fastRankHistogramBlock start length).2.1 =
      fastRadicalRangeCount start length 16 := by
  exact rankHistogramList_second _

theorem fastRankHistogramBlock_third (start length : Nat) :
    (fastRankHistogramBlock start length).2.2.1 =
      fastRadicalRangeCount start length 4 := by
  exact rankHistogramList_third _

theorem fastRankHistogramBlock_fourth (start length : Nat) :
    (fastRankHistogramBlock start length).2.2.2 =
      fastRadicalRangeCount start length 1 := by
  exact rankHistogramList_fourth _

def addRankHistogram (a b : RankHistogram) : RankHistogram :=
  (a.1 + b.1, a.2.1 + b.2.1,
    a.2.2.1 + b.2.2.1, a.2.2.2 + b.2.2.2)

theorem fastRankHistogramBlock_add (start m n : Nat) :
    fastRankHistogramBlock start (m + n) =
      addRankHistogram (fastRankHistogramBlock start m)
        (fastRankHistogramBlock (start + m) n) := by
  apply Prod.ext
  · change (fastRankHistogramBlock start (m + n)).1 =
      (fastRankHistogramBlock start m).1 +
        (fastRankHistogramBlock (start + m) n).1
    rw [fastRankHistogramBlock_first, fastRankHistogramBlock_first,
      fastRankHistogramBlock_first, fastRadicalRangeCount_add]
  · apply Prod.ext
    · change (fastRankHistogramBlock start (m + n)).2.1 =
        (fastRankHistogramBlock start m).2.1 +
          (fastRankHistogramBlock (start + m) n).2.1
      rw [fastRankHistogramBlock_second, fastRankHistogramBlock_second,
        fastRankHistogramBlock_second, fastRadicalRangeCount_add]
    · apply Prod.ext
      · change (fastRankHistogramBlock start (m + n)).2.2.1 =
          (fastRankHistogramBlock start m).2.2.1 +
            (fastRankHistogramBlock (start + m) n).2.2.1
        rw [fastRankHistogramBlock_third, fastRankHistogramBlock_third,
          fastRankHistogramBlock_third, fastRadicalRangeCount_add]
      · change (fastRankHistogramBlock start (m + n)).2.2.2 =
          (fastRankHistogramBlock start m).2.2.2 +
            (fastRankHistogramBlock (start + m) n).2.2.2
        rw [fastRankHistogramBlock_fourth, fastRankHistogramBlock_fourth,
          fastRankHistogramBlock_fourth, fastRadicalRangeCount_add]

set_option maxRecDepth 100000 in
set_option maxHeartbeats 0 in
theorem first256_fast_histogram :
    fastRankHistogramBlock 0 256 = (1, 59, 196, 0) := by
  decide +kernel

set_option maxRecDepth 100000 in
set_option maxHeartbeats 0 in
theorem first512_fast_histogram :
    fastRankHistogramBlock 0 512 = (1, 91, 420, 0) := by
  decide +kernel

end Spin.UniformMixer
