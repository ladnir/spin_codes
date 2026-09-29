import SpinCodes.Majorant.Segment0
import SpinCodes.Majorant.Segment1
import SpinCodes.Majorant.Segment2
import SpinCodes.Majorant.Segment3
import SpinCodes.Majorant.Segment4
import SpinCodes.Majorant.Segment5
import SpinCodes.Majorant.Segment6
import SpinCodes.Majorant.Segment7
import SpinCodes.Majorant.Segment8
import SpinCodes.Majorant.Segment9
import SpinCodes.Majorant.Segment10
import SpinCodes.Majorant.Segment11
import SpinCodes.Majorant.Segment12
import SpinCodes.Majorant.Segment13
import SpinCodes.Majorant.Segment14
import SpinCodes.Majorant.Segment15
import SpinCodes.Majorant.Segment16
import SpinCodes.Majorant.Segment17
import SpinCodes.Majorant.Segment18
import SpinCodes.Majorant.RefinedData
import SpinCodes.Numeric.BACentral
namespace Spin.Majorant
open Spin.Numeric Spin.Numeric.Fix
theorem segment0_real {a b w : ℝ}
    (ha : a ∈ Set.Icc (0:ℝ) 1) (hb : b ∈ Set.Icc (0:ℝ) 1)
    (hw : w ∈ Set.Icc ((13 / 125):ℝ) (18296026121 / 120000000000)) :
    MajorantClaim (833 / 500) (-43311 / 250000) a b w := by
  apply S0.cover (s := (833 / 500)) (i := (-43311 / 250000))
  · norm_num [Mem, scale]
  · norm_num [Mem, scale]
  · simpa [Mem, scale] using ha
  · simpa [Mem, scale] using hb
  · constructor <;> norm_num [scale] <;> linarith [hw.1, hw.2]
theorem segment1_real {a b w : ℝ}
    (ha : a ∈ Set.Icc (0:ℝ) 1) (hb : b ∈ Set.Icc (0:ℝ) 1)
    (hw : w ∈ Set.Icc ((18296026121 / 120000000000):ℝ) (791599208623 / 5000000000000)) :
    MajorantClaim (83 / 50) (-3446583973879 / 20000000000000) a b w := by
  apply S1.cover (s := (83 / 50)) (i := (-3446583973879 / 20000000000000))
  · norm_num [Mem, scale]
  · norm_num [Mem, scale]
  · simpa [Mem, scale] using ha
  · simpa [Mem, scale] using hb
  · constructor <;> norm_num [scale] <;> linarith [hw.1, hw.2]
theorem segment2_real {a b w : ℝ}
    (ha : a ∈ Set.Icc (0:ℝ) 1) (hb : b ∈ Set.Icc (0:ℝ) 1)
    (hw : w ∈ Set.Icc ((791599208623 / 5000000000000):ℝ) (1609032935677 / 10000000000000)) :
    MajorantClaim (33 / 20) (-5335812508647 / 31250000000000) a b w := by
  apply S2.cover (s := (33 / 20)) (i := (-5335812508647 / 31250000000000))
  · norm_num [Mem, scale]
  · norm_num [Mem, scale]
  · simpa [Mem, scale] using ha
  · simpa [Mem, scale] using hb
  · constructor <;> norm_num [scale] <;> linarith [hw.1, hw.2]
theorem segment3_real {a b w : ℝ}
    (ha : a ∈ Set.Icc (0:ℝ) 1) (hb : b ∈ Set.Icc (0:ℝ) 1)
    (hw : w ∈ Set.Icc ((1609032935677 / 10000000000000):ℝ) (2469788513329 / 15000000000000)) :
    MajorantClaim (163 / 100) (-3350558688107 / 20000000000000) a b w := by
  apply S3.cover (s := (163 / 100)) (i := (-3350558688107 / 20000000000000))
  · norm_num [Mem, scale]
  · norm_num [Mem, scale]
  · simpa [Mem, scale] using ha
  · simpa [Mem, scale] using hb
  · constructor <;> norm_num [scale] <;> linarith [hw.1, hw.2]
theorem segment4_real {a b w : ℝ}
    (ha : a ∈ Set.Icc (0:ℝ) 1) (hb : b ∈ Set.Icc (0:ℝ) 1)
    (hw : w ∈ Set.Icc ((2469788513329 / 15000000000000):ℝ) (4271577070219 / 25000000000000)) :
    MajorantClaim (8 / 5) (-40647089344673 / 250000000000000) a b w := by
  apply S4.cover (s := (8 / 5)) (i := (-40647089344673 / 250000000000000))
  · norm_num [Mem, scale]
  · norm_num [Mem, scale]
  · simpa [Mem, scale] using ha
  · simpa [Mem, scale] using hb
  · constructor <;> norm_num [scale] <;> linarith [hw.1, hw.2]
theorem segment5_real {a b w : ℝ}
    (ha : a ∈ Set.Icc (0:ℝ) 1) (hb : b ∈ Set.Icc (0:ℝ) 1)
    (hw : w ∈ Set.Icc ((4271577070219 / 25000000000000):ℝ) (13956016256281 / 75000000000000)) :
    MajorantClaim (31 / 20) (-77022601619127 / 500000000000000) a b w := by
  apply S5.cover (s := (31 / 20)) (i := (-77022601619127 / 500000000000000))
  · norm_num [Mem, scale]
  · norm_num [Mem, scale]
  · simpa [Mem, scale] using ha
  · simpa [Mem, scale] using hb
  · constructor <;> norm_num [scale] <;> linarith [hw.1, hw.2]
theorem segment6_real {a b w : ℝ}
    (ha : a ∈ Set.Icc (0:ℝ) 1) (hb : b ∈ Set.Icc (0:ℝ) 1)
    (hw : w ∈ Set.Icc ((13956016256281 / 75000000000000):ℝ) (428613422794653 / 2000000000000000)) :
    MajorantClaim (7 / 5) (-31533292681423 / 250000000000000) a b w := by
  apply S6.cover (s := (7 / 5)) (i := (-31533292681423 / 250000000000000))
  · norm_num [Mem, scale]
  · norm_num [Mem, scale]
  · simpa [Mem, scale] using ha
  · simpa [Mem, scale] using hb
  · constructor <;> norm_num [scale] <;> linarith [hw.1, hw.2]
theorem segment7_real {a b w : ℝ}
    (ha : a ∈ Set.Icc (0:ℝ) 1) (hb : b ∈ Set.Icc (0:ℝ) 1)
    (hw : w ∈ Set.Icc ((428613422794653 / 2000000000000000):ℝ) (499776708257287 / 2000000000000000)) :
    MajorantClaim (6 / 5) (-832718284462267 / 10000000000000000) a b w := by
  apply S7.cover (s := (6 / 5)) (i := (-832718284462267 / 10000000000000000))
  · norm_num [Mem, scale]
  · norm_num [Mem, scale]
  · simpa [Mem, scale] using ha
  · simpa [Mem, scale] using hb
  · constructor <;> norm_num [scale] <;> linarith [hw.1, hw.2]
theorem segment8_real {a b w : ℝ}
    (ha : a ∈ Set.Icc (0:ℝ) 1) (hb : b ∈ Set.Icc (0:ℝ) 1)
    (hw : w ∈ Set.Icc ((499776708257287 / 2000000000000000):ℝ) (144603079418107 / 500000000000000)) :
    MajorantClaim (1 / 1) (-16647078810249 / 500000000000000) a b w := by
  apply S8.cover (s := (1 / 1)) (i := (-16647078810249 / 500000000000000))
  · norm_num [Mem, scale]
  · norm_num [Mem, scale]
  · simpa [Mem, scale] using ha
  · simpa [Mem, scale] using hb
  · constructor <;> norm_num [scale] <;> linarith [hw.1, hw.2]
theorem segment9_real {a b w : ℝ}
    (ha : a ∈ Set.Icc (0:ℝ) 1) (hb : b ∈ Set.Icc (0:ℝ) 1)
    (hw : w ∈ Set.Icc ((144603079418107 / 500000000000000):ℝ) (663870278625033 / 2000000000000000)) :
    MajorantClaim (4 / 5) (30683842683431 / 1250000000000000) a b w := by
  apply S9.cover (s := (4 / 5)) (i := (30683842683431 / 1250000000000000))
  · norm_num [Mem, scale]
  · norm_num [Mem, scale]
  · simpa [Mem, scale] using ha
  · simpa [Mem, scale] using hb
  · constructor <;> norm_num [scale] <;> linarith [hw.1, hw.2]
theorem segment10_real {a b w : ℝ}
    (ha : a ∈ Set.Icc (0:ℝ) 1) (hb : b ∈ Set.Icc (0:ℝ) 1)
    (hw : w ∈ Set.Icc ((663870278625033 / 2000000000000000):ℝ) (755265762055919 / 2000000000000000)) :
    MajorantClaim (3 / 5) (909341020092481 / 10000000000000000) a b w := by
  apply S10.cover (s := (3 / 5)) (i := (909341020092481 / 10000000000000000))
  · norm_num [Mem, scale]
  · norm_num [Mem, scale]
  · simpa [Mem, scale] using ha
  · simpa [Mem, scale] using hb
  · constructor <;> norm_num [scale] <;> linarith [hw.1, hw.2]
theorem segment11_real {a b w : ℝ}
    (ha : a ∈ Set.Icc (0:ℝ) 1) (hb : b ∈ Set.Icc (0:ℝ) 1)
    (hw : w ∈ Set.Icc ((755265762055919 / 2000000000000000):ℝ) (10640568152347 / 25000000000000)) :
    MajorantClaim (2 / 5) (4161516955371 / 25000000000000) a b w := by
  apply S11.cover (s := (2 / 5)) (i := (4161516955371 / 25000000000000))
  · norm_num [Mem, scale]
  · norm_num [Mem, scale]
  · simpa [Mem, scale] using ha
  · simpa [Mem, scale] using hb
  · constructor <;> norm_num [scale] <;> linarith [hw.1, hw.2]
theorem segment12_real {a b w : ℝ}
    (ha : a ∈ Set.Icc (0:ℝ) 1) (hb : b ∈ Set.Icc (0:ℝ) 1)
    (hw : w ∈ Set.Icc ((10640568152347 / 25000000000000):ℝ) (46255924070167 / 100000000000000)) :
    MajorantClaim (1 / 5) (15724076464601 / 62500000000000) a b w := by
  apply S12.cover (s := (1 / 5)) (i := (15724076464601 / 62500000000000))
  · norm_num [Mem, scale]
  · norm_num [Mem, scale]
  · simpa [Mem, scale] using ha
  · simpa [Mem, scale] using hb
  · constructor <;> norm_num [scale] <;> linarith [hw.1, hw.2]
theorem segment13_real {a b w : ℝ}
    (ha : a ∈ Set.Icc (0:ℝ) 1) (hb : b ∈ Set.Icc (0:ℝ) 1)
    (hw : w ∈ Set.Icc ((46255924070167 / 100000000000000):ℝ) (49374602183332977403672480068036683319333956789 / 103845937170696552570609926584401920000000000000)) :
    MajorantClaim (1 / 10) (297841147503783 / 1000000000000000) a b w := by
  apply S13.cover (s := (1 / 10)) (i := (297841147503783 / 1000000000000000))
  · norm_num [Mem, scale]
  · norm_num [Mem, scale]
  · simpa [Mem, scale] using ha
  · simpa [Mem, scale] using hb
  · constructor <;> norm_num [scale] <;> linarith [hw.1, hw.2]
theorem segment14_real {a b w : ℝ}
    (ha : a ∈ Set.Icc (0:ℝ) 1) (hb : b ∈ Set.Icc (0:ℝ) 1)
    (hw : w ∈ Set.Icc ((49374602183332977403672480068036683319333956789 / 103845937170696552570609926584401920000000000000):ℝ) (81669610035739215919018503711223740925 / 170141183460469231731687303715884105728)) :
    MajorantClaim (9 / 100) (102967997603379096932443767333991244233 / 340282366920938463463374607431768211456) a b w := by
  apply S14.cover (s := (9 / 100)) (i := (102967997603379096932443767333991244233 / 340282366920938463463374607431768211456))
  · norm_num [Mem, scale]
  · norm_num [Mem, scale]
  · simpa [Mem, scale] using ha
  · simpa [Mem, scale] using hb
  · constructor <;> norm_num [scale] <;> linarith [hw.1, hw.2]
theorem segment15_real {a b w : ℝ}
    (ha : a ∈ Set.Icc (0:ℝ) 1) (hb : b ∈ Set.Icc (0:ℝ) 1)
    (hw : w ∈ Set.Icc ((81669610035739215919018503711223740925 / 170141183460469231731687303715884105728):ℝ) (82519260580058937111451454383981630225 / 170141183460469231731687303715884105728)) :
    MajorantClaim (7 / 100) (53117391002404332784602253741220096935 / 170141183460469231731687303715884105728) a b w := by
  apply S15.cover (s := (7 / 100)) (i := (53117391002404332784602253741220096935 / 170141183460469231731687303715884105728))
  · norm_num [Mem, scale]
  · norm_num [Mem, scale]
  · simpa [Mem, scale] using ha
  · simpa [Mem, scale] using hb
  · constructor <;> norm_num [scale] <;> linarith [hw.1, hw.2]
theorem segment16_real {a b w : ℝ}
    (ha : a ∈ Set.Icc (0:ℝ) 1) (hb : b ∈ Set.Icc (0:ℝ) 1)
    (hw : w ∈ Set.Icc ((82519260580058937111451454383981630225 / 170141183460469231731687303715884105728):ℝ) (41684710442490683359023934633467251675 / 85070591730234615865843651857942052864)) :
    MajorantClaim (1 / 20) (109535552428011023053662565657799459079 / 340282366920938463463374607431768211456) a b w := by
  apply S16.cover (s := (1 / 20)) (i := (109535552428011023053662565657799459079 / 340282366920938463463374607431768211456))
  · norm_num [Mem, scale]
  · norm_num [Mem, scale]
  · simpa [Mem, scale] using ha
  · simpa [Mem, scale] using hb
  · constructor <;> norm_num [scale] <;> linarith [hw.1, hw.2]
theorem segment17_real {a b w : ℝ}
    (ha : a ∈ Set.Icc (0:ℝ) 1) (hb : b ∈ Set.Icc (0:ℝ) 1)
    (hw : w ∈ Set.Icc ((41684710442490683359023934633467251675 / 85070591730234615865843651857942052864):ℝ) (84219921256861908840067417679371902625 / 170141183460469231731687303715884105728)) :
    MajorantClaim (3 / 100) (112870329263410277722384480428476839213 / 340282366920938463463374607431768211456) a b w := by
  apply S17.cover (s := (3 / 100)) (i := (112870329263410277722384480428476839213 / 340282366920938463463374607431768211456))
  · norm_num [Mem, scale]
  · norm_num [Mem, scale]
  · simpa [Mem, scale] using ha
  · simpa [Mem, scale] using hb
  · constructor <;> norm_num [scale] <;> linarith [hw.1, hw.2]
theorem segment18_real {a b w : ℝ}
    (ha : a ∈ Set.Icc (0:ℝ) 1) (hb : b ∈ Set.Icc (0:ℝ) 1)
    (hw : w ∈ Set.Icc ((84219921256861908840067417679371902625 / 170141183460469231731687303715884105728):ℝ) (129223289418938324344194200630890939964085310021 / 259614842926741381426524816461004800000000000000)) :
    MajorantClaim (1 / 100) (58119563056842377037993588567825857659 / 170141183460469231731687303715884105728) a b w := by
  apply S18.cover (s := (1 / 100)) (i := (58119563056842377037993588567825857659 / 170141183460469231731687303715884105728))
  · norm_num [Mem, scale]
  · norm_num [Mem, scale]
  · simpa [Mem, scale] using ha
  · simpa [Mem, scale] using hb
  · constructor <;> norm_num [scale] <;> linarith [hw.1, hw.2]

theorem majorant_left {a b w : ℝ}
    (ha : a ∈ Set.Icc (0:ℝ) 1) (hb : b ∈ Set.Icc (0:ℝ) 1)
    (hw : w ∈ Set.Icc ((13:ℝ)/125) (1/2))
    (hf : a/2 ≤ b ∧ b ≤ 1-a/2 ∧ b/2 ≤ w ∧ w ≤ 1-b/2 ∧ 0 ≤ a ∧ a ≤ 1) :
    gBA a + piBA a b + piBA b w ≤ refined.toFun w := by
  by_cases h0 : w ≤ ((18296026121 / 120000000000):ℝ)
  · have hwj : w ∈ Set.Icc ((13 / 125):ℝ) (18296026121 / 120000000000) := by constructor <;> linarith [hw.1]
    exact le_trans (segment0_real ha hb hwj hf) (active_real0 hwj)
  by_cases h1 : w ≤ ((791599208623 / 5000000000000):ℝ)
  · have hwj : w ∈ Set.Icc ((18296026121 / 120000000000):ℝ) (791599208623 / 5000000000000) := by constructor <;> linarith [hw.1]
    exact le_trans (segment1_real ha hb hwj hf) (active_real1 hwj)
  by_cases h2 : w ≤ ((1609032935677 / 10000000000000):ℝ)
  · have hwj : w ∈ Set.Icc ((791599208623 / 5000000000000):ℝ) (1609032935677 / 10000000000000) := by constructor <;> linarith [hw.1]
    exact le_trans (segment2_real ha hb hwj hf) (active_real2 hwj)
  by_cases h3 : w ≤ ((2469788513329 / 15000000000000):ℝ)
  · have hwj : w ∈ Set.Icc ((1609032935677 / 10000000000000):ℝ) (2469788513329 / 15000000000000) := by constructor <;> linarith [hw.1]
    exact le_trans (segment3_real ha hb hwj hf) (active_real3 hwj)
  by_cases h4 : w ≤ ((4271577070219 / 25000000000000):ℝ)
  · have hwj : w ∈ Set.Icc ((2469788513329 / 15000000000000):ℝ) (4271577070219 / 25000000000000) := by constructor <;> linarith [hw.1]
    exact le_trans (segment4_real ha hb hwj hf) (active_real4 hwj)
  by_cases h5 : w ≤ ((13956016256281 / 75000000000000):ℝ)
  · have hwj : w ∈ Set.Icc ((4271577070219 / 25000000000000):ℝ) (13956016256281 / 75000000000000) := by constructor <;> linarith [hw.1]
    exact le_trans (segment5_real ha hb hwj hf) (active_real5 hwj)
  by_cases h6 : w ≤ ((428613422794653 / 2000000000000000):ℝ)
  · have hwj : w ∈ Set.Icc ((13956016256281 / 75000000000000):ℝ) (428613422794653 / 2000000000000000) := by constructor <;> linarith [hw.1]
    exact le_trans (segment6_real ha hb hwj hf) (active_real6 hwj)
  by_cases h7 : w ≤ ((499776708257287 / 2000000000000000):ℝ)
  · have hwj : w ∈ Set.Icc ((428613422794653 / 2000000000000000):ℝ) (499776708257287 / 2000000000000000) := by constructor <;> linarith [hw.1]
    exact le_trans (segment7_real ha hb hwj hf) (active_real7 hwj)
  by_cases h8 : w ≤ ((144603079418107 / 500000000000000):ℝ)
  · have hwj : w ∈ Set.Icc ((499776708257287 / 2000000000000000):ℝ) (144603079418107 / 500000000000000) := by constructor <;> linarith [hw.1]
    exact le_trans (segment8_real ha hb hwj hf) (active_real8 hwj)
  by_cases h9 : w ≤ ((663870278625033 / 2000000000000000):ℝ)
  · have hwj : w ∈ Set.Icc ((144603079418107 / 500000000000000):ℝ) (663870278625033 / 2000000000000000) := by constructor <;> linarith [hw.1]
    exact le_trans (segment9_real ha hb hwj hf) (active_real9 hwj)
  by_cases h10 : w ≤ ((755265762055919 / 2000000000000000):ℝ)
  · have hwj : w ∈ Set.Icc ((663870278625033 / 2000000000000000):ℝ) (755265762055919 / 2000000000000000) := by constructor <;> linarith [hw.1]
    exact le_trans (segment10_real ha hb hwj hf) (active_real10 hwj)
  by_cases h11 : w ≤ ((10640568152347 / 25000000000000):ℝ)
  · have hwj : w ∈ Set.Icc ((755265762055919 / 2000000000000000):ℝ) (10640568152347 / 25000000000000) := by constructor <;> linarith [hw.1]
    exact le_trans (segment11_real ha hb hwj hf) (active_real11 hwj)
  by_cases h12 : w ≤ ((46255924070167 / 100000000000000):ℝ)
  · have hwj : w ∈ Set.Icc ((10640568152347 / 25000000000000):ℝ) (46255924070167 / 100000000000000) := by constructor <;> linarith [hw.1]
    exact le_trans (segment12_real ha hb hwj hf) (active_real12 hwj)
  by_cases h13 : w ≤ ((49374602183332977403672480068036683319333956789 / 103845937170696552570609926584401920000000000000):ℝ)
  · have hwj : w ∈ Set.Icc ((46255924070167 / 100000000000000):ℝ) (49374602183332977403672480068036683319333956789 / 103845937170696552570609926584401920000000000000) := by constructor <;> linarith [hw.1]
    exact le_trans (segment13_real ha hb hwj hf) (active_real13 hwj)
  by_cases h14 : w ≤ ((81669610035739215919018503711223740925 / 170141183460469231731687303715884105728):ℝ)
  · have hwj : w ∈ Set.Icc ((49374602183332977403672480068036683319333956789 / 103845937170696552570609926584401920000000000000):ℝ) (81669610035739215919018503711223740925 / 170141183460469231731687303715884105728) := by constructor <;> linarith [hw.1]
    exact le_trans (segment14_real ha hb hwj hf) (active_real14 hwj)
  by_cases h15 : w ≤ ((82519260580058937111451454383981630225 / 170141183460469231731687303715884105728):ℝ)
  · have hwj : w ∈ Set.Icc ((81669610035739215919018503711223740925 / 170141183460469231731687303715884105728):ℝ) (82519260580058937111451454383981630225 / 170141183460469231731687303715884105728) := by constructor <;> linarith [hw.1]
    exact le_trans (segment15_real ha hb hwj hf) (active_real15 hwj)
  by_cases h16 : w ≤ ((41684710442490683359023934633467251675 / 85070591730234615865843651857942052864):ℝ)
  · have hwj : w ∈ Set.Icc ((82519260580058937111451454383981630225 / 170141183460469231731687303715884105728):ℝ) (41684710442490683359023934633467251675 / 85070591730234615865843651857942052864) := by constructor <;> linarith [hw.1]
    exact le_trans (segment16_real ha hb hwj hf) (active_real16 hwj)
  by_cases h17 : w ≤ ((84219921256861908840067417679371902625 / 170141183460469231731687303715884105728):ℝ)
  · have hwj : w ∈ Set.Icc ((41684710442490683359023934633467251675 / 85070591730234615865843651857942052864):ℝ) (84219921256861908840067417679371902625 / 170141183460469231731687303715884105728) := by constructor <;> linarith [hw.1]
    exact le_trans (segment17_real ha hb hwj hf) (active_real17 hwj)
  by_cases h18 : w ≤ ((129223289418938324344194200630890939964085310021 / 259614842926741381426524816461004800000000000000):ℝ)
  · have hwj : w ∈ Set.Icc ((84219921256861908840067417679371902625 / 170141183460469231731687303715884105728):ℝ) (129223289418938324344194200630890939964085310021 / 259614842926741381426524816461004800000000000000) := by constructor <;> linarith [hw.1]
    exact le_trans (segment18_real ha hb hwj hf) (active_real18 hwj)
  have hwj : w ∈ Set.Icc ((129223289418938324344194200630890939964085310021 / 259614842926741381426524816461004800000000000000):ℝ) (1 / 2) := by constructor <;> linarith [hw.2]
  have hc := le_trans (ba_central ha hf.1 hf.2.1 hf.2.2.1 hf.2.2.2.1) central_rounding
  have hactive := active_real19 hwj
  norm_num at hactive hc
  exact le_trans hc hactive

theorem majorant_pointwise {a b w : ℝ}
    (ha : a ∈ Set.Icc (0:ℝ) 1) (hb : b ∈ Set.Icc (0:ℝ) 1)
    (hw : w ∈ Set.Icc ((13:ℝ)/125) (112/125))
    (hf : a/2 ≤ b ∧ b ≤ 1-a/2 ∧ b/2 ≤ w ∧ w ≤ 1-b/2) :
    gBA a + piBA a b + piBA b w ≤ refined.toFun w := by
  by_cases h : w ≤ (1:ℝ)/2
  · exact majorant_left ha hb ⟨hw.1, h⟩ ⟨hf.1, hf.2.1, hf.2.2.1, hf.2.2.2, ha.1, ha.2⟩
  · have hwr : 1-w ∈ Set.Icc ((13:ℝ)/125) (1/2) := by constructor <;> linarith [hw.2]
    have hfr : a/2 ≤ b ∧ b ≤ 1-a/2 ∧ b/2 ≤ 1-w ∧ 1-w ≤ 1-b/2 ∧ 0 ≤ a ∧ a ≤ 1 := by
      exact ⟨hf.1, hf.2.1, by linarith [hf.2.2.2], by linarith [hf.2.2.1], ha.1, ha.2⟩
    have hr := majorant_left ha hb hwr hfr
    rwa [piBA_reflect, refined_reflect] at hr

/-- The paper's variational exponent. Restricting to feasible pairs implements
the paper's value -∞ for infeasible paths. -/
noncomputable def baExponent (w : ℝ) : ℝ := sSup
  ((fun p : ℝ × ℝ => gBA p.1 + piBA p.1 p.2 + piBA p.2 w) ''
    {p | p.1 ∈ Set.Icc (0:ℝ) 1 ∧ p.2 ∈ Set.Icc (0:ℝ) 1 ∧
      p.1/2 ≤ p.2 ∧ p.2 ≤ 1-p.1/2 ∧ p.2/2 ≤ w ∧ w ≤ 1-p.2/2})

/-- Equation ba-spectrum-majorant, with every affine support checked by the
kernel, including all feasible boundary points. -/
theorem baExponent_le_refined {w : ℝ} (hw : w ∈ Set.Icc ((13:ℝ)/125) (112/125)) :
    baExponent w ≤ refined.toFun w := by
  apply csSup_le
  · refine ⟨_, ⟨(0,0), ?_, rfl⟩⟩
    norm_num
    constructor <;> linarith [hw.1, hw.2]
  · rintro y ⟨⟨a,b⟩, ⟨ha,hb,hf⟩, rfl⟩
    exact majorant_pointwise ha hb hw hf

end Spin.Majorant
