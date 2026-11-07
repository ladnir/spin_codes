import SpinCodes.UniformMixer.ScalableDenseLeaf
import SpinCodes.UniformMixer.ScalableDenseCoverCore


/-! Generated numerical resource pilot; Lean checks every proposition. -/

namespace Spin.UniformMixer.ScalableDenseChunkPilot

open Spin.Numeric
open Spin.Structured.DenseOccupationFixed
open Spin.Structured.DenseGeometry

set_option maxRecDepth 8192
set_option maxHeartbeats 1000000

private def slope : QInput := ⟨833, 500⟩
private def intercept : QInput := ⟨-43311, 250000⟩

private def w11P : QInput := ⟨3418172182599063, 2500000000000000000⟩
private def w11Y : QInput := ⟨3449198973084403, 10000000000000000⟩
private def w11Z : QInput := ⟨249764041923714602568781432389, 250000000000000000000000000000⟩
private def w11Radius : QInput := ⟨38810780444344446253964500401, 40000000000000000000000000000⟩
private def w11Vector : Fin 3 → QInput := ![⟨1, 1⟩, ⟨261491572407949, 1000000000000000000⟩, ⟨1287880737392121, 5000000000000000000⟩]
private theorem w11Transfer :
    transferRatCheck
      (transferQInputRat w11P * transferQInputRat w11Y)
      (transferQInputRat w11Z) (transferQInputRat w11Radius)
      (fun i => transferQInputRat (w11Vector i)) = true := by
  decide +kernel

private def w12P : QInput := ⟨15648434568125509, 10000000000000000000⟩
private def w12Y : QInput := ⟨46350850090699891, 100000000000000000⟩
private def w12Z : QInput := ⟨499274187147198594219610214687, 500000000000000000000000000000⟩
private def w12Radius : QInput := ⟨954636400210748300021870908063, 1000000000000000000000000000000⟩
private def w12Vector : Fin 3 → QInput := ![⟨1, 1⟩, ⟨2626727268736291, 10000000000000000000⟩, ⟨160406931671781, 625000000000000000⟩]
private theorem w12Transfer :
    transferRatCheck
      (transferQInputRat w12P * transferQInputRat w12Y)
      (transferQInputRat w12Z) (transferQInputRat w12Radius)
      (fun i => transferQInputRat (w12Vector i)) = true := by
  decide +kernel

private def w121P : QInput := ⟨12238011055893281, 10000000000000000000⟩
private def w121Y : QInput := ⟨2183047231417087, 3125000000000000⟩
private def w121Z : QInput := ⟨998288993131874559843639944831, 1000000000000000000000000000000⟩
private def w121Radius : QInput := ⟨473373537047097221953545663561, 500000000000000000000000000000⟩
private def w121Vector : Fin 3 → QInput := ![⟨1, 1⟩, ⟨105309830869797, 400000000000000000⟩, ⟨512354494801741, 2000000000000000000⟩]
private theorem w121Transfer :
    transferRatCheck
      (transferQInputRat w121P * transferQInputRat w121Y)
      (transferQInputRat w121Z) (transferQInputRat w121Radius)
      (fun i => transferQInputRat (w121Vector i)) = true := by
  decide +kernel

private theorem vertex0 :
    scalableVertexCheck slope intercept w11P w11Y
      w11Radius w11Z ⟨202751, 327680000⟩ ⟨13, 125⟩ 50
      (-100000000000000000000000) = true := by
  decide +kernel

private theorem vertex1 :
    scalableVertexCheck slope intercept w11P w11Y
      w11Radius w11Z ⟨202751, 327680000⟩ ⟨55736026121, 480000000000⟩ 50
      (-100000000000000000000000) = true := by
  decide +kernel

private theorem vertex2 :
    scalableVertexCheck slope intercept w11P w11Y
      w11Radius w11Z ⟨34093, 40960000⟩ ⟨13, 125⟩ 50
      (-100000000000000000000000) = true := by
  decide +kernel

private theorem vertex3 :
    scalableVertexCheck slope intercept w11P w11Y
      w11Radius w11Z ⟨34093, 40960000⟩ ⟨30776026121, 240000000000⟩ 50
      (-100000000000000000000000) = true := by
  decide +kernel

private theorem vertex4 :
    scalableVertexCheck slope intercept w11P w11Y
      w11Radius w11Z ⟨58187, 81920000⟩ ⟨13, 125⟩ 50
      (-100000000000000000000000) = true := by
  decide +kernel

private theorem vertex5 :
    scalableVertexCheck slope intercept w11P w11Y
      w11Radius w11Z ⟨58187, 81920000⟩ ⟨30776026121, 240000000000⟩ 50
      (-100000000000000000000000) = true := by
  decide +kernel

private theorem vertex6 :
    scalableVertexCheck slope intercept w11P w11Y
      w11Radius w11Z ⟨58187, 81920000⟩ ⟨55736026121, 480000000000⟩ 50
      (-100000000000000000000000) = true := by
  decide +kernel

private theorem vertex7 :
    scalableVertexCheck slope intercept w11P w11Y
      w11Radius w11Z ⟨851, 1310720⟩ ⟨13, 125⟩ 50
      (-100000000000000000000000) = true := by
  decide +kernel

private theorem vertex8 :
    scalableVertexCheck slope intercept w11P w11Y
      w11Radius w11Z ⟨851, 1310720⟩ ⟨30776026121, 240000000000⟩ 50
      (-100000000000000000000000) = true := by
  decide +kernel

private theorem vertex9 :
    scalableVertexCheck slope intercept w11P w11Y
      w11Radius w11Z ⟨851, 1310720⟩ ⟨55736026121, 480000000000⟩ 50
      (-100000000000000000000000) = true := by
  decide +kernel

private theorem vertex10 :
    scalableVertexCheck slope intercept w12P w12Y
      w12Radius w12Z ⟨34093, 40960000⟩ ⟨18296026121, 120000000000⟩ 50
      (-100000000000000000000000) = true := by
  decide +kernel

private theorem vertex11 :
    scalableVertexCheck slope intercept w12P w12Y
      w12Radius w12Z ⟨34093, 40960000⟩ ⟨30776026121, 240000000000⟩ 50
      (-100000000000000000000000) = true := by
  decide +kernel

private theorem vertex12 :
    scalableVertexCheck slope intercept w12P w12Y
      w12Radius w12Z ⟨58187, 81920000⟩ ⟨18296026121, 120000000000⟩ 50
      (-100000000000000000000000) = true := by
  decide +kernel

private theorem vertex13 :
    scalableVertexCheck slope intercept w12P w12Y
      w12Radius w12Z ⟨58187, 81920000⟩ ⟨30776026121, 240000000000⟩ 50
      (-100000000000000000000000) = true := by
  decide +kernel

private theorem vertex14 :
    scalableVertexCheck slope intercept w121P w121Y
      w121Radius w121Z ⟨12047, 20480000⟩ ⟨13, 125⟩ 50
      (-100000000000000000000000) = true := by
  decide +kernel

private theorem vertex15 :
    scalableVertexCheck slope intercept w121P w121Y
      w121Radius w121Z ⟨12047, 20480000⟩ ⟨18296026121, 120000000000⟩ 50
      (-100000000000000000000000) = true := by
  decide +kernel

private theorem vertex16 :
    scalableVertexCheck slope intercept w121P w121Y
      w121Radius w121Z ⟨12047, 20480000⟩ ⟨30776026121, 240000000000⟩ 50
      (-100000000000000000000000) = true := by
  decide +kernel

private theorem vertex17 :
    scalableVertexCheck slope intercept w121P w121Y
      w121Radius w121Z ⟨12047, 20480000⟩ ⟨55736026121, 480000000000⟩ 50
      (-100000000000000000000000) = true := by
  decide +kernel

private theorem vertex18 :
    scalableVertexCheck slope intercept w121P w121Y
      w121Radius w121Z ⟨202751, 327680000⟩ ⟨13, 125⟩ 50
      (-100000000000000000000000) = true := by
  decide +kernel

private theorem vertex19 :
    scalableVertexCheck slope intercept w121P w121Y
      w121Radius w121Z ⟨202751, 327680000⟩ ⟨55736026121, 480000000000⟩ 50
      (-100000000000000000000000) = true := by
  decide +kernel

private theorem vertex20 :
    scalableVertexCheck slope intercept w121P w121Y
      w121Radius w121Z ⟨58187, 81920000⟩ ⟨18296026121, 120000000000⟩ 50
      (-100000000000000000000000) = true := by
  decide +kernel

private theorem vertex21 :
    scalableVertexCheck slope intercept w121P w121Y
      w121Radius w121Z ⟨58187, 81920000⟩ ⟨30776026121, 240000000000⟩ 50
      (-100000000000000000000000) = true := by
  decide +kernel

private theorem vertex22 :
    scalableVertexCheck slope intercept w121P w121Y
      w121Radius w121Z ⟨851, 1310720⟩ ⟨30776026121, 240000000000⟩ 50
      (-100000000000000000000000) = true := by
  decide +kernel

private theorem vertex23 :
    scalableVertexCheck slope intercept w121P w121Y
      w121Radius w121Z ⟨851, 1310720⟩ ⟨55736026121, 480000000000⟩ 50
      (-100000000000000000000000) = true := by
  decide +kernel

private def supportLine : ℚ × ℚ := ((833 / 500 : ℚ), (-43311 / 250000 : ℚ))
private theorem support_mem :
    supportLine ∈ Spin.Majorant.refined.supports := by
  decide +kernel

private abbrev target : ℝ → ℝ → Prop :=
  DenseTarget (1 / 10000000) 10000000

private def leaf0A0 : QInput := ⟨12047, 20480000⟩
private def leaf0A1 : QInput := ⟨202751, 327680000⟩
private def leaf0X0 : QInput := ⟨13, 125⟩
private def leaf0X1 : QInput := ⟨55736026121, 480000000000⟩
private theorem leaf0_rate {α x : ℝ}
    (hα : α ∈ Set.Icc leaf0A0.real leaf0A1.real)
    (hx : x ∈ Set.Icc leaf0X0.real leaf0X1.real) :
    ScalableDensePointRate α x slope.real intercept.real
      (1 / 10000000) 10000000 := by
  have h := scalable_dense_checked_leaf slope intercept
    w121P w121Y w121Radius w121Z
    leaf0A0 leaf0A1 leaf0X0 leaf0X1
    w121Vector (10000000 : ℝ)
    (by decide) (by decide) (by decide) (by decide)
    (by decide) (by decide)
    (by decide) (by decide) (by decide) (by decide)
    (by norm_num [w121P, QInput.real])
    (by norm_num [w121P, QInput.real])
    (by norm_num [w121Y, QInput.real])
    (by norm_num [w121Y, QInput.real])
    (by norm_num [w121Radius, QInput.real])
    (by norm_num [w121Z, QInput.real])
    (by norm_num [w121Z, QInput.real])
    (by norm_num [leaf0A0, QInput.real])
    (by norm_num [leaf0A1, QInput.real])
    (by norm_num [leaf0X0, QInput.real])
    (by norm_num [leaf0X1, QInput.real])
    (by norm_num : (0 : ℝ) ≤ 10000000)
    (by intro j; fin_cases j <;> norm_num [w121Vector, QInput.real])
    w121Transfer vertex14 vertex17 vertex18 vertex19 hα hx
  simpa [w121Vector, QInput.real] using h

private theorem leaf0_claim : DenseClaim target
    ⟨(12047 / 20480000 : ℚ), (202751 / 327680000 : ℚ), (13 / 125 : ℚ), (55736026121 / 480000000000 : ℚ)⟩ := by
  intro α x hr
  have hα : α ∈ Set.Icc leaf0A0.real leaf0A1.real := by
    norm_num [Rect.Contains, leaf0A0, leaf0A1, QInput.real] at hr ⊢
    exact ⟨hr.1, hr.2.1⟩
  have hx : x ∈ Set.Icc leaf0X0.real leaf0X1.real := by
    norm_num [Rect.Contains, leaf0X0, leaf0X1, QInput.real] at hr ⊢
    exact ⟨hr.2.2.1, hr.2.2.2⟩
  refine ⟨supportLine, support_mem, ?_⟩
  simpa [target, supportLine, slope, intercept, QInput.real] using
    (leaf0_rate hα hx)

private def leaf1A0 : QInput := ⟨202751, 327680000⟩
private def leaf1A1 : QInput := ⟨851, 1310720⟩
private def leaf1X0 : QInput := ⟨13, 125⟩
private def leaf1X1 : QInput := ⟨55736026121, 480000000000⟩
private theorem leaf1_rate {α x : ℝ}
    (hα : α ∈ Set.Icc leaf1A0.real leaf1A1.real)
    (hx : x ∈ Set.Icc leaf1X0.real leaf1X1.real) :
    ScalableDensePointRate α x slope.real intercept.real
      (1 / 10000000) 10000000 := by
  have h := scalable_dense_checked_leaf slope intercept
    w11P w11Y w11Radius w11Z
    leaf1A0 leaf1A1 leaf1X0 leaf1X1
    w11Vector (10000000 : ℝ)
    (by decide) (by decide) (by decide) (by decide)
    (by decide) (by decide)
    (by decide) (by decide) (by decide) (by decide)
    (by norm_num [w11P, QInput.real])
    (by norm_num [w11P, QInput.real])
    (by norm_num [w11Y, QInput.real])
    (by norm_num [w11Y, QInput.real])
    (by norm_num [w11Radius, QInput.real])
    (by norm_num [w11Z, QInput.real])
    (by norm_num [w11Z, QInput.real])
    (by norm_num [leaf1A0, QInput.real])
    (by norm_num [leaf1A1, QInput.real])
    (by norm_num [leaf1X0, QInput.real])
    (by norm_num [leaf1X1, QInput.real])
    (by norm_num : (0 : ℝ) ≤ 10000000)
    (by intro j; fin_cases j <;> norm_num [w11Vector, QInput.real])
    w11Transfer vertex0 vertex1 vertex7 vertex9 hα hx
  simpa [w11Vector, QInput.real] using h

private theorem leaf1_claim : DenseClaim target
    ⟨(202751 / 327680000 : ℚ), (851 / 1310720 : ℚ), (13 / 125 : ℚ), (55736026121 / 480000000000 : ℚ)⟩ := by
  intro α x hr
  have hα : α ∈ Set.Icc leaf1A0.real leaf1A1.real := by
    norm_num [Rect.Contains, leaf1A0, leaf1A1, QInput.real] at hr ⊢
    exact ⟨hr.1, hr.2.1⟩
  have hx : x ∈ Set.Icc leaf1X0.real leaf1X1.real := by
    norm_num [Rect.Contains, leaf1X0, leaf1X1, QInput.real] at hr ⊢
    exact ⟨hr.2.2.1, hr.2.2.2⟩
  refine ⟨supportLine, support_mem, ?_⟩
  simpa [target, supportLine, slope, intercept, QInput.real] using
    (leaf1_rate hα hx)

private def leaf2A0 : QInput := ⟨851, 1310720⟩
private def leaf2A1 : QInput := ⟨58187, 81920000⟩
private def leaf2X0 : QInput := ⟨13, 125⟩
private def leaf2X1 : QInput := ⟨55736026121, 480000000000⟩
private theorem leaf2_rate {α x : ℝ}
    (hα : α ∈ Set.Icc leaf2A0.real leaf2A1.real)
    (hx : x ∈ Set.Icc leaf2X0.real leaf2X1.real) :
    ScalableDensePointRate α x slope.real intercept.real
      (1 / 10000000) 10000000 := by
  have h := scalable_dense_checked_leaf slope intercept
    w11P w11Y w11Radius w11Z
    leaf2A0 leaf2A1 leaf2X0 leaf2X1
    w11Vector (10000000 : ℝ)
    (by decide) (by decide) (by decide) (by decide)
    (by decide) (by decide)
    (by decide) (by decide) (by decide) (by decide)
    (by norm_num [w11P, QInput.real])
    (by norm_num [w11P, QInput.real])
    (by norm_num [w11Y, QInput.real])
    (by norm_num [w11Y, QInput.real])
    (by norm_num [w11Radius, QInput.real])
    (by norm_num [w11Z, QInput.real])
    (by norm_num [w11Z, QInput.real])
    (by norm_num [leaf2A0, QInput.real])
    (by norm_num [leaf2A1, QInput.real])
    (by norm_num [leaf2X0, QInput.real])
    (by norm_num [leaf2X1, QInput.real])
    (by norm_num : (0 : ℝ) ≤ 10000000)
    (by intro j; fin_cases j <;> norm_num [w11Vector, QInput.real])
    w11Transfer vertex7 vertex9 vertex4 vertex6 hα hx
  simpa [w11Vector, QInput.real] using h

private theorem leaf2_claim : DenseClaim target
    ⟨(851 / 1310720 : ℚ), (58187 / 81920000 : ℚ), (13 / 125 : ℚ), (55736026121 / 480000000000 : ℚ)⟩ := by
  intro α x hr
  have hα : α ∈ Set.Icc leaf2A0.real leaf2A1.real := by
    norm_num [Rect.Contains, leaf2A0, leaf2A1, QInput.real] at hr ⊢
    exact ⟨hr.1, hr.2.1⟩
  have hx : x ∈ Set.Icc leaf2X0.real leaf2X1.real := by
    norm_num [Rect.Contains, leaf2X0, leaf2X1, QInput.real] at hr ⊢
    exact ⟨hr.2.2.1, hr.2.2.2⟩
  refine ⟨supportLine, support_mem, ?_⟩
  simpa [target, supportLine, slope, intercept, QInput.real] using
    (leaf2_rate hα hx)

private def leaf3A0 : QInput := ⟨12047, 20480000⟩
private def leaf3A1 : QInput := ⟨851, 1310720⟩
private def leaf3X0 : QInput := ⟨55736026121, 480000000000⟩
private def leaf3X1 : QInput := ⟨30776026121, 240000000000⟩
private theorem leaf3_rate {α x : ℝ}
    (hα : α ∈ Set.Icc leaf3A0.real leaf3A1.real)
    (hx : x ∈ Set.Icc leaf3X0.real leaf3X1.real) :
    ScalableDensePointRate α x slope.real intercept.real
      (1 / 10000000) 10000000 := by
  have h := scalable_dense_checked_leaf slope intercept
    w121P w121Y w121Radius w121Z
    leaf3A0 leaf3A1 leaf3X0 leaf3X1
    w121Vector (10000000 : ℝ)
    (by decide) (by decide) (by decide) (by decide)
    (by decide) (by decide)
    (by decide) (by decide) (by decide) (by decide)
    (by norm_num [w121P, QInput.real])
    (by norm_num [w121P, QInput.real])
    (by norm_num [w121Y, QInput.real])
    (by norm_num [w121Y, QInput.real])
    (by norm_num [w121Radius, QInput.real])
    (by norm_num [w121Z, QInput.real])
    (by norm_num [w121Z, QInput.real])
    (by norm_num [leaf3A0, QInput.real])
    (by norm_num [leaf3A1, QInput.real])
    (by norm_num [leaf3X0, QInput.real])
    (by norm_num [leaf3X1, QInput.real])
    (by norm_num : (0 : ℝ) ≤ 10000000)
    (by intro j; fin_cases j <;> norm_num [w121Vector, QInput.real])
    w121Transfer vertex17 vertex16 vertex23 vertex22 hα hx
  simpa [w121Vector, QInput.real] using h

private theorem leaf3_claim : DenseClaim target
    ⟨(12047 / 20480000 : ℚ), (851 / 1310720 : ℚ), (55736026121 / 480000000000 : ℚ), (30776026121 / 240000000000 : ℚ)⟩ := by
  intro α x hr
  have hα : α ∈ Set.Icc leaf3A0.real leaf3A1.real := by
    norm_num [Rect.Contains, leaf3A0, leaf3A1, QInput.real] at hr ⊢
    exact ⟨hr.1, hr.2.1⟩
  have hx : x ∈ Set.Icc leaf3X0.real leaf3X1.real := by
    norm_num [Rect.Contains, leaf3X0, leaf3X1, QInput.real] at hr ⊢
    exact ⟨hr.2.2.1, hr.2.2.2⟩
  refine ⟨supportLine, support_mem, ?_⟩
  simpa [target, supportLine, slope, intercept, QInput.real] using
    (leaf3_rate hα hx)

private def leaf4A0 : QInput := ⟨851, 1310720⟩
private def leaf4A1 : QInput := ⟨58187, 81920000⟩
private def leaf4X0 : QInput := ⟨55736026121, 480000000000⟩
private def leaf4X1 : QInput := ⟨30776026121, 240000000000⟩
private theorem leaf4_rate {α x : ℝ}
    (hα : α ∈ Set.Icc leaf4A0.real leaf4A1.real)
    (hx : x ∈ Set.Icc leaf4X0.real leaf4X1.real) :
    ScalableDensePointRate α x slope.real intercept.real
      (1 / 10000000) 10000000 := by
  have h := scalable_dense_checked_leaf slope intercept
    w11P w11Y w11Radius w11Z
    leaf4A0 leaf4A1 leaf4X0 leaf4X1
    w11Vector (10000000 : ℝ)
    (by decide) (by decide) (by decide) (by decide)
    (by decide) (by decide)
    (by decide) (by decide) (by decide) (by decide)
    (by norm_num [w11P, QInput.real])
    (by norm_num [w11P, QInput.real])
    (by norm_num [w11Y, QInput.real])
    (by norm_num [w11Y, QInput.real])
    (by norm_num [w11Radius, QInput.real])
    (by norm_num [w11Z, QInput.real])
    (by norm_num [w11Z, QInput.real])
    (by norm_num [leaf4A0, QInput.real])
    (by norm_num [leaf4A1, QInput.real])
    (by norm_num [leaf4X0, QInput.real])
    (by norm_num [leaf4X1, QInput.real])
    (by norm_num : (0 : ℝ) ≤ 10000000)
    (by intro j; fin_cases j <;> norm_num [w11Vector, QInput.real])
    w11Transfer vertex9 vertex8 vertex6 vertex5 hα hx
  simpa [w11Vector, QInput.real] using h

private theorem leaf4_claim : DenseClaim target
    ⟨(851 / 1310720 : ℚ), (58187 / 81920000 : ℚ), (55736026121 / 480000000000 : ℚ), (30776026121 / 240000000000 : ℚ)⟩ := by
  intro α x hr
  have hα : α ∈ Set.Icc leaf4A0.real leaf4A1.real := by
    norm_num [Rect.Contains, leaf4A0, leaf4A1, QInput.real] at hr ⊢
    exact ⟨hr.1, hr.2.1⟩
  have hx : x ∈ Set.Icc leaf4X0.real leaf4X1.real := by
    norm_num [Rect.Contains, leaf4X0, leaf4X1, QInput.real] at hr ⊢
    exact ⟨hr.2.2.1, hr.2.2.2⟩
  refine ⟨supportLine, support_mem, ?_⟩
  simpa [target, supportLine, slope, intercept, QInput.real] using
    (leaf4_rate hα hx)

private def leaf5A0 : QInput := ⟨58187, 81920000⟩
private def leaf5A1 : QInput := ⟨34093, 40960000⟩
private def leaf5X0 : QInput := ⟨13, 125⟩
private def leaf5X1 : QInput := ⟨30776026121, 240000000000⟩
private theorem leaf5_rate {α x : ℝ}
    (hα : α ∈ Set.Icc leaf5A0.real leaf5A1.real)
    (hx : x ∈ Set.Icc leaf5X0.real leaf5X1.real) :
    ScalableDensePointRate α x slope.real intercept.real
      (1 / 10000000) 10000000 := by
  have h := scalable_dense_checked_leaf slope intercept
    w11P w11Y w11Radius w11Z
    leaf5A0 leaf5A1 leaf5X0 leaf5X1
    w11Vector (10000000 : ℝ)
    (by decide) (by decide) (by decide) (by decide)
    (by decide) (by decide)
    (by decide) (by decide) (by decide) (by decide)
    (by norm_num [w11P, QInput.real])
    (by norm_num [w11P, QInput.real])
    (by norm_num [w11Y, QInput.real])
    (by norm_num [w11Y, QInput.real])
    (by norm_num [w11Radius, QInput.real])
    (by norm_num [w11Z, QInput.real])
    (by norm_num [w11Z, QInput.real])
    (by norm_num [leaf5A0, QInput.real])
    (by norm_num [leaf5A1, QInput.real])
    (by norm_num [leaf5X0, QInput.real])
    (by norm_num [leaf5X1, QInput.real])
    (by norm_num : (0 : ℝ) ≤ 10000000)
    (by intro j; fin_cases j <;> norm_num [w11Vector, QInput.real])
    w11Transfer vertex4 vertex5 vertex2 vertex3 hα hx
  simpa [w11Vector, QInput.real] using h

private theorem leaf5_claim : DenseClaim target
    ⟨(58187 / 81920000 : ℚ), (34093 / 40960000 : ℚ), (13 / 125 : ℚ), (30776026121 / 240000000000 : ℚ)⟩ := by
  intro α x hr
  have hα : α ∈ Set.Icc leaf5A0.real leaf5A1.real := by
    norm_num [Rect.Contains, leaf5A0, leaf5A1, QInput.real] at hr ⊢
    exact ⟨hr.1, hr.2.1⟩
  have hx : x ∈ Set.Icc leaf5X0.real leaf5X1.real := by
    norm_num [Rect.Contains, leaf5X0, leaf5X1, QInput.real] at hr ⊢
    exact ⟨hr.2.2.1, hr.2.2.2⟩
  refine ⟨supportLine, support_mem, ?_⟩
  simpa [target, supportLine, slope, intercept, QInput.real] using
    (leaf5_rate hα hx)

private def leaf6A0 : QInput := ⟨12047, 20480000⟩
private def leaf6A1 : QInput := ⟨58187, 81920000⟩
private def leaf6X0 : QInput := ⟨30776026121, 240000000000⟩
private def leaf6X1 : QInput := ⟨18296026121, 120000000000⟩
private theorem leaf6_rate {α x : ℝ}
    (hα : α ∈ Set.Icc leaf6A0.real leaf6A1.real)
    (hx : x ∈ Set.Icc leaf6X0.real leaf6X1.real) :
    ScalableDensePointRate α x slope.real intercept.real
      (1 / 10000000) 10000000 := by
  have h := scalable_dense_checked_leaf slope intercept
    w121P w121Y w121Radius w121Z
    leaf6A0 leaf6A1 leaf6X0 leaf6X1
    w121Vector (10000000 : ℝ)
    (by decide) (by decide) (by decide) (by decide)
    (by decide) (by decide)
    (by decide) (by decide) (by decide) (by decide)
    (by norm_num [w121P, QInput.real])
    (by norm_num [w121P, QInput.real])
    (by norm_num [w121Y, QInput.real])
    (by norm_num [w121Y, QInput.real])
    (by norm_num [w121Radius, QInput.real])
    (by norm_num [w121Z, QInput.real])
    (by norm_num [w121Z, QInput.real])
    (by norm_num [leaf6A0, QInput.real])
    (by norm_num [leaf6A1, QInput.real])
    (by norm_num [leaf6X0, QInput.real])
    (by norm_num [leaf6X1, QInput.real])
    (by norm_num : (0 : ℝ) ≤ 10000000)
    (by intro j; fin_cases j <;> norm_num [w121Vector, QInput.real])
    w121Transfer vertex16 vertex15 vertex21 vertex20 hα hx
  simpa [w121Vector, QInput.real] using h

private theorem leaf6_claim : DenseClaim target
    ⟨(12047 / 20480000 : ℚ), (58187 / 81920000 : ℚ), (30776026121 / 240000000000 : ℚ), (18296026121 / 120000000000 : ℚ)⟩ := by
  intro α x hr
  have hα : α ∈ Set.Icc leaf6A0.real leaf6A1.real := by
    norm_num [Rect.Contains, leaf6A0, leaf6A1, QInput.real] at hr ⊢
    exact ⟨hr.1, hr.2.1⟩
  have hx : x ∈ Set.Icc leaf6X0.real leaf6X1.real := by
    norm_num [Rect.Contains, leaf6X0, leaf6X1, QInput.real] at hr ⊢
    exact ⟨hr.2.2.1, hr.2.2.2⟩
  refine ⟨supportLine, support_mem, ?_⟩
  simpa [target, supportLine, slope, intercept, QInput.real] using
    (leaf6_rate hα hx)

private def leaf7A0 : QInput := ⟨58187, 81920000⟩
private def leaf7A1 : QInput := ⟨34093, 40960000⟩
private def leaf7X0 : QInput := ⟨30776026121, 240000000000⟩
private def leaf7X1 : QInput := ⟨18296026121, 120000000000⟩
private theorem leaf7_rate {α x : ℝ}
    (hα : α ∈ Set.Icc leaf7A0.real leaf7A1.real)
    (hx : x ∈ Set.Icc leaf7X0.real leaf7X1.real) :
    ScalableDensePointRate α x slope.real intercept.real
      (1 / 10000000) 10000000 := by
  have h := scalable_dense_checked_leaf slope intercept
    w12P w12Y w12Radius w12Z
    leaf7A0 leaf7A1 leaf7X0 leaf7X1
    w12Vector (10000000 : ℝ)
    (by decide) (by decide) (by decide) (by decide)
    (by decide) (by decide)
    (by decide) (by decide) (by decide) (by decide)
    (by norm_num [w12P, QInput.real])
    (by norm_num [w12P, QInput.real])
    (by norm_num [w12Y, QInput.real])
    (by norm_num [w12Y, QInput.real])
    (by norm_num [w12Radius, QInput.real])
    (by norm_num [w12Z, QInput.real])
    (by norm_num [w12Z, QInput.real])
    (by norm_num [leaf7A0, QInput.real])
    (by norm_num [leaf7A1, QInput.real])
    (by norm_num [leaf7X0, QInput.real])
    (by norm_num [leaf7X1, QInput.real])
    (by norm_num : (0 : ℝ) ≤ 10000000)
    (by intro j; fin_cases j <;> norm_num [w12Vector, QInput.real])
    w12Transfer vertex13 vertex12 vertex11 vertex10 hα hx
  simpa [w12Vector, QInput.real] using h

private theorem leaf7_claim : DenseClaim target
    ⟨(58187 / 81920000 : ℚ), (34093 / 40960000 : ℚ), (30776026121 / 240000000000 : ℚ), (18296026121 / 120000000000 : ℚ)⟩ := by
  intro α x hr
  have hα : α ∈ Set.Icc leaf7A0.real leaf7A1.real := by
    norm_num [Rect.Contains, leaf7A0, leaf7A1, QInput.real] at hr ⊢
    exact ⟨hr.1, hr.2.1⟩
  have hx : x ∈ Set.Icc leaf7X0.real leaf7X1.real := by
    norm_num [Rect.Contains, leaf7X0, leaf7X1, QInput.real] at hr ⊢
    exact ⟨hr.2.2.1, hr.2.2.2⟩
  refine ⟨supportLine, support_mem, ?_⟩
  simpa [target, supportLine, slope, intercept, QInput.real] using
    (leaf7_rate hα hx)

def root : Rect := ⟨(12047 / 20480000 : ℚ), (34093 / 40960000 : ℚ), (13 / 125 : ℚ), (18296026121 / 120000000000 : ℚ)⟩

theorem root_claim : DenseClaim target root := by
  exact (denseClaim_density target ⟨(12047 / 20480000 : ℚ), (34093 / 40960000 : ℚ), (13 / 125 : ℚ), (18296026121 / 120000000000 : ℚ)⟩ (30776026121 / 240000000000 : ℚ) ((denseClaim_alpha target ⟨(12047 / 20480000 : ℚ), (34093 / 40960000 : ℚ), (13 / 125 : ℚ), (30776026121 / 240000000000 : ℚ)⟩ (58187 / 81920000 : ℚ) ((denseClaim_density target ⟨(12047 / 20480000 : ℚ), (58187 / 81920000 : ℚ), (13 / 125 : ℚ), (30776026121 / 240000000000 : ℚ)⟩ (55736026121 / 480000000000 : ℚ) ((denseClaim_alpha target ⟨(12047 / 20480000 : ℚ), (58187 / 81920000 : ℚ), (13 / 125 : ℚ), (55736026121 / 480000000000 : ℚ)⟩ (851 / 1310720 : ℚ) ((denseClaim_alpha target ⟨(12047 / 20480000 : ℚ), (851 / 1310720 : ℚ), (13 / 125 : ℚ), (55736026121 / 480000000000 : ℚ)⟩ (202751 / 327680000 : ℚ) (leaf0_claim) (leaf1_claim))) (leaf2_claim))) ((denseClaim_alpha target ⟨(12047 / 20480000 : ℚ), (58187 / 81920000 : ℚ), (55736026121 / 480000000000 : ℚ), (30776026121 / 240000000000 : ℚ)⟩ (851 / 1310720 : ℚ) (leaf3_claim) (leaf4_claim))))) (leaf5_claim))) ((denseClaim_alpha target ⟨(12047 / 20480000 : ℚ), (34093 / 40960000 : ℚ), (30776026121 / 240000000000 : ℚ), (18296026121 / 120000000000 : ℚ)⟩ (58187 / 81920000 : ℚ) (leaf6_claim) (leaf7_claim))))

theorem root_claim_global :
    DenseClaim (DenseTarget (1 / 10000000) 10000000) root := by
  exact root_claim

#print axioms root_claim

#print axioms w11Transfer
#print axioms vertex0
#print axioms vertex23

end Spin.UniformMixer.ScalableDenseChunkPilot
