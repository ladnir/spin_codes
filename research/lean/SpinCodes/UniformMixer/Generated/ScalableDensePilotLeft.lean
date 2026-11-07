import SpinCodes.UniformMixer.ScalableDenseLeaf
import SpinCodes.UniformMixer.ScalableDenseCoverCore


/-! Generated numerical resource pilot; Lean checks every proposition. -/

namespace Spin.UniformMixer.ScalableDensePilotLeft

open Spin.Numeric
open Spin.Structured.DenseOccupationFixed
open Spin.Structured.DenseGeometry

set_option maxRecDepth 8192
set_option maxHeartbeats 1000000

private def slope : QInput := ⟨9, 100⟩
private def intercept : QInput := ⟨102967997603379096932443767333991244233, 340282366920938463463374607431768211456⟩

private def w74P : QInput := ⟨2475000089226973, 2500000000000000⟩
private def w74Y : QInput := ⟨12125046787559153, 25000000000000000⟩
private def w74Z : QInput := ⟨123595556909896105236687524711, 1000000000000000000000000000000⟩
private def w74Radius : QInput := ⟨46999411517631, 500000000000000000000000000000⟩
private def w74Vector : Fin 3 → QInput := ![⟨1, 1⟩, ⟨17650026616931, 125000000000000⟩, ⟨675507456001233, 12500000000000000⟩]
private theorem w74Transfer :
    transferRatCheck
      (transferQInputRat w74P * transferQInputRat w74Y)
      (transferQInputRat w74Z) (transferQInputRat w74Radius)
      (fun i => transferQInputRat (w74Vector i)) = true := by
  decide +kernel

private def w82P : QInput := ⟨99900000169545189, 100000000000000000⟩
private def w82Y : QInput := ⟨48500091705249537, 100000000000000000⟩
private def w82Z : QInput := ⟨2471910491445593750277965467, 20000000000000000000000000000⟩
private def w82Radius : QInput := ⟨93998509190313, 1000000000000000000000000000000⟩
private def w82Vector : Fin 3 → QInput := ![⟨1, 1⟩, ⟨1081874976172229, 5000000000000000⟩, ⟨1017555032111227, 10000000000000000⟩]
private theorem w82Transfer :
    transferRatCheck
      (transferQInputRat w82P * transferQInputRat w82Y)
      (transferQInputRat w82Z) (transferQInputRat w82Radius)
      (fun i => transferQInputRat (w82Vector i)) = true := by
  decide +kernel

private def w111P : QInput := ⟨49863101813168781, 50000000000000000⟩
private def w111Y : QInput := ⟨50426842715046771, 100000000000000000⟩
private def w111Z : QInput := ⟨1544943816139595642593979887, 12500000000000000000000000000⟩
private def w111Radius : QInput := ⟨93998294429031, 1000000000000000000000000000000⟩
private def w111Vector : Fin 3 → QInput := ![⟨1, 1⟩, ⟨667611956984541, 500000000000000⟩, ⟨1001433752214203, 1000000000000000⟩]
private theorem w111Transfer :
    transferRatCheck
      (transferQInputRat w111P * transferQInputRat w111Y)
      (transferQInputRat w111Z) (transferQInputRat w111Radius)
      (fun i => transferQInputRat (w111Vector i)) = true := by
  decide +kernel

private theorem vertex0 :
    scalableVertexCheck slope intercept w74P w74Y
      w74Radius w74Z ⟨1270001, 1280000⟩ ⟨316186482737629508933514007290843083540271318859, 664613997892457936451903530140172288000000000000⟩ 50
      (-100000000000000000000000) = true := by
  decide +kernel

private theorem vertex1 :
    scalableVertexCheck slope intercept w74P w74Y
      w74Radius w74Z ⟨1270001, 1280000⟩ ⟨49374602183332977403672480068036683319333956789, 103845937170696552570609926584401920000000000000⟩ 50
      (-100000000000000000000000) = true := by
  decide +kernel

private theorem vertex2 :
    scalableVertexCheck slope intercept w74P w74Y
      w74Radius w74Z ⟨1270001, 1280000⟩ ⟨790938778754819906208810355365628484592013285671, 1661534994731144841129758825350430720000000000000⟩ 50
      (-100000000000000000000000) = true := by
  decide +kernel

private theorem vertex3 :
    scalableVertexCheck slope intercept w74P w74Y
      w74Radius w74Z ⟨5090003, 5120000⟩ ⟨316186482737629508933514007290843083540271318859, 664613997892457936451903530140172288000000000000⟩ 50
      (-100000000000000000000000) = true := by
  decide +kernel

private theorem vertex4 :
    scalableVertexCheck slope intercept w74P w74Y
      w74Radius w74Z ⟨5090003, 5120000⟩ ⟨49374602183332977403672480068036683319333956789, 103845937170696552570609926584401920000000000000⟩ 50
      (-100000000000000000000000) = true := by
  decide +kernel

private theorem vertex5 :
    scalableVertexCheck slope intercept w74P w74Y
      w74Radius w74Z ⟨5090003, 5120000⟩ ⟨790938778754819906208810355365628484592013285671, 1661534994731144841129758825350430720000000000000⟩ 50
      (-100000000000000000000000) = true := by
  decide +kernel

private theorem vertex6 :
    scalableVertexCheck slope intercept w82P w82Y
      w82Radius w82Z ⟨1, 1⟩ ⟨316186482737629508933514007290843083540271318859, 664613997892457936451903530140172288000000000000⟩ 50
      (-100000000000000000000000) = true := by
  decide +kernel

private theorem vertex7 :
    scalableVertexCheck slope intercept w82P w82Y
      w82Radius w82Z ⟨1, 1⟩ ⟨49374602183332977403672480068036683319333956789, 103845937170696552570609926584401920000000000000⟩ 50
      (-100000000000000000000000) = true := by
  decide +kernel

private theorem vertex8 :
    scalableVertexCheck slope intercept w82P w82Y
      w82Radius w82Z ⟨1, 1⟩ ⟨790938778754819906208810355365628484592013285671, 1661534994731144841129758825350430720000000000000⟩ 50
      (-100000000000000000000000) = true := by
  decide +kernel

private theorem vertex9 :
    scalableVertexCheck slope intercept w82P w82Y
      w82Radius w82Z ⟨5110001, 5120000⟩ ⟨316186482737629508933514007290843083540271318859, 664613997892457936451903530140172288000000000000⟩ 50
      (-100000000000000000000000) = true := by
  decide +kernel

private theorem vertex10 :
    scalableVertexCheck slope intercept w82P w82Y
      w82Radius w82Z ⟨5110001, 5120000⟩ ⟨49374602183332977403672480068036683319333956789, 103845937170696552570609926584401920000000000000⟩ 50
      (-100000000000000000000000) = true := by
  decide +kernel

private theorem vertex11 :
    scalableVertexCheck slope intercept w82P w82Y
      w82Radius w82Z ⟨5110001, 5120000⟩ ⟨790938778754819906208810355365628484592013285671, 1661534994731144841129758825350430720000000000000⟩ 50
      (-100000000000000000000000) = true := by
  decide +kernel

private theorem vertex12 :
    scalableVertexCheck slope intercept w111P w111Y
      w111Radius w111Z ⟨2550001, 2560000⟩ ⟨316186482737629508933514007290843083540271318859, 664613997892457936451903530140172288000000000000⟩ 50
      (-100000000000000000000000) = true := by
  decide +kernel

private theorem vertex13 :
    scalableVertexCheck slope intercept w111P w111Y
      w111Radius w111Z ⟨2550001, 2560000⟩ ⟨49374602183332977403672480068036683319333956789, 103845937170696552570609926584401920000000000000⟩ 50
      (-100000000000000000000000) = true := by
  decide +kernel

private theorem vertex14 :
    scalableVertexCheck slope intercept w111P w111Y
      w111Radius w111Z ⟨2550001, 2560000⟩ ⟨790938778754819906208810355365628484592013285671, 1661534994731144841129758825350430720000000000000⟩ 50
      (-100000000000000000000000) = true := by
  decide +kernel

private theorem vertex15 :
    scalableVertexCheck slope intercept w111P w111Y
      w111Radius w111Z ⟨5090003, 5120000⟩ ⟨316186482737629508933514007290843083540271318859, 664613997892457936451903530140172288000000000000⟩ 50
      (-100000000000000000000000) = true := by
  decide +kernel

private theorem vertex16 :
    scalableVertexCheck slope intercept w111P w111Y
      w111Radius w111Z ⟨5090003, 5120000⟩ ⟨49374602183332977403672480068036683319333956789, 103845937170696552570609926584401920000000000000⟩ 50
      (-100000000000000000000000) = true := by
  decide +kernel

private theorem vertex17 :
    scalableVertexCheck slope intercept w111P w111Y
      w111Radius w111Z ⟨5090003, 5120000⟩ ⟨790938778754819906208810355365628484592013285671, 1661534994731144841129758825350430720000000000000⟩ 50
      (-100000000000000000000000) = true := by
  decide +kernel

private theorem vertex18 :
    scalableVertexCheck slope intercept w111P w111Y
      w111Radius w111Z ⟨5110001, 5120000⟩ ⟨316186482737629508933514007290843083540271318859, 664613997892457936451903530140172288000000000000⟩ 50
      (-100000000000000000000000) = true := by
  decide +kernel

private theorem vertex19 :
    scalableVertexCheck slope intercept w111P w111Y
      w111Radius w111Z ⟨5110001, 5120000⟩ ⟨49374602183332977403672480068036683319333956789, 103845937170696552570609926584401920000000000000⟩ 50
      (-100000000000000000000000) = true := by
  decide +kernel

private theorem vertex20 :
    scalableVertexCheck slope intercept w111P w111Y
      w111Radius w111Z ⟨5110001, 5120000⟩ ⟨790938778754819906208810355365628484592013285671, 1661534994731144841129758825350430720000000000000⟩ 50
      (-100000000000000000000000) = true := by
  decide +kernel

private def supportLine : ℚ × ℚ := ((9 / 100 : ℚ), (102967997603379096932443767333991244233 / 340282366920938463463374607431768211456 : ℚ))
private theorem support_mem :
    supportLine ∈ Spin.Majorant.refined.supports := by
  decide +kernel

private abbrev target : ℝ → ℝ → Prop :=
  DenseTarget (1 / 10000000) 10000000

private def leaf0A0 : QInput := ⟨1270001, 1280000⟩
private def leaf0A1 : QInput := ⟨5090003, 5120000⟩
private def leaf0X0 : QInput := ⟨49374602183332977403672480068036683319333956789, 103845937170696552570609926584401920000000000000⟩
private def leaf0X1 : QInput := ⟨316186482737629508933514007290843083540271318859, 664613997892457936451903530140172288000000000000⟩
private theorem leaf0_rate {α x : ℝ}
    (hα : α ∈ Set.Icc leaf0A0.real leaf0A1.real)
    (hx : x ∈ Set.Icc leaf0X0.real leaf0X1.real) :
    ScalableDensePointRate α x slope.real intercept.real
      (1 / 10000000) 10000000 := by
  have h := scalable_dense_checked_leaf slope intercept
    w74P w74Y w74Radius w74Z
    leaf0A0 leaf0A1 leaf0X0 leaf0X1
    w74Vector (10000000 : ℝ)
    (by decide) (by decide) (by decide) (by decide)
    (by decide) (by decide)
    (by decide) (by decide) (by decide) (by decide)
    (by norm_num [w74P, QInput.real])
    (by norm_num [w74P, QInput.real])
    (by norm_num [w74Y, QInput.real])
    (by norm_num [w74Y, QInput.real])
    (by norm_num [w74Radius, QInput.real])
    (by norm_num [w74Z, QInput.real])
    (by norm_num [w74Z, QInput.real])
    (by norm_num [leaf0A0, QInput.real])
    (by norm_num [leaf0A1, QInput.real])
    (by norm_num [leaf0X0, QInput.real])
    (by norm_num [leaf0X1, QInput.real])
    (by norm_num : (0 : ℝ) ≤ 10000000)
    (by intro j; fin_cases j <;> norm_num [w74Vector, QInput.real])
    w74Transfer vertex1 vertex0 vertex4 vertex3 hα hx
  simpa [w74Vector, QInput.real] using h

private theorem leaf0_claim : DenseClaim target
    ⟨(1270001 / 1280000 : ℚ), (5090003 / 5120000 : ℚ), (49374602183332977403672480068036683319333956789 / 103845937170696552570609926584401920000000000000 : ℚ), (316186482737629508933514007290843083540271318859 / 664613997892457936451903530140172288000000000000 : ℚ)⟩ := by
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

private def leaf1A0 : QInput := ⟨5090003, 5120000⟩
private def leaf1A1 : QInput := ⟨2550001, 2560000⟩
private def leaf1X0 : QInput := ⟨49374602183332977403672480068036683319333956789, 103845937170696552570609926584401920000000000000⟩
private def leaf1X1 : QInput := ⟨316186482737629508933514007290843083540271318859, 664613997892457936451903530140172288000000000000⟩
private theorem leaf1_rate {α x : ℝ}
    (hα : α ∈ Set.Icc leaf1A0.real leaf1A1.real)
    (hx : x ∈ Set.Icc leaf1X0.real leaf1X1.real) :
    ScalableDensePointRate α x slope.real intercept.real
      (1 / 10000000) 10000000 := by
  have h := scalable_dense_checked_leaf slope intercept
    w111P w111Y w111Radius w111Z
    leaf1A0 leaf1A1 leaf1X0 leaf1X1
    w111Vector (10000000 : ℝ)
    (by decide) (by decide) (by decide) (by decide)
    (by decide) (by decide)
    (by decide) (by decide) (by decide) (by decide)
    (by norm_num [w111P, QInput.real])
    (by norm_num [w111P, QInput.real])
    (by norm_num [w111Y, QInput.real])
    (by norm_num [w111Y, QInput.real])
    (by norm_num [w111Radius, QInput.real])
    (by norm_num [w111Z, QInput.real])
    (by norm_num [w111Z, QInput.real])
    (by norm_num [leaf1A0, QInput.real])
    (by norm_num [leaf1A1, QInput.real])
    (by norm_num [leaf1X0, QInput.real])
    (by norm_num [leaf1X1, QInput.real])
    (by norm_num : (0 : ℝ) ≤ 10000000)
    (by intro j; fin_cases j <;> norm_num [w111Vector, QInput.real])
    w111Transfer vertex16 vertex15 vertex13 vertex12 hα hx
  simpa [w111Vector, QInput.real] using h

private theorem leaf1_claim : DenseClaim target
    ⟨(5090003 / 5120000 : ℚ), (2550001 / 2560000 : ℚ), (49374602183332977403672480068036683319333956789 / 103845937170696552570609926584401920000000000000 : ℚ), (316186482737629508933514007290843083540271318859 / 664613997892457936451903530140172288000000000000 : ℚ)⟩ := by
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

private def leaf2A0 : QInput := ⟨1270001, 1280000⟩
private def leaf2A1 : QInput := ⟨5090003, 5120000⟩
private def leaf2X0 : QInput := ⟨316186482737629508933514007290843083540271318859, 664613997892457936451903530140172288000000000000⟩
private def leaf2X1 : QInput := ⟨790938778754819906208810355365628484592013285671, 1661534994731144841129758825350430720000000000000⟩
private theorem leaf2_rate {α x : ℝ}
    (hα : α ∈ Set.Icc leaf2A0.real leaf2A1.real)
    (hx : x ∈ Set.Icc leaf2X0.real leaf2X1.real) :
    ScalableDensePointRate α x slope.real intercept.real
      (1 / 10000000) 10000000 := by
  have h := scalable_dense_checked_leaf slope intercept
    w74P w74Y w74Radius w74Z
    leaf2A0 leaf2A1 leaf2X0 leaf2X1
    w74Vector (10000000 : ℝ)
    (by decide) (by decide) (by decide) (by decide)
    (by decide) (by decide)
    (by decide) (by decide) (by decide) (by decide)
    (by norm_num [w74P, QInput.real])
    (by norm_num [w74P, QInput.real])
    (by norm_num [w74Y, QInput.real])
    (by norm_num [w74Y, QInput.real])
    (by norm_num [w74Radius, QInput.real])
    (by norm_num [w74Z, QInput.real])
    (by norm_num [w74Z, QInput.real])
    (by norm_num [leaf2A0, QInput.real])
    (by norm_num [leaf2A1, QInput.real])
    (by norm_num [leaf2X0, QInput.real])
    (by norm_num [leaf2X1, QInput.real])
    (by norm_num : (0 : ℝ) ≤ 10000000)
    (by intro j; fin_cases j <;> norm_num [w74Vector, QInput.real])
    w74Transfer vertex0 vertex2 vertex3 vertex5 hα hx
  simpa [w74Vector, QInput.real] using h

private theorem leaf2_claim : DenseClaim target
    ⟨(1270001 / 1280000 : ℚ), (5090003 / 5120000 : ℚ), (316186482737629508933514007290843083540271318859 / 664613997892457936451903530140172288000000000000 : ℚ), (790938778754819906208810355365628484592013285671 / 1661534994731144841129758825350430720000000000000 : ℚ)⟩ := by
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

private def leaf3A0 : QInput := ⟨5090003, 5120000⟩
private def leaf3A1 : QInput := ⟨2550001, 2560000⟩
private def leaf3X0 : QInput := ⟨316186482737629508933514007290843083540271318859, 664613997892457936451903530140172288000000000000⟩
private def leaf3X1 : QInput := ⟨790938778754819906208810355365628484592013285671, 1661534994731144841129758825350430720000000000000⟩
private theorem leaf3_rate {α x : ℝ}
    (hα : α ∈ Set.Icc leaf3A0.real leaf3A1.real)
    (hx : x ∈ Set.Icc leaf3X0.real leaf3X1.real) :
    ScalableDensePointRate α x slope.real intercept.real
      (1 / 10000000) 10000000 := by
  have h := scalable_dense_checked_leaf slope intercept
    w111P w111Y w111Radius w111Z
    leaf3A0 leaf3A1 leaf3X0 leaf3X1
    w111Vector (10000000 : ℝ)
    (by decide) (by decide) (by decide) (by decide)
    (by decide) (by decide)
    (by decide) (by decide) (by decide) (by decide)
    (by norm_num [w111P, QInput.real])
    (by norm_num [w111P, QInput.real])
    (by norm_num [w111Y, QInput.real])
    (by norm_num [w111Y, QInput.real])
    (by norm_num [w111Radius, QInput.real])
    (by norm_num [w111Z, QInput.real])
    (by norm_num [w111Z, QInput.real])
    (by norm_num [leaf3A0, QInput.real])
    (by norm_num [leaf3A1, QInput.real])
    (by norm_num [leaf3X0, QInput.real])
    (by norm_num [leaf3X1, QInput.real])
    (by norm_num : (0 : ℝ) ≤ 10000000)
    (by intro j; fin_cases j <;> norm_num [w111Vector, QInput.real])
    w111Transfer vertex15 vertex17 vertex12 vertex14 hα hx
  simpa [w111Vector, QInput.real] using h

private theorem leaf3_claim : DenseClaim target
    ⟨(5090003 / 5120000 : ℚ), (2550001 / 2560000 : ℚ), (316186482737629508933514007290843083540271318859 / 664613997892457936451903530140172288000000000000 : ℚ), (790938778754819906208810355365628484592013285671 / 1661534994731144841129758825350430720000000000000 : ℚ)⟩ := by
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

private def leaf4A0 : QInput := ⟨2550001, 2560000⟩
private def leaf4A1 : QInput := ⟨5110001, 5120000⟩
private def leaf4X0 : QInput := ⟨49374602183332977403672480068036683319333956789, 103845937170696552570609926584401920000000000000⟩
private def leaf4X1 : QInput := ⟨316186482737629508933514007290843083540271318859, 664613997892457936451903530140172288000000000000⟩
private theorem leaf4_rate {α x : ℝ}
    (hα : α ∈ Set.Icc leaf4A0.real leaf4A1.real)
    (hx : x ∈ Set.Icc leaf4X0.real leaf4X1.real) :
    ScalableDensePointRate α x slope.real intercept.real
      (1 / 10000000) 10000000 := by
  have h := scalable_dense_checked_leaf slope intercept
    w111P w111Y w111Radius w111Z
    leaf4A0 leaf4A1 leaf4X0 leaf4X1
    w111Vector (10000000 : ℝ)
    (by decide) (by decide) (by decide) (by decide)
    (by decide) (by decide)
    (by decide) (by decide) (by decide) (by decide)
    (by norm_num [w111P, QInput.real])
    (by norm_num [w111P, QInput.real])
    (by norm_num [w111Y, QInput.real])
    (by norm_num [w111Y, QInput.real])
    (by norm_num [w111Radius, QInput.real])
    (by norm_num [w111Z, QInput.real])
    (by norm_num [w111Z, QInput.real])
    (by norm_num [leaf4A0, QInput.real])
    (by norm_num [leaf4A1, QInput.real])
    (by norm_num [leaf4X0, QInput.real])
    (by norm_num [leaf4X1, QInput.real])
    (by norm_num : (0 : ℝ) ≤ 10000000)
    (by intro j; fin_cases j <;> norm_num [w111Vector, QInput.real])
    w111Transfer vertex13 vertex12 vertex19 vertex18 hα hx
  simpa [w111Vector, QInput.real] using h

private theorem leaf4_claim : DenseClaim target
    ⟨(2550001 / 2560000 : ℚ), (5110001 / 5120000 : ℚ), (49374602183332977403672480068036683319333956789 / 103845937170696552570609926584401920000000000000 : ℚ), (316186482737629508933514007290843083540271318859 / 664613997892457936451903530140172288000000000000 : ℚ)⟩ := by
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

private def leaf5A0 : QInput := ⟨5110001, 5120000⟩
private def leaf5A1 : QInput := ⟨1, 1⟩
private def leaf5X0 : QInput := ⟨49374602183332977403672480068036683319333956789, 103845937170696552570609926584401920000000000000⟩
private def leaf5X1 : QInput := ⟨316186482737629508933514007290843083540271318859, 664613997892457936451903530140172288000000000000⟩
private theorem leaf5_rate {α x : ℝ}
    (hα : α ∈ Set.Icc leaf5A0.real leaf5A1.real)
    (hx : x ∈ Set.Icc leaf5X0.real leaf5X1.real) :
    ScalableDensePointRate α x slope.real intercept.real
      (1 / 10000000) 10000000 := by
  have h := scalable_dense_checked_leaf slope intercept
    w82P w82Y w82Radius w82Z
    leaf5A0 leaf5A1 leaf5X0 leaf5X1
    w82Vector (10000000 : ℝ)
    (by decide) (by decide) (by decide) (by decide)
    (by decide) (by decide)
    (by decide) (by decide) (by decide) (by decide)
    (by norm_num [w82P, QInput.real])
    (by norm_num [w82P, QInput.real])
    (by norm_num [w82Y, QInput.real])
    (by norm_num [w82Y, QInput.real])
    (by norm_num [w82Radius, QInput.real])
    (by norm_num [w82Z, QInput.real])
    (by norm_num [w82Z, QInput.real])
    (by norm_num [leaf5A0, QInput.real])
    (by norm_num [leaf5A1, QInput.real])
    (by norm_num [leaf5X0, QInput.real])
    (by norm_num [leaf5X1, QInput.real])
    (by norm_num : (0 : ℝ) ≤ 10000000)
    (by intro j; fin_cases j <;> norm_num [w82Vector, QInput.real])
    w82Transfer vertex10 vertex9 vertex7 vertex6 hα hx
  simpa [w82Vector, QInput.real] using h

private theorem leaf5_claim : DenseClaim target
    ⟨(5110001 / 5120000 : ℚ), (1 / 1 : ℚ), (49374602183332977403672480068036683319333956789 / 103845937170696552570609926584401920000000000000 : ℚ), (316186482737629508933514007290843083540271318859 / 664613997892457936451903530140172288000000000000 : ℚ)⟩ := by
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

private def leaf6A0 : QInput := ⟨2550001, 2560000⟩
private def leaf6A1 : QInput := ⟨5110001, 5120000⟩
private def leaf6X0 : QInput := ⟨316186482737629508933514007290843083540271318859, 664613997892457936451903530140172288000000000000⟩
private def leaf6X1 : QInput := ⟨790938778754819906208810355365628484592013285671, 1661534994731144841129758825350430720000000000000⟩
private theorem leaf6_rate {α x : ℝ}
    (hα : α ∈ Set.Icc leaf6A0.real leaf6A1.real)
    (hx : x ∈ Set.Icc leaf6X0.real leaf6X1.real) :
    ScalableDensePointRate α x slope.real intercept.real
      (1 / 10000000) 10000000 := by
  have h := scalable_dense_checked_leaf slope intercept
    w111P w111Y w111Radius w111Z
    leaf6A0 leaf6A1 leaf6X0 leaf6X1
    w111Vector (10000000 : ℝ)
    (by decide) (by decide) (by decide) (by decide)
    (by decide) (by decide)
    (by decide) (by decide) (by decide) (by decide)
    (by norm_num [w111P, QInput.real])
    (by norm_num [w111P, QInput.real])
    (by norm_num [w111Y, QInput.real])
    (by norm_num [w111Y, QInput.real])
    (by norm_num [w111Radius, QInput.real])
    (by norm_num [w111Z, QInput.real])
    (by norm_num [w111Z, QInput.real])
    (by norm_num [leaf6A0, QInput.real])
    (by norm_num [leaf6A1, QInput.real])
    (by norm_num [leaf6X0, QInput.real])
    (by norm_num [leaf6X1, QInput.real])
    (by norm_num : (0 : ℝ) ≤ 10000000)
    (by intro j; fin_cases j <;> norm_num [w111Vector, QInput.real])
    w111Transfer vertex12 vertex14 vertex18 vertex20 hα hx
  simpa [w111Vector, QInput.real] using h

private theorem leaf6_claim : DenseClaim target
    ⟨(2550001 / 2560000 : ℚ), (5110001 / 5120000 : ℚ), (316186482737629508933514007290843083540271318859 / 664613997892457936451903530140172288000000000000 : ℚ), (790938778754819906208810355365628484592013285671 / 1661534994731144841129758825350430720000000000000 : ℚ)⟩ := by
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

private def leaf7A0 : QInput := ⟨5110001, 5120000⟩
private def leaf7A1 : QInput := ⟨1, 1⟩
private def leaf7X0 : QInput := ⟨316186482737629508933514007290843083540271318859, 664613997892457936451903530140172288000000000000⟩
private def leaf7X1 : QInput := ⟨790938778754819906208810355365628484592013285671, 1661534994731144841129758825350430720000000000000⟩
private theorem leaf7_rate {α x : ℝ}
    (hα : α ∈ Set.Icc leaf7A0.real leaf7A1.real)
    (hx : x ∈ Set.Icc leaf7X0.real leaf7X1.real) :
    ScalableDensePointRate α x slope.real intercept.real
      (1 / 10000000) 10000000 := by
  have h := scalable_dense_checked_leaf slope intercept
    w82P w82Y w82Radius w82Z
    leaf7A0 leaf7A1 leaf7X0 leaf7X1
    w82Vector (10000000 : ℝ)
    (by decide) (by decide) (by decide) (by decide)
    (by decide) (by decide)
    (by decide) (by decide) (by decide) (by decide)
    (by norm_num [w82P, QInput.real])
    (by norm_num [w82P, QInput.real])
    (by norm_num [w82Y, QInput.real])
    (by norm_num [w82Y, QInput.real])
    (by norm_num [w82Radius, QInput.real])
    (by norm_num [w82Z, QInput.real])
    (by norm_num [w82Z, QInput.real])
    (by norm_num [leaf7A0, QInput.real])
    (by norm_num [leaf7A1, QInput.real])
    (by norm_num [leaf7X0, QInput.real])
    (by norm_num [leaf7X1, QInput.real])
    (by norm_num : (0 : ℝ) ≤ 10000000)
    (by intro j; fin_cases j <;> norm_num [w82Vector, QInput.real])
    w82Transfer vertex9 vertex11 vertex6 vertex8 hα hx
  simpa [w82Vector, QInput.real] using h

private theorem leaf7_claim : DenseClaim target
    ⟨(5110001 / 5120000 : ℚ), (1 / 1 : ℚ), (316186482737629508933514007290843083540271318859 / 664613997892457936451903530140172288000000000000 : ℚ), (790938778754819906208810355365628484592013285671 / 1661534994731144841129758825350430720000000000000 : ℚ)⟩ := by
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

def root : Rect := ⟨(1270001 / 1280000 : ℚ), (1 / 1 : ℚ), (49374602183332977403672480068036683319333956789 / 103845937170696552570609926584401920000000000000 : ℚ), (790938778754819906208810355365628484592013285671 / 1661534994731144841129758825350430720000000000000 : ℚ)⟩

theorem root_claim : DenseClaim target root := by
  exact (denseClaim_alpha target ⟨(1270001 / 1280000 : ℚ), (1 / 1 : ℚ), (49374602183332977403672480068036683319333956789 / 103845937170696552570609926584401920000000000000 : ℚ), (790938778754819906208810355365628484592013285671 / 1661534994731144841129758825350430720000000000000 : ℚ)⟩ (2550001 / 2560000 : ℚ) ((denseClaim_density target ⟨(1270001 / 1280000 : ℚ), (2550001 / 2560000 : ℚ), (49374602183332977403672480068036683319333956789 / 103845937170696552570609926584401920000000000000 : ℚ), (790938778754819906208810355365628484592013285671 / 1661534994731144841129758825350430720000000000000 : ℚ)⟩ (316186482737629508933514007290843083540271318859 / 664613997892457936451903530140172288000000000000 : ℚ) ((denseClaim_alpha target ⟨(1270001 / 1280000 : ℚ), (2550001 / 2560000 : ℚ), (49374602183332977403672480068036683319333956789 / 103845937170696552570609926584401920000000000000 : ℚ), (316186482737629508933514007290843083540271318859 / 664613997892457936451903530140172288000000000000 : ℚ)⟩ (5090003 / 5120000 : ℚ) (leaf0_claim) (leaf1_claim))) ((denseClaim_alpha target ⟨(1270001 / 1280000 : ℚ), (2550001 / 2560000 : ℚ), (316186482737629508933514007290843083540271318859 / 664613997892457936451903530140172288000000000000 : ℚ), (790938778754819906208810355365628484592013285671 / 1661534994731144841129758825350430720000000000000 : ℚ)⟩ (5090003 / 5120000 : ℚ) (leaf2_claim) (leaf3_claim))))) ((denseClaim_density target ⟨(2550001 / 2560000 : ℚ), (1 / 1 : ℚ), (49374602183332977403672480068036683319333956789 / 103845937170696552570609926584401920000000000000 : ℚ), (790938778754819906208810355365628484592013285671 / 1661534994731144841129758825350430720000000000000 : ℚ)⟩ (316186482737629508933514007290843083540271318859 / 664613997892457936451903530140172288000000000000 : ℚ) ((denseClaim_alpha target ⟨(2550001 / 2560000 : ℚ), (1 / 1 : ℚ), (49374602183332977403672480068036683319333956789 / 103845937170696552570609926584401920000000000000 : ℚ), (316186482737629508933514007290843083540271318859 / 664613997892457936451903530140172288000000000000 : ℚ)⟩ (5110001 / 5120000 : ℚ) (leaf4_claim) (leaf5_claim))) ((denseClaim_alpha target ⟨(2550001 / 2560000 : ℚ), (1 / 1 : ℚ), (316186482737629508933514007290843083540271318859 / 664613997892457936451903530140172288000000000000 : ℚ), (790938778754819906208810355365628484592013285671 / 1661534994731144841129758825350430720000000000000 : ℚ)⟩ (5110001 / 5120000 : ℚ) (leaf6_claim) (leaf7_claim))))))

theorem root_claim_global :
    DenseClaim (DenseTarget (1 / 10000000) 10000000) root := by
  exact root_claim

#print axioms root_claim

#print axioms w74Transfer
#print axioms vertex0
#print axioms vertex20

end Spin.UniformMixer.ScalableDensePilotLeft
