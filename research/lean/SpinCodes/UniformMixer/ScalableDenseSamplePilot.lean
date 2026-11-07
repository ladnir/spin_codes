import SpinCodes.UniformMixer.ScalableDenseLeaf

/-! Resource pilot: two adjacent boxes sharing one transfer witness. -/

noncomputable section
set_option maxRecDepth 8192
set_option maxHeartbeats 1000000

namespace Spin.UniformMixer

open Spin.Numeric Spin.Structured.DenseOccupationFixed

private def sampleSlope : QInput := ⟨-833, 500⟩
private def sampleIntercept : QInput := ⟨373189, 250000⟩
private def sampleP : QInput := ⟨618749999097849, 625000000000000⟩
private def sampleY : QInput := ⟨12874999158079689, 25000000000000000⟩
private def sampleZ : QInput :=
  ⟨61797754152116686151078193289, 500000000000000000000000000000⟩
private def sampleRadius : QInput :=
  ⟨93998318722704, 1000000000000000000000000000000⟩
private def sampleVector : Fin 3 → QInput :=
  ![⟨1, 1⟩, ⟨2694812531397891, 1000000000000000⟩,
    ⟨126938283539277, 125000000000000⟩]

private theorem sample_transfer :
    transferCheck (Fix.mul sampleP.interval sampleY.interval)
      sampleZ.interval sampleRadius.interval
      (fun i => (sampleVector i).interval) = true := by
  decide +kernel

private def sampleA0 : QInput := ⟨310001, 320000⟩
private def sampleA1 : QInput := ⟨1, 1⟩
private def sampleX0 : QInput := ⟨854343973879, 960000000000⟩
private def sampleX1 : QInput := ⟨112, 125⟩

private theorem sample00 :
    scalableVertexCheck sampleSlope sampleIntercept sampleP sampleY
      sampleRadius sampleZ sampleA0 sampleX0 50
      (-100000000000000000000000) = true := by decide +kernel
private theorem sample01 :
    scalableVertexCheck sampleSlope sampleIntercept sampleP sampleY
      sampleRadius sampleZ sampleA0 sampleX1 50
      (-100000000000000000000000) = true := by decide +kernel
private theorem sample10 :
    scalableVertexCheck sampleSlope sampleIntercept sampleP sampleY
      sampleRadius sampleZ sampleA1 sampleX0 50
      (-100000000000000000000000) = true := by decide +kernel
private theorem sample11 :
    scalableVertexCheck sampleSlope sampleIntercept sampleP sampleY
      sampleRadius sampleZ sampleA1 sampleX1 50
      (-100000000000000000000000) = true := by decide +kernel

private def sampleB0 : QInput := ⟨2470009, 2560000⟩
private def sampleB1 : QInput := ⟨310001, 320000⟩
private def sampleY0 : QInput := ⟨6875463973879, 7680000000000⟩

private theorem sample20 :
    scalableVertexCheck sampleSlope sampleIntercept sampleP sampleY
      sampleRadius sampleZ sampleB0 sampleY0 50
      (-100000000000000000000000) = true := by decide +kernel
private theorem sample21 :
    scalableVertexCheck sampleSlope sampleIntercept sampleP sampleY
      sampleRadius sampleZ sampleB0 sampleX1 50
      (-100000000000000000000000) = true := by decide +kernel
private theorem sample30 :
    scalableVertexCheck sampleSlope sampleIntercept sampleP sampleY
      sampleRadius sampleZ sampleB1 sampleY0 50
      (-100000000000000000000000) = true := by decide +kernel
private theorem sample31 :
    scalableVertexCheck sampleSlope sampleIntercept sampleP sampleY
      sampleRadius sampleZ sampleB1 sampleX1 50
      (-100000000000000000000000) = true := by decide +kernel

#print axioms sample_transfer
#print axioms sample00
#print axioms sample31

end Spin.UniformMixer
