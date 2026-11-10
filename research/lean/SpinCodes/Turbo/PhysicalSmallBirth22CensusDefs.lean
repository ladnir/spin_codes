import SpinCodes.Turbo.PhysicalSmallBirth22Data
import SpinCodes.Turbo.PhysicalSmallBirth22NumericCensusDefs

namespace Spin.Turbo
open Spin.UniformMixer Spin.Structured

theorem smallBirthOneCode22_actual (a : SmallBirthOne22) :
    smallBirthOneCode22 a =
      9 * (physicalBirthIndex22 (stepFeedback (packetStepInput22 (smallBirthOneWord22 a)))).val +
        stepWeight (packetStepInput22 (smallBirthOneWord22 a)) 0 := by
  unfold smallBirthOneCode22
  rw [← smallBirthSingleImages22_checked a, MapSpectrum.fastWeight_eq,
    smallBirthPackedExpansionWeight22, smallBirthOnePacked22_input,
    smallBirthClassWeight22_actual, smallBirthOneWord22_weight, smallBirthLabel22_weight]

theorem smallBirthTwoCode22_actual (a : SmallBirthTwo22) :
    smallBirthTwoCode22 a.1.val.1 a.1.val.2 a.2.1 a.2.2 =
      9 * (physicalBirthIndex22 (stepFeedback (packetStepInput22 (smallBirthTwoWord22 a)))).val +
        stepWeight (packetStepInput22 (smallBirthTwoWord22 a)) 0 := by
  unfold smallBirthTwoCode22
  rw [← smallBirthSingleImages22_checked (a.1.val.1,a.2.1),
    ← smallBirthSingleImages22_checked (a.1.val.2,a.2.2),
    ← PackedMap.eval_xor, ← PackedMap.eval_xor]
  change 9 * (smallBirthClassWeight22 (MapSpectrum.fastWeight 64
    (PackedMap.eval productionPackedRows (PackedMap.eval smallBirthFeedbackRows22
      (smallBirthTwoPacked22 a))))).val + _ = _
  rw [MapSpectrum.fastWeight_eq, smallBirthPackedExpansionWeight22,
    smallBirthTwoPacked22_input, smallBirthClassWeight22_actual,
    smallBirthTwoWord22_weight, smallBirthLabel22_weight, smallBirthLabel22_weight]

end Spin.Turbo


