import SpinCodes.Turbo.RS8SelectedValues22
import SpinCodes.Turbo.PhysicalThinning22

/-!
# The regional value law at a physical 16-slot step

The generic selected-value law uses exactly the same complete four-bit
experiment as the previously checked physical-step thinning module.
This identity pins the regional route bridge to the existing operator
input law; the conditional distribution of selected sets across steps
is a separate obligation.
-/

noncomputable section

namespace Spin.Turbo

theorem regionalPacketSupport22_step (x : PacketStep16) :
    regionalPacketSupport22 x = packetStepSupport22 x := by
  rfl

theorem regionalSelectedGiven22_step
    (T : Finset (Fin 16)) :
    regionalSelectedGiven22 T = selectedGiven22 T := by
  rfl

#print axioms regionalSelectedGiven22_step

end Spin.Turbo
