import SpinCodes.UniformMixer.ScalableDenseCoverCore
import SpinCodes.UniformMixer.Generated.ScalableDensePilotLeft
import SpinCodes.UniformMixer.Generated.ScalableDensePilotRight

/-! Two independently checked numerical chunks compose across one exact cut. -/

namespace Spin.UniformMixer.ScalableDenseAssemblyPilot

open Spin.Structured.DenseGeometry

def root : Rect := ⟨(1270001 / 1280000 : ℚ), (1 / 1 : ℚ), (49374602183332977403672480068036683319333956789 / 103845937170696552570609926584401920000000000000 : ℚ), (395941961288156086979430514821335018037341631359 / 830767497365572420564879412675215360000000000000 : ℚ)⟩

theorem root_claim :
    DenseClaim (DenseTarget (1 / 10000000) 10000000) root := by
  apply denseClaim_density (DenseTarget (1 / 10000000) 10000000) root (790938778754819906208810355365628484592013285671 / 1661534994731144841129758825350430720000000000000 : ℚ)
  · simpa [root, Spin.UniformMixer.ScalableDensePilotLeft.root] using Spin.UniformMixer.ScalableDensePilotLeft.root_claim_global
  · simpa [root, Spin.UniformMixer.ScalableDensePilotRight.root] using Spin.UniformMixer.ScalableDensePilotRight.root_claim_global

#print axioms root_claim

end Spin.UniformMixer.ScalableDenseAssemblyPilot
