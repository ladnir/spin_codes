/-
Shared parameters of the interval covers.
-/
import SpinCodes.Numeric.BAEvalDefs

namespace Spin.Cover

/-- `-7.68·10⁻⁸`, the paper's constant, at `scale = 10^30`. -/
def thr : Int := -76800000000000000000000

/-- Series terms per logarithm.  Ten gives a worst-case enclosure width of
`3·5⁻²¹ ≈ 3.2·10⁻¹³`, against a margin to the threshold of about `6·10⁻¹¹`. -/
def terms : Nat := 10

end Spin.Cover
