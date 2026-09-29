import SpinCodes.Native

/-!
# Audit of the public native result

Compile this module to print the full types and recursive axiom dependencies of
the four public results. The declarations are accessed through `SpinCodes.Native`,
so this also checks the documented public import. The original explicit statement
pins remain in `Structured.ConcreteNativeTheoremPin` and
`Structured.ConcreteNativeDistanceConsequencesFinalPin`.

Successful import checking reuses built dependencies. It is distinct from a clean
source replay; see `FINAL_REPRODUCTION.md` for both workflows.
-/

open Spin.Structured.ConcreteNativeFamily

-- Failure probability for the original outer/inner product law.
#check @native_minimum_distance_failure_tendsto
#print axioms native_minimum_distance_failure_tendsto

-- Exact rate at every native index, for every setup realization.
#check @realizedCode_rate
#print axioms realizedCode_rate

-- Strict relative-distance success under the same law.
#check @relative_distance_success_tendsto
#print axioms relative_distance_success_tendsto

-- Eventual existence includes positive mass under the actual setup law.
#check @eventually_exists_rate_half_distance_gt_eleven_percent
#print axioms eventually_exists_rate_half_distance_gt_eleven_percent
