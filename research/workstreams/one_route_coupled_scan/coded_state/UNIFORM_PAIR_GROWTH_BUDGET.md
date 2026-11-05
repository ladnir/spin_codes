# Couple terminal growth to the early-episode budget

The unchanged three-shear model now closes every no-full-reset history
with at most eight defective incoming batches at 46% distance, with
failure probability at most 2^-78. The probability is over setup and
quantifies over nonzero J-valued messages, J=GF(2^32). This does not cover
all eight-defect histories: the full-reset branch still has 69 residual
cells. The earlier complete seven-defect certificate is preserved.

The construction, state size, RS outer, and regional route are unchanged.
This continuation refines the [localized-growth proof](UNIFORM_PAIR_LOCALIZED_GROWTH.md).

## Why the first coupling was insufficient

Fix h,a,c,T as before: h bounds defective incoming batches, a indicates
an active first state, c counts live intervals, and T counts nonzero
incoming states. There are e=h-c+a internal defects after paying for
activations from zero. In the residual no-full-reset profiles, the proper-
state bound forces a full state, so the final live interval ends at Q_(L-1).
All earlier live intervals are proper and closed.

Let E be the total length of those early intervals, u the minimum number
of internal defects their lengths require, and R_E their retained regional
cover bound. The previous per-batch cap S implies that reaching 512 state
coordinates requires at least g defects, where g is the smallest integer
with (2^(g+1)-1)S>=512. Thus u+g<=e is necessary.

This constraint alone does not improve the worst eight-defect profile:
a=1,c=7,T=1955 retains exponent +14000. Its maximizing early profile already
has u=0. The issue is not just allocation of the defect budget. The cap S
was maximized over input batches unrelated to the final interval where
growth must occur.

## Retain the final interval's complete regions

For a fixed endpoint tuple, the terminal input interval has

    B_T=T-E+1=16F+r batches,  0<=r<16.

It ends at the final regional boundary. Therefore it contains F completely
included regions and, if r>0, one initial region with at least r/16 of
its positions in the input support I. The complete input support has mass
mu=(T+c)/16 and at most J_0=2c-a-1 partial regions. It touches at most
R_E+F+1[r>0] regions.

In the retained row moment calculation, restrict the count of complete
regions to at least F and the total touched-region count to the bound
above. After removing a chosen region, adjust these bounds accordingly.
The verifier still maximizes over all admissible complete/partial profiles;
it does not assume balanced occupancies or independent row events.

For a batch in a complete terminal region, remove one guaranteed complete
region. The eligibility test needs at least 128 hits among the remaining
regions. Multiply its joint-row moment by (1/16)^s for membership of s
specified rows in the chosen batch, exactly as in the localized proof.

For a batch in the initial terminal region, remove a region whose occupancy
in I is a/16 for some a in {r,...,16}. All F guaranteed complete terminal
regions remain. This remains true if earlier episodes fill the initial
terminal region completely in I. Maximize over these possible a values.

The resulting caps are P for initial-region batches and U for complete-
region batches. They bound actual input weights simultaneously on the
new good-route event. This separation is important: activation occurs in
the initial region, not at whichever batch maximizes the global cap.

## Respect the order of the two input regions

Only e'=e-u internal defects remain available in the terminal live interval.
When r>0, activation consumes one of the r initial-region input batches.
At most r-1 subsequent defective updates can use inputs from that region.
Every such update precedes every update whose input lies in a complete
terminal region.

For k initial-region defective updates followed by e'-k complete-region
defective updates, the deterministic state-growth bound is

    W_k = 2^(e'-k) (2^(k+1)-1) P + (2^(e'-k)-1) U.

Maximize over 0<=k<=min(e',r-1). Clean updates cannot increase support,
and allowing all e' defects is conservative. If r=0, use
(2^(e'+1)-1)U instead. A maximum below512 excludes the endpoint/early-cost
profile on the good route. Otherwise retain its previous rank and volume
bound; reaching512 is necessary, not sufficient for a bad word.

## Keep the fixed endpoints in the route union

The generic route union used up to h+2 arbitrary intervals. It is larger
than needed here. The c live intervals have a fixed final end, and their
first start is fixed when a=1. There are at most

    (L+1)^(2c-a-1) < 2^(13(2c-a-1))

endpoint tuples for a given a,c. The tuple determines E, the early lengths,
their minimum defect cost u, and the regional-cover bound used in the cap.
These quantities need no additional independent union.

There are at most 2(h+1) choices of a,c, at most L chosen batches, and two
cap types. Request

    13(2c-a-1) + ceil(log2(2(h+1))) + 12 + 1 + 80

bits in each cap. Summing over all tuples, batches, and types leaves a
common exception at most2^-80. This retains the 80-bit route budget;
it does not relax the target to obtain the result.

The endpoint tuple is fixed before routing is sampled. Simultaneity then
allows the actual tuple and defect placement to depend on the sampled
code and message. Restricting the old inner count to this good-route event
requires no independence between old and new route events.

## Finite replay and intermediate diagnostics

The previous no-full-reset h=8 verifier left45 profiles after its first
localized exclusion. Budget coupling alone leaves the worst unchanged.
Retaining terminal-region geometry resolves28 and leaves17, all with
a=1,c=6 and three internal defects.

For the latter's T=1955 profile, the split-cap count with the generic
endpoint union gives P=33,U=36,r=15. Its state envelope can still reach516.
The sharper endpoint union resolves this remaining loss without changing
the construction or the route exception budget.

All45 retained profiles then pass. The worst new exponent is -1909 at
a=1,c=6,T=1959. Its retained early cell has E=39,u=0,R_E=9 and caps
P=41,U=29,r=1. Its growth envelope exceeds512, so this cell is handled by
the old rank/volume bound rather than falsely declared impossible.

The complete replay recomputes the full range instead of assuming the
saved residual list:

- 69551 profile cells checked;
- 706 old residuals;
- 661 excluded by the preceding localized bound;
- the remaining45 covered by this refinement;
- no residual, with overall good-route exponent at most -230.

The old regional-rank event, the preceding localized event, and the new
split-cap event each have exception probability at most2^-80. Therefore

    2^-230 + 3*2^-80 <= 2^-78.

This proves the stated no-full-reset result at W=floor(.46N). It is not
a whole-code certificate and does not replace the extension-field analysis.

## Full-reset diagnostic and next step

The unmodified full-reset refinement was also run at h=8. It reduces341
old residuals to69. The remaining groups, with noncontiguous T values
allowed between the listed endpoints, are:

| a,c | Count | First T | Last T | Largest remaining exponent |
| --- | ---: | ---: | ---: | ---: |
| 0,6 | 5 | 1965 | 1969 | +3295 |
| 0,7 | 2 | 1967 | 1968 | +11956 |
| 1,5 | 5 | 1983 | 1987 | +7960 |
| 1,6 | 14 | 1970 | 1986 | +13662 |
| 1,7 | 24 | 1959 | 1985 | +21279 |
| 1,8 | 19 | 1963 | 1984 | +37396 |

These are failed upper bounds, not exhibited low-weight outputs. A full
state followed by zero can lie in an earlier interval, so the terminal-end
assumption above must not be imported into this branch. Next retain the
actual interval that grows to full support and then resets. First check
its endpoint count and local input geometry; then share its growth and
decay defect budgets with the other intervals. Higher defect counts remain
open. No structural requirement for redesign has been established.

## Source and checks

`diagnose_pair_growth_budget.py --test` passes3071 terminal-geometry checks,
1375 ordered-growth checks, and375 restricted one-row moment checks.
The joint-row justification is the retained positive-moment argument,
not an inference of independence from these tests.

`--frontier --split` checks the selected45 profiles; `--cover` recomputes
the entire h=8 no-full-reset range. `--terminal` reproduces the intermediate
single-cap test. `--split --coarse-endpoints` preserves the split-cap
comparison with the old endpoint union. The full cover took about48.3s
locally; its cap is60s. Selected/frontier passes have30s/45s caps.

The h8 full-reset diagnostic uses the preserved
`diagnose_pair_full_reset_frontier.py --integrate --defects 8` and took41.1s
under its45s cap. All jobs were serial, used Python-B and bounded caches,
and wrote no numerical files. These are verifier times, not encoder timings.
The construction, certified control, and frozen earlier proofs are preserved.
