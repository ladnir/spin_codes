# A complete sparse-cohort proposal at seventeen active groups

The fixed-r17 sparse family now has a complete numerical proposal below
−50 bits. The new terms are FLOAT calculations and still require directed
replay. They do not establish the whole XOR2 distance target.

The construction remains
`xor2-gf256-h4-alternating-swapped-pairs-k20-rate-half`, with forbidden
weight 734003. Write \((n_1,n_2,n_H)\) for the counts of active groups
having one, two, or at least three nonzero RS constituents. This note
covers \(n_H=0\) and \(n_1+n_2=17\). It retains every indexed source
support and every fixed constituent-row assignment. The
[coverage roadmap](HD_XOR_COVERAGE_PLAN.md) specifies the remaining
whole-code obligations and the preceding failed controls.

## Bounded joint fit

The two weakest new cells had ordered source matching counts \((a,b)=(0,7)\):
cohort \((2,15,0)\) with multiplicity \(\mu=15\), and cohort \((3,14,0)\)
with \(9\le\mu\le14\). The objective was the log2 sum of their class
bounds, including all source, cohort, constituent, projective, and potential
factors. A trial used one common Chernoff mark \(z\) and potential
\(v=(1,2^p)\) for both cells.

The fit started from the preceding literal witness and allowed exactly
32 objective evaluations. Its domain was \(0.99\le z\le0.9998\) and
\(-100\le p\le-20\). Each trial built one regional bank and one
`MultiplicityCapsFloat` evaluator, then reused its \(\mu=14,15\) caps.
There was no source-class sweep during fitting.

The best trial was number 24:

\[
z=\frac{1122141909660549}{1125899906842624},\qquad
p=-\frac{7932457}{131072}.
\]

The objective improved from −35.49648830521795 to −249.23496677842917,
a gain of 213.73847847321122 bits. The two cell bounds became
−250.9893055561747 and −249.74215787146204. The bounded fit took
12.39 seconds. All 32 trials, including unsuccessful ones, are retained.

## One frozen-bank screen

The gain exceeded the predeclared 10-bit gate, so the winning literal
witness was screened once across all 81 source classes and all 18 sparse
cohorts. No further fitting occurred. The run built one raw regional bank
and 81 source-conditioned evaluators. It used the exact canonical-cap
reuse identity established by the earlier tests: truncated multiplicity
bins require only upper caps 8, 14, 15, 16, and 17.

The positive support polynomials and class counts are unchanged. For a
cohort \((n_1,n_2,0)\), the aggregate RS support floor is
\(17(n_1+2n_2)\), and the normalization exponent is
\(16(n_1+2n_2)\). Each class pays its exact source count,
\(\binom{17}{n_1}\), \(8^{n_1}\), and fixed-pair assignment count once.
Its row sets remain fixed across all 32 columns. The source counts sum
to \(\binom{8192}{17}\), and the complete partition has 2673 cells.

The screen took 99.64 seconds. It retained its raw new-bank values, then
selected the smallest value per cell from four sources of evidence:

- the frozen two-bank FLOAT screen;
- this one new fixed-bank FLOAT screen;
- the saved sibling \((1,16,0)\) FLOAT schedule;
- matching 256- and 384-bit directed homogeneous cell ceilings.

The input artifacts are read-only and SHA-256 pinned. Selection is a
deterministic choice of proof witness for each message class, not a choice
based on sampled setup. It therefore incurs no witness-bank multiplicity.
The resulting numerical union is explicitly marked `MIXED-PRELIMINARY`:
using some directed inputs does not certify its new FLOAT terms.

| Complete-source family | Selected log2 proposal |
|---|---:|
| \((0,17,0)\), with cellwise selection | −50.89537379226055 |
| \((1,16,0)\) | −88.58262308674374 |
| \((2,15,0)\) | −157.1945677323168 |
| \((3,14,0)\) | −158.78720114244052 |
| All 17 nonhomogeneous sparse cohorts | −88.58262308674374 |
| All 18 sparse cohorts | −50.895373792254034 |

The original homogeneous directed aggregate can instead remain untouched.
Combining its ceiling −50.891481584655 with the selected nonhomogeneous
FLOAT sum gives −50.891481584648496. Thus the proposal does not depend on
improving the already certified homogeneous family.

The new bank alone still fails: its full sparse union is +145.5717852715084,
and its nonhomogeneous union is +18.644286257894294. These raw controls
remain in the artifact. The selected schedule uses 1861 old two-bank cells,
532 new-bank cells, 207 directed homogeneous cells, and 73 sibling FLOAT
cells. Every selected literal and its source are recorded.

The complete proposal is now dominated by directed homogeneous terms.
Its largest cell is source \((4,0)\), cohort \((0,17,0)\), \(\mu=17\),
at −51.614684956866. The largest selected nonhomogeneous cell is source
\((8,6)\), cohort \((1,16,0)\), \(\mu=16\), at −89.42241976617515.
No additional witness search followed the full screen.

## Directed replay without rebuilding every cell

A replay should accept the saved literal schedule, not saved FLOAT values.
For a selected cohort range, validate the complete set of keys
\((a,b,n_1,\text{multiplicity bin})\) and reconstruct all exact counts.
Use one declared Arb precision context for the entire replay.

First build `arb_pair_regions(z,17)` once per distinct rational \(z\).
For each ordered source class and literal pair \((z,p)\), construct the
labelled regional subset tables and one `MultiplicityCapsArb` evaluator.
Cache each requested canonical upper cap. Reuse those caps in
`sparse_count_arb(caps,n1,n2)` from the independently reviewed
`hd_xor_sparse_mixed_certificate.py`. That function builds the positive
column polynomial, takes its 32nd power, retains the aggregate RS tail,
and applies the exact normalization. Class multiplicities and the single
projective/potential endpoint factors are applied outside that polynomial.

The selected 17-cohort nonhomogeneous schedule has 2268 cells, but needs
only six distinct raw regional banks, 309 source/witness evaluators, and
514 canonical cap requests. It does not require 2268 fresh evaluator builds.
No symmetry between the two matching orders is assumed.

A cheaper preferred replay also preserves the complete \((1,16,0)\)
calculation already running separately. Replay only \(n_1\ge2\):

| Fresh replay object | Count |
|---|---:|
| Class cells | 1944 |
| Distinct rational \(z\) values | 3 |
| Ordered-source/witness evaluators | 242 |
| Requested canonical caps | 362 |

Those 16 cohorts have FLOAT union −156.78144282140883. Combining that
proposal with the preserved homogeneous ceiling and the sibling FLOAT
aggregate −70.35414804665 gives −50.89147958788594. This also passes the
−50 proposal target. It avoids replaying the improved sibling selection
merely to reduce an already negligible contribution. The sibling term
must remain labelled FLOAT until its separate directed replay completes.

At 256 and 384 bits, a new adapter should independently recompute every
selected remaining class, compare exact rational cell ceilings, and sum
the rounded masses in a fresh higher-precision context. The adapter must
identify its cohort range so that the preserved families are charged once.
No such new directed replay was launched for this note. Other group counts
and cohorts with dense groups remain outside this fixed-r17 sparse result.

## Reproduction and pins

Seven new tests cover the bounded fit, shared trial bank, cap identity,
scope, input pins, and complete saved-artifact recombination. Together with
the eight shared-cohort tests, all 15 tests pass. No test rebuilds a full
operator bank; missing optional artifact fixtures are skipped explicitly.

From `theory`, use fresh output paths:

```powershell
C:/Python314/python.exe -B -m unittest hd_xor_sparse_cohort_refine_tests hd_xor_sparse_cohorts_tests
C:/Python314/python.exe -B hd_xor_sparse_cohort_refine.py --fit-output refine_fit_new.json --screen-output refine_screen_new.json
```

The second command saves its fit first. It performs the full screen only
if that bounded fit clears the 10-bit gain gate. Existing output paths
are rejected.

| File in `theory` | SHA-256 |
|---|---|
| `hd_xor_sparse_cohort_refine.py` | `8fd7dbf51bfdc732c8e3a5bccbe217fb8b617083349b9f4113b2e0c18d223fc4` |
| `hd_xor_sparse_cohort_refine_tests.py` | `219e8c29afbbf9c3e0b899d05f099a305496dcd468c260761d5f4315772e5570` |
| `hd_xor_sparse_cohort_refine32_0_7_fit.json` | `e06525fc59ecf60ca27d502c58fcd9626a56ec6744342a93f18e5350eb16a90f` |
| `hd_xor_sparse_cohort_refine32_all_sources_float.json` | `116383e032d3813361dcdaffcc8a32e53f5c06876536914d3ded25decc8013e1` |

Next: finish the separate sibling receipt, then replay the remaining
\(n_1\ge2\) cohort range using shared directed cap banks. After that
fixed-r17 milestone, extend the coverage method across group counts and
establish a multi-dense-group interface; do not substitute this sparse
family for the whole-code goal.
