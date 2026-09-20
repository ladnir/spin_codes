# Removing slack from the K=2^16 IMT bound

This study keeps the selected half-rate BCH-256 outer, IMT maps, routing law,
and (t,s)=(128,19) fixed. It changes only the analysis. The accepted paper
certificate remains `FULL_M16_VERIFIED.json`; new records are separate overlays.
No encoder, performance result, or data release is changed.

## Results

Margins are minus log2 of rational upper bounds, not estimates of true failure.

| Analysis | Q1 margin | Dense margin | Combined full margin |
|---|---:|---:|---:|
| Accepted certificate | 43.5928335021 | 42.3171086801 | 41.8183267960 |
| Expanded tilt search | 43.5952368529 | unchanged | 41.8190288676 |
| Retain injection shells | 44.3853417276 | unchanged | 42.0085686753 |
| Injection shells plus eight dense replacements | 44.3853417276 | 42.7724739304 | 42.3643587462 |

The sparse Q=2..63 contribution is retained exactly in every full sum.
The final gain is 0.5460319502 bits. All new numerical bounds passed 512-bit
replay. The Q1 replays use linear region products instead of repeated squaring.
The original 746-file integration authentication also passes.
These refinements are not yet inputs to the paper's evidence adapter.

The bounded tilt search evaluates 39 witnesses, extending the six original
ones. Its small gain shows that this particular tuning grid was not the
main source of slack. It does not prove global optimality of the tilts.

## Retaining the state shell after activation

Use the recurrence Y=X+Aq, q'=Mq+CX. For each fixed nonzero q, the fresh
transvection has marginal law

    Mq ~ (1/2) delta_q + (1/2) Uniform(nonzero states).

The existing representation has zero mass Z, arbitrary nonzero mass D, and
uniform-shell majorants C_v. Here v=wt(Aq), and a_v counts nonzero states
in that shell. The selected shells are 48, 56, 64, 72, and 80.
The notation C_v denotes a coordinate, not the feedback map C.

The old zero-to-singleton transition assigned all mass to D. It therefore
forgot which image weights were possible immediately after activation.
Direct evaluation of the 128 feedback columns gives:

| Image weight wt(A C e_i) | 48 | 56 | 64 | 72 | 80 |
|---|---:|---:|---:|---:|---:|
| Number of positions i | 3 | 32 | 59 | 32 | 2 |

Add coordinates H_v. Each H_v bounds the total mass of an arbitrary measure
supported on shell v; it does not assert uniformity within that shell.
The refined representation bounds a weighted measure by

    Z delta_0 + nu_D + sum_v C_v U_v + sum_v nu_v,

where nu_D has nonzero support and mass at most D, U_v is uniform on shell v,
and nu_v has support in shell v and mass at most H_v.

Fix z=exp(-lambda), lambda>0, and m=2^19-1. All transfer entries include the
factor z^wt(Y). The original rows remain unchanged except singleton activation
from zero. If h_v counts feedback positions in shell v, that row becomes

    Z -> H_v: z h_v / 128.

This is exact: the output is the singleton input and the next state is C e_i.
The feedback columns are nonzero and distinct.

For zero input, a state in H_v emits weight v. The lazy branch stays in its
shell, and the fresh branch is uniform over nonzero states. Thus

    H_v -> H_v: z^v / 2,
    H_v -> C_w: z^v a_w / (2m).

For a uniform singleton input, its emitted moment at any fixed state in shell v is

    f_v = (v z^(v-1) + (128-v) z^(v+1)) / 128.

Let c_v be the maximum of z^wt(e_i+A C e_i)/128 over feedback positions
whose image weight is v, or zero if there are none. Distinct feedback columns
imply that at most one singleton cancels a given state on the lazy branch.
The following row is therefore valid:

    H_v -> Z:   c_v/2 + f_v/(2m),
    H_v -> D:   f_v/2,
    H_v -> C_w: f_v a_w/(2m).

To justify the fresh branch, condition on the singleton input. Translating
a uniform nonzero state gives probability at most 1/m at every target.
Average the emitted weight over singleton inputs to obtain f_v/m.
The lazy nonzero mass is at most f_v/2 and may safely be assigned to D.
Its separate cancellation allowance can double-count mass, which is conservative.

Every row bounds the appropriate weighted measures, so positive products
preserve domination. Q1 has at most one input one per region. The existing
region products and outer-weight coefficient recurrence therefore apply to
these twelve coordinates. Terminal summation bounds the total weighted mass.
The same BCH spectrum inequalities then bound the Q1 union.

`q1_injection_shells.py` implements this extension. Exact rational tests
enumerate all source states and singleton inputs in a small independent-map
example. They check domination by the represented measures, not entrywise
comparison between two different representations.

## Dense replacements

The accepted dense cover has 1,096 leaves. Its eight largest power-of-two
ceilings account for 97.3305673540% of the dense upper bound.
`dense_slack.py` reuses their witnesses and the existing coupled-input bound.
For each selected rectangle, it evaluates:

1. the same parent bound with a finer dyadic ceiling;
2. the sum of two bounds obtained by bisecting its density interval.

The replacement is the minimum of these two complete bounds and the old
ceiling. Both children cover the parent; possible shared boundaries only
overcount. Untouched leaves keep their old ceilings. Exact rational arithmetic
replaces the selected terms in the original dense sum and combines it with
the new Q1 term and original sparse term.

Bisection helps two of these eight rectangles and worsens the isolated bound
for the others. Keeping the old parent as a candidate prevents regression.
The overlay authenticates its source records and checks the original complete
partition; it does not infer full coverage from the eight selected regions.

## Reproduction

Run from the permute_conv repository root, with python-flint and the retained
original evidence available. Use fresh output names when generating new records.
The following commands replay the local records without overwriting producers:

```text
python -B workstreams/inner_design/finite_migration/q1_slack.py --output workstreams/inner_design/finite_migration/Q1_SLACK_v1.json --verify
python -B workstreams/inner_design/finite_migration/q1_injection_shells.py --output workstreams/inner_design/finite_migration/Q1_INJECTION_SHELLS_v1.json --verify
python -B workstreams/inner_design/finite_migration/dense_slack.py --output workstreams/inner_design/finite_migration/DENSE_SLACK_v1.json --verify
python -B -m unittest discover -s workstreams/inner_design/finite_migration -p test_q1_slack.py
python -B paper/check_finite_integration.py
```

Replay writes a separate `_replay.json` file and refuses to overwrite it.
To regenerate, omit `--verify` and supply fresh output paths. The dense overlay
currently reads the named injection-shell record above and its replay receipt.
No numerical data from this study should be committed with the source changes.

## Interpretation and next step

Some of the short-length loss is demonstrably analysis slack: the encoder did
not change when these bounds improved. This does not determine the true failure
probability or eliminate the gap to the previous inner's certificate.

Next, retune the witnesses of the few dominant dense rectangles rather than
replaying or subdividing the entire cover. Any further Q1 refinement should
retain more information after subsequent singleton inputs, where H_v still
falls back to D. Check the complete union before investing in either track;
improving only one contribution eventually leaves the other dominant.
