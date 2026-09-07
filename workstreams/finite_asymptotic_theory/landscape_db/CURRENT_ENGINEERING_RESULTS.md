# Current engineering results: BCH, RM, and random constituents

Snapshot: 2026-09-07. This account covers exact-spectrum BCH through length
128, RM, and random constituents. BCH-256 is excluded by the requested scope.
This is the entry point for the current parameter study, its complete BCH
follow-up, and the separate fixed-RM finite reference certificate.
The numerical appendix collects [all 130 BCH geometries and matched family
references](CURRENT_RESULTS_TABLES.md). The broader historical research archive
is indexed separately; an old ledger's empty certificate view is not a claim
that the repository contains no certificates.

The main result is a regime-dependent engineering picture. BCH-128 retains
negligible higher-occupation corrections in the useful tested region, whereas
BCH-64 develops a visible correction as the message grows or the state shrinks.
Larger epochs can create a first-moment obstruction that is almost invisible
in Q1. RM and random have useful Q1 scaling explanations and selected complete
results, but lack the same refined grid-wide occupation evidence.

## What the numbers mean

K is message dimension in bits, B is the rate-half constituent block length,
T is the input size of an inner epoch, and S is the inner state size. There
are L=K/(B/2) outer rows and N=2K output bits. Except where stated otherwise,
the distance target is 10%, the state starts at zero, output precedes update,
state carries across regions, and there is no final flush.

The fixed-constituent setup uses one outer code repeatedly, independent row
and region permutations, and fresh independent nonzero field multipliers.
The inner maps are fixed. Their precise generators and kernel spectra are
part of the instance; matching B,T,S does not establish matching encoders.

Q counts nonzero outer rows. If U_Q bounds that occupation's expected bad-word
count, define U_full=sum_Q U_Q and the full margin M=-log2(U_full). Markov's
inequality turns a small full first moment into a failure-probability bound.
A Q1 margin concerns only one contribution. Write R=sum_(Q>=2) U_Q/U_1;
adding the other bounds costs Delta=log2(1+R) margin bits. These ratios compare
selected upper bounds, not the ratios of true bad-word contributions.

| Evidence | What it supports |
| --- | --- |
| Complete outward certificate | The stated finite probability bound for the pinned construction and setup |
| Replayed full binary64 bound | An engineering margin with complete occupation coverage; outward rounding remains separate |
| Q1 or other partial bound | Only the named occupations; a passing margin is not a complete distance claim |
| Positive first-moment lower obstruction | A small first-moment bound is impossible there; actual code failure is not established |
| Weak full upper bound | Coverage is complete, but the selected bound is inconclusive |
| Engineering model | Conditional scaling guidance; neither a uniform theorem nor a certificate |

## Evidence inventory

| Study | Numerical scope | Current interpretation |
| --- | --- | --- |
| Earlier broad finite grid | 3,108 native tuples; 263,262 observations; four positive full diagnostics | All occupations evaluated, but most dense bounds are weak |
| Matched family surfaces | 1,292 family evaluations over 387 geometries | Q1-only comparison of exact BCH/RM, a random mean, and random spectrum caps |
| Refined BCH-64/128 audit | 130 geometries: 77 useful full bounds, seven weak full bounds, 46 first-moment obstructions | No missing geometries; 68 useful bounds have R<=0.1; 76 have R<1 |
| Fixed RM(4,9) reference | K=2^16, B=512, T=64, S=14 | Recorded outward certificate: 42.57798 margin bits |

These counts describe overlapping studies and must not be added as if they
were disjoint experiments. The RM certificate is a separate pinned
construction. It is not automatically represented in the broad database view,
and this consolidation does not rerun its expensive proof producers.

The BCH audit authenticates 602 files and reconciles aggregation within
2.54e-15 bits. Its twelve fixed-calibration model comparisons pass independent
90-digit arithmetic checks within 5.35e-15 bits. The 46 positive lower witnesses
were also replayed at 90 digits. This checks numerical consistency and
provenance; it does not measure the slack of the underlying inequalities.

## BCH: message size and the 40-bit question

For T=64,S=20, the complete margins at the three requested message sizes are:

| K | BCH-64 full margin | BCH-128 full margin | BCH-128 shortfall from 40 |
| --- | ---: | ---: | ---: |
| 2^16 | 14.27122 | 36.56453 | 3.43547 |
| 2^18 | 12.36048 | 34.72901 | 5.27099 |
| 2^20 | 10.37192 | 32.76960 | 7.23040 |

No BCH-64/128 point in the refined grid reaches 40 bits. The best sampled
BCH-128 point is K=2^13,T=64,S=20, with 38.14275 full margin bits. This is a
statement about our present bounds and tested maps, not impossibility for
BCH-128 or a proof that its actual distance guarantee fails.

Increasing S reaches a plateau because a low-weight outer word can activate
the inner state late. The row-count factor then supplies the large-K slope:
when the per-row transfer stabilizes, doubling K costs approximately one
Q1 margin bit. Small-K finite-epoch behavior does not follow that slope
uniformly, which explains why the best sampled point is not the smallest K.

There is identifiable Q1 proof slack. In the simplified persistent-state onset
model, BCH-128 has a 50.099-bit exact-probability intercept versus 46.804 bits
for its Chernoff bound, before the outer row-count factor. The roughly
3.295-bit gap is inside that idealization. It motivates a tighter distributional
calculation, but does not guarantee recovering 3.295 bits for the real encoder.
K=2^16 is therefore the closest requested BCH-128 target worth tightening;
the larger-K deficits exceed this particular modeled gap.

Higher occupations produce an additional K-dependent effect:

| Constituent | K | T,S | Full margin | Loss relative to Q1 |
| --- | --- | --- | ---: | ---: |
| BCH-64 | 2^20 | 64,20 | 10.37192 | 0.0150234 |
| BCH-64 | 2^24 | 64,20 | 6.17126 | 0.222935 |
| BCH-64 | 2^26 | 64,20 | 3.57442 | 0.820138 |
| BCH-128 | 2^20 | 64,20 | 32.76960 | 5.77436e-9 |
| BCH-128 | 2^24 | 64,20 | 28.78254 | 9.10589e-8 |
| BCH-128 | 2^26 | 64,20 | 26.78319 | 3.63971e-7 |

Writing U_Q=choose(L,Q) H_Q(L) separates row counting from transfer. If H2/H1
stabilizes, U2/U1=(L-1)H2/(2H1) grows approximately linearly with K. Its
normalized coefficient is about twenty bits smaller for BCH-128 than BCH-64
at the large-K anchors. This explains the visibly different aggregation costs.
A two-term model calibrated only at K=2^20 is checked at every integer exponent
21..26. At the largest endpoint it is optimistic by 0.07283 bits for BCH-64
and conservative by 0.01360 bits for BCH-128. No result beyond exponent 26,
or uniform interpolation between sampled points, follows from that fit.

## BCH: state size and epoch size

At K=2^20,T=64, both constituents have useful full bounds for every tested
S=10..20. Reducing S from 20 to 16 costs only 0.15638 bits for BCH-64 and
0.27469 bits for BCH-128 at this K. At K=2^24 those costs are 0.27004 and
0.27442 bits: the BCH-64 state penalty increases as higher occupations grow.

At S=10,K=2^20, BCH-64 has 6.67426 full margin bits and a 0.58438-bit
aggregation loss. BCH-128 has 27.04102 full margin bits and only 2.14e-6 bits
of aggregation loss. At BCH-64 S=10,K=2^22 the full margin is 3.57541 bits
but the loss reaches 1.69492 bits. Q1 is no longer the majority of that
selected upper bound. Thus a small-state approximation validated at one K
cannot be extended freely to larger messages.

T=128 is attractive once enough state is available. At K=2^20:

| Constituent | T | S | Full margin | Higher-occupation loss |
| --- | ---: | ---: | ---: | ---: |
| BCH-64 | 128 | 19 | 10.35330 | 0.0158728 |
| BCH-64 | 128 | 20 | 10.36538 | 0.0150961 |
| BCH-128 | 128 | 18 | 32.69995 | 6.91928e-9 |
| BCH-128 | 128 | 19 | 32.73912 | 6.19266e-9 |
| BCH-128 | 128 | 20 | 32.75921 | 5.81929e-9 |

At S=20, changing T=64 to T=128 costs only 0.00654 margin bits for BCH-64
and 0.01039 for BCH-128. But the supported small-state range changes sharply:

| T at K=2^20 | First-moment obstructions, both blocks | Weak bounds | Useful full bounds |
| ---: | --- | --- | --- |
| 64 | Tested S=7..8 | S=9, both blocks | S=10..20, both blocks |
| 128 | Tested S=8..16 | S=17, both; also BCH-64 S=18 | BCH-64 S=19..20; BCH-128 S=18..20 |
| 256 | Every tested S=9..20 | None needed to establish the obstruction | None in this slice |

Multi-bit inputs in the syndrome kernel can keep the state zero, leaving
output equal to input. Collective counts of these trajectories can obstruct
a small first moment even while Q1 barely changes with T. For example,
BCH-128 T256/S20 has a lower exponent near 102,883.94 on its expected bad-word
count despite a Q1 margin near 32.75 bits. This is why the T tradeoff requires
all-occupation evidence. Fewer epochs means fewer state updates, but total
implementation cost also depends on the maps and arithmetic; no speed ranking
is established by these bounds.

The refined full T=128 comparison above is at K=2^20. Full checks at K=2^16
and 2^18 have not been supplied by this BCH grid. No current full bound in
this exact-spectrum BCH collection reaches 40 bits at T=128. Changing T
alone has not shown a meaningful Q1 gain for BCH-128; the smaller-K Q1
calculation remains a more relevant tightening target.

## The fixed-RM finite reference certificate

The [fixed RM certificate](../small_k_replay/RM2SUB_RM49_FINITE_CERTIFICATE.md)
uses RM(4,9) [512,256,32], K=2^16, T64/S14, and its pinned inner map. Its
outward bound gives at least 42.57798 margin bits over all Q=1..256. This is
a complete 40-bit reference and should not be confused with historical RM
entries still marked under review in a different database ledger.

This statement concerns its fresh-independent-multiplier setup. The
consolidation makes no PRG-replacement, decoder, or end-to-end implementation
claim and performs no fresh replay of its certificate producers.

## RM and random: the same coordinates, different constituent trends

At K=2^20,T=64,S=20, matched Q1 references are:

| Family | B | Q1 margin | Dominant weight |
| --- | ---: | ---: | ---: |
| BCH | 64 | 10.387 | 12 |
| BCH | 128 | 32.770 | 22 |
| RM | 128 | 10.619 | 16 |
| RM | 512 | 39.201 | 32 |
| Random mean | 128 | 21.196 | 13 |
| Random mean | 256 | 60.087 | 26 |
| Random mean | 512 | 137.351 | 52 |

These are partial bounds at matched maps. RM-512 is overwhelmingly dominated
by its weight-32 shell at this anchor; the random mean is spread across many
shells. The largest random-512 shell supplies only about 7.8% of its Q1 bound.
The minimum-shell onset explanation predicts a square-root leading block-size
scale along the rate-half RM family, whereas the random mean onset contribution
is approximately 0.3B at 10% distance. These are explanatory model statements;
finite transfer effects and higher occupations remain necessary for full bounds.
The two exact BCH constituents in the plateau regime do not determine a robust
asymptotic block-size law. This account keeps its BCH conclusions within the
exact-spectrum range.

For RM-512, the T64/S20 Q1 margins at K=2^16,2^18,2^20 are 43.04255,
41.16950,39.20062 bits. Increasing S from 16 to 20 at K=2^20 buys only about
0.116 bits. At K=2^18 the Q1 value motivates a full-bound investigation; it
is not an existing 40-bit certificate. At K=2^20 the current Q1 bound itself
falls below 40. The separate S14 certificate establishes one smaller-K anchor,
not the full RM engineering surface.

The random mean and the conditional random caps must stay separate:

| Random B, at K=2^20,T64/S20 | Ensemble Q1 margin | Conditional Q1 margin | Q1 plus spectrum-event failure |
| --- | ---: | ---: | ---: |
| 256 | 60.08666 | 7.75122 | 7.75122 |
| 512 | 137.35061 | 77.57586 | 59.99999 |

The mean averages one uniformly sampled constituent reused across rows. The
caps bound a good-spectrum event and charge its 2^-60 failure budget once.
Neither is the measured spectrum of a chosen random code. The earlier broad
ledger does have full random-[512,256] diagnostics at K=2^16,T64,S=12,18,20,
all approximately 60 bits after that shared charge. A displayed 60.0 is rounded;
it does not certify a strict greater-than-60-bit result. Outward replay remains
separate. This is useful full evidence at three points, not an audited random
surface comparable to the BCH study.

The repository also contains an [asymptotic random-outer theorem](../RANDOM_OUTER_RM2SUB_CERTIFICATE.md)
with fixed T128/S19, B approximately 9 log2 N, and relative-distance target
11%. Its outer injections are sampled independently between rows. That is a
different ensemble from the reused-constituent surface above, and its o(1)
failure statement does not supply a finite 40-bit parameter choice.

## Planning around 40 bits within the current scope

| K | Established finite reference | Other promising work |
| --- | --- | --- |
| 2^16 | RM-512 T64/S14, outward margin 42.57798 | Tighten BCH-128 Q1; outward-check the random-512 T64/S12 full diagnostic |
| 2^18 | No complete 40-bit reference collected here at this exact K | Investigate the full RM-512 bound beyond its 41.16950-bit Q1 screen; tighten BCH-128 Q1 |
| 2^20 | No complete 40-bit reference collected here at this exact K | Investigate full random-512 bounds; current BCH-128 and RM-512 Q1 bounds are below 40 |

For specifically T=128, no full 40-bit exact-BCH result is established here.
The useful BCH-64/128 full comparisons at K=2^20 use S=19..20 and S=18..20,
respectively, but their margins remain below 40. For smaller K, a passing Q1
screen would still need complete occupation evidence. RM and random remain
separate candidate families; their Q1 surfaces cannot settle this question.

The next priorities are tighter BCH-128 Q1 evaluation at K=2^16, the seven
weak BCH grid points below, and selected higher-occupation checks for RM and
random before extending their engineering conclusions. No new search is part
of this consolidation.

Seven BCH-64/128 points remain weak: T64/S9 for both blocks at K=2^20;
BCH-64 T128/S17 and S18; BCH-128 T128/S17; and BCH-64 T64/S10 and S12 at
K=2^24. The [detailed audit](BCH_DOMINANCE_ANALYSIS.md) identifies their
problematic intervals. At the last S12 point, the full deficit is only
1.644 bits and Q5..256 is responsible. At S10, even the selected Q1..4
union is weak. A dense-only refinement cannot resolve that case.

## Where to find the evidence

- [Numerical appendix](CURRENT_RESULTS_TABLES.md): all 130 BCH geometries,
  all matched family references at the three requested K values, the earlier
  positive full diagnostics, and input fingerprints.
- [BCH occupation analysis](BCH_DOMINANCE_ANALYSIS.md): complete-bound method,
  cancellation obstructions, Q1 ratios, model checks, and exact weak-point list.
- [BCH growth explanation](BCH_GROWTH_ANALYSIS.md): row counting, onset model,
  and the distinction between transfer effects and constituent spectrum effects.
- [Family surfaces](CONSTITUENT_ENGINEERING_SURFACES.md): RM/random trends,
  model assumptions, and primary spectrum coverage.
- [Broad grid findings](GRID_FINDINGS.md): the earlier 3,108-tuple ledger,
  including weak complete covers and four positive diagnostics.

Local generated figures include bch_full_bound_k_scaling_v8.png,
bch_full_bound_k_state_scaling_v8.png, bch_full_bound_state_scaling_v8.png,
bch_full_bound_coverage_v8.png, and bch_q1_vs_dense_tradeoff_v8.png. The family
figures include constituent_engineering_comparison.png and
engineering_surface_slices.png. They remain local and are not published in Git.

To check this numerical snapshot from this directory, with its existing local
research inputs, run the following sequentially:

~~~powershell
python audit_bch_evidence_grid_v1.py
python collect_current_results_v1.py --check
~~~

The collector checks input fingerprints and exact coverage of the 130-point
BCH table. It reads existing evidence and performs no parameter search. To
refresh the appendix intentionally, run it without --check and review the
resulting documentation diff. A source-only checkout must first reproduce
or restore the ignored research inputs. Producers and proof receipts used by
existing results remain frozen; this consolidation changes no numerical bound.
