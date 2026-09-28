# Complete support coverage through 58 active groups

The independent-row construction now has complete support coverage for
every occupancy q = 1,...,58. Their combined bad-event contribution is
below 6.786362e-14, giving more than 43.744354 bits of margin for this
restricted message class. Occupancies 59 through 2048 are not covered by
these runs, so this is not a full-code minimum-distance certificate.

The construction is BCH[256,128] with K = 2^20, N = 2^21, four rows per
group, independent row shuffles, and two-update IMT(128,19). An occupancy
q means exactly q of the 2048 groups have a nonzero message. The bad event
is an output of weight at most 209715 from any message in that class.
The probability is over ideal setup sampled once and fixed for all
messages. [BRIDGE.md](BRIDGE.md) gives the end-to-end counting argument.

Each occupancy result includes every support vector in {38,...,256}^q,
every ordered message tuple, and all binomial(2048,q) group locations.
All-one penalties are coupled to their reciprocal outer weights. These
runs use input-weight tilt one and the CDF cover, not the later
shell/prefix-rank cover. The original runs through occupancy 35 use the
base inner envelope. The extension to 36 through 40 additionally uses
full-feedback, window-histogram, and joint-cancellation refinements. The
hill-climb replay below further tightens occupancy 35 with density bounds.

## One active group

The one-group run uses 384-bit outward arithmetic and a lower tilt grid.
Its complete cover has eight retained leaves after 25 splits. The
computed outward upper begins 6.782959315129205655e-14, with reported margin
43.7450784890161848 bits. We use the conservative decimal upper

    U_1 = 6.782960e-14.

Reproduce from the repository root:

```sh
python -B research/workstreams/permutation_locality/independent_rows/verify.py --groups 1 --precision 384 --tilts .00016 .00025 .00032 .0004 .0005 .00064 --penalties 1 --max-splits 200 --target-bits 39
```

The requested 39-bit search target is a stopping parameter, not the
resulting claim. The returned outward bound gives the stronger margin
listed above.

## Occupancies two through thirty-two

The batch uses 192-bit outward arithmetic, tilts .001, .0032, .008, .016,
and .024, and penalties .75 and 1. Every requested occupancy returns a
complete support cover and an outward certificate. No interpolation in
q is used.

```powershell
$proofOccupancies = 2..32
python -B research/workstreams/permutation_locality/independent_rows/verify.py --occupancies @proofOccupancies --tilts .001 .0032 .008 .016 .024 --penalties .75 1 --precision 192 --target-bits 52 --max-splits 100
```

The table rounds each reported margin down to six decimal places. For
each row q, the constructed bad-event upper is strictly below 2^(-m_q),
where m_q is the displayed number. The q = 1 row comes from the separate
384-bit run; rows 2 through 32 come from this batch. Rows 33 through 35
come from the extensions below.

| Active groups q | Conservative margin m_q, bits |
|---|---:|
| 1 | 43.745078 |
| 2 | 57.971149 |
| 3 | 144.747640 |
| 4 | 191.755087 |
| 5 | 54.867014 |
| 6 | 172.146445 |
| 7 | 269.796734 |
| 8 | 351.046767 |
| 9 | 421.644605 |
| 10 | 484.111527 |
| 11 | 539.972470 |
| 12 | 590.214974 |
| 13 | 635.517315 |
| 14 | 676.226097 |
| 15 | 713.203288 |
| 16 | 746.252858 |
| 17 | 775.786131 |
| 18 | 801.978492 |
| 19 | 824.702415 |
| 20 | 838.818195 |
| 21 | 810.434479 |
| 22 | 767.325970 |
| 23 | 677.031042 |
| 24 | 610.644982 |
| 25 | 548.679912 |
| 26 | 534.767718 |
| 27 | 526.417752 |
| 28 | 497.829334 |
| 29 | 441.972777 |
| 30 | 374.936521 |
| 31 | 348.285581 |
| 32 | 282.426186 |
| 33 | 166.347539 |
| 34 | 80.464846 |
| 35 | 43.538276 |

These are margins for upper bounds, not estimates of the true failure
probabilities. Their variation reflects the chosen witnesses and stopping
rule as well as the occupancy. In particular, the q = 5 search stops
without splitting its initial support box once it meets its budget.

The batch's outward sum of certificate uppers is enclosed by

    [3.39751895646005887657879428449780047276698572150323766481e-17
       +/- 7.48e-75].

The corresponding reported margin is enclosed by

    [54.7082960125517650330240653625605371379731140086396180539
       +/- 4.70e-56].

These enclosures concern the computed upper bound and its logarithm.
They do not give a lower bound on the actual failure probability.

## Occupancies thirty-three and thirty-four

The extension retains the smaller tilts needed for sparse supports and
adds tilt .032. At 192-bit precision, both occupancies produce complete
outward support covers after 25 splits. Their outward uppers are below
8.402353e-51 and 5.993321e-25, respectively. The covers have 514 and 533
retained leaves.

```sh
python -B research/workstreams/permutation_locality/independent_rows/verify.py --occupancies 33 34 --tilts .001 .0032 .008 .016 .024 .032 --penalties .75 1 --precision 192 --target-bits 52 --max-splits 100
```

The original extension request also included occupancies 35 and 36. Its
q = 35 screening result did not meet the requested 52-bit budget. The
separate replay below now certifies that occupancy with a different
allocation. At q = 36, the original run ended with binary64 log2 upper
+5.909 after 100 splits and was not replayed. The refined run below now
certifies that occupancy.

## Occupancy thirty-five

The same witness grid gives a complete outward cover at q = 35 after
100 splits, with 1745 retained leaves. At 192-bit precision, the computed
upper and its margin are enclosed by

    [7.82839396033940951908298352465338472248267975268279202695e-14
       +/- 1.22e-71],
    [43.5382769676644893569953018780502154822241555955477930435
       +/- 5.95e-56] bits.

We use the conservative decimal upper

    U_35 = 7.828394e-14 < 2^-43.

```sh
python -B research/workstreams/permutation_locality/independent_rows/verify.py --groups 35 --tilts .001 .0032 .008 .016 .024 .032 --penalties .75 1 --precision 192 --target-bits 41 --max-splits 100
```

The command's 41-bit target is a search parameter. The verifier requires
a binary64 proposal two bits stronger before starting outward replay, so
this command uses a 43-bit proposal gate. A requested target of 43 would
instead require a 45-bit proposal and did not replay this cover. The
reported outward upper was independently checked to be below 2^-43;
the certificate does not rely on the binary64 score. No verifier
arithmetic was changed for this run.

We allocate 2^-43 to occupancy 35 rather than requiring the earlier
uniform 2^-52 budget. The aggregate calculation below shows that this
larger allocation remains compatible with a 40-bit full-proof target.

## Occupancies thirty-six through forty

The refined envelope gives full support covers for all five occupancies
with 192-bit outward arithmetic. Each cover uses 25 splits, includes
every support vector and active-group location, and meets a 52-bit
contribution budget. The table rounds the verified margins downward.

| Active groups q | Margin, bits | Retained leaves |
|---|---:|---:|
| 36 | 414.265557 | 443 |
| 37 | 372.215581 | 446 |
| 38 | 329.327462 | 495 |
| 39 | 303.436167 | 555 |
| 40 | 266.547892 | 600 |

The outward sum for these five occupancies is below 5.768846e-81.
These runs do not use the newer density-column or row-conditioning
refinements; they preserve the existing unrestricted CDF cover.

```sh
python -B research/workstreams/permutation_locality/independent_rows/verify.py --occupancies 36 37 38 39 40 --tilts .016 .024 .032 .04 --penalties .75 1 --precision 192 --target-bits 52 --max-splits 100 --full-feedback 6 --window-histogram 8 --joint-cancellation
```

The run also retained eight exact-dyadic operator memo files in
`tmp/independent-row-q36-40-operators`. The reproduction command above
does not require them and rebuilds its operators from the checked inputs.

## Historical aggregate before the density refinements

Summing 2^(-m_q) over the rounded table for q = 2,...,32, with 384-bit
outward arithmetic, gives an upper beginning

    3.39752076065073062983814924178424862039815198011881438279e-17.

Adding the conservative U_1 above gives an upper beginning

    6.78635752076065073062983814924178424862039815198011881438e-14,

whose margin is greater than 43.744355 bits. An upper for that value is
6.786357521e-14. Adding the outward decimal uppers for q = 33 and 34
still gives less than 6.786358e-14 for occupancies 1 through 34.

Adding the conservative q = 35 upper gives the exact decimal sum

    6.786358e-14 + 7.828394e-14 = 1.4614752e-13.

An independent 384-bit outward calculation gives a margin greater than
42.637639885500592382 bits for this conservative sum.

Adding the q = 36,...,40 contribution and rounding upward gives

    U_(1..40) < 1.4614753e-13 < 2^(-42.637639).

The decimal comparison was checked with exact rational arithmetic; an
independent 384-bit calculation gives margin 42.6376397867856... bits.

There are 2008 remaining occupancies, q = 41,...,2048. Requiring an upper
below 2^-52 for each would still suffice for the full 40-bit target:

    1.4614753e-13 + 2008 * 2^-52
      = 5.920130966894628666341304779052734375e-13
      < 2^-40.

The margin of this prospective aggregate is more than 40.619436 bits.
This budget calculation is not a certificate for the missing occupancies.
The next step is to cover q = 41,...,2048, either individually or through
rigorous range bounds, while keeping their outward aggregate within budget.

The shared-row-shuffle route and its earlier certificates are unchanged.

## Hill-climb replay: occupancy thirty-five

The complete q = 35 cover now uses the stronger inner refinements and
both sparse-side and dense-side tilts. It passes 192-bit outward replay
with 25 splits and 410 retained leaves:

    U_35 < 9.617964e-131,
    -log2(U_35) > 431.906849 bits.

This covers every support vector, every compatible outer word, and every
active-group location. It replaces the earlier 7.828394e-14 upper for
this occupancy. The improvement combines earlier inner refinements with
the new density bounds; it is not attributed to the density bounds alone.

For tilts .016 and .024, the replay uses penalties .75 and 1, full feedback
through six windows, output histograms through eight windows, joint
cancellation, and conditioned-column density. For tilts .032, .04, and
.048, it additionally uses the exact six-window feedback-density bound,
with penalties .75, .9, and 1. All region operators have degree 64; the
cover uses their exact prefix through its own occupancy. The prefix
placement identity is checked separately by the verifier.

Replacing only q = 35 in the preceding aggregate gives

    U_(1..40) < 6.786359e-14,
    -log2(U_(1..40)) > 43.744355 bits.

The decimal comparison was checked with exact rational arithmetic.
At 384-bit precision, the conservative upper's margin is enclosed by
43.744355577445772825... bits. The remaining 2008 occupancies would
still fit a 52-bit allocation apiece, leaving a prospective aggregate
below 5.137292e-13. This remains a budget, not proof of those occupancies.

The following Python driver reproduces the two-grid calculation from the
repository root. It deliberately omits local operator memoization. Each
grid bounds the same construction; combining their witnesses does not
combine different encoders or setup distributions.

```python
import sys
from types import SimpleNamespace
sys.path.insert(0, 'research/workstreams/permutation_locality/independent_rows')
from verify import (build_operators, authenticated_caps, weighted_cdf_upper,
                    integer_cdf, Q, cover, TAIL_TERMINAL)

args = SimpleNamespace(
    groups=64, precision=192, operator_cache=None,
    tilts=['.016', '.024'], penalties=['.75', '1'],
    full_feedback=6, window_histogram=8, joint_cancellation=True,
    column_density=True, feedback_density=0, weight_tilt='1',
    probe_supports=[], probe_vector=[], retain_parents=True,
    joint_witness=True, joint_top=2, target_bits=52, max_splits=100,
    screen_only=False)
operators = build_operators(args)
args.tilts = ['.032', '.04', '.048']
args.penalties = ['.75', '.9', '1']
args.feedback_density = 6
operators.update(build_operators(args))
caps = authenticated_caps()
counts = {rho: integer_cdf(weighted_cdf_upper(caps, 1 << 128,
                                             full_weight=1/Q(rho)))
          for rho in args.penalties}
for q in [35, 40, 41, 42, 43, 44, 45, 46, 47, 48]:
    args.groups = q
    cover(args, operators, counts, TAIL_TERMINAL)
```

Each call has its own complete support cover and reports whether its
outward replay succeeds. The requested occupancy list alone is not a
claim that all its entries passed.

## Complete extension through occupancy forty-eight

Every occupancy from 41 through 48 passes the two-grid driver above.
Each cover uses 25 splits and includes all support vectors, compatible
outer words, and active-group locations. The margins below are rounded
down from the 192-bit outward replay.

| Active groups q | Margin, bits | Retained leaves |
|---|---:|---:|
| 41 | 313.183280 | 583 |
| 42 | 290.478863 | 561 |
| 43 | 257.037183 | 506 |
| 44 | 220.409541 | 521 |
| 45 | 180.969960 | 525 |
| 46 | 142.732415 | 567 |
| 47 | 112.625029 | 606 |
| 48 | 78.281562 | 675 |

Their outward sum is below 2.722085e-24, with margin greater than
78.281562 bits. The same run independently replays q = 40 with 610
leaves and upper below 4.473997e-99, improving its margin to more than
326.709317 bits. The aggregate below conservatively retains the older
q = 36,...,40 sum instead of relying on that improvement.

Combining the earlier q = 1,...,34 upper, the improved q = 35 upper,
the earlier q = 36,...,40 sum, and this new batch gives

    U_(1..48) < 6.786358e-14 + 9.617964e-131
                + 5.768846e-81 + 2.722085e-24
              < 6.786359e-14.

Both comparisons were checked using exact rational arithmetic. The
conservative upper has margin greater than 43.744355 bits. There are
2000 remaining occupancies. Allocating 2^-52 to each would give

    6.786359e-14 + 2000 * 2^-52 < 2^-40.

That last inequality is only the remaining budget. It does not certify
any occupancy above 48; those require their own complete covers or
rigorous range bounds.

## Occupancy forty-nine

Adding tilts .052, .056, and .06 at rho = .9, with both density
refinements, gives a complete q = 49 cover. Its 192-bit outward replay
uses 50 splits and 1429 retained leaves:

    U_49 < 1.388120e-43,
    -log2(U_49) > 142.369776 bits.

To reproduce it, extend the Python driver above with these statements:

```python
args.groups = 64
args.tilts = ['.052', '.056', '.06']
args.penalties = ['.9']
operators.update(build_operators(args))
args.groups = 49
cover(args, operators, counts, TAIL_TERMINAL)
```

Adding U_49 to the preceding unrounded decimal aggregate still gives
U_(1..49) < 6.786359e-14, with more than 43.744355 bits of margin.
The larger margin at q = 49 than at q = 48 reflects different witness
grids and cover refinement depths; it is not a claim that the code's
actual failure probability decreases between those occupancies.

## Extension through occupancy fifty-six

The same sixteen-witness grid passes the next seven complete CDF covers.
Each result includes every compatible support vector, outer word, and
active-group location. The table rounds uppers upward and margins downward
from the 192-bit outward replay.

| Active groups q | Bad-event upper | Margin, bits | Splits | Retained leaves |
|---|---:|---:|---:|---:|
| 50 | 4.705469e-31 | 100.745432 | 50 | 1475 |
| 51 | 5.042458e-22 | 70.748291 | 50 | 1553 |
| 52 | 9.114171e-44 | 142.976724 | 75 | 2557 |
| 53 | 1.229410e-36 | 119.291446 | 75 | 2590 |
| 54 | 1.876797e-30 | 98.749570 | 75 | 2632 |
| 55 | 3.680288e-24 | 77.846455 | 75 | 2656 |
| 56 | 3.300553e-42 | 137.798272 | 100 | 3843 |

To reproduce these results, use the preceding extended driver with
`args.groups` set successively to 50 through 56. Keep `max_splits=100`
and the full sixteen-witness operator set. Adding these seven decimal
uppers to the preceding unrounded aggregate still gives

    U_(1..56) < 6.786359e-14.

The outward sum for q = 49,...,56 is below 5.079261e-22, with margin
greater than 70.737799 bits. The complete aggregate comparison was checked
with exact rationals. There are 1992
remaining occupancies; allocating 2^-52 to each still fits the full
2^-40 budget. Neither this budget calculation nor the selected dense
support-class results cover those remaining occupancies.

## Extension through occupancy fifty-eight

The same sixteen-witness grid passes complete CDF covers for q = 57
and 58. Both use 100 splits and 192-bit outward replay. Their bounds
include every support vector, compatible outer word, and active-group
location. Decimal uppers are rounded upward; margins are rounded down.

| Active groups q | Bad-event upper | Margin, bits | Retained leaves |
|---|---:|---:|---:|
| 57 | 1.164709e-32 | 106.081729 | 3893 |
| 58 | 2.794438e-20 | 64.956003 | 3941 |

Reproduce with the preceding extended driver, setting `args.groups` to
57 and then 58 and `max_splits=175`. The search stops at 100 splits
once its proposal passes the gate for outward replay.

Using the conservative unrounded aggregate components above gives

    U_(1..58) < 6.786358e-14 + 9.617964e-131 + 5.768846e-81
                + 2.722085e-24 + 5.079261e-22
                + 1.164709e-32 + 2.794438e-20
              < 6.786362e-14.

The second comparison was checked with exact rationals. At 384-bit
precision, the last upper has margin 43.7443549396834... bits.
Allocating 2^-52 to each of the 1990 remaining occupancies still gives
less than 2^-40 in total. This is a budget check only; q = 59,...,2048
remain unproved.

The q = 59 attempt with the same grid ended at 175 splits and 6449
leaves. Its binary64 log2 upper was -48.463562, above the -54 proposal
gate for a 52-bit target. The driver therefore did not replay it with
outward arithmetic. That score adds no certificate to the aggregate.
The q = 60 attempt also ended without outward replay: 175 splits,
6677 leaves, and binary64 log2 upper -15.829523. Stronger local
refinements are being tested separately before repeating these covers.
