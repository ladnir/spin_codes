# Complete two-bit covers and the dense gap

The active objective is a full K=2^20 proof for the two-bit ensemble in
[README.md](README.md). The original target was 10% distance and aggregate
setup-failure probability below 2^-40. The first complete proof may use
a lower distance, a weaker margin, or adjusted inner parameters. Establish
one complete operating point before optimizing its guarantees.
The computations below have not changed the construction or its randomness.

**First closure obtained and improved:** the entire occupancy range now
closes through 9.25% distance with margin above 49.11 bits. See
[FIRST_CLOSURE.md](FIRST_CLOSURE.md) for the whole-code statement and
verification scope. The unresolved results below concern higher distances.

## Verified coverage

At the original 10% target, complete support covers now close
**every occupancy q=1,...,400**.
Here q counts active pairs of BCH rows; either one or both rows in a pair
can be nonzero. Each cover includes unequal support sizes, all nonzero
pair messages, and all binomial(4096,q) possible sets of active pairs.

The q=1 bound is the independently replayed 49.11287738250-bit bound from
the first study. The outward sum for q=2,...,32 is at most

    4.5171974135816698842599845084094819236552611219596324318e-29,

or more than 94.16048678695 bits. The additional q=33,...,64 sum is at most
1.4427910980635686894753287614370816879986095029831553822e-170,
or more than 564.19891370392 bits. The q=65,...,128 sum is at most
2.2682358966865575500112642043722153103766594852269342476e-108,
or more than 357.58666355937 bits. The q=129,...,168 sum is at most
4.24665451017929054540387112105021711997561435276509589593e-22,
or more than 70.99609134561 bits. An additional complete-cover batch
certified q in {169,176,192,256}, with total at most
4.52259019369536235669788420327178398013796945683408024130e-50,
or more than 163.91925546966 bits. The next completed batch certified
q in {170,200,224,272,288,320}, with total at most
4.46441832510727888000770175218229095099873035762571401440e-64,
or more than 210.44492585582 bits. Finally, the complete q=171,...,320
batch has aggregate at most
4.02297732405715215987657330260697225712006296491508682850e-27,
or more than 87.68379495587 bits. Hence the union of all verified
occupancies is below
2^-49.11. This is a restricted-class result, not a full-code certificate.

Selected individual margins, rounded downward:

| q | Complete-cover margin |
|---:|---:|
| 1 | 49.11287738250 |
| 2 | 94.16048678695 |
| 8 | 392.47436262634 |
| 16 | 550.23276981365 |
| 24 | 580.61214420342 |
| 32 | 1013.83702843384 |
| 37 | 564.19927210731 |
| 48 | 955.79314453594 |
| 64 | 958.17237870067 |
| 65 | 964.05450536929 |
| 128 | 357.58666356269 |
| 168 | 70.99609160949 |
| 169 | 163.91925546966 |
| 170 | 210.44492585582 |
| 176 | 2461.07712292270 |
| 192 | 3215.82475702593 |
| 200 | 3407.46111818353 |
| 224 | 3154.65295381186 |
| 256 | 2782.28092494307 |
| 272 | 3504.13970961511 |
| 288 | 3760.35297010680 |
| 320 | 3647.68550877863 |
| 321 | 3623.11986125303 |
| 352 | 2965.11328563007 |
| 384 | 1367.81571146121 |
| 400 | 348.15607948548 |

Every q=2,...,170 passed using a single box {38,...,256}^q. Covers through
400 require at most one symmetric split, giving at most q+1 leaves.
Occupancies above 400 remain open. The nonmonotone
margins reflect a finite tilt
grid and upper-bound slack; they are not estimates of the true distance.

`full_cover.py` reconstructs the local and regional operators. `cover.py`
selects rational support witnesses with binary64 proposals, then replays
the entire cover outward. It reuses tested interval folds and symmetric
box splitting without modifying the four-bit drivers. The location factor
is explicitly binomial(4096,q). Failed numerical covers return no bound.

```powershell
$twoBitOccupancies = 2..32
python -B research/workstreams/permutation_locality/two_bit/full_cover.py --occupancies $twoBitOccupancies --tilts .00064 .001 .0016 .0025 .004 .0064 .01 .016 .025 --max-splits 75
$twoBitOccupancies = 33..64
python -B research/workstreams/permutation_locality/two_bit/full_cover.py --occupancies $twoBitOccupancies --tilts .01 .016 .02 .025 .032 .04 .048 .064 --max-splits 40
$twoBitOccupancies = 65..128
python -B research/workstreams/permutation_locality/two_bit/full_cover.py --occupancies $twoBitOccupancies --tilts .016 .02 .025 .028 .032 .036 .04 .048 .056 .064 .08 .096 --max-splits 0 --certificate-dir tmp/two-bit-covers
```

The optional `--certificate-dir` saves the rational witnesses and support
partition as small JSON files. `--replay` reconstructs the operators and
checks the exact partition and label multiplicities before recomputing
every bound; stored upper bounds are never trusted. Generated witnesses
belong outside version control.

The save/replay CLI was checked end to end at q=2 with a degree-two
census: it independently reconstructed the 94.16048678695-bit bound.
The q=65,...,128 witnesses were also saved. An initial tilt grid beginning
at .04 missed q=65 despite its large actual certified margin; subdivision
did not repair that witness-choice problem. Extending the grid downward
closed the whole batch without a single support split.

The newer [profile-measure bound](PROFILE_MEASURE.md) improves the counting
interface and supplies a complementary dense branch. It is not needed
for the complete q<=128 results above. At q=4096 it certifies the full box
{160,...,192}^4096, including unequal supports, with over 53578 bits of
margin. That restricted box must not be reported as a complete q=4096
certificate. Mixed light/central/heavy supports remain open.

The density-fold save/replay smoke test independently reproduced
94.20498869799 bits at q=2. The completed q=129,...,256 batch certified
q=129,...,168 using this fold. The remaining occupancies in that batch
did not meet its 55-bit budget and remain open. For example, the q=169
proposal was -46.28338 in log2 units, whereas q=170 was -19.2843.
These proposals were not outward certificates. Reproduce the batch with:

```powershell
$twoBitOccupancies = 129..256
python -B research/workstreams/permutation_locality/two_bit/full_cover.py --occupancies $twoBitOccupancies --tilts .04 .044 .048 .052 .056 .064 .072 .08 .096 .112 .128 --max-splits 0 --density-fold --certificate-dir tmp/two-bit-density-covers
```

The stronger output-moment calculation then closed q=169 and the three
additional occupancies listed above. It evaluated exact output polynomials
through all 64 local packet occupancies. These witnesses were saved under
`tmp/two-bit-bridge-covers`. Reproduce their run with:

```sh
python -B research/workstreams/permutation_locality/two_bit/full_cover.py --occupancies 169 176 192 256 512 --tilts .064 .072 .096 .128 .192 --max-splits 1 --density-fold --output-degree 64 --certificate-dir tmp/two-bit-bridge-covers
```

The q=512 cover did not pass. The first bridge witnesses predate storage
of the output-degree setting; pass `--output-degree 64` when replaying
those files. New witnesses record this setting. The subsequent
conditional-subset refinement can make reruns stronger than this table.

With that refinement, the unsplit q=169 margin improves to
235.54684763874 bits. The next split-cover batch additionally certified
q in {170,200,224,272,288,320}; q=384 remained unresolved at binary64
log2 upper +260.04647. Its witnesses include their operator settings:

```sh
python -B research/workstreams/permutation_locality/two_bit/full_cover.py --occupancies 170 200 224 272 288 320 384 --tilts .064 .08 .096 .128 .16 --max-splits 1 --density-fold --output-degree 64 --certificate-dir tmp/two-bit-bridge-subset-split
```

The complete q=171,...,320 batch finished with all 150 occupancies
outward-verified. Its witnesses are saved under the directory below:

```powershell
$twoBitOccupancies = 171..320
python -B research/workstreams/permutation_locality/two_bit/full_cover.py --occupancies $twoBitOccupancies --tilts .064 .08 .096 .128 .16 --max-splits 1 --density-fold --output-degree 64 --certificate-dir tmp/two-bit-bridge-complete
```

A finer tilt grid then certified q=321,352,384,400. Their outward sum is
at most 1.56522545481862163048716512008292130743158109397876464830e-105,
or more than 348.15607948548 bits. Adding these occupancies preserves the
2^-49.11 bound on all verified classes. The same run left q=448 unresolved
at binary64 log2 upper +6008.56764:

```sh
python -B research/workstreams/permutation_locality/two_bit/full_cover.py --occupancies 321 352 384 400 448 --tilts .088 .096 .104 .112 .12 --max-splits 1 --density-fold --output-degree 64 --certificate-dir tmp/two-bit-bridge-retuned
```

The earlier q=321,...,384 grid started at .128 and failed throughout.
This was a witness-search limitation: the successful q=321 cover uses
.096, .104, and .112. Larger tilts do not subsume smaller ones.
A complete q=322,...,400 batch filled every intervening occupancy using
.096, .104, .112, and .12. All 79 occupancies passed outward verification.
Their aggregate is at most
1.56522545481862163054107294480565428012263387995753002568e-105,
giving more than 348.15607948548 bits. Its witnesses are under
`tmp/two-bit-bridge-322-400`. Together with q=321, this completes the
contiguous sparse range through 400 and preserves its 2^-49.11 bound.

```powershell
$twoBitOccupancies = 322..400
python -B research/workstreams/permutation_locality/two_bit/full_cover.py --occupancies $twoBitOccupancies --tilts .096 .104 .112 .12 --max-splits 1 --density-fold --output-degree 64 --certificate-dir tmp/two-bit-bridge-322-400
```

The [row-mixture method](MIXTURE_BOUND.md) now combines light, central,
and heavy classes without enumerating their full count simplex. Its
two-dimensional cover is a separate dense proof branch. Complete closure
still requires the intermediate occupancies and an outward aggregate.
The [lifted refinement](LIFTED_MIXTURE.md) retains the central-component
fraction and contracts cells using integer component counts. Full-domain
attempts at q_min=1024 and q_min=512 are not yet certificates. The latter
uses theta=2/5; its light-only endpoint at q=512 has an outward margin
above 1209 bits, but that endpoint does not bound the remaining mixtures.
The q_min=1024 run stopped at its 10,000-cell budget with 15 pending
subtrees, not a mathematical failure. Resuming the q_min=512 search
reverified its earlier leaves but still required very fine cells. The
[affine moment refinement](AFFINE_MOMENT.md) now keeps the outer count
coupled to the inner moment. It repaired one whole rectangle from log2
upper +333.69829 to below -293.04967. Its full-domain attempt remains open.

## First closure with adjustable parameters

The 10% sparse result is reusable: for any smaller threshold d, the event
that a nonzero message has output weight at most d is a subset of the
already bounded event. Thus changing only the claimed distance costs no
encoding time and requires no new sparse proof.

The dense drivers accept `--threshold`, the inclusive bad-output-weight
cutoff. New records store it explicitly; records without that field mean
209715. Replay inherits the recorded cutoff and rejects an override that
changes it. Resume may lower the cutoff and recomputes all accepted bounds.
It may not silently raise a previously verified cutoff.

The completed affine cover handles every q>=401 at d=167772, the floor of
0.08*N. It uses the same BCH outer, two-bit route, and IMT(128,19) inner
with two updates. All 540 terminal leaves passed or were proved empty;
the dense aggregate has more than 86.79627 bits of margin. Independent
256-bit replay and aggregation with the sparse lemmas confirm a full-code
margin above 49.11 bits.

```sh
python -B research/workstreams/permutation_locality/two_bit/affine_mixture.py --minimum-groups 401 --threshold 167772 --max-cells 2000 --max-depth 48 --output tmp/two-bit-affine-dense-401-d08.json
python -B research/workstreams/permutation_locality/two_bit/affine_mixture.py --replay tmp/two-bit-affine-dense-401-d08.json
```

After this first closure, increase the target distance and compare
additional mixing or larger state sizes when needed.
The original parameters are candidates, not constraints on this search.

That pivot produced complete dense covers at 5% and 6% without any inner
change. The 5% run needed 205 visited cells (103 leaves); the 6% run needed
399 cells (200 leaves). Their dense margins exceed 312.52 and 294.49 bits,
respectively. These cover all q>=401 and join the preceding sparse lemmas.
`closure.py` checks the occupancy coverage and recomputes the aggregate
while replaying the dense component; its sparse inputs are the documented
prior replay bounds, not freshly rerun sparse calculations.

The 7% lifted cover also completed, with 1451 visited cells, 726 leaves,
and dense margin above 80.04077802883 bits. The stronger affine method
closed 8% in 1079 visited cells. Both passed independent 256-bit replay;
their combined margins remain above 49.11 bits. These are sufficient
proof bounds, not estimates of the actual code's distance.

For hill climbing, `--retarget` permits an explicitly changed cutoff.
It reconstructs the saved partition, recomputes every old leaf at the new
cutoff, and requeues every bound that misses the per-cell budget. A failed
old witness is never silently retained as proved. Ordinary replay still
rejects changed cutoffs, and ordinary resume cannot raise them. The first
retarget run seeks 9% using the 8% partition.

An intermediate 8.5% retarget completed after 854 new cells, with 820
terminal leaves and no unresolved cells. Independent 256-bit replay gives
dense margin 74.6299285691175 bits and full aggregate margin greater than
49.11287660853 bits. The witness is
`tmp/two-bit-affine-dense-401-d085.json`; the original 8% witness is preserved.
The 9% run completed after 2256 new cells, with 1453 terminal leaves and
no gaps. Its independent 256-bit replay gives dense margin above
73.84716062499 bits and whole-code margin above 49.11287658688 bits.
The 9.25% run completed after 2302 new cells with 2235 terminal leaves.
Independent 256-bit replay gives dense margin above 71.23535577341 bits
and whole-code margin above 49.11287632260 bits. The 9.5% attempt stopped
with unresolved cells, even after increasing the subdivision depth.
An updated proof split assigned q=401,...,408 to the sparse method.
All eight full support covers passed at the 9.5% cutoff, with aggregate
margin above 1491.07046549378 bits. The q>=409 dense attempt remained
incomplete. The extension through q=424 completed; its last sixteen
saved support covers give aggregate margin above 206.24224259823 bits.
The q>=425 dense run stopped at its budget with 2733 accepted leaves
and 199 unresolved cells. Inspection on 2026-09-28 checked the sparse
records' scopes and partitions, but did not replay their numerical bounds.
The 9.5% effort is now paused in favor of joint route--inner tuning.
Its saved witnesses remain local and do not constitute a full certificate.
[HILL_CLIMB.md](HILL_CLIMB.md)
records those diagnostics and the three-update attempt, which is not yet
a full-code certificate.

## An interval method, and why it has not closed the dense range

Let R_r be the regional operator when exactly r packets are active. With
q active pairs and a common support witness p, define

    R_q(p) = sum_{r=0}^q binomial(q,r) p^r (1-p)^(q-r) R_r.

Introduce an auxiliary marking probability theta, independent of the
construction. Put G=4096 and

    beta_q(theta) = binomial(G,q) theta^q (1-theta)^(G-q),
    H(v) = sum_{j=0}^64 binomial(64,j) v^j (1-v)^(64-j) T_j.

Binomial thinning and independent virtual slot occupancies give the exact
matrix identity

    sum_q beta_q(theta) R_q(p) = H(theta*p)^64.

All matrices are nonnegative. Retaining one term therefore gives the
entrywise bound R_q(p) <= H(theta*p)^64 / beta_q(theta). Applying it to
each of 256 regions is valid even though the actual active-pair labels
are shared across regions. This is a coefficient inequality, not a change
to the setup distribution.

Let F(p) bound the entire nonzero-pair outer sum

    sum_{u=38}^{256} a(u) / (binomial(256,u) p^u (1-p)^(256-u)),

using the exact shell/CDF interfaces. For fixed lambda, theta, and p, the
resulting first-moment upper bound is

    B(q) = exp(lambda*209715) e_zero H(theta*p)^16384 tau
           * binomial(G,q) F(p)^q / beta_q(theta)^256.

For 0<theta<1 its successive ratio is

    B(q+1)/B(q)
      = F(p) ((q+1)/(G-q))^255 ((1-theta)/theta)^256.

The ratio increases with q, so B is log-convex on integer occupancies.
Consequently one frozen witness at both endpoints bounds an entire
integer interval. Its union bound is the interval length times the larger
endpoint bound. Different optimized witnesses at the endpoints do not
justify this conclusion. The boundary theta=1 permits only q=G.

`dense_cover.py` implements this method. Exact small noncommuting-matrix
tests check the thinning identity. Other tests check log-convexity and the
interval-length factor. **The current dense experiments do not close.**
In particular, this broad relaxation is positive even at q=32, where the
direct regional method already closes. Its extra conditioning loss is
therefore material, and failure of this bound is not evidence of a bad code.

The strongest tested broad bound at q=4096 still has log2 upper about
+815016, even after the mass-based refinement below. It cannot contribute
to the final failure budget. Do not replace the verified direct covers
with this relaxation.

## Local refinements and diagnostics

The new `conditional_atom` bound extends the exact feedback census beyond
its cutoff. Fix a shape with j packets. Leave d packets of a censused
subshape unexposed and expose the other j-d slots and lane masks. The
remaining packets are uniformly placed subject to avoiding those slots.
The avoidance probability is binomial(64-j+d,d)/binomial(64,d). If A_d
bounds every feedback atom of the unconditioned subshape, then

    A_d * binomial(64,d)/binomial(64-j+d,d)

bounds every atom after this conditioning. Include the zero atom in A_d,
because exposed feedback translates the target syndrome. The driver takes
the minimum over feasible censused subshapes and the trivial bound 1.
With the production injectivity check, d=1 recovers the old bound. Exact
small-map tests also include noninjective one-packet maps; the refinement
does not assume that a generic map has the production map's injectivity.

`--output-degree 64` evaluates the exact-law positive output polynomial
for all local packet shapes, not just those in the feedback census.
Pointwise output bounds now use |v-W|, where v is the expanded-state
weight and W is the input weight. This is stronger than max(0,v-W) when
the input is heavy. Neither change alters the encoder.

`dense_cover.py --mass-columns` uses a different complete bound for two
coupled transfer columns. Let M bound mature mass, let h bound its
pointwise output tilt, and let m bound its average output tilt for every
fixed nonzero entering state. Let a_* bound every feedback atom and let
a_+ bound nonzero feedback atoms. On the lazy branch, zero-return mass is
at most alpha*M*min(m,h*a_+), while every nonzero output-state atom is at
most alpha*M*min(m,h*a_*), with alpha=1/4. Replace the full density-based
columns by these mass-based columns; do not independently minimize their
entries. Existing refresh contributions remain unchanged.

The zero-feedback atom must be included in a_* for a density bound.
The production audit through eight occupied slots found no shape whose
zero atom exceeds its largest nonzero atom, so making this distinction
explicit did not change the previously reported numerical results.
Exhaustive small tests also cover maps where the zero atom does dominate.

This local refinement removes one source of artificial state persistence,
but does not resolve the dense gap. At large occupancies the proof still
maximizes packet shapes in the inner independently of the outer's profile
multiplicities. Retaining the pair profile (single-bit versus double-bit
packets) is the next substantial refinement to investigate.

The coarse selected-event screen at q=512 first exposed a binary64
underflow before matrix squaring. The proposal routine now normalizes
the initial matrix before squaring and rejects lost numerical mass as
evidence. Outward verification never used that failed screen. A regression
test checks a matrix scaled by 2^-900.

```sh
python -B research/workstreams/permutation_locality/two_bit/dense_cover.py --groups 32 64 128 256 512 1024 2048 4095 4096 --tilts .032 .064 .096 .128 .192 .256 --output-degree 64 --outward
python -B research/workstreams/permutation_locality/two_bit/dense_cover.py --groups 32 64 128 256 512 1024 2048 4095 4096 --tilts .016 .032 .064 .128 .256 .512 1 2 --cutoff 4 --output-degree 64 --mass-columns --outward
python -B -m unittest discover -s research/workstreams/permutation_locality/two_bit -p 'test_*.py'
```

Forty-nine tests pass. Exact convolution preparation is shared across output
tilts. A positive hypergeometric recurrence accelerates the binary64
placement screen; the outward regional polynomial remains unchanged.
The outward Poisson-binomial coefficient calculation now groups repeated
probabilities and uses Arb polynomial products. Exact rational recurrence
tests cover up to 320 factors, including unequal probabilities and tiny
tail coefficients. This removes the Python quadratic loop without changing
the bound or trusting binary64 coefficients.
No four-bit proof source, production encoder, or production parameter
selection has been changed.
