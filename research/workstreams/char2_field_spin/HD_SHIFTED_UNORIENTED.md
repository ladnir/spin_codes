# Removing singleton orientation from the coverage domain

The alternating route admits a regional bound depending only on total
support weight r and the full-pair counts a,b. It covers every split
between active even and odd groups. This reduces the coverage domain
from (E,O,a,b) to (r,a,b), and coefficient transforms from three
variables to two.

The directed grid-128 replay for (r,a,b)=(992,450,420) gives the exact
log2 upper ceiling

    -63581397015111 / 125000000000,

approximately -508.651176121. Replays at 256-bit and 384-bit Arb
precision returned the same ceiling. This receipt covers all 93
feasible choices 450<=E<=542, with O=992-E. It does not cover
other triples or establish a whole-code bound.

## A coefficientwise orientation envelope

Fix the weight marker z and exponent 0<alpha<=1. Let
C[u,v,f] be the unmarked 2-by-2 coefficient of one step's established
oriented polynomial. The step contains 32 source pairs; u,v,f count
types 10,01,11. These coefficients include the local geometry counts
and the entrywise alpha power of each state-operator cap.

For an unsigned singleton count s and full count f, define the matrix

    B[s,f][i,j] = max_(0<=u<=s) C[u,s-u,f][i,j] / binom(s,u).

Let B(X,y) be the resulting two-variable matrix polynomial. Expanding
(x10+x01)^s gives

    C(x10,x01,y) <= B(x10+x01,y)

coefficientwise and entrywise. The binomial divisor is necessary: the
substitution itself supplies binom(s,u) orientations of each monomial.

Partition each step into good arrangements, with at most one singleton,
and the remaining arrangements. Write G,D for their oriented polynomials.
The same coefficientwise envelope holds for G and D separately.
For the good part it is an equality:

    G(x10,x01,y) = G_unsigned(x10+x01,y).

As established in `HD_SHIFTED_GOOD_PAIRS.md`, the good-good two-step
polynomial P already depends on singleton orientation only through
X=x10+x01. Therefore define

    T_unsigned = P + G_unsigned D_unsigned
                   + D_unsigned G_unsigned + D_unsigned^2.

The oriented two-step polynomial is coefficientwise bounded by
T_unsigned(x10+x01,y). Positive matrix products preserve this order,
including through the 64th power. All products retain their matrix
order; commutation is not assumed.

For regional totals U,V,F and S=U+V, coefficient extraction gives

    [x10^U x01^V y^F] T_oriented^64
      <= binom(S,U) [X^S y^F] T_unsigned^64.

The oriented and unsigned regional normalizers satisfy

    multinomial(M;n0,U,V,F)
      = multinomial(M;n0,S,F) binom(S,U).

The extra binomial factor cancels exactly. Hence a valid regional
operator for every orientation split is

    R_unsigned(r,f)
      = [X^(r-2f) y^f] T_unsigned^64
        / multinomial(M;M-r+f,r-2f,f),

where M=4096. This proof uses the established oriented regional bound.
It introduces no independence assumption about singleton labels.

Positive markers for X and y cancel through the unsigned normalizer.
There is no additional factor 2^S in the regional bound. If original
markers were t10,t01,t11, the initial unsigned proposal uses
(t10+t01,t11). Other positive marker choices are also valid witnesses.

## Counting all even/odd splits

Let N3(r,a,b) count indexed cyclic supports of length 2M with the
specified weight and full-pair counts. This is the sum of the exact
four-parameter counts over E+O=r. For a nonconstant support, put
c=r-a-b. The existing function

    hd_shifted_pairs_count.profile_count_upper(M,r,a,b)

returns a sound upper bound using the minimum of the cyclic-run count

    B(2M,r,c) = (2M)/c * binom(r-1,c-1) binom(2M-r-1,c-1)

and the two one-matching counts

    binom(M,f) binom(M-f,r-2f) 2^(r-2f),    f in {a,b}.

If r=2a or r=2b, contraction gives the exact count B(M,r/2,c).
The function handles constant and invalid classes separately.
For the displayed triple its log2 upper is 1413.500493891. This
already counts every even/odd split; it must not be multiplied by
the number of splits again.

## Replay and numerical checks

From the theory directory, the directed replays are

    python hd_shifted_unoriented_certificate.py --grid 128 --precision 256
    python hd_shifted_unoriented_certificate.py --grid 128 --precision 384

The executable uses the saved asymmetric checkpoint, exact rational
marker sums, rounded rational packet caps, and outward interval FFT
arithmetic. Both runs returned the rational ceiling above and marked
`whole_code_certificate=False`.

The separate FLOAT implementation is `hd_shifted_unoriented_screen.py`.
At the same checkpoint it obtained these all-orientation bounds:

| Grid | log2 upper screen |
| --- | ---: |
| 64 | -288.252943826 |
| 128 | -508.651173482 |
| 256 | -526.209292585 |

The five tests in `hd_shifted_unoriented_screen_tests.py` check the
coefficientwise envelope, binomial cancellation for a constant kernel,
the unchanged good coefficients, all 93 histogram splits, degree-64
folding, and the global convex bound below. Run them with

    python -m unittest hd_shifted_unoriented_screen_tests -v

## Coarse boxes without coefficient transforms

A global positive-marker bound evaluates T_unsigned at its markers
instead of extracting coefficients. For a fixed witness, its regional
numerator is independent of the support profile. Only the outer density,
support count, marker powers, and multinomial normalizers vary.

For this purpose, use the one-matching count

    N3(r,a,b) <= M0 * 2^s0,
    Hj = (M-r+fj,r-2fj,fj),
    Mj = multinomial(M;Hj),     s_j=r-2f_j,
    f0=a,                      f1=b.

Let ell_sj and ell_fj be the log2 markers for region type j. Absorb
the fixed matrix product and fixed outer terms into K. The log2 bound
for one triple is

    F(r,a,b) = K + alpha*r*(2048-128*log2(255)) + s0
               -127*log2(M0) -128*log2(M1)
               -128*sum_j(s_j*ell_sj + f_j*ell_fj).

The function -log(multinomial) is a sum of convex log-factorials,
up to a constant. Every histogram entry is affine in (r,a,b).
Thus F is convex on the domain of nonnegative histograms. Its maximum
on a rectangular box occurs at a corner.

If all eight corners have valid histograms, evaluate only those corners.
Add log2 of the number of integer triples to bound the union over
the whole box. Infeasible triples only make this count looser. The
argument covers every even/odd split of each feasible triple. It uses
neither sampled interior points nor independent maximization of regional
normalizers. The helper `global_box_bound` implements this FLOAT screen.

Two directed global-box receipts illustrate the scale. They bound the
union over their entire domains, including every even/odd split:

| Domain | Number of triples | Directed log2 upper |
| --- | ---: | ---: |
| 4032<=r<=4160, 1536<=a,b<=1664 | 2146689 | -337.920794262092 |
| 6016<=r<=6272, 2372<=a,b<=2628 | 16974593 | -30008.472335239249 |

The exact ceilings are -84480198565523/250000000000 for the first
box and -30008472335239249/1000000000000 for the second. Both
boxes returned the same exact ceiling at 256-bit and 384-bit precision.
The witnesses use the common exact exponent

    alpha = 4097437497366089/36028797018963968.

The first uses

    z = 7722805725266083/144115188075855872,
    singleton marker = 244777405/305135501,
    full marker      = 1189188654/621493969.

The second uses

    z = 7151016493358955/576460752303423488,
    singleton marker = 973542133/268519977,
    full marker      = 9096903675/859563932.

Both regional operators use the same marker pair within each witness.
These exact rational choices are obtained by the executable's fixed
conversion of the original floating proposals. Reproduce the directed
receipts with

    python -B hd_shifted_unoriented_certificate.py --coarse middle --precision 256
    python -B hd_shifted_unoriented_certificate.py --coarse middle --precision 384
    python -B hd_shifted_unoriented_certificate.py --coarse dense --precision 256
    python -B hd_shifted_unoriented_certificate.py --coarse dense --precision 384

At the lower-weight triple (1024,128,128), a short search of this
unsigned global bound remained positive. These global boxes therefore
complement the regional alias method; they do not replace it everywhere.

## A route-independent terminal range

The existing inverse-image argument handles sufficiently dense supports
without distinguishing matching histograms. Let N=2097152 and
D=734003. A q-ary Hamming-ball bound at positive mark t gives

    ball(t) = N*log2(1+255*t) - D*log2(t).

Fix a route and inner setup. Invertibility gives each candidate output
one inverse input, whose active-group set is determined. The established
outer pointwise bound and projective scaling then bound all supports
with at least r0 active groups by

    ball(t) - (128*r0+1)*log2(255).

There is no subset-count factor and no additional sum over r. The
optimal rational mark for this expression is

    t = D / ((N-D)*255) = 734003/347602995.

A 256-bit and 384-bit Arb replay using
`hd_randomized_outer_proof.ball_log2` and `dense_class_log2` gives,
for all r>=7649, the same exact ceiling

    log2 upper <= -146842460599147/500000000000,

approximately -293.684921198294. The evaluation at r=7648 is
positive, so this particular ball envelope does not certify
that earlier endpoint. The optimized mark improves the existing fixed
mark's terminal range, which started at r=7665.

The underlying enclosure can be reproduced from the theory directory:

```python
from fractions import Fraction
from flint import ctx
from hd_randomized_outer_proof import Parameters, ball_log2, dense_class_log2
p = Parameters(field=256, state=16, limit=734003)
for precision in (256, 384):
    with ctx.workprec(precision):
        ball = ball_log2(p, mark=Fraction(734003, 347602995))
        print(dense_class_log2(7649, p, ball=ball))
```

The next coverage pass can therefore remove the terminal range, accept
the directed global boxes, and use the two-variable alias
atlas for unresolved boxes. It must maintain disjoint coverage and sum
the accepted domain bounds. No result in this note certifies the entire
remaining support domain.

## A bounded FLOAT bulk-coverage plan

The separate driver `theory/hd_shifted_bulk_search.py` constructs a
disjoint box partition for 4096<=r<=7648. The saved plan is
`theory/hd_shifted_bulk_plan.json`. Its accepted leaves cover
1,538,866,458 feasible unsigned triples using 251 boxes and 140
distinct witnesses. Each triple includes every even/odd support split.
The search produces FLOAT proposals; their later directed replays are
recorded below.

For each fixed witness, the driver computes the positive regional
matrices once. It then evaluates the profile-dependent affine and
log-factorial terms on each box. It considers both one-matching
support counts, taking the smaller of the two complete box envelopes.
It does not choose a different witness at each vertex.

A box may cross a histogram boundary. The driver intersects it with
the domain

    max(0,r-M) <= a,b <= r/2,

and uses the rational vertices returned by
`hd_shifted_cover_backend.histogram_vertices`. Convexity of each complete
one-matching envelope still bounds its maximum by these vertices.
The union factor is the exact feasible integer-triple count supplied
by the same backend. For nonconstant supports, the all-full point
a=b=r/2 is excluded. Half-integer vertices are retained when needed.

The saved pass spent exactly 400 witness proposals. It used breadth-first
subdivision, with denser weight slabs visited first at each depth.
Every witness has alpha=.11372673628845778; its weight marker is
.23149023644229288 raised to the slab midpoint divided by 2048.
Singleton/full marker proposals use the box's mean clipped vertex
histograms, followed by a full-box check. There are 12 distinct weight
markers among the used witnesses. Witness quality is heuristic;
acceptance tests the entire clipped box, not a selection of profiles.

The accepted leaves have aggregate FLOAT log2 upper -210.83675373360637.
Their distribution is:

| Weight range | Accepted triples | Feasible triples in range |
| --- | ---: | ---: |
| 4096..5119 | 184224255 | 3314374656 |
| 5120..6143 | 746000380 | 1702713344 |
| 6144..7167 | 532127935 | 627922944 |
| 7168..7648 | 76513888 | 76513888 |

Inspection of the retained leaves gives the stronger endpoint 6656:
the plan screens the entire unsigned domain for 6656<=r<=7648.
Its 494 unresolved leaves all lie below 6656. They are retained in the
artifact; failure to accept them within this proposal budget does not
show that the global bound fails there. The range 2048<=r<4096 was not
searched in this pass.

Every accepted leaf records its inclusive ranges, a specific witness
index, exact integer-triple count, and FLOAT union bound. The witness
records contain the weight marker, exponent, and regional singleton/full
log2 markers. Directed replay must convert these to exact positive
markers and recompute the matrix and scalar bounds; it must not reuse
the FLOAT log bounds as certificates. Existing certified domains may
overlap this plan. The coverage ledger subtracts overlaps for exact
disjoint accounting. Summing overlapping event bounds is also a valid
union bound, but can be looser and does not give an exact domain count.

Reproduce the proposed partition and run its sixteen tests from the theory
directory:

    python -B hd_shifted_bulk_search.py --first 4096 --last 7648 --proposals 400 --compact
    python -B -m unittest hd_shifted_bulk_search_tests -v

The tests verify cached evaluation against the earlier global helper,
clipped-domain domination, exact disjoint partition counts, saved
witness replays, and complete coverage of the displayed dense band.
The JSON files are local replay artifacts excluded by the repository's
existing ignore rule. Their six integration tests skip with an explicit
reason when the artifacts are absent; the other ten tests remain
lightweight. Plan regeneration uses the explicit bounded commands above.
The immediate next step is directed replay of the saved leaves,
prioritizing 6656..7648 to extend the contiguous terminal range. After
that, refine the unresolved lower slabs or use the regional alias
envelope where the global positive bound remains too loose.

A bounded continuation reuses this bank and refines only the ten
unresolved leaves in 6144..7167. It closes that entire range with 17
additional witness proposals. The separate artifact
`theory/hd_shifted_bulk_plan2.json` contains 95 accepted leaves and no
unresolved leaves. Of these, 76 were inherited from the first plan;
19 new leaves cover the remaining 95,795,009 triples with aggregate
FLOAT log2 upper -1228.683970413118. Their field `inherited` is false.
Witness indices refer to the second artifact's own compact bank.

The two plans together therefore screen every feasible unsigned triple
for 6144<=r<=7648, including every even/odd split. The second plan's
full domain has 627,922,944 triples. Its partition test checks exact
counts, disjointness, and containment of each new leaf in a formerly
unresolved leaf. The first plan was not modified. Reproduce the
continuation with

    python -B hd_shifted_bulk_search.py --first 6144 --last 7167 --proposals 400 --resume hd_shifted_bulk_plan.json --compact

Directed replay need only process the 19 new leaves in addition to
the first plan. No new result in this continuation covers r<6144.

A third bounded pass refines only the first plan's unresolved leaves
in 5120..6143. It uses exactly 400 additional proposals. The local
artifact `theory/hd_shifted_bulk_plan3.json` records 383 accepted leaves,
of which 104 were inherited and 279 are new. They cover 1,458,395,260
of the range's 1,702,713,344 triples. The new leaves add 712,394,880
triples. Its 156 unresolved leaves contain 244,318,084 triples and all
end at or below r=6015. Hence the three plans cover every unsigned
triple for 6016<=r<=7648. Regenerate this third proposal with

    python -B hd_shifted_bulk_search.py --first 5120 --last 6143 --proposals 400 --resume hd_shifted_bulk_plan.json --compact

The first plan has now passed directed replay of all 251 leaves at
256-bit and 384-bit precision. Both give aggregate log2 ceiling

    -42167350746281/200000000000,

approximately -210.836753731405. The second plan has passed all 95
leaves at both precisions with the same aggregate ceiling. The third
plan has passed all 383 leaves at both precisions, giving

    -52265110262917/500000000000,

approximately -104.530220525834. These aggregates overlap through
inherited leaves. Their sum remains a valid upper bound; ledger
reconciliation is needed for exact disjoint counts, not for validity
of that union bound. Restricting the receipts to the separate bands
6016..6143, 6144..7167, 7168..7648, and r>=7649 gives the directed
combined ceiling -104530220525833/1000000000000. The domain contains
860,977,921 unsigned triples. It does not include the remaining
lower-weight domain; see `HD_SHIFTED_DOMAIN_COVERAGE.md` for the ledger.

## Why the first witness bank stalled at low empty-pair counts

The first bank fixed alpha=.11372673628845778 while changing the weight
and singleton/full markers. This restriction is not required by the
proof: every 0<alpha<=1 gives an admissible fractional-moment bound.
Its best value can depend strongly on the support histogram.

At (r,a,b)=(5120,1280,1280), refining only the two common regional
marker odds at z=.023538031086279037 and the bank's fixed alpha leaves
the global bound at +21052.213179. A bounded 100-iteration refinement
of z, alpha, and those two marker odds instead gives -94546.545970.
The new witness is

    z = .011103584869424775,
    alpha = .5357886737301887,
    (singleton odds, full odds) = (6.22715071206585, 8.584298238117446).

Both regions use the same marker pair. This comparison is for one
profile, so box width cannot explain the fixed-alpha failure. The
zero-empty endpoint (5120,1024,1024) also closes: at the new z,alpha,
marker odds (12.333868598163452,13.999999999253784) give FLOAT
-5788.208157. No alias extraction or additional orientation information
is needed for these profile witnesses.

Fixed-profile markers need not be good box markers. In particular, the
first new point witness still gives +55121.294422 on the original
difficult box. Refining the marker odds for the complete clipped box
closes it without changing z or alpha. Three directed replays give:

| Box (r; a; b) | Marker odds, both regions | Triples | Directed log2 upper |
| --- | --- | ---: | ---: |
| 5120..5247; 1024..1278; 1024..1231 | (7.905242031272326, 9.62627537577853) | 3716736 | -244.602115362210 |
| 5120..5375; 1280..1535; 1280..1535 | (6.283482967970975, 8.617029057088011) | 16777216 | -52742.002215632226 |
| 5120..5247; 1152..1407; 1152..1407 | (6.907696045558403, 8.905543123712771) | 8388608 | -31487.133480863329 |

The 256-bit and 384-bit replays return the same exact ceilings, respectively
-24460211536221/100000000000,
-26371001107816113/500000000000, and
-31487133480863329/1000000000000. Their exact weight marker and exponent are
1600195221773379/144115188075855872 and
1206488835680305/2251799813685248. Each positive marker is the Python
binary float 2**odds converted to a fraction with denominator at most
10^9. These boxes overlap; their individual receipts do not imply that
all profiles below 6016 are covered.

Six fixed FLOAT checkpoints, including these three, are stored in
`hd_shifted_bulk_search.py`. No witness optimization runs during replay:

    python -B hd_shifted_bulk_search.py --alpha-checkpoints
    python -B hd_shifted_bulk_search.py --alpha-checkpoints --precision 256
    python -B hd_shifted_bulk_search.py --alpha-checkpoints --precision 384

The structural diagnosis is therefore twofold: the bank's fixed alpha
caused genuine pointwise slack, and a witness optimized at one profile
can incur substantial extra slack across a box. Neither observation
requires a change to the construction or a stronger regional operator.
The next coverage pass should add the admissible larger-alpha family
and optimize complete-box marker envelopes before attributing remaining
failures to loss of orientation or spectral information.

## An optional adaptive-alpha continuation

The optional `AdaptiveAlphaProposer` implements that correction without
changing the earlier search's default behavior. It first adds the six
validated larger-alpha witnesses to the available bank. A box that still
fails receives four independently adjustable regional marker odds.
The proposer first tries alpha=.5357886737301887, followed by
.25, .5, .75, and 1, all at z=.011103584869424775. A short refinement
of z is allowed only if those marker refinements fail.

Every objective evaluation uses the complete clipped-box envelope.
Finding one sufficient witness ends refinement for that box. The work
limits are explicit: one global objective-call budget, one marker-call
limit per alpha, one weight-refinement call limit, and the outer search's
box-proposal limit. Rejected numerical proposals also consume the
objective-call budget. Reconstruction of inherited witnesses is a
separate fixed initialization cost.

The saved fourth plan, `theory/hd_shifted_bulk_plan4.json`, resumes only
the third plan's unresolved leaves. It FLOAT-covers all 1,702,713,344
feasible triples in 5120<=r<=6143, including every even/odd split.
Its 542 accepted leaves comprise 383 inherited leaves and 159 new leaves;
no unresolved leaf remains. The new leaves contain 244,318,084 triples.
The initial 156 unresolved leaves required only three net additional
subdivisions. Thus this improvement comes primarily from better
witnesses, not from replacing the domain by tiny boxes.

The run used 46 box proposals and 9,153 of 12,000 allowed objective
evaluations. It performed 87 marker refinements and three short weight
refinements. Among the accepted new witnesses, only two parameter pairs
were needed: the common z above with alpha=.5357886737301887 or .25.
The 159 new leaves use 47 marker witnesses. Incremental directed replay
therefore needs 159 leaf checks, 47 regional marker pairs, and two
cached unmarked coefficient constructions.

The weakest new leaf has FLOAT log2 bound -163.27786414114024. The
complete plan's FLOAT aggregate remains -104.53022052533925 because inherited
leaves dominate it. All 542 leaves now pass directed replay at 256 and
384 bits, with exact log2 ceiling -52265110262917/500000000000.
The first three artifacts remain
unchanged. Regenerate the fourth proposal explicitly with

    python -B hd_shifted_bulk_search.py --first 5120 --last 6143 --proposals 156 --resume hd_shifted_bulk_plan3.json --adaptive-alpha --max-evaluations 12000 --marker-maxfev 200 --weight-maxfev 32 --compact

The [domain ledger](HD_SHIFTED_DOMAIN_COVERAGE.md) joins this band with
the earlier bands and dense tail. Every r>=5120 is covered; smaller
supports remain unresolved. This is not a whole-code certificate.

Two bounded probes below this band did not justify another bulk sweep.
At (4096,0,0) and (4608,512,512), the same five-alpha family with limited
weight refinement remained positive, approximately +169704.5 and
+81515.8. Each probe used 1,032 objective calls. These are not globally
optimized values and do not establish a limitation of the regional
operator. Broader weight-marker refinement or a structural endpoint
argument should precede a lower-band partition search.

For example, (4096,0,0) contains exactly the two alternating supports.
They force E=0 or E=4096. The unsigned local envelope nevertheless
maximizes over all local singleton orientation counts. Comparing this
envelope with the forced oriented operator can isolate one possible
source of slack. No such comparison is assumed in the fourth plan.

## Run-aligned boxes and a coupled support sum

The isolated helper `theory/hd_shifted_run_geometry.py` uses coordinates
(r,c,d), where c=r-a-b and d=a-b. The inverse map is

    a=(r-c+d)/2,     b=(r-c-d)/2.

It accepts a rectangular run box only when c>=1 and all eight mapped
corners have nonnegative matching histograms. This excludes constant
supports and avoids clipping. Mapped vertices may be half-integral;
the scalar envelope uses the same log-Gamma extension as the existing
backend. Integer profiles satisfy r-c+d even.

Let n_r,n_c,n_d be the three inclusive interval lengths. For each
interval, let e be its number of even integers minus its number of odd
integers. The exact number of feasible integer profiles in an accepted
box is

    (n_r*n_c*n_d + e_r*e_c*e_d)/2.

Each e is zero for an even-length interval. Otherwise it is +1 or -1
according to the parity of the lower endpoint. This profile count is
ledger metadata. It is not an additional factor in the weighted support
sum below.

Fix positive rational run markers x,y and let

    G(x,y)=trace([[1,x*y],[1,x]]^(2*M)).

For every indexed binary support U, this generating function includes
the weight x^r(U)*y^c(U). Let B be an accepted run box. Suppose a common
regional numerator gives a per-support bound K*2^S(r,a,b) throughout B.
Here S is the backend's scalar with no support-count term. Then

    sum_(U in B) K*2^S(U)
      <= K*G(x,y)*max_(p in B) 2^(S(p)-r*log2(x)-c*log2(y)).

The inequality follows by factoring x^r*y^c from each summand and
bounding the remaining factor uniformly. The generating function may
include supports outside B; their positive weights only enlarge the
bound. This argument counts all supports once and needs no extra
profile-count multiplier.

The tilted scalar S-r*log2(x)-c*log2(y) is convex because the subtracted
term is affine. Its maximum is therefore bounded by the eight mapped
vertices. Combining these terms before maximization retains information
lost by multiplying a separate support-count maximum by a scalar
maximum. The helper takes the smaller of the coupled bound and the
earlier separated bound.

The regional alias targets follow directly from the mapped vertices:

    region 0: singleton=c-d, full=(r-c+d)/2;
    region 1: singleton=c+d, full=(r-c-d)/2.

Their extrema are rounded outward before querying the existing alias
tables. `run_box_receipt` reuses a `CachedBound` object's regional
numerators and scalar; it does not modify the shared backend.

At center (r,c,d)=(1024,224,0), the saved grid-128 witness and run
markers x=25/31,y=7/775 give the following directed 256-bit results:

| Radius (r,c,d) | Integer profiles | Coupled log2 upper |
| --- | ---: | ---: |
| (32,4,4) | 2633 | +130.430907181772 |
| (64,8,4) | 9869 | +729.679709196220 |
| (64,8,8) | 18641 | +776.353300782787 |

The first box's separated bound was +184.756487806604. Coupling the
run tilt improves it by about 54.33 bits, but does not close it. The
requested box therefore covers 2,633 parity-valid profiles, not ten
thousand. The larger tests also remain positive. These computations
validate the geometry and expose the remaining slack; they do not
establish a large-box distance certificate.

Run the six geometry tests and the first displayed replay with

    python -B -m unittest hd_shifted_run_geometry_tests -v
    python -B hd_shifted_run_geometry.py --grid 128 --precision 256

The tests exhaust all admissible run boxes for M=1,2,3,4 against actual
binary cycles. They also check half-integral vertices, exact parity,
the weighted-sum inequality, alias target ranges, and absence of an
extra profile-count factor. Further coverage work was stopped pending
the separate structural review of the alternating-support endpoint.
All previous route receipts and frozen plans remain unchanged.
