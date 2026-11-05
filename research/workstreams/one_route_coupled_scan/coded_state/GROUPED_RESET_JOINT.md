# Joint lengths, inactive groups, and a paired-reset alternative

## Follow-up: relaxing the exploratory target to 45%

Update: these proposals now have complete directed certificates, including
the d=1 running-sum variant. See [the certificate and feedback reduction](GROUPED_45_CERTIFICATE.md).
The paragraphs below retain the initial FLOAT-stage record.

At the user's suggestion, keep the group32 construction fixed and lower
the bad-output threshold to W=floor(.45N). The simpler two-dimensional
bound in `striped_reset_probe.py` now covers the entire class range for
both profiles. It does not need the joint length/activity refinement.

| Outer / total state | Accepted rectangles | Worst log2 class bound | Required per-class bits |
|---|---:|---:|---:|
| RS[512,256] / 256 | 300 | -69.79359 | 65 |
| RS[256,128] / 512 | 309 | -65.77481 | 63 |

These are complete FLOAT coverage proposals at a 40-bit whole-code target,
not directed certificates. Independent coverage checks and directed
arithmetic replay remain required. The construction, reset-only timings,
and uniform-feedback assumption are unchanged. In particular, the roughly
2.8ms group32 timing excludes outer, routing, and feedback.

Reproduce serially with `python -B coded_state/striped_reset_probe.py
--outer 512 --group 32 --cover --distance-permille 450`, then replace
`--outer 512` by `--outer 256`. Both runs retain the existing bounded cover
and generate no files. The 49.1% diagnostics below remain preserved.

The next step is to verify these simpler complete covers, then reduce
feedback cost using 45% as a certified control before climbing toward49%.

The target remains K=2^20, rate1/2 over GF(2^128), distance >49.1%, and
setup failure below2^-40. Both outer/route profiles and full uniform
feedback are unchanged. All setup randomness is sampled once and fixed
for every message. This note contains new counting bounds and local
checks, not a completed certificate or new native timings.

The [group32 candidate](RESET_TIMING_AND_GROUPS.md) saves reset work but
permits a closed excursion to use only a few reset groups. The previous
bound charged all possible excursion lengths while allowing the output
capacity of the shortest excursions. Those advantages cannot always occur
together. The new count retains their common parameter.

## Inactive groups cost probability or output capacity

Write b=Gd for the batch and state size, with G independent reset groups
of size d. Fix c ordered excursions occupying v batches. The first
r=c-1 excursions are closed; the last is allowed to remain open. For
closed excursion i, let its length be 2+e_i and let a_i be the number of
groups containing any nonzero output. Fix these groups, not just their
number. Put

    A = sum_i a_i,
    H = sum_i d(G-a_i)e_i,
    L = v-2r-1.

The final excursion has length1+e_last, and sum_i e_i+e_last=L.
The parameter H counts inactive output coordinates in the nonterminal,
noninitial batches of closed excursions. It excludes terminal coordinates:
the required terminal reset already fixes the entire output batch.

Fix the route and a nonzero outer input before exposing feedback. At each
live batch, the independent uniform feedback makes the entire output
uniform on F^b. Requiring inactive coordinates to vanish costs
q^(-d(G-a_i)). Each prescribed terminal reset costs q^-b. Sequential
exposure, discarding histories with the wrong state pattern, gives

    Pr_inner[the fixed pattern] <= q^(-br-H).

At each activation the state is zero, so X=Y. Thus the input is supported
in the same activation-restricted set S as before:

    |S| = b(v-r)+dA.

This forward argument still counts projective outer inputs in S. It is
not a claim that an adversarial message was chosen before setup.

For the inverse argument, inactive groups remove coordinates from every
batch of their closed excursions. The exact allowed output capacity is

    E-H, where E=b(v-2r)+2dA.

Closed reset equations have rank dA. The randomized inverse-input
coordinates still number a=b(v-c). Therefore output-volume counting and
the outer moment use

    q^(-dA) Volume_q(E-H,W) * q^(D(A_coords)-a),
    W=floor(.491N).

Here A_coords is the set of noninitial excursion batches; it is distinct
from the integer A. Inactive output coordinates do not remove inverse
randomness at a live batch under full uniform feedback.

## Sum masks and lengths together

For a route marker tau=q^nu, 0<nu<=1, the forward argument weights a
pattern by tau^-H. For any alpha>0 define

    F_alpha(x,z) = sum_{j=1}^G binom(G,j) x^j
                              / (1-z*tau^(-alpha*d(G-j))).

For fixed v,c,A, the exact weighted count of group masks and lengths is

    C_alpha = [x^A z^L] F_alpha(x,z)^r / (1-z).

The final factor counts the unrestricted final length. The remaining
placement count is binom(n-v+c,c), where n=N/b. Every coefficient is
nonnegative, so any x>0 and0<z<1 gives the valid bound

    C_alpha <= x^-A z^-L F_alpha(x,z)^r / (1-z).

Use alpha=1 for the forward bound. Optimizing x and z affects tightness,
not validity. The FLOAT implementation also retains the old unweighted
combinatorial bound as an independent alternative.

For the inverse bound, use a linear upper bound

    log_q Volume_q(E-H,W) <= C + alpha(E-H).

Available bounds include the full-space bound (C,alpha)=(0,1), the support
count (W,1/128), and entropy tangents. For an anchor u>W, the latter use

    alpha = log(u/(u-W))/(128 log2),
    C = W + Ent(u,W)/(128 log2) - alpha*u.

The chosen anchors have alpha<=1. Concavity proves the bound for capacity
above W. Below W, the tangent remains above the full-space bound because
its value at W is at least W and its slope is at most1. The same
coefficient C_alpha then sums the inverse count. This retains H without
enumerating every possible total H or treating it as independent of length.

The existing regional moment bounds are unchanged. Take the smaller
complete forward or inverse class bound and then union over classes.
Do not multiply the two bounds or multiply independently averaged events.

## Bounded diagnostics

`group32_joint.py` examines the same selected classes as the previous
screen. It is FLOAT-only and does not establish a setup-failure margin.

| Profile | Previous group32 worst | Joint forward only | Joint forward and inverse |
|---|---:|---:|---:|
| RS[512,256], state256,105 classes | +851.80 | +303.70 | +27.01 |
| RS[256,128], state512,100 classes | +1443.94 | +772.22 | +772.22 |

Entries are log2 upper bounds before capping probabilities at1. Positive
entries are failed upper bounds, not low-weight witnesses. Finer rational
route markers reduce the selected larger-outer worst class to about+5.81;
this does not certify the class or complete the range.

The inverse improvement matters. At larger-outer v4224,c128,A508, the bound
changes from+303.70 to-329.37 using the original marker set. The remaining
larger-outer example is v4224,c128,A254. No statement about the maximum over
all classes follows from this selected screen.

A separate variance cap uses actual short-excursion geometry. If f is
the number of fully active closed excursions, then

    f <= floor((A-r)/(G-1)).

There are m=16 batches per region in both profiles. For activation support,
the sum of Bernoulli variances for one outer block is at most

    min(B/4, c/2,
        1/2 + (r+A/G)/m + f*(1/2-2/m) + H/(md)).

The inverse-coordinate version replaces (r+A/G)/m by r/m and2/m by1/m.
Disjoint interval contributions may be split before bounding variance;
combining them in one region only decreases that sum. A Bennett-tail
envelope is valid but did not improve the tested critical bounds. This
route refinement is retained as a diagnostic, not claimed as a win.

## Cheap mixing: reject a weak option, retain paired resets

A permutation of whole reset groups is not the needed mixing. In this
uniform-feedback model it can be moved to time-dependent permutations of
input and output groups. Those input permutations remain inside regions,
where the route is already uniform. Output permutations preserve weight.
Thus this change does not improve the ensemble's distance distribution.

A sparse companion matrix is less trivial. Over GF256 with modulus0x11b,
the exact checks find irreducible polynomials

    x^8+x^4+x^3+x+2,
    x^16+x^4+x^3+x+15.

Their companions have no return of a cyclic basis vector within8192 steps.
Each update costs one group scaling, three group XORs, and a cyclic buffer
rotation. Consecutive powers applied to a physical group span all G group
directions. These rank identities over GF256 persist over F. They do not
say every F-valued state is cyclic after extending the field.

However, these mixers still leave an exact short-excursion weakness.
For most physical groups, P e_g=e_(g+1). A two-batch output can follow
that shift. Its unrolled reset columns share one group direction, so any
d+1 distinct local columns have rank d. The same low-weight reset survives.
`group32_mixer.py` checks this negative control alongside the favorable
consecutive-rank examples. A long period alone is not sufficient.

The more direct continuation couples each pair of32-symbol reset groups
into one64-symbol MDS reset group. This keeps total state b, full uniform
feedback, and the outer/route. It uses the existing degree63 novel-moment
construction, with the required64-point transform and translations. It is
not merely an extra XOR or a free coordinate relabeling. Its cost must be
measured; the group32 timing cannot be assigned to it.

Every active paired group now needs at least65 output symbols to close
and contributes rank64. The same joint proof applies with d=64. The
bounded screens pass all84 selected larger-outer classes, with worst
log2 bound below-521.84, and all100 selected original-outer classes, with
worst below-185.54. These are encouraging screens, not full certificates.

## Reproduction and next gate

Run serially with Python -B:

- `coded_state/group32_joint.py` for the group32 joint screen;
- add `--group 64` for paired resets;
- add `--selected --fine` for the finer critical-class markers;
- `coded_state/test_group32_joint.py` for local counting and algebra checks;
- `coded_state/group32_mixer.py` for the sparse-mixer checks and counterexample;
- `coded_state/group_reset_cover.py --outer 512 --group 64` for a bounded
  three-dimensional FLOAT coverage proposal.

The coverage proposal bounds each rectangle uniformly in v,c,A, rather
than checking only its center. It uses bounded positive generating-function
markers, linear maxima, and the retained cardinality/partial-region moment
envelopes. Its output is not a directed certificate, even if every rectangle
passes. A separate directed replay and independent coverage check are still
required before a whole-code claim.

The next gate is full class coverage for paired resets, followed by directed
verification if the proposal closes. Keep the original group32 candidate
open; its improved count remains useful, but do not spend an unbounded
search on the same failing intermediate classes. No production or paper
changes, coefficient dumps, or new benchmark results are part of this work.

### Coverage status and a separate enumerator attempt

The larger-outer, group64 three-dimensional cover has not completed. A
150-second bounded proposal accepted2234 rectangles after4476 checks,
then expired with9 frontier rectangles. Refining linear maximization to
the feasible polytope r<=A<=Gr, v>=2r+1 changed the subdivision path; a
180-second run also expired, after3583 checks and1787 accepted rectangles.
These are incomplete proposals, not evidence that the unseen classes pass.
No directed replay was attempted. Small rounded-parameter caches only
propose x,z; every bound reevaluates the actual rectangle and marker.

Additional fine-marker spot checks cover intermediate activity values
missing from the original screen. For example, larger-outer d64 at
v4224,c128,A190 gives about-317.45, and original-outer d64 at
v2048,c32,A155 gives about-208.99. These remain individual FLOAT bounds.

There is also an alternative output enumerator, preserved in
`group_reset_enumerator.py`. For0<=z<=1, the zero coset of a linear code
maximizes its weight enumerator among cosets. To see this, expand the
coset indicator in additive characters. The Fourier factors of the
coordinate weight function are1+(q-1)z at zero and1-z otherwise. They
are nonnegative; every character phase is at most1, with equality for
the zero coset. This proves the claimed maximum without assuming the
syndrome is random.

For one closed group of length l>=2, choose2d output coordinates and fix
the others. The remaining reset equations define a coset of a[2d,d]
MDS code. Writing P_2d(z) for its zero-coset enumerator and
S(z)=1+(q-1)z gives

    P_dl(z) <= P_2d(z) S(z)^(d(l-2)),
    P_2d(z) <= 1 + sum_{w=d+1}^{2d} binom(2d,w) q^(w-d) z^w.

Consequently all output group choices can be summed before taking the
fractional moment, using

    Volume <= inf_{0<z<=1} z^-W P_2d(z)^(Gr) S(z)^(b(v-2r)).

The bound is uniform in the lengths but loses the syndrome/length coupling.
At larger-outer group32, v4224,c128, it gives about+731.09 after routing and
placement counting; v4352,c256 gives about+1300.43. It does not supersede
the joint activity/length argument. Preserve it as an independently valid
alternative, not a successful closure.
