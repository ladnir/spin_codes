# Cheap coded-state inners: local facts and global proof target

Status: research candidates, not whole-code distance certificates. The goal is
at least 49% symbol distance with 40-bit setup-failure margin at K=2^20.
The certified >49.1% window construction remains unchanged. The new
[ideal reset ladder](IDEAL_RESET_LADDER.md) adds whole-code proof controls:
state512 with the original outer, or state256 with RS[512,256], both with
MDS resets and uniform feedback. Neither certifies the cheap candidates
below or the256-state running sum. No new native timing is claimed.

The [reset-coordinate reduction](CHEAP_RESET_COORDINATES.md) now realizes
both ideal controls with structured reset transforms and unchanged uniform
feedback. It preserves their certificates algebraically. Native reset and
end-to-end timing comparisons remain open; it is not a certificate for the
older structured-feedback candidates below.

Update2026-10-22: the [random-feedback bridge](RANDOM_FEEDBACK_BRIDGE.md)
gives a genuine reset-supported obstruction to the64-state target under
the current route distribution. Its setup failure is at least1.88% at
49.1%, even with ideal feedback randomness. The partial results below are
preserved, not promoted to a whole-code claim. Continue from a larger-memory
ideal-randomness bridge; do not keep tuning this excluded ensemble's proof.

The next [structural alternative](POST_DIAGONAL.md) moves the diagonal
after A. At the same full-field multiplication count as two taps, it gives
uniform state refresh from every nonzero state and a coordinate inverse law.
It remains a research candidate, with reset and hole enumeration unresolved.
Its [sparse-source certificate](SPARSE_PROOF.md) now closes all messages in
1..31 outer blocks, including all resets and state defects, at >49.1%/50 bits.
The [regional extension](ROUTE_EPISODES.md) closes32..64 blocks at47 bits;
combined failure is below2^-46. Larger source-support classes remain open.
An inner-only [column-scaled variant](QUIET_BRANCH.md) adds one GF256
scaling per output and retains those certificates. It also closes every
entirely quiet trajectory at >49.1%/66 bits, with free terminal state.
Mixed live/quiet trajectories remain unresolved; no native timing is claimed.

The [alternating-map variant](ALTERNATING_MAPS.md) now alternates two disjoint
evaluation halves. It bounds holes across an unchanged-state pair by63,
rather than126, at the same full-field product count. The established
partial certificates transfer; the whole-code mixed-trajectory sum remains
open. A joint two-observation rank bound and exact toy checks isolate the
next sparse-update proof question without asserting generic full rank.
The [one-half update theorem](HALF_UPDATES.md) now proves feedback-pair
weight at least 150 for weight-one live updates in the actual 64-state
maps. Its restricted global classes pass, but full-length mixed trajectories
still require a sharper outer-shortening or output-rank count.
The [actual-support refinement](FEEDBACK_SUPPORT.md) now supplies both
and closes the full-length live class with alternating zero-output batches
at >49.1%/61 bits. Other trajectory classes remain unresolved.
The [changing-pair extension](CHANGING_PAIRS.md) permits up to32 arbitrary
nonzero intervening batches at the same restricted-class bound. It also
derives fresh two-state ranks and a within-coset MDS property of A_1 B.
Independent pair counting still loses cross-boundary constraints; no
whole-code result or construction change is claimed.

The [partial proof](PARTIAL_PROOF.md) now gives exact partial-state output
ranks, a separated state/kernel sampling description, and all-message bounds
for bounded-deficit trajectory classes at49% and49.1%. Larger-deficit
trajectories remain unresolved; these are not full distance certificates.

## Fixed construction and probability space

Keep F=GF(2^128), its E=GF256 subfield, the fixed RS[256,128] outer,
and the existing independent uniform permutations within 256 regions.
There are N=2K output symbols. Divide the routed input into batches X_t
of b=128 symbols. State Q_t has s=64 symbols, with Q_0=0 and no final-state
constraint. All diagonal coefficients below are sampled independently of
the route, once during setup, and fixed for every message. Distance must
hold for every nonzero message, including messages selected after setup.

## Build the maps from one polynomial interface

Use the existing Cantor ordering of E. Let Y be its first 128 points and
Z its next 64 points. For q in F^64, interpolate the degree-<64 polynomial
u with evaluations q on Z, and define A q as its evaluations on Y.
For y in F^128, interpolate the degree-<128 polynomial f on Y, and define
B y as its evaluations on Z. These are fixed E-linear maps, hence F-linear.

Polynomial uniqueness gives BA=I. Root counting gives the following facts
over F, not merely over E:

- im(A) is an [128,64,65] MDS code.
- ker(B) is an [128,64,65] generalized RS code.
- For q != 0, wt(q)+wt(Aq) >=129.
- For y != 0, wt(y)+wt(By) >=65.

Indeed, the first combined vector evaluates a degree-<64 polynomial on
192 distinct points. The second evaluates a degree-<128 polynomial there.
Every square minor of B is nonzero: it is the systematic parity map of
that [192,128,65] code. For fixed supports S of y and I of By, the space
of permitted y has dimension max(0, |S|+|I|-64). This follows by shortening
the same MDS code to S union I.

These maps need not use dense matrices. Write Delta(x)=x^2+x and
v_i=Delta composed i times. In the novel polynomial basis, split
f=u+v_6 v, with deg(u),deg(v)<64. On Z, v_6 is the constant
alpha=beta_1, which is neither zero nor one. Therefore:

- A: inverse 64-point transform on Z; zero-pad; forward 128-point transform.
- B: inverse 128-point transform; form u+alpha*v; forward 64-point transform on Z.

The existing TransformCantor routines support both cosets. This is an
algorithmic cost argument, not a native timing measurement.

## Iteration 1: retain state across zero-output batches

The initial proposal was Y_t=X_t+A D_t Q_{t-1}, Q_t=B Y_t.
A zero output batch immediately erases the state. This is a local failure
mode for a proof based on long active intervals, not an exhibited outer codeword.

Instead, keep the same maps and use

    Y_t = X_t + A D_t Q_{t-1}
    Q_t = Q_{t-1} + B Y_t.

Sample every diagonal entry of D_t from F minus {0,1}. The additional
computation is one XOR per state symbol. For empty input,
Q_t=(I+D_t)Q_{t-1}, so its support cannot shrink. For empty output,
Q_t=Q_{t-1}. For a fixed candidate output sequence, every Q_t is determined
without knowing D: Q_t=sum_{j<=t} B Y_j. Consequently the inverse-input
distribution remains accessible to a first-moment proof.

If r=wt(Q_{t-1})>0 and X_t=0, then wt(Y_t)>=129-r.
If a live state shuts down, B Y_t=Q_{t-1} and
B X_t=(I+D_t)Q_{t-1}. Hence wt(X_t)+r and wt(Y_t)+r are each at least65.
Also wt(X_t)+wt(Y_t)>=129-r. Neither input nor output individually must
be heavy when the state is heavy. From an empty state, an input in ker(B)
bypasses the state, but has at least65 nonzero symbols.

At K20, this uses 0.5 full-F multiplications per output and 16 MiB of raw
diagonal coefficients. The retained certified inner uses nine products per
output and 288 MiB of prepared keys. Small-field transforms and routing
still cost work; multiplication counts are not runtime predictions.

## Obstruction: active-state rank is not full-batch randomness

Fix the route and an output y. Put I_t=supp(Q_{t-1}), r_t=|I_t|, and
R=sum_t r_t. The inverse input lies in an affine translate of
V_I=direct_sum_t im(A[:,I_t]). Its random coordinates exclude both zero
and the corresponding Q entry. Each allowed inverse input has probability
(|F|-2)^(-R). If C is the routed outer code, a valid upper bound is

    Pr_D[inverse(y) in C] <= |F|^dim(C intersect V_I) / (|F|-2)^R.

An empty affine intersection has probability zero. Replacing the
intersection dimension by min(K,R) throws away the needed compatibility
and does not close the bound. Likewise, a separate union bound over
state-support patterns can spend roughly one extra bit per state symbol.
Near half distance there is little entropy slack. These are failures of
coarse bounds, not evidence of an actual low-distance codeword.

The required step is to bound the joint outer/trajectory constraints.
Equivalently, with outer parity check H, study the rank of the concatenated
columns H_t A[:,I_t], and retain affine compatibility when that rank drops.
The old certificate assumed uniform feedback in F^b and cannot be reused.

## Iteration 2: buy back full rank with two diagonal taps

Rather than insist on the smallest multiplication count, consider one
product per output. Let K0 q evaluate (v_6+alpha)u on Y, where u interpolates
q on Z. Define A2=A+K0. Then B A2=I, and [A,A2] is invertible.
To see this, im(K0)=ker(B), while BA=I; the two spaces are complementary.
Both A and A2 separately generate [128,64,65] MDS codes, because the
multipliers 1+alpha and alpha on the two halves of Y are nonzero.

Sample each pair (d,e) uniformly from {(d,e) in F^2: d+e !=1}.
Use independent pairs across state coordinates and batches, and set

    Y_t = X_t + A D_t Q_{t-1} + A2 E_t Q_{t-1}
    Q_t = Q_{t-1} + B Y_t.

On empty input, Q_t=(I+D_t+E_t)Q_{t-1}; support is preserved.
On empty output, state is unchanged. The maps can share their large
transform: interpolate Dq and Eq on Z to u and v; evaluate
u+(1+alpha)v+v_6*v on Y. Relative to the one-tap feedback, this adds one
64-point inverse transform, 64 small-field scalings, and 64 full-F products.
It does not require two independent 128-point forward transforms.

For fixed state with r nonzero coordinates, [A,A2] maps the 2r random
coordinates injectively. The feedback has exactly (|F|(|F|-1))^r possible
values, all equiprobable. When r=64, its distribution is uniform in F^128
conditioned on 64 affine inequalities. Its maximum density relative to
uniform is (|F|/(|F|-1))^64. This restores the old local uniformity argument
for fully supported states, with negligible density overhead, but does not
establish it for partially supported states. The inverse state is still
determined by y independently of both diagonal tapes.

For the two-tap candidate, replace V_I by the direct sum of spans of
A[:,I_t] and A2[:,I_t] in the inverse membership bound, and replace its
denominator by (|F|(|F|-1))^R. The unresolved cases now have partial state
support. A proof must count the constraints needed to create and sustain
those states, including empty and terminal intervals. Merely asserting
that typical states are dense would not prove all-message distance.

Raw setup is32 MiB at K20. No timing, finite certificate, or promotion
has been claimed for either new inner.

## Checks and next decision

`check_maps.py` compares novel-basis transforms with independent Lagrange
interpolation. It exhausts the 65,535 nonzero states and kernel words in
the GF16 analogue (b=8,s=4), checks all494 square minors of B, and checks
the actual GF256 b=128,s=64 maps. It also checks inverse recurrences,
the local zero-output example, and the complementary two-tap map.
These finite checks support implementation correctness; the root-counting
arguments above establish the stated local properties over F.

Next discriminating test: derive a joint support/rank enumerator for the
two-tap candidate, using the outer constraints when the state is partial.
Keep the one-tap design as the lower-cost comparison. If the proof requires
too much state or extra full-field arithmetic, reassess the recurrence
rather than tune a parameter sweep. Native timing should follow a viable
global bound, not substitute for one.
