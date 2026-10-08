# Regional proof for a vector-state observation inner

The vector-state inner has an exact two-state support law and a validated
regional rank envelope. Neither bound closes the K20 half-rate, 40%-distance
target in the tested classes. This note preserves the proof components and
the bounded failed screens; it is not a distance certificate or a timing claim.

Unlike the earlier [two-state candidate](HD_ADDITIVE_TWO_STATE.md), this
construction has no prefix before routing. Its state dimension is variable,
and its observation uses a fresh uniform vector. The earlier
[mixed-interval cut](HD_TWO_STATE_CUTS.md) concerns a different construction;
its numerical rejection cannot be transferred without rechecking the geometry.

## Construction and probability space

Fix a field F of size `q`, a GRS[b,k] constituent, and `R` outer rows.
The message dimension is `K=kR` and output length is `N=bR`. Encode the rows
independently with the fixed constituent. Collect each outer coordinate into
one region of length `R`. Independently sample a uniform permutation of the
row labels in each region, then concatenate the regions.

The inner state `S` belongs to F^r and starts at zero. At every position,
independently sample `d,f` from all of F^r and `c` from F excluding zero.
For input `x`, output and update are

\[
 y=cx+f\mathbin{\cdot}S,\qquad S'=S+dx.
\]

Keep every output and discard the terminal state. State carries cross all
region boundaries. All coefficients and permutations are sampled once and
then fixed for every message. The input/output map is triangular with
nonzero diagonal `c`, so it is invertible. No determinant conditioning is
imposed on a local input/state matrix.

The evaluated candidate uses `q=65536`, `b=128`, `k=64`, `R=16384`, and
`r=8`. Hence `K=2^20`, `N=2^21`, and rate is exactly one half. The forbidden
output weights for 40% distance are 0 through `ceil(2N/5)-1=838860`.
Eight parallel F lanes occupy each 128-bit record. A scalar-code distance
guarantee would imply the same nonzero-record distance, but the screens below
do not establish such a guarantee.

## Exact fixed-message law

Fix a message and the routed input values. A state flag records whether
the entire vector `S` is zero or nonzero. Let `0<z<=1` mark a nonzero output,
and define

\[
 a=q^{-r},\qquad g(z)=\frac{1+(q-1)z}{q}.
\]

The following matrices give exact output-weight moments. Rows and columns
index the old and new state flags, in order zero, live:

\[
 T_0=\begin{pmatrix}1&0\\0&g(z)\end{pmatrix},\qquad
 T_1=\begin{pmatrix}az&(1-a)z\\ag(z)&(1-a)g(z)\end{pmatrix}.
\]

The subscript indicates whether the input is zero or nonzero. A zero input
does not change the state. With a nonzero input, `S+dx` is uniform in F^r
and independent of the output, since `d` is independent of `c,f`.
When `S` is live, its inner product with uniform `f` is uniform in F.
When `S=0` and `x` is nonzero, `cx` is nonzero.

For a region containing exactly `h` live inputs, define its exact kernel

\[
 K_h(z)=\binom Rh^{-1}[t^h](T_0+tT_1)^R.
\]

The coefficient retains the order of noncommuting matrix products. A uniform
regional permutation produces a uniform `h`-subset of live positions.
Conditional on the occupancy profile `(h_1,...,h_b)`, the exact moment is

\[
 (1,0)K_{h_1}(z)\cdots K_{h_b}(z)\binom11.
\]

This formula carries state across regions. Replacing the product by separate
zero-initialized regional moments would describe a different encoder.

## Source counts and the regional majorant

Let `p=b-k` and `d_min=p+1`. For a fixed coordinate subset of size `w`, an
MDS code has the same exact-support count `a_w`, independent of the subset:

\[
 a_0=1,\qquad a_w=0\quad(1\le w<d_{\min}),\qquad
 a_w=(q-1)\sum_{j=0}^{w-d_{\min}}(-1)^j
       \binom{w-1}{j}q^{w-d_{\min}-j}\quad(w\ge d_{\min}).
\]

If `e_w` denotes the elementary symmetric polynomial in `b` variables, the
all-message occupancy polynomial is

\[
 P(t_1,\ldots,t_b)^R,\qquad
 P=\sum_{w=0}^b a_we_w(t_1,\ldots,t_b).
\]

Its coefficients count messages with each regional occupancy profile.
An independent coordinate shuffle inside every row is therefore unnecessary
for these all-message bounds. Such a shuffle makes a fixed message's support
uniform conditional on weight, but the MDS support counts already have the
required symmetry after summing messages.

For a computationally smaller bound, fix the number `A` of active outer
rows. On every feasible nonzero support, `a_w<=q^(w-p)`. Extending the
support sum to all subsets yields the following upper bound on the
probability that some message in this class has weight at most `D`:

\[
 B_A(z)=\frac{\binom RA q^{-pA}}{(q-1)z^D}
 (1,0)\left(\sum_{h=0}^A\binom Ah q^h K_h(z)\right)^b\binom11.
\]

The divisor `q-1` identifies nonzero scalar multiples, which always have
the same support under a fixed linear encoder. The binomial normalization
inside `K_h` remains regional; there is no global `binom(N,W)` factor.

The support extension includes infeasible low-weight rows. One valid
refinement marks total raw weight by `lambda>=1`, replaces `q^h` by
`(q*lambda)^h`, and multiplies by `lambda^(-d_min*A)`. Every actual class-A
word has at least `d_min*A` live coordinates. This refinement does not assume
independent row supports after conditioning on regional occupancies.

## A separate rank envelope

Counting all field-valued messages can overcount a rare surviving subspace.
The rank envelope instead fixes each outer support and tracks its surviving
linear source space. Initially, a feasible nonzero row support of size `w`
has dimension `w-p`. Introduce a real mark `1<=tau<=q`.

The conceptual history retains every scalar input form declared live and
every earlier whole-state map declared live. Discard a branch if a required
scalar or whole-state map becomes identically zero on the surviving space.
The two matrix indices do not replace this history. Keeping impossible
histories that pass these necessary checks can only increase the bound.

Output weight is expanded using

\[
 z^{\mathbf1[y\ne0]}=z+(1-z)\mathbf1[y=0].
\]

The first term imposes no output constraint. It must not impose output
liveness. The second term adds the zero-output linear equation.

For zero input and live old state, let `h>=1` be the rank of the state map
on the current source space. Uniform `f` gives a uniform functional in its
`h`-dimensional span. If `h=1`, every nonzero functional kills the required
old-state liveness; only probability `1/q` survives. If `h>=2`, a zero
functional has probability at most `q^-2`; every other functional drops
source rank by one. Consequently its expected rank factor is at most

\[
 G=\tau^{-1}+(1-\tau^{-1})q^{-2}.
\]

Charging exactly `tau^-1` here would be invalid. For independent state forms
`(u,v)`, the actual factor is `q^-2+(1-q^-2)/tau`, which is strictly larger.

For live input `X`, the observation admits the sharper charge `tau^-1`.
Condition on `c` and write `W=span(S_1,...,S_r)`. If `X` is outside `W`,
every output equation has rank one. If `X` lies in `W`, the equation is
uniform in `W`. Its `q-1` nonzero multiples of `X` kill required input
liveness. With `h=dim(W)`, the remaining expected factor is at most

\[
 q^{-h}+\frac{1-q^{1-h}}\tau\le\tau^{-1}.
\]

A return to zero imposes `S_j+d_jX=0` for all `r` coordinates. Expose the
fresh `d_j` sequentially. If `S_j` is independent of live `X`, the equation
drops rank one. If it is proportional to `X`, only one value of `d_j`
preserves input liveness. Each step therefore costs at most `tau^-1`.
The joint return cost is at most `tau^-r`.

For the return-and-zero-output entry, impose the return equations first.
Every surviving old-live history then has `S=-dX` with `d!=0`. Fresh `f`
preserves the live input under the output-zero equation exactly when
`f dot d=c`, with probability `1/q`. This order gives a stronger joint
bound. Constraint intersections and their final rank do not depend on order;
a discarded live history cannot be restored by further constraints.

Define `H0=z+(1-z)G`, `H1=z+(1-z)/tau`, and `Hr=z+(1-z)/q`.
The resulting positive rank envelope is

\[
 U_0=\begin{pmatrix}1&0\\0&H0\end{pmatrix},\qquad
 U_1=\begin{pmatrix}az&(1-a)z\\\tau^{-r}Hr&H1\end{pmatrix}.
\]

Its birth row is exact: when the old state is identically zero, `dX` is
zero exactly when `d=0`, without any source-rank restriction. The other
entries are upper bounds, not transition probabilities.

Replace `T0,T1` by `U0,U1` in the regional coefficient formula. The class-A
rank bound then replaces `q` by `tau` in the source factors and replaces
the projective divisor `q-1` by `tau`. A surviving nonzero witness leaves
source dimension at least one, so `tau^(remaining dimension)/tau>=1`.
This is a different counting argument from the fixed-message moment.
The same minimum-weight mark is valid.

## Bounded K20 screens

The scripts evaluated GF65536, state dimension eight, half rate, and 40%
distance. The table reports floating upper-bound margins; negative values
are failed proof bounds, not low-weight codewords.

| Outer | Active rows | Fixed-message envelope | Rank envelope |
|---|---:|---:|---:|
| [128,64] | 1 | 187.341 bits | 187.341 bits |
| [128,64] | 16 | -7,298.473 bits | -582.277 bits |
| [128,64] | 64 | -47,711.269 bits | -5,256.326 bits |
| [256,128] | 64 | not evaluated | -9,826.929 bits |

The bounded optimizers converged for the earlier, slightly weaker envelope.
The displayed rank values replay those fixed witnesses through the final
local tightening; no further optimization was run. At active counts 16 and
64, optimizing the minimum-weight mark returned `lambda=1`. The larger
outer control also returned `lambda=1` and remained negative.

The exact-message A=64 source, output, and regional terms are approximately
`-64952`, `12636`, and `100028` bits. Rank relaxation changes them to
`-9467`, `8620`, and `6102` bits. The large cancellation illustrates why
neither the full-field source count nor an isolated local return probability
is a sufficient success criterion.

The chosen marks are stored in `SCREEN_WITNESSES` in
[`hd_vector_region.py`](theory/hd_vector_region.py). Five exact tests check
the GF4 local law, regional chronology and normalization, the complete
positive operator with previous live histories, joint return charges, and
a small source majorant. An independent agent audited the rank argument.

```text
python -B research/workstreams/char2_field_spin/theory/hd_vector_region_tests.py
python -B research/workstreams/char2_field_spin/theory/hd_vector_region.py --replay
```

No directed complete class cover, native implementation, or timing was
performed for this candidate. The new proof tools are preserved, but these
failed gates do not justify native work on this route. The next design
comparison should change the source/route structure or provide a materially
tighter joint rank description, rather than repeat the same local search.
