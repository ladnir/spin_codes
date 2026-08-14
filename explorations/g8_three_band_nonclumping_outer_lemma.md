# A nonclumping outer lemma for the fixed three-band layout

## Purpose and scope

The current full-support diagnostic loses the incidence pattern of the outer
message blocks.  Its outer columns depend only on the nine packet-weight
counts.  The fixed three-band layout contains more information: every band
tile contains 64 blocks, a block pair meets in at most one band, and the layout
contains no Pasch configuration.  This note gives a conditional-Finner lemma
that retains the first two facts.  It also identifies the additional local
majorant needed to use the four-block facts.

The lemma concerns the unpunctured data layout.  A final outer witness must add
a multivariate correction for the 128 graph replacements.  One valid choice
adds
\[
128\max_{|j-k|\leq1}\log_2(t_k/t_j)
\]
to the base-two outer constant.  Scalar support corrections from the star
certificates are not interchangeable with this packet-class correction.
Binary64 evaluation of the lemma is a discovery diagnostic until an outward
implementation replays every frozen parameter.

## Fixed incidence structure

Let

\[
I:=\mathbb Z_{256}\times\{0,\ldots,63\}
\]

index the 16,384 data-message blocks.  For
\(b\in\{0,1,2\}\), define slopes

\[
(s_0,s_1,s_2):=(0,9,20)
\]

and tile maps

\[
\tau_b(t,\ell):=t+s_b\ell\pmod {256}.
\]

For a band tile \((b,u)\), define

\[
E_{b,u}:=\{i\in I:\tau_b(i)=u\}.
\]

The committed certificate establishes the following properties.

1. Every set \(E_{b,u}\) has size 64.
2. For distinct bands \(b,c\), every intersection
   \(E_{b,u}\cap E_{c,v}\) has size at most one.
3. No four blocks form two disjoint colliding pairs in every band.
4. The numbers of five-edge four-block motifs are

   \[
   (M_0,M_1,M_2)=(79616,122880,79360),
   \]

   where \(b\) identifies the band containing the single colliding pair.

The script `scripts/certify_three_band_tile_map.py` verifies these statements
by exhaustive integer computation.

## Local tile moments

Fix positive packet fugacities
\(t=(t_0,\ldots,t_8)\), with \(t_0=1\).  Define

\[
R_w(t):=
\frac{[z^w]\left(\sum_{j=0}^8 {8\choose j}t_jz^j\right)^8}
{{64\choose w}},
\qquad 0\leq w\leq64.
\]

This is the exact expected packet monomial after a uniform permutation of 64
binary lanes whose total weight is \(w\).  For each band tile \((b,u)\), let

\[
F_{b,u,t}:(\mathbb F_2^{64})^{E_{b,u}}\longrightarrow\mathbb R_{\geq0}
\]

be its packet-fugacity factor after averaging the committed within-band
coordinate permutations.  The unpunctured outer moment is

\[
Z_\tau(t):=
\sum_{x\in(\mathbb F_2^{64})^I}
\prod_{b=0}^2\prod_{u\in\mathbb Z_{256}}F_{b,u,t}(x|_{E_{b,u}}).
\]

This equation defines the interface needed below.  A diagnostic evaluator may
construct the factors from the fixed-band projection spectra.  A certificate
must reconstruct them from the committed code and layout.

For \(S\subseteq I\), condition on

\[
x_i\ne0\quad\Longleftrightarrow\quad i\in S.
\]

The remaining variables \((x_i)_{i\in S}\) are independent and uniform in
\(\mathbb F_2^{64}\setminus\{0\}\).  Define the tile occupancy

\[
r_{b,u}(S):=|S\cap E_{b,u}|.
\]

For \(0\leq r\leq64\), define the uniform local envelope

\[
m_{b,r}(t):=
\max_{\substack{u\in\mathbb Z_{256},\ A\subseteq E_{b,u}\\|A|=r}}
\mathbb E\!\left[
F_{b,u,t}(X_A,0_{E_{b,u}\setminus A})^3
\right].
\]

The expectation is over independent uniform nonzero values \((X_i)_{i\in A}\).
The maximum makes the definition valid even when the averaged local factor is
not block-symmetric.

The committed projection spectra do not directly evaluate this cube moment.
Surjectivity nevertheless gives an explicit upper envelope.  Let
\(c_0:=42\) and \(c_1=c_2:=43\).  Let \(R_w(t)\) denote the exact packet
moment of one 64-lane column of weight \(w\), and define

\[
A_d(t):=\sum_{w=0}^d {d\choose w}R_w(t)^3.
\]

For \(c\in\{42,43\}\), define

\[
\overline m_{c,r}(t):=
\frac{1}{(2^{64}-1)^r}
\sum_{d=0}^r(-1)^{r-d}{r\choose d}
2^{(64-c)d}A_d(t)^c.
\]

To derive this identity, fix the coordinate permutations before taking the
cube.  Apply inclusion--exclusion to the condition that all \(r\) inputs are
nonzero.  If \(d\) inputs remain unrestricted, surjectivity gives
\(2^{(64-c)d}\) fibers, and the \(c\) projected coordinates are independent.
Thus, \(\overline m_{c,r}\) is the exact cube moment of the fixed-permutation
tile monomial.  Jensen's inequality gives

\[
m_{b,r}(t)\leq \overline m_{c_b,r}(t).
\]

The alternating sum requires high-precision arithmetic.  At unit fugacity,
\(A_d=2^d\), and inclusion--exclusion gives
\(\overline m_{c,r}=1\) exactly.

## Conditional-Finner lemma

**Lemma 1 (nonclumping conditional-Finner bound).**  For every positive
fugacity vector \(t\),

\[
Z_\tau(t)
\leq
\sum_{S\subseteq I}(2^{64}-1)^{|S|}
\prod_{b=0}^2\prod_{u\in\mathbb Z_{256}}
\overline m_{c_b,r_{b,u}(S)}(t)^{1/3}.
\tag{1}
\]

At \(t=(1,\ldots,1)\), both sides equal \(2^{64|I|}=2^K\).

**Proof.**  Fix \(S\subseteq I\).  Conditional on the displayed activity
event, each active block variable occurs in exactly one factor from each band.
Finner's inequality with coefficient \(1/3\) therefore gives

\[
\mathbb E\!\left[
\prod_{b,u}F_{b,u,t}(X|_{E_{b,u}})
\,\middle|\,
X_i\ne0\Longleftrightarrow i\in S
\right]
\leq
\prod_{b,u}\overline m_{c_b,r_{b,u}(S)}(t)^{1/3}.
\]

Exactly \((2^{64}-1)^{|S|}\) message assignments have active set \(S\).
Summing the conditional inequality over all \(S\) proves (1).  When all
fugacities equal one, every local factor and every local envelope equals one.
The binomial theorem then gives the normalization identity.  \(\square\)

The ordinary three-band Finner bound applies Finner before exposing \(S\).
Lemma 1 permits an incidence-aware summation after the local cube moments are
known.

## Exact occupancy state

For a set \(S\), define

\[
h_{b,r}(S):=
|\{u\in\mathbb Z_{256}:r_{b,u}(S)=r\}|,
\qquad 0\leq r\leq64.
\]

For every band \(b\),

\[
\sum_rh_{b,r}=256,
\qquad
\sum_rrh_{b,r}=|S|.
\tag{2}
\]

Define the number of colliding block pairs in band \(b\) by

\[
T_b(S):=\sum_r {r\choose2}h_{b,r}(S).
\]

Pair capacity gives the exact constraint

\[
T_0(S)+T_1(S)+T_2(S)\leq {|S|\choose2}.
\tag{3}
\]

For one fixed band and one histogram \(h_b\), the number of sets with that
histogram is

\[
N_b(h_b)=
\frac{256!}{\prod_rh_{b,r}!}
\prod_r{64\choose r}^{h_{b,r}}.
\tag{4}
\]

Consequently, the number of sets with a specified histogram triple is at most
\(\min_bN_b(h_b)\).  Equations (1)--(4) give a finite exact relaxation with
195 nonnegative integer histogram coordinates.  The relaxation uses exact
band balance and pair capacity without assuming random tile labels.

## A coefficient form of the pair-capacity refinement

The histogram sum can be collapsed further.  Fix \(s\in\{0,\ldots,|I|\}\)
and \(\eta\geq0\).  For each band define

\[
P_{b,\eta}(y):=
\sum_{r=0}^{64}{64\choose r}
\overline m_{c_b,r}(t)e^{-3\eta {r\choose2}}y^r.
\]

**Corollary 2 (pair-capacity coefficient bound).**  The contribution to (1)
from active sets of size \(s\) is at most

\[
(2^{64}-1)^s
\inf_{\eta\geq0}
e^{\eta {s\choose2}}
\prod_{b=0}^2
\left([y^s]P_{b,\eta}(y)^{256}\right)^{1/3}.
\tag{5}
\]

**Proof.**  For every active set of size \(s\), equation (3) implies

\[
1\leq e^{\eta({s\choose2}-\sum_bT_b(S))}.
\]

Multiply its summand in (1) by the right-hand side.  Apply Hölder across the
three bands.  In one band, the 256 tile sets partition \(I\), so the resulting
sum is the coefficient of \(P_{b,\eta}^{256}\).  Taking the infimum over
\(\eta\) proves (5).  \(\square\)

For a cheaper upper bound, replace each coefficient by

\[
[y^s]P(y)^{256}\leq y_0^{-s}P(y_0)^{256},
\qquad y_0>0.
\]

This replacement reduces one evaluation to three 65-term log-sum-exp
computations.

## Four-block state

Pair occupancy does not encode the no-Pasch certificate.  Define
\(P(S)\) as the number of four-block subsets that split into two colliding
pairs in every band.  Define \(Q_b(S)\) as the number of five-edge four-block
subsets whose single-pair band is \(b\).  The committed layout gives

\[
P(S)=0,
\qquad
0\leq Q_b(S)\leq M_b.
\tag{6}
\]

Thus a motif-aware exact state is

\[
\sigma(S):=
\left(
|S|,(h_{b,r})_{b,r},P(S),Q_0(S),Q_1(S),Q_2(S)
\right).
\]

Equation (1) may be summed over these states, with (2), (3), and (6) imposed.
This statement uses all committed incidence facts.  It does not yet yield a
smaller computable bound because the local Finner weight depends only on the
histograms.  A numerical use of (6) requires a motif-sensitive local
majorant; omitting that requirement would silently discard the main counting
problem.

One suitable interface is a pair-interaction majorant.  For each band, seek
\(a_b\geq0\) and \(\rho_b\geq0\) such that

\[
\overline m_{c_b,r}(t)^{1/3}
\leq
\overline m_{c_b,0}(t)^{1/3}a_b^r(1+\rho_b)^{{r\choose2}}
\quad\text{for }0\leq r\leq64.
\tag{7}

Expanding the right-hand side produces a gas of colored collision edges.
Pair capacity makes the edge colors disjoint.  The six-edge Pasch coefficient
vanishes.  The five-edge four-vertex contribution is bounded using
\((M_0,M_1,M_2)\).  A theorem using this expansion must also bound every
connected term with at least seven collision edges.  A tree-graph or polymer
remainder is appropriate only when the fitted \(\rho_b\) values satisfy its
convergence condition.

## Expected leverage at the current worst cell

The 72 vertices of `h2:073` have physical weights between 720,910 and
1,245,174.  The graph word contributes at most 128 nonzero coordinates.
Every realizing outer word therefore has at least 5,632 to 9,727 nonzero data
blocks, even when every active block has weight 128.

This density has two consequences.

1. Pair capacity alone is unlikely to remove hundreds of thousands of bits.
   For \(s>190\), the tile-capacity maximum
   \(\sum_bT_b(S)\leq94.5s\) is already below \({s\choose2}\).
   Equation (3) can still help tilted states, but it is not a dense-regime
   obstruction.
2. The no-Pasch and five-edge counts can matter only when the fitted local
   interaction weights are large enough to favor simultaneous clumping in
   several bands.

A reasonable discovery expectation is 0--30,000 bits from (5).  A convergent
motif correction may plausibly add 10,000--100,000 bits.  Neither estimate is
a theorem.  A gain near the current 580,000-bit residual would require a much
stronger local outer lemma or additional construction randomness.

## Bounded diagnostic plan

The first diagnostic should not optimize a new outer witness.  It should use
the current worst-leaf dual barycenters and their frozen packet fugacities.

1. Reconstruct \(R_w(t)\) and verify the committed band-surjectivity
   certificates.  The current projection tables do not directly serialize
   the required nine-fugacity cube moments.
2. Compute the 195 values \(\overline m_{c_b,r}(t)\) from the displayed
   inclusion--exclusion formula.
3. Check \(\overline m_{c_b,r}(1)=1\) and the global normalization
   \(Z_\tau(1)=2^K\).
4. Sum (5) over every \(s\in\{0,\ldots,16384\}\).  A packet profile does not
   determine \(s\).  Use at most 64 fixed \(\eta\)-values and coefficient
   Chernoff bounds.  Start with one leading leaf, then extend to all eight if
   the measured gain justifies the work.
5. Fit the least pair majorants (7).  Stop the motif branch if a certified
   polymer convergence condition fails.
6. If convergence holds, enumerate connected collision motifs through six
   edges.  Insert the exact zero Pasch count and the three exact five-edge
   counts.
7. Compare the resulting frozen affine outer prices with the conditioned-row,
   total-spectrum, and BL2 prices.  Subtract the packet-profile normalization
   once in the component replay.

The state tables contain only `3 * 65` local rows.  A Chernoff coefficient
evaluation costs `O(3 * 65)` operations per `(s,eta)` pair.  Scanning every
\(s\leq16384\) and 64 fixed tilts requires about 200 million elementary
log-domain updates per leaf.  The first probe may evaluate one leaf before all
eight.  It must not restrict the active-set sizes unless a separate tail bound
covers every omitted size.

The diagnostic artifact should bind the manifest, geometry checkpoint,
component catalogues, projection spectra, tile-map certificate, and all local
moment tables.  It must label (5) separately from any motif correction and
record whether the polymer convergence gate passed.

## Bounded diagnostic result

The complete active-set scan was run once for h2:073.  It used every
\(s\in\{0,\ldots,16384\}\) and eight fixed collision tilts.  The optimized
log-sum equals the zero-tilt log-sum:

\[
\log_2 Z_{\mathrm{pair}}=
\log_2 Z_{\eta=0}=390432.47233220655.
\]

Thus, pair capacity saves zero bits in the complete sum.  The support-size-one
term saves about \(0.0022\) bits, but dense support sizes dominate.  The
artifact is out/g8_pair_capacity_nonclumping_h2_073.json, with SHA-256
151ae23c9ec820d873137eac72e621fa3e4ed962a777308052edbe0a91a5e553.

The singleton-matched pair majorants have
\(\rho_{42}\approx0.01818\) and \(\rho_{43}\approx0.01862\).  Simple
uncentered and centered polymer criteria both fail.  The corresponding gate
values are approximately \(9.6\) and \(19.3\), while each criterion requires
a value at most one.  Therefore, the present proof does not certify a
six-edge motif remainder.  The current nonclumping branch stops here.
