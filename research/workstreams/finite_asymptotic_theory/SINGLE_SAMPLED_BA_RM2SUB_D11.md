# Linear-time distance certificate for one sampled BA outer with RM2Sub-S19

## Certified statement

This note defines a scalable Structured SPIN variant. It is not a theorem
about the frozen \(B=256\) ParityFanout construction.

For each positive integer \(m\), set \(L_m:=128m\). Let \(B_m\) be the least
positive multiple of \(24\) that satisfies

\[
  B_m\ge \frac{39}{4}\log_2(L_mB_m),
  \qquad N_m:=L_mB_m.
  \tag{1}
\]

Then

\[
  B_m=\frac{39}{4}\log_2N_m+O(1),
  \qquad
  \frac{N_{m+1}-N_m}{N_m}=o(1).
  \tag{2}
\]

At outer length \(B=B_m\), take \(B/24\) copies of the extended binary
Golay \([24,12,8]\) code. Sample two independent uniform permutations of
\([B]\). Place one length \(B\) accumulator after each permutation. Denote
the resulting rate-half BA-3 code by \(\mathcal O_B\).

The setup samples \(\mathcal O_B\) once and reuses it at all \(L=L_m\)
outer-block positions. It then samples these independent objects:

1. one uniform coordinate permutation for each outer-block position;
2. one uniform permutation of \([L]\) for each of the \(B\) regions; and
3. one uniform nonzero \(\mathbb F_{2^{19}}\) multiplier for each RM2Sub
   epoch.

The inner uses the fixed RM2Sub-S19 maps with step length \(128\). Let
\(\mathcal C_m\) be the resulting binary linear code. The certificate proves

\[
  \Pr\!\left[
    d_{\min}(\mathcal C_m)\le\lfloor0.11N_m\rfloor
  \right]=o(1).
  \tag{3}
\]

The probability in (3) includes the two BA interleavers, all route
permutations, and all RM2Sub multipliers. The code has rate exactly one half.
Its ordinary and transposed encoders both use \(O(N_m)\) bit operations.
Setup sampling and storage also use \(O(N_m)\) operations and words. Since
the failure probability is \(o(1)\), a deterministic setup with all stated
properties exists for every sufficiently large admissible length.

## One sampled outer

Let \(A_{\mathcal O_B}(w)\) be the realized BA spectrum, and define

\[
  \overline A_B(w):=\mathbb E[A_{\mathcal O_B}(w)].
  \tag{4}
\]

The expectation in (4) is only over the two BA interleavers. The exact Golay
polynomial and the exact accumulator input--output enumerator determine
\(\overline A_B\).

The linear-time outer certificate proves an event \(\mathcal G_B\) with

\[
  \Pr[\mathcal G_B]=1-o(1).
  \tag{5}
\]

On \(\mathcal G_B\), every nonzero outer word has relative weight in

\[
  \mathcal W:=[0.104,0.896],
  \tag{6}
\]

and

\[
  A_{\mathcal O_B}(w)\le B^2\overline A_B(w)
  \quad\text{for every }w.
  \tag{7}
\]

The setup pays the failure probability in (5) once. Reusing the realized
code does not multiply that probability by \(L\).

## Certified BA spectrum majorant

Write

\[
  h(x):=-x\ln x-(1-x)\ln(1-x).
\]

Let

\[
  G(u):=1+759u^8+2576u^{12}+759u^{16}+u^{24},
\]

\[
  g(a):=\inf_{u>0}
  \left\{\frac1{24}\ln G(u)-a\ln u\right\},
\]

and, on the accumulator feasibility region, let

\[
  \pi(a,b)
  :=b h\!\left(\frac{a}{2b}\right)
   +(1-b)h\!\left(\frac{a}{2(1-b)}\right)-h(a).
\]

Set \(\pi(a,b):=-\infty\) outside that region. The limiting expected BA
spectrum exponent is bounded by

\[
  a_{\rm BA}(x)
  :=\sup_{0\le a,b\le1}[g(a)+\pi(a,b)+\pi(b,x)].
  \tag{8}
\]

`certify_golay_ba_concave_majorant.py` constructs a piecewise-affine
function \(\widehat a_{\rm BA}\) on \(\mathcal W\). It proves

\[
  a_{\rm BA}(x)\le\widehat a_{\rm BA}(x)
  \quad(x\in\mathcal W).
  \tag{9}
\]

The function is the lower envelope of affine supports. Therefore it is
concave. It uses 15 supports on each side of \(1/2\). The central support is

\[
  \widehat a_{\rm BA}(x)\le\frac{\ln2}{2},
  \tag{10}
\]

which follows exactly because the code has at most \(2^{B/2}\) words.
The outward verifier exhausted 101,910 boxes. Its largest accepted upper
endpoint was

\[
  -2.2000492744786192\times10^{-9}.
  \tag{11}
\]

Uniform coefficient and Stirling bounds applied to the exact finite BA
enumerator give

\[
  \frac1B\ln\overline A_B(w)
  \le
  \widehat a_{\rm BA}(w/B)
  +O\!\left(\frac{\ln B}{B}\right)
  \tag{12}
\]

uniformly for \(w/B\in\mathcal W\). The two finite weight sums contribute
only polynomial factors. The \(B^2\) factor in (7) is absorbed by the same
remainder.

## Why mixed row weights reduce to their mean

Fix \(Q\) active outer-block positions and set \(\alpha:=Q/L\). Suppose
their relative output weights are \(x_1,\ldots,x_Q\), and define

\[
  x:=\frac1Q\sum_{j=1}^Qx_j.
  \tag{13}
\]

The reduction needs more than Jensen's inequality. The region permutations
erase the row labels at subexponential cost.

For each active row \(j\), introduce independent Bernoulli-\(x_j\) bits and
condition their length \(B\) sum to equal \(Bx_j\). Since \(Bx_j\) is a
mode, the inverse conditioning probability is at most \(B+1\).

Now fix one region. Before its region permutation, the \(L\) bits are
independent Bernoulli variables with average success probability

\[
  q:=\alpha x.
  \tag{14}
\]

Let \(S\) be their sum. After a uniform region permutation, every vector of
weight \(k\) has probability \(\Pr[S=k]/\binom Lk\). Chernoff's inequality
and the method-of-types lower bound give

\[
  \Pr[S=k]\le e^{-L D_{\rm KL}(k/L\Vert q)},
\]

\[
  \Pr[\operatorname{Bin}(L,q)=k]
  \ge\frac1{L+1}e^{-L D_{\rm KL}(k/L\Vert q)}.
\]

Therefore the permuted region law is pointwise at most \(L+1\) times the
iid Bernoulli-\(q\) law. Across all \(B\) regions, the cost is

\[
  (L+1)^B=\exp(o(N)).
  \tag{15}
\]

The row-conditioning cost is at most

\[
  (B+1)^Q=\exp(o(N))
  \tag{16}
\]

whenever \(\alpha\) has a positive limit. Equations (15)--(16), (9), and
Jensen's inequality reduce all row-weight mixtures to
\(\widehat a_{\rm BA}(x)\).

Every routed array in this class has exactly \(N\alpha x\) one-bits. On that
support, the pointwise likelihood ratio between iid Bernoulli-\(q\) bits and
iid Bernoulli-(r) bits is

\[
  \exp\!\left(ND_{\rm KL}(\alpha x\Vert r)\right).
\]

Thus the RM2Sub transfer may use any reference probability \(r\), provided
the exponent pays this relative-entropy term.

## Positive occupation

Fix reference probabilities \(p,y\in(0,1)\). The reference bit is one with
probability \(r:=py\). The relative-entropy chain rule gives

\[
  D_{\rm KL}(\alpha x\Vert py)
  \le D_{\rm KL}(\alpha\Vert p)
  +\alpha D_{\rm KL}(x\Vert y).
  \tag{17}
\]

Let \(T_3(r,z)\) be the proved three-state RM2Sub transfer for iid
Bernoulli-(r) input. For fixed \((\alpha,x)\), it is enough to find
\(p,y,z\in(0,1)\) for which

\[
\begin{split}
  \Phi(\alpha,x;p,y,z)
  :={}&\alpha\widehat a_{\rm BA}(x)
   +D_{\rm KL}(\alpha\Vert p)
   +\alpha D_{\rm KL}(x\Vert y)\\
  &+\frac1{128}\ln\rho(T_3(py,z))
   -0.11\ln z<0.
  \tag{18}
\end{split}
\]

`certify_golay_ba_rm2sub_joint_interval.py` proves (18) for every

\[
  10^{-4}\le\alpha\le1,
  \qquad 0.104\le x\le0.896.
  \tag{19}
\]

Each box receives one rational \(p,y,z\) witness and one rational positive
Collatz vector. For a fixed affine spectrum support, the exponent is convex
in \((\alpha,q=\alpha x)\). Four vertex checks therefore certify each box.
The verifier accepted 621 boxes. Its largest outward upper endpoint was

\[
  -3.194048334767875\times10^{-7}.
  \tag{20}
\]

The quantifier order is \(\sup_{\alpha,x}\inf_{p,y,z}\). The witness is
fixed before any vertex of its box is checked.

## Fixed and growing sparse occupation

Define the BA likelihood excess relative to a random rate-half outer by

\[
  \ell(x)
  :=\widehat a_{\rm BA}(x)-h(x)+\frac{\ln2}{2}.
  \tag{21}
\]

On each affine segment, \(\ell\) is convex. It is therefore enough to check
the segment endpoints. Outward interval arithmetic proves

\[
  \sup_{x\in\mathcal W}\ell(x)<0.01281.
  \tag{22}
\]

For each active row, the symmetrized selected-outer counting measure is
therefore at most \(\exp((0.01281+o(1))B)\) times the random rate-half
counting measure. For occupation \(Q\), the product comparison costs at most
\(\exp((0.01281+o(1))BQ)\). This comparison uses the same route and RM2Sub
randomness as the random-outer certificates.

The uniform comparison in (22) is sufficient for \(Q=1,2\). For fixed
\(Q\ge3\), the weight-coupled transfer in
`WEIGHT_COUPLED_FIXED_OCCUPATION.md` retains a separate positive fugacity for
each BA weight segment. Its 62 outward endpoint checks prove a natural decay
coefficient greater than

\[
  0.07196093457399053
  \tag{23}
\]

per active row before choosing its position. With block constant \(39/4\),
the remaining fixed-occupation decay is greater than
\(0.0008689160550218037QB\). The certificates for \(Q=1,2\) retain larger
margins after the uniform transfer (22).

For \(Q\to\infty\) and \(Q/L\le10^{-4}\), the exact four-state random-outer
certificate gives natural decay coefficient

\[
  \kappa_0:=\frac{780897}{6665600}.
\]

Equation (22) leaves

\[
  \kappa_{\rm BA}
  :=\kappa_0-0.01281
  =\frac{86938833}{833200000}
  >0.10434329.
  \tag{24}
\]

The support choice costs at most \((4\ln2/39+o(1))QB\). Thus the remaining
natural margin is greater than

\[
  \kappa_{\rm BA}-\frac{4\ln2}{39}>0.0325484.
  \tag{25}
\]

The \(B^2\) selection factor per active row and the binomial local-limit
remainders contribute \(o(QB)\) when \(Q\to\infty\). For fixed \(Q\), they
contribute \(o_Q(B)\). `certify_golay_ba_rm2sub_sparse.py` records the exact
rational comparisons behind (22)--(25).

## Completion of the distance proof

Condition on \(\mathcal G_B\). Partition nonzero messages by their
occupation \(Q\).

- The fixed-occupation bounds cover every constant \(Q\).
- The four-state bound covers \(Q\to\infty\) with \(Q/L\le10^{-4}\).
- The joint interval certificate covers \(Q/L\ge10^{-4}\).

Every occupation sequence has a subsequence in one of these classes. The
fixed range is finite. The sparse range has a negative exponent on the
scale \(QB\). The positive range has a uniform negative exponent on the
scale \(N\). Summing over at most \(L\) occupation values preserves \(o(1)\).
Markov's inequality proves the conditional version of (3). Equation (5)
then proves (3) in the complete setup probability space.

## Linear-time integration

The companion linear-time certificate uses the same Golay--BA-3 outer
pipeline and the same fixed RM2Sub-S19 constituent. The present theorem
changes the admissible block schedule and the distance analysis. It does not
insert a dense linear map.

Each outer encoder applies constant-size Golay circuits, stored coordinate
permutations, and accumulators. Its transpose applies the transposed circuits
in reverse order. The route consists only of stored permutations. Each
RM2Sub epoch uses a fixed-size linear circuit and one field multiplication.
Thus the companion circuit decomposition applies without a new cost
assumption.

## Rate, cost, and arbitrary lengths

The BA-3 outer has dimension \(B/2\). The RM2Sub recurrence is invertible for
every fixed multiplier schedule. Hence \(\mathcal C_m\) has dimension
\(N_m/2\).

One BA encoder consists of constant-size Golay circuits, two permutations,
and two accumulators. It and its transpose cost \(O(B)\). Applying the same
encoder at all \(L\) positions costs \(O(N)\). The route and fixed RM2Sub-S19
step circuits also cost \(O(N)\) in either direction. Generating and storing
the two BA interleavers, the route permutations, and the field multipliers
uses \(O(N)\) setup work and storage.

For an arbitrary requested output length \(N\), take the next admissible
length \(N^+\). Restrict its input to \(\lfloor N/2\rfloor\) coordinates,
and puncture a fixed set of \(N^+-N\) outputs. Equation (2) gives
\(N^+-N=o(N)\). The wrapper therefore has rate \(1/2-o(1)\), relative
distance at least \(0.11-o(1)\), and \(O(N)\) ordinary and transposed work.
The wrapper adds no randomness.

## Scope

The proved logarithmic constant is \(39/4=9.75\). The fixed-occupation
margin is positive but small, and the common norm bound is tightest at
\(Q=3\). A smaller constant first requires a sharper \(Q=3\) transfer.
After that step, a weight-coupled four-state sparse transfer may become
limiting.

This theorem concerns the scalable Golay--BA-3/RM2Sub-S19 variant. It does
not authenticate the modeled spectrum of the frozen ParityFanout outer.
