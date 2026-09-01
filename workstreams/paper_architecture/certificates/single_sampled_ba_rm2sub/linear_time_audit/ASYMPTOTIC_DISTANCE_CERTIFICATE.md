# Asymptotic linear-time distance certificate

## Certified statement

This note gives a complete asymptotic certificate for one Structured SPIN
family. The family has rate one half, relative distance at least \(0.101\),
and linear ordinary and transposed encoding work.
`CERTIFICATE_MANIFEST.json` binds the proof receipts, verifiers, and frozen
dependencies by SHA-256 hash.

Write

\[
  h(x):=-x\ln x-(1-x)\ln(1-x)
\]

for binary entropy in natural units, with the usual continuous values at the
endpoints.

For each positive integer \(m\), define

\[
  B_m:=24m,
\]

\[
  L_m:=128
  \left\lceil
    \frac{\exp(\sqrt{B_m})}{128}
  \right\rceil,
  \qquad
  n_m:=B_mL_m.
\tag{1}
\]

Then

\[
  B_m\longrightarrow\infty,
  \qquad
  \frac{B_m}{\ln n_m}\longrightarrow\infty,
  \qquad
  B_m=o(n_m).
\tag{2}
\]

The divisibility conditions in (1) align the length-24 outer constituent and
the length-128 inner epoch. They are not distance assumptions.

For each sufficiently large \(m\), the setup procedure below produces a
binary linear code

\[
  \mathcal C_m(\mathsf{Setup}_m):
  \mathbb F_2^{n_m/2}\to\mathbb F_2^{n_m}
\]

with the following properties:

\[
  \Pr_{\mathsf{Setup}_m}
  [d_{\min}(\mathcal C_m(\mathsf{Setup}_m))<0.101n_m]
  =o(1),
\tag{3}
\]

\[
  C_{\mathcal C_m}^{\rm enc}
  =O(n_m),
  \qquad
  C_{\mathcal C_m}^{\mathsf T}
  =O(n_m).
\tag{4}
\]

The probability in (3) is over the route and inner randomness defined below;
the two local BA interleavers have already been fixed by the deterministic
selection argument. In particular, (3) implies that, for every sufficiently
large \(m\), at least one deterministic setup realization has the same rate,
distance, and cost bounds.

## Selected linear-time outer

Let \(\mathcal G_{24}\) be the extended binary Golay \([24,12,8]\) code. For
block length \(B=B_m\), start with the direct sum of \(B/24\) copies of
\(\mathcal G_{24}\). Apply two accumulator stages, each preceded by a
selected interleaver on \([B]\). Denote the resulting BA-3 code by
\(\mathcal O_B\).

`LINEAR_OUTER_CERTIFICATE.md` proves that the interleavers can be selected so
that every nonzero outer word has relative weight in

\[
  \mathcal W=[0.104,0.896]
\tag{5}
\]

and the spectrum \(A_B(w)\) satisfies

\[
  \max_{w:A_B(w)>0}
  \frac1B
  \ln\frac{A_B(w)2^B}{\binom Bw}
  <0.36+o(1).
\tag{6}
\]

The proof uses the exact Golay polynomial and the exact accumulator
input--output enumerator. It treats sparse active Golay blocks separately
from linear-density active blocks. The outward-rounded results are recorded
in `OUTER_INTERVAL_RECEIPT.json`.

The ordinary circuit for \(\mathcal O_B\) consists of constant-size Golay
circuits, two permutations, and two accumulators. Its transpose reverses
these circuits. Therefore

\[
  C_{\mathcal O_B}^{\rm enc},
  C_{\mathcal O_B}^{\mathsf T}=O(B).
\tag{7}
\]

The complete outer is the direct sum of \(L=L_m\) copies of the same selected
code \(\mathcal O_B\). It has dimension \(LB/2=n_m/2\).

## Route and structured inner

The setup algorithm samples the following independent objects:

1. For each of the \(L\) outer-block positions, a uniform permutation of its
   \(B\) output coordinates.
2. For each of the \(B\) transposed regions, a uniform permutation of its
   \(L\) coordinates.
3. For each inner epoch, a scalar
   \(\alpha_i\gets\mathbb F_{2^{19}}^*\).

The first permutations transpose the \(L\)-by-\(B\) outer array into \(B\)
regions of length \(L\). The second permutations randomize positions inside
each region.

The inner uses the frozen maps

\[
  P:\mathbb F_2^{128}\to\mathbb F_{2^{19}},
  \qquad
  A:\mathbb F_{2^{19}}\to\mathbb F_2^{128}.
\]

For epoch \(i\), it computes

\[
  Y_i:=X_i+A(Q_i),
  \qquad
  Q_{i+1}:=\alpha_iQ_i+P(X_i),
  \qquad
  Q_0:=0.
\tag{8}
\]

The frozen receipts prove

\[
  P\circ A=0,
  \qquad
  \min_{q\ne0}\operatorname{wt}(A(q))=48,
  \qquad
  d_{\min}(\ker P)=6.
\tag{9}
\]

They also give the complete weight spectra of \(\ker P\) and
\(\operatorname{im}A\).

For fixed setup randomness, (8) is triangular in the epochs: knowing the
past inputs determines \(Q_i\), and then \(X_i=Y_i+A(Q_i)\). Thus the inner
is invertible. The complete code therefore retains the outer dimension
\(n/2\).

## Distance proof

Fix a nonzero message. Let \(a\) be the number of active outer-block
positions, and set

\[
  \tau:=a/L.
\]

The proof partitions messages by the asymptotic behavior of \(a\) and
\(\tau\).

### Bounded and vanishing occupation

Equation (6) implies the outer spectrum envelope

\[
  a_{\rm out}(\beta)
  \le h(\beta)-\ln2+0.36
  \qquad(\beta\in\mathcal W).
\tag{10}
\]

Choose the live-time threshold

\[
  \rho:=0.27>
  0.101\cdot\frac{128}{48}.
\]

The one-active exponent in `ONE_ACTIVE_UPPER_BOUND.md` satisfies

\[
  \chi(0.27,2^{19}-1)<-1.30.
\tag{11}
\]

Combining (10) and (11) gives the certified margin

\[
  \kappa
  \ge
  \ln2-0.36-e^{-1.30}
  >0.0606.
\tag{12}
\]

`DISTANCE_CONSTANTS_RECEIPT.json` records the outward-rounded verification of
(11)--(12). Since \(B/\ln n\to\infty\), the one-active union bound is
\(o(1)\). `FIXED_OCCUPATION_UPPER_BOUND.md` extends the same conclusion to
every bounded occupation.

If \(a\to\infty\) and \(\tau\to0\), then

\[
  at/L=t\tau\longrightarrow0.
\]

`SPARSE_GROWING_OCCUPATION.md` applies in this regime. Equation (12) and
\(B/\ln n\to\infty\) make its total first moment \(o(1)\).

### Positive occupation

`BULK_OCCUPATION_TRANSFER.md` defines the exact full-state transfer and a
rigorous two-state envelope. The envelope does not assume a uniform live
state. It uses the ambient intersection numbers

\[
  J_{d,k,\ell}
  =
  \binom d{(k+d-\ell)/2}
  \binom{128-d}{(k-d+\ell)/2}
\]

and the frozen kernel spectrum to bound every kernel-shell intersection.
Thus the frozen marginal receipts and \(P\circ A=0\) suffice; no joint
kernel-shift receipt is assumed.

For \(0<\tau\le10^{-3}\), choose

\[
  z(\tau):=1-\frac{11}{5}\tau,
  \qquad
  \xi(\tau,v):=\eta(v)\tau,
\tag{13}
\]

where \(v=x/\tau\). `verify_small_tau_interval.py` constructs a positive
two-state norm witness for each certified \(v\)-interval. It proves that the
bulk exponent divided by \(\tau\) is at most

\[
  -0.00826.
\tag{14}
\]

For \(10^{-3}\le\tau\le0.985\),
`verify_middle_interval.py` fixes one rational \(z\)-witness on each
occupation slab. It then covers the coefficient-density range with rational
boxes and permits one \(\xi\)-witness per box. This quantifier order matches
the bulk gate:

\[
  \inf_z\sup_x\inf_\xi.
\]

Every box has a strictly negative outward-rounded upper endpoint.
`MIDDLE_INTERVAL_RECEIPT.json` records all four covered ranges and the hashes
of both verifiers.

### High occupation

Let

\[
  d_{1/2}:=\limsup_{B\to\infty}\frac{D_B(1/2)}B.
\]

Equation (6) gives \(d_{1/2}<0.36\). The direct comparison with a uniform
input closes every occupation satisfying

\[
  \tau>
  \frac{h(0.101)}{\ln2-0.36}.
\tag{15}
\]

Outward rounding gives

\[
  \frac{h(0.101)}{\ln2-0.36}
  <0.982373<0.985.
\tag{16}
\]

Thus (16) overlaps the interval proof. At full occupation, the negative
exponent margin is greater than \(0.00587\).

### Union over occupations

The first moment is partitioned by the integer occupation \(a\). Every
asymptotic sequence of occupation values has a subsequence in one of the
following classes:

- bounded \(a\);
- \(a\to\infty\) and \(\tau\to0\);
- \(\liminf\tau>0\).

The cited fixed- and sparse-occupation bounds already include the choice of
active outer positions. Their negative exponent is on the \(aB\) scale and
is uniform once \(B/\ln n\to\infty\). The positive-occupation verifier has a
strict uniform margin on each of its four ranges; its exponent is on the
\(n\) scale. Hence summing the at most \(L+1\) occupation levels preserves an
\(o(1)\) total first moment. Equivalently, a failure of this conclusion would
select an occupation sequence and then one of the three subsequences above,
contradicting its corresponding uniform bound. Markov's inequality proves
(3).

## Cost proof

There are \(L=n/B\) local outer encoders. Equations (7) and (2) give

\[
  L\,O(B)=O(n)
\]

ordinary and transposed outer work. The coordinate route and its transpose
apply \(O(n)\) stored indices.

The inner parameters \((t,s)=(128,19)\) are fixed. The frozen ordinary and
transposed step circuits therefore have constant cost. Repeating them
\(n/128\) times costs \(O(n)\). Setup stores \(O(n)\) permutation indices and
\(O(n/128)\) field scalars. This proves (4).

## Scope

The certified distance is \(0.101\). The rate-one-half GV value is
approximately \(0.110027\). The current conservative transfer does not
certify \(0.102\): its middle-density and high-density bounds leave a gap near
full occupation. This limitation concerns the proof envelope, not an upper
bound on the true distance of Structured SPIN.
