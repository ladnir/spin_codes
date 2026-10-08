# Blockwise field-SPIN: a cheap exact refresh law

## Status and construction interface

This note supplies a local random-state law for the blockwise restart in
[HD_BLOCKWISE_RESTART.md](HD_BLOCKWISE_RESTART.md). It does not certify distance
or report an implementation timing. The primary coefficient field is
\(\mathbb F=\mathrm{GF}(2^{16})\), of size \(q=65536\), not the 16-element field.
Distance is measured in nonzero field symbols. The same coefficients can act
on eight field lanes in each 128-bit record.

After the outer encoder and the single regional permutation, divide the input
into packets \(X_i\in\mathbb F^t\). Fix maps
\(A:\mathbb F^s\to\mathbb F^t\) and
\(C:\mathbb F^t\to\mathbb F^s\). Starting with \(S_0=0\), use

\[
 Y_i=X_i+AS_i,\qquad S_{i+1}=M_iS_i+CX_i.
\]

This is the blockwise SPIN order: emit from the old state, refresh that state,
and add the feedback of the **raw input**. There is no extra preaccumulator.
The state carries across regional boundaries, and the final state is free.
The maps \(A,C\) are fixed; independently sampled refresh maps \(M_i\) are
fixed for all messages after setup. They are independent of the regional
permutations. In particular, the local law below does not resample a code for
each encoded message.

The full input-to-output map is triangular with identity diagonal blocks.
Given the state and \(Y_i\), recover \(X_i=Y_i+AS_i\), then update the state.
Thus the full map is invertible for every choice of the \(M_i\), even without
assuming their invertibility. We use invertible refreshes to prevent a live
state from disappearing on a zero-input step. Neither \(CA=0\) nor
\(C=A^T\) is required.

For an output dual packet \(U_i\), the transpose recurrence is

\[
 V_i=U_i+C^Tr_{i+1},\qquad
 r_i=A^TU_i+M_i^Tr_{i+1},\qquad r_{\rm final}=0.
\]

If a refresh contains several rank-one factors, the transpose applies their
transposes in reverse order. This fixes the implementation interface before
choosing fast circuits for \(A\) and \(C\).

## An invertible rank-one refresh

Sample \(u\) uniformly from \(\mathbb F^s\setminus\{0\}\). Conditional on
\(u\), sample \(v\) uniformly subject to \(v^Tu\ne1\), and set

\[
 M=I+uv^T.
\]

The matrix determinant lemma gives \(\det M=1+v^Tu\ne0\). Its inverse is
\(I+uv^T/(1+v^Tu)\), since the field has characteristic two. This is an
invertible rank-one update; it is not generally a strict transvection, whose
pairing would be restricted to zero.

Fix a nonzero state \(z\). Put \(D=q^s-1\). Direct counting gives

\[
 \Pr[Mz=y]=
 \begin{cases}
 q^{s-1}/D,&y=z,\\
 1/D,&y\in\langle z\rangle\setminus\{0,z\},\\
 (1-1/q)/D,&y\notin\langle z\rangle,\\
 0,&y=0.
 \end{cases}
\]

For a target outside the line of \(z\), \(u\) must be one of the \(q-1\)
nonzero multiples of \(y-z\). The two forms \(v^Tu\) and \(v^Tz\) are
independent. Restricting the first to values other than one leaves the required
second value with probability \(1/q\). This gives \((q-1)/(qD)\).
For \(y=cz\), \(c\ne0,1\), the required pairing is \(c-1\), an allowed
value; its probability is \(1/(q-1)\), giving \(1/D\). The remaining
probability is the atom at \(z\).

There is a useful positive decomposition. On nonzero states, let \(I\) leave
the state unchanged, \(L\) replace it uniformly within its nonzero scalar
line, and \(J\) replace it uniformly among all nonzero vectors. Then

\[
 P=aI+(q^{-1}-a)L+(1-q^{-1})J,\qquad
 a=\frac{q^{s-1}-1}{q^s-1}.
\]

All coefficients are nonnegative. The three operators commute and satisfy
\(L^2=L\), \(J^2=J\), and \(LJ=J\). Consequently, \(r\) independent
refresh rounds, with no feedback addition between them, have the exact law

\[
 P^r=a^rI+(q^{-r}-a^r)L+(1-q^{-r})J.
\]

The total nonuniform mass is \(q^{-r}\). For the primary field it is
\(2^{-16}\) after one round and \(2^{-32}\) after two. These are mixture
weights, not setup-failure probabilities or whole-code security margins.
Over the binary field, each nonzero scalar line consists of one vector, so
\(I=L\) and the formula recovers the binary IMT lazy-refresh law.

### Why strict transvections need a different formula

If instead \(v\) is uniform in \(u^\perp\), then
\(T=I+uv^T\) is a strict transvection. For \(q>2\), \(Tz\) cannot be a
different nonzero scalar multiple of \(z\). Hence the binary expression
\(q^{-1}\delta_z+(1-q^{-1})J\) is false for this sampling rule.

Multiplying the strict transvection by an independent uniform
\(\lambda\in\mathbb F^*\) repairs the symmetry and gives
\(q^{-1}L+(1-q^{-1})J\). This costs an extra scalar-state multiplication.
The invertible rank-one update above already has a positive exact law without
that multiplication, so it is the cheaper candidate to investigate first.

For example, at \(q=4,s=2,z=(1,0)\), the invertible update assigns
probabilities \(4/15,1/15,1/20\) to \((1,0),(2,0),(0,1)\), respectively.
These unequal values explicitly rule out a fresh-uniform-state interpretation.

## Feedback and weighted output histories

Condition on the route, all earlier refreshes, a fixed message, and its old
state \(z\). Write \(c=CX_i\). The current emitted weight is the deterministic
quantity \(\operatorname{wt}(X_i+Az)\). The fresh random state is
\(Mz+c\), or the translated \(P^r\) law for several rounds.

If \(z=0\), the successor is deterministically \(c\). If \(z\ne0\), the
probability of returning to zero is exactly

\[
 a^r\mathbf1[c=z]
 +\frac{q^{-r}-a^r}{q-1}
   \mathbf1[c\in\langle z\rangle\setminus\{0\}]
 +\frac{1-q^{-r}}{q^s-1}\mathbf1[c\ne0].
\]

In particular, zero input cannot turn a live state off. A nonzero input can
cancel the refresh only through its actual feedback syndrome \(CX_i\).
The identity and line terms retain information about that syndrome and the
incoming state; they must not be replaced by an independent uniform state.

For a generating-function variable \(0<\theta\le1\), the exact local
operator multiplies by
\(\theta^{\operatorname{wt}(X_i+Az)}\), then applies the translated
refresh law. Earlier output weights generally bias the incoming state
measure. The fixed-state conditional law remains valid for that weighted
measure, but does not make the measure uniform.

The next proof task is therefore an emission-and-syndrome envelope for the
chosen fixed maps: bound weighted fibers of \(CX_i\), cancellations in
\(X_i+Az\), and the two nonuniform branches. The binary SPIN proof provides
the architecture, not a ready-made q-ary distance certificate.

### A fixed-vector law is not a quotient-rank law

One rank-one update cannot supply the same rank contraction as a fresh dense
uniform matrix. Take state forms equal to the coordinate basis and require
the new state to equal the old one. Their difference is \(uv^T\), of rank
at most one for every setup. Its rank is zero exactly when \(v=0\), with
probability

\[
 p_0=\frac1{q^{s-1}(q-1)}.
\]

Thus its rank moment is \(p_0+(1-p_0)/\tau\), not \(\tau^{-s}\).
The tests retain this counterexample. Any proof using source-space rank
charges needs a separate composable lemma; the near-uniform marginal law
alone does not justify such a lemma.

## Sampling, cost, and initial parameter choices

For exact setup sampling, fix an index \(j\) with \(u_j\ne0\). Sample the
other \(s-1\) coordinates of \(v\) uniformly and choose the pairing
\(d\) uniformly from \(\mathbb F\setminus\{1\}\). Solve the one linear
equation \(v^Tu=d\) for \(v_j\). This samples each permitted \(v\) once.
Alternatively, reject a uniform \(v\) when its pairing equals one; the
rejection probability is exactly \(1/q\). All inversions belong to setup.

Applying \(M\) requires a dot product and an axpy: \(2s\) field
coefficient-times-record products and \(2s-1\) additions per packet.
Applying a dense \(s\)-by-\(s\) matrix requires \(s^2\) products.
A strict transvection followed by independent scalar multiplication requires
\(3s\) products. These are arithmetic counts, not native instruction counts
or measured times; each record carries eight \(\mathrm{GF}(2^{16})\) lanes.

If \(A\) generates an MDS \([t,s]\) code and \(C\) is an MDS parity-check
map, their useful deterministic guarantees are as follows.

| Packet \(t\) | State \(s\) | Minimum nonzero \(Az\) weight | Minimum nonzero \(X\in\ker C\) weight | Refresh products per raw symbol |
|---:|---:|---:|---:|---:|
| 32 | 16 | 17 | 17 | 1 |
| 64 | 16 | 49 | 17 | 1/2 |
| 64 | 32 | 33 | 33 | 1 |

The square choices \(t=s\) have no feedback kernel, but their general
emission distance is only one. They should not be assumed superior merely
because their state is larger. Conversely, dense implementations of both
fixed rectangular maps cost \(2ts\) products per packet, overwhelming the
refresh cost. The intended experiment requires fixed fast MDS circuits,
with their actual spectra and costs, rather than dense maps justified only
by their minimum distances.

Maps over an embedded subfield can retain their rank and MDS minors over
\(\mathbb F\). This can make the fixed circuits cheaper. However, sampling
the refresh coefficients only from that subfield does **not** give the law
proved above on arbitrary \(\mathbb F\)-states. The stated law requires
full-field \(u,v\).

An explicit periodic schedule may refresh only once every \(\ell\) packets,
using \(M_i=I\) otherwise. That reduces refresh work, but the law applies
only at the refresh boundary. Intervening emissions and feedback additions
must be retained in the period operator. They cannot be replaced by
\(P^\ell\). Start with one refresh per packet; consider sparser refreshes
only after the complete bound exposes the available slack.

## Validation and recommendation

The independent small-field helper is
[theory/hd_blockwise_inner.py](theory/hd_blockwise_inner.py). Its
[eight tests](theory/hd_blockwise_inner_tests.py) pass, covering:

- all permitted-map inverses over the selected small fields;
- exact fixed-state image distributions, binary recovery, and q-ary
  counterexamples to a fresh-uniform claim;
- scalar-symmetrized strict transvections and two-round kernel composition;
- exact pivot-sampler multiplicities and feedback-return probabilities;
- literal streaming inverse and 243 transpose basis checks, including a
  prescribed identity-refresh packet and multi-round refresh;
- the failure of the unsupported full-rank-contraction interpretation.

Recommended next step: compare the blockwise \((t,s)=(32,16)\) and
\((64,32)\) interfaces, with \((64,16)\) as a cheaper-state control. First
derive their q-ary emission/syndrome envelope and apply the existing regional
geometry gates. Keep fixed-Cauchy maps with independent diagonal randomness
as a competing construction: this exact rank-one law is an additional
candidate, not evidence that IMT must beat the existing diagonal design.
No native implementation or full certificate is claimed by this checkpoint.
