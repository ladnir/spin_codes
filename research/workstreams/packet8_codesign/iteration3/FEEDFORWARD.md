# A same-state feedforward option

The H2-weighted 16-bit comparison has about 122 zero-to-birth transitions.
Strengthening those emissions is a distinct option from enlarging the
state or tightening the return coefficient. The construction below keeps
the original 16-bit byte-native A and C. No new numerical screen or encoder
implementation is claimed here.

## Construction and identities

For the eight-byte maps A(a,b)_h=a+hb and
C(X)=(sum X_h,sum hX_h), define F=I+AC. The verified identity CA=0 gives

\[
 F^2=I,\qquad FA=A,\qquad CF=C.
\]

The candidate step emits Y=FX+Aa and updates a'=Ma+CX. The update uses the
same fresh independent transitive randomizer as before. The state starts
at zero, remains continuous, and is not flushed.

The diagonal input map F is invertible, so the finite causal transform
remains invertible before the outer code is applied. This is a new
construction, not a coordinate relabeling of the original routed code:
F changes which bit weights are emitted for fixed sparse inputs.

For example, a nonzero input byte v at h=0 produces seven copies of v
under F, with the h=0 output byte zero. Its weight becomes 7 wt(v), rather
than wt(v). This example motivates a birth analysis; it is not a bound
for arbitrary inputs or source states.

## The transpose uses the existing two maps

Fix reverse input w and future adjoint state b. Let p=A^T w. The reverse
step is

\[
 x^*=w+C^T(b+p),\qquad b^*=p+M^Tb.
\]

Thus it still evaluates A^T once and C^T once. The old state b, not b+p,
must enter M^T. Using M^T(b+p) would implement a different recurrence.

In `spin/experiments/packet8_codesign/Fast.cpp`, the existing feedback
computation supplies p, and the existing expansion can consume b+p.
For two state bytes with two payload halves, forming b+p adds four ZMM
XORs. The algebra adds no GFNI application and changes no routed store.

The scheduling cost is real: current `pairStep` calls emit packets while
accumulating feedback. Feedforward needs the complete p before emitting
the first packet. One fixed-width implementation could retain all eight
packed packets, occupying 16 ZMM registers, then compute p and scatter.
Together with state and expansion temporaries, this approaches the
register budget. A reread/repack fallback adds 16 VBMI and 16 GFNI operations
per physical step. Neither schedule has been implemented or benchmarked.

The retained baseline uses 16 GFNI operations for packing, four for C^T
expansion, four for A^T feedback, and eight for the GL(2,GF256) update per
ordinary physical step. The proposed fold preserves these arithmetic
counts if it can retain the packed inputs without spills.
At the final reverse step, the old kernel omits feedback because no earlier
state needs it. Feedforward still needs p for that step's output; only this
boundary adds a feedback evaluation that the baseline omitted.

## What must change in the local proof

Fix occupancy j, a uniformly selected j-byte support, and independent
uniform nonzero input labels. For 0<z<=1, define

\[
 W^F_j(s)=E[z^{wt(FX)}1_{CX=s}],\qquad
 M^F_j(a)=E[z^{wt(FX+Aa)}].
\]

The old weighted feedback measure cannot be reused: feedforward changes
the weight attached to each birth. However, CF=C and the definition of F
give the exact zero-syndrome identity

\[
 W^F_j(0)=E[z^{wt(X)}1_{CX=0}]=W_j(0).
\]

Let L=65535, let U_j be the old emission moment averaged over nonzero
states, and let V_j=E[z^{wt(X)}]. Write B_j^F=sum_s W_j^F(s). Translating
the sum over all states by CX yields

\[
 L U^F_j=L U_j+V_j-B^F_j.
\]

This identity computes the new uniform-source row without enumerating
M^F_j at every state. If only a lower bound on B_j^F is available, it gives
an upper bound on U_j^F. In particular, dropping B_j^F is safe. Subtracting
an upper bound on B_j^F would reverse the required inequality.

For the remaining birth-family rows, define

\[
 J_j(s,b)=E[z^{wt(X+Ab)}1_{CX=s}].
\]

Then W_j^F(s)=J_j(s,s) and M_j^F(a)=sum_s J_j(s,a+s). These joint terms
retain the correlation between feedback and emitted weight. Neither
product of marginals nor a 16-bit denominator change supplies them.

## Cheap exact pieces

The one-byte input set has 2040 elements. Direct enumeration gives the
complete W_1^F measure and its total weight polynomial. Two-byte inputs
have 28*255^2=1820700 possibilities. Since every two-byte restriction of C
has rank 16, their feedback syndromes can also be accumulated directly.
Both censuses can stream into 65536 state bins; no 24-bit census is needed.

The exact one-byte birth family is supported on 2040 nonzero states.
Its transition into another one-byte step requires 2040^2=4161600 input
pairs. For a pair (X,Y), the relevant exponent is

\[
 wt(FX)+wt(FY+A CX).
\]

A positive histogram of these exponents, divided by 2040^2 and normalized
by the birth mass, gives the exact birth-to-emission moment for every z.
Together with the uniform-row identity, this supplies the three-coordinate
q1 operator without an all-state feedforward emission census.

Exact W_2^F is not the same as an exact q2 gate. The latter also needs
moments of dense two-byte birth families entering later occupied steps.
Naively pairing every two-byte input with every other input is not a
bounded strategy.

## A rigorous majorant for higher occupancies

An existing joint method can be adapted to eight-bit packets; the frozen
four-bit implementation must not be used by changing metadata alone.
Let P_j(s,w)=Pr[CX=s,wt(X)=w] for the literal byte-native feedback map.
Its exact character polynomial is

\[
 \widehat P_j(t,w)=
 \frac{[u^j y^w]}{\binom8j255^j}
 \prod_{h=0}^7
 \left[1+u\left((1+y)^{8-r_h}(1-y)^{r_h}-1\right)\right],
\]

where r_h=wt(C_h^Tt). Inverse Walsh transformation recovers P_j.
This distribution is independent of z and can be prepared once.

For a byte-weight profile r of a 64-bit word, let d(r,j,w) be the minimum
weight after adding an input with exactly j nonzero bytes and bit weight w.
A packet dynamic program computes d exactly. At a byte of weight r_h,
an inactive input costs r_h. An active input of bit weight ell in 1,...,8
costs |r_h-ell|. Track both total occupied bytes and total input weight.

The following are pointwise upper bounds because 0<z<=1:

\[
 W_j^F(s)\le\sum_wP_j(s,w)z^{d(profile(As),j,w)},
\]

\[
 M_j^F(a)\le\sum_{s,w}P_j(s,w)
 z^{d(profile(A(a+s)),j,w)}.
\]

The second expression is a sum of XOR convolutions over the 16-bit state.
Walsh transforms compute all entering-state bounds without enumerating
all input words. There are 261 feasible (j,w) pairs for j=0,...,8. Streaming
one occupancy at a time limits the largest state-by-weight array to
65536*57 entries. Integer character coefficients or outward arithmetic
are required for a certificate; floating clipping does not suffice.

These inequalities can retain exact one- and two-byte births and bound
higher births pointwise. Normalizing a pointwise upper birth measure by
its own total mass gives a valid dominated birth coordinate, injected
with that mass. A separately proved total cap must not be substituted by
rescaling the pointwise upper measure downward.

## Recommended bounded gate

First generate the exact one-byte and birth-pair histograms, verify the
transpose on tiny instances, and run the fresh q1 operator. Next combine
exact low-occupancy births with the joint packet-distance upper bounds
under the already valid capped-route geometry. This targets the observed
birth contribution without asserting that the full tail will improve.
The positive actual24 H3 gate remains the stronger current evidence.
