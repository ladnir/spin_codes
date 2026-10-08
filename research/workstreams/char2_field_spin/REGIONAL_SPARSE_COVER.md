# A complete sparse support cover and its missing rank transfer

The GRS[32,22] candidate still has no complete distance certificate. This note gives a cover of every sparse bad word, a conditional diagonal-rank bound, and an exact packet interface for strengthening that bound. The evaluated cover is too loose. Its failure is not a construction lower bound.

The dense-input theorem in [REGIONAL_GLOBAL_ROUTE.md](REGIONAL_GLOBAL_ROUTE.md) is unchanged. The remaining class has nonzero inner input weight below \(0.24N\) and output weight at most \(e=\lceil N/10\rceil-1\).

The current ensemble uses independent regional permutations, independent nonzero GF256 input diagonals, and a four-symbol state. A regional orientation is sampled once per region, not per packet. The archived inner and the anti-diagonal repair are distinguished below.

## A support cover that includes every word

Write the packet recurrence as

\[
 Y_t=X_t+A_tQ_t,\qquad Q_{t+1}=Q_t+C_tX_t,
 \qquad C_tA_t=0.
\]

Each packet has sixteen coordinates. Every nonzero state emits at least thirteen nonzero coordinates. Every nonzero vector in \(\ker C_t\) has weight at least five. These statements hold for both the archived inner and the repaired inner.

For one word, let \(I_t,J_t\) be its input and output supports. Let \(b_t\) indicate whether its incoming state is nonzero. Then every packet satisfies the following necessary conditions:

* If \(b_t=0\), then \(X_t=Y_t\), so \(I_t=J_t\).
* If \(b_t=1\), then \(|I_t\cup J_t|\ge13\).
* A nonzero packet with \(b_t=b_{t+1}=0\) has \(|I_t|=|J_t|\ge5\).
* Both a birth \(0\to1\) and a return \(1\to0\) have nonempty input and output supports.
* A zero-output packet preserves the state: \(Q_{t+1}=Q_t\).

The last identity follows from \(X_t=A_tQ_t\) and \(C_tA_t=0\). Equivalently, the state recurrence can use outputs:

\[
 Q_{t+1}=Q_t+C_tY_t.
\]

Starting from zero, the sequence \((b_t)\) consists of complete state excursions and at most one unfinished final excursion. If there are \(c\) complete excursions, \(u\in\{0,1\}\) unfinished excursions, and \(z\) nonzero zero-state-to-zero-state packets, then

\[
 2c+u+5z\le\operatorname{wt}(Y).
\]

If the input and output weights are \(W,E\), the number of packets with nonzero incoming state is at most \((W+E)/13\). Every actual word defines such a cover. Unlike a rectangular-window family, this cover permits arbitrary row supports, multiple excursions per region, and excursions across boundaries.

## Exact conditional rank charge from input diagonals

Fix a route, the regional orientations, and an exact input support \(I\). The shortened outer space has dimension

\[
 D(I)=\sum_{\text{outer rows }r}\max\{0,|I\cap r|-10\}.
\]

If a row meets \(I\) in between one and ten coordinates, the exact support is impossible and is discarded. Otherwise write \(\ell_i\) for each required nonzero coordinate form on the shortened outer space.

Maintain one common span \(\mathcal W\) of the zero equations revealed so far. A branch is dead if any required \(\ell_i\) lies in \(\mathcal W\). This is a conservative test: passing it need not guarantee simultaneous nonzero values over a finite field.

At a coordinate \(i\in I\setminus J\), the zero-output equation is

\[
 d_i\ell_i+c_i=0.
\]

Here \(d_i\) is uniform in \(\mathbb F_q^*\), independently of the current history. The form \(c_i\) is fixed before \(d_i\) is revealed because packet emission uses the pre-update state.

On an alive branch, \(\ell_i\) is nonzero in the quotient by \(\mathcal W\). If \(c_i\) is independent of \(\ell_i\), every diagonal value adds one rank. If \(c_i=\beta\ell_i\), every value except possibly \(d_i=\beta\) kills the required form. The sole possible survivor has probability at most \(1/(q-1)\).

Therefore, for \(1\le\tau\le q-1\), the alive potential \(\tau^{-\operatorname{rank}\mathcal W}\) contracts by at most \(1/\tau\) per selected coordinate. A surviving nonzero witness has equation rank at most \(D(I)-1\). Thus

\[
 \Pr_D[\text{a word realizes }I,J]
 \le \min\{1,\tau^{D(I)-1-|I\setminus J|}\}.
\]

The common equation span must be retained across packets. This is a rank/existence argument, not a union over \(q^{D(I)}\) words. An independent audit checked the two quotient cases and the reveal order. Exhaustive GF4 tests cover 720 local cases.

For a complete finite upper bound, sum over the support covers above and average the clipped conditional bound over routes. The clipping belongs inside that average. The numerical relaxation below drops this clipping and several state constraints; it cannot be promoted to a certificate.

## Evaluation of the complete but relaxed cover

Put \(g=\tau^{-1}\). Mark input and output weights by \(x,y\). The coefficient of a support pair includes \(g^{|I\setminus J|}\). Define

\[
 A(x,y)=\sum_{u=13}^{16}\binom{16}{u}(gx+y+xy)^u.
\]

A positive two-state transfer that covers every packet is

\[
 T_{00}=1+\sum_{w=5}^{16}\binom{16}{w}(xy)^w,
 \quad T_{01}=(1+xy)^{16}-1,
\]
\[
 T_{11}=A(x,y),\qquad
 T_{10}=A(x,y)-\sum_{u=13}^{16}\binom{16}{u}[(gx)^u+y^u].
\]

The final expression is implemented as a positive sum of support types, not by numerical subtraction. It excludes empty input or empty output at a return.

Route averaging uses the shortened dimension and the minimum outer support size eleven. A tilted Bernoulli comparison gives the valid relaxation

\[
 \mathbb E_P[\mathbf1_{\rm feasible}\tau^{D(I)}]
 \le (L+1)^{32}\lambda^{-\alpha N}
 \left[(1-\alpha)^{32}+
 \tau^{-10}\eta^{-11}(1-\alpha+\tau\lambda\eta\alpha)^{32}\right]^L,
\]

where \(\alpha=|I|/N\), \(\lambda>0\), and \(\eta\ge1\). The extra tilt \(\eta\) upper-bounds the indicator that a nonempty row support has size at least eleven. The same column-conditioning and arithmetic--geometric-mean argument as in the dense theorem makes this uniform over regional profiles.

Combining this relaxation with the transfer and optimizing numerically gives the following positive normalized logarithms:

| Inverse-input density | Packet-only cover | Constant-state empty-input runs |
|---:|---:|---:|
| .001 | +.07487 | +.024995 |
| .01 | +.09808 | +.09663 |
| .10 | +.16855 | +.16855 |
| .2399 | +.02645 | +.02645 |

These are exploratory upper-envelope values, with finite boundary factors omitted. Positive values do not establish decay. They neither prove that every parameter choice fails nor imply that the construction is bad. The large slack remains even after the complete support cover and the stronger diagonal rank charge are used.

## Constant-state runs: archived and repaired emission laws

For the archived maps, \(A_\beta=A_0T_\beta\). Write \(v(u)=(1,u,u^2,u^3)^T\). Direct calculation gives

\[
 v(u)^TT_\beta=(1+\beta u)^3
 v\!\left(\frac{u}{1+\beta u}\right)^T.
\]

Across four phases, zero occurs four times as an evaluation point, while sixty nonzero evaluation points are distinct. A nonzero cubic can vanish at zero and at two other available points, producing six zero coordinates. The exact minimum over a full empty-input cycle is therefore 58, not 61. The replay checks all 1770 choices of the two nonzero roots attaining 58.

The run screen retains one cubic root mask throughout each empty-input run. Its series sums masks over the sixty-one distinct points, respecting repeated phases. It permits a new mask at a region boundary. This improves the very sparse envelope but does not close the middle class.

The separate anti-diagonal repair uses the four-by-four reversal matrix \(R\):

\[
 C_\beta=T_\beta C_0,\qquad A_\beta=A_0RT_\beta.
\]

Since \(T_\beta^TRT_\beta=R\), emission rows become \(v(u+\beta)^TR\). Its sixty-four points are distinct, so a full unchanged-state cycle emits at least 61 symbols. This is a changed inner, not a new interpretation of the archived map. It requires the corresponding adjoint change because the local map is no longer symmetric.

## An output-first packet interface for the repaired inner

The remaining diagonals must stay fresh when the state update is analyzed. The following factorization enforces that chronology.

Fix the prior history and a packet support pair \(I,J\). At label \(z_j\), put

\[
 c_j=v(z_j),\qquad a_j=v(z_j)^TR.
\]

Let \(S\) be the incoming vector of state forms in the common quotient. First impose all zero-output equations. At \(i\notin I\cup J\), impose \(a_iS=0\). At \(i\in I\setminus J\), reveal \(d_i\) and impose \(d_i\ell_i+a_iS=0\). This phase reveals none of the diagonals indexed by \(I\cap J\).

In the resulting quotient, the next state is exactly

\[
 \boxed{S'=F_JS+\sum_{j\in I\cap J}c_jd_j\ell_j,
 \qquad F_J=I_4+\sum_{j\in J}c_ja_j.}
\]

Indeed, \(S'=S+\sum_jc_jY_j\), and the already imposed equations remove all terms outside \(J\). Substituting \(Y_j=d_j\ell_j+a_jS\) gives the formula.

Conditioning on survival of the first phase preserves independent uniform diagonals on \(I\cap J\), provided survival tests use only the revealed history. In particular, the tests must not condition on unrevealed output-nonzero forms. Those requirements may be tested jointly when their diagonals are revealed, or omitted conservatively. They cannot justify resampling the same diagonal.

Only \(|I\cap J|\le|J|\) fresh state updates remain. All current forms belong to a space of dimension at most \(4+|I\cap J|\). With at most 1.6 nonzero outputs per packet on average, the small-output packet cases are the natural first finite kernels.

For \(J=\varnothing\), \(F_J=I_4\): the state is unchanged. For \(|J|=1\), \(F_J\) is an invertible transvection because \(a_jc_j=0\). Invertibility is not true for arbitrary \(J\). For two labels, with \(g=(z_1+z_2)^3\),

\[
 \det F_J=1+g^2,
\]

which can vanish. A transfer must retain the actual \(F_J\), not replace it by an invertible map. The replay checks the packet identity in 128 exact GF256 cases, including zero-output packets.

### A direct rank charge from the missing state constraints

Suppose the active incoming state has form-rank \(r\). Put

\[
 w=|I|,\quad v=|J|,\quad a=|I\cap J|,
 \quad z=16-w-v+a.
\]

If \(z\ge4\), four independent emission rows force the state to zero, so the active-state branch is dead. Otherwise the \(z\) zero-input, zero-output equations add rank at least \(\max\{0,r+z-4\}\). This follows from the rank inequality for the product of their emission matrix with the state-form matrix.

After these deterministic equations, the \(w-a\) fresh diagonal equations each contract the alive rank potential by \(1/\tau\). The complete zero-output phase therefore contracts by at most

\[
 \boxed{\tau^{-[(w-a)+\max\{0,r+z-4\}]}
       =\tau^{-\max\{w-a,r+12-v\}}.}
\]

For \(r=4\), this is \(\tau^{-(16-v)}\): **every zero output contributes a rank charge**, including positions where the input is zero. The earlier cover charged only \(w-a\) positions.

The surviving state-form rank is at most \(\min\{r,4-z\}\) before applying \(F_J\) and the remaining fresh updates. Thus the new charge and the subsequent state-rank transition must be kept together; charging full rank repeatedly without updating the state rank would be invalid.

### Zero-output and one-output rank transitions

Zero-output packets leave the quotient state unchanged after their equations have been imposed. For one output at label \(w\), put \(F=I_4+c_wa_w\), and let \(U\) be the image of the state-form map after the zero-output phase. Write \(r=\dim U\).

The identity \(a_wc_w=0\) gives \(F^2=I_4\) and \(Fc_w=c_w\). Therefore:

* If \(w\notin I\), the update is \(FS\), so its rank is exactly \(r\).
* If \(\ell_w\) is independent of the state row span, the new image is \(\operatorname{span}(FU,c_w)\). Its rank is \(r+1\) unless \(c_w\in U\), in which case it is \(r\).
* If \(\ell_w\) belongs to the state row span, the rank remains \(r\), except possibly for one nonzero diagonal value. That exceptional value can reduce rank by at most one.

For the last claim, write \(\ell_w=uS\). The update restricts \(F+d_wc_wu\) to \(U\). A rank drop requires \(c_w\in FU\); since \(Fc_w=c_w\), this is equivalent to \(c_w\in U\). The remaining scalar cancellation equation has at most one solution for \(d_w\).

When \(r<4\), at most \(r\) distinct injection labels have \(c_w\in U\), because any four distinct Vandermonde columns are independent. This is a count of exceptional labels, not a probability claim unless their sampling law is supplied.

These statements use the diagonal's unconditional fresh law after the zero-output phase. Requiring the surviving output form to be nonzero only removes outcomes; it is not imposed before invoking that law. The replay exhausts 4032 dependent-form configurations over GF4, including every nonzero diagonal, and also checks the independent-form and absent-input cases.

### Fresh input forms supplied by the shortened outer

There is deterministic freshness that need not be guessed by the state operator. Fix an exact row support of size \(r_i>10\). The shortened GRS row has dimension \(d_i=r_i-10\), and any \(d_i\) coordinate forms are independent.

Designate its first \(d_i\) active coordinates in stream order as fresh. Each such form is independent of all earlier input forms: earlier coordinates of the same row are independent, while other rows occupy direct-sum message spaces. The state and all prior zero equations involve only earlier input forms.

A region contains at most one coordinate from each source row. Thus the designated fresh forms in one packet are also mutually independent, and remain independent of the other current-row forms. Their exact uses are:

* A fresh input forced to zero output contributes exactly one equation-rank pivot, unless the required nonzero input makes the branch impossible. This relation does not reduce the rank of the old state forms.
* A fresh input at an allowed nonzero output remains independent of the state after the zero-output phase. It supplies the independent-form update above.
* If at least two such fresh allowed-output inputs survive in a packet, the following formal emission polynomial has no evaluation root. A candidate root would have to annihilate both distinct injection columns, requiring its label to equal both input labels.

The last statement uses independence of the fresh input forms: their coefficients cannot cancel against the old state or against each other.

Freshness is part of the row-support pattern. It can be counted by scanning outer coordinates backward: the last ten active coordinates of a nonzero row are dependent, and all earlier active coordinates are fresh. Use twelve backward states: counts zero through ten, followed by a state recording that a fresh coordinate has occurred. Accept only the empty-row state and the fresh-coordinate state. This supplies a finite support automaton; it is not permission to choose the fresh/dependent pattern independently in each packet.

### Exact impulse geometry for the remaining fresh updates

For the repaired map,

\[
 v(z)^TRv(w)=(z+w)^3.
\]

Writing the state polynomial as \(P(z)=v(z)^TRS\), a fresh update at label \(w\) changes it to

\[
 P'(z)=P(z)+d\ell(z+w)^3.
\]

This gives finite transition restrictions stronger than an arbitrary three-root mask:

* A one-symbol birth has exactly the distinct root \(w\).
* If \(\ell\) is independent of the current coefficient span, every formal root except possibly \(w\) disappears.
* If \(\ell\) belongs to that span and its rank is \(r\ge2\), the polynomial modulo \(\ell\) has rank \(r-1\). At most \(5-r\) evaluation labels can become new formal roots. Each such label, other than \(w\), corresponds to at most one value of the fresh nonzero diagonal.
* At rank one, write \(P=f\ell\). The transition is the explicit scalar cubic pencil \(f+d(z+w)^3\); a total cancellation produces zero state and must be classified separately.

The bound \(5-r\) uses the MDS property: a subspace of dimension \(5-r\) cannot contain more than \(5-r\) of these distinct evaluation rows. The unique-diagonal claim uses the nonzero coefficient \((z+w)^3\) for \(z\ne w\).

## Exact next proof obligation

For each packet output support \(J\), apply the zero-output phase, its rank charge, the actual matrix \(F_J\), and the remaining fresh cubic updates in that order. The transfer must retain the common zero-equation span through its state rank and formal-root information. It must keep the dependence between rank loss, new roots, and the same sampled diagonal.

The precise missing step is a finite upper operator for those joint transitions, coupled to the outer freshness automaton and uniform over the remaining dependent forms. Rank alone or root count alone is not yet proved sufficient state. The identities above reduce the low-output kernels to spaces of dimension at most \(4+|J|\), and give a direct charge absent from the failed envelope. They do not yet evaluate the resulting whole-code probability.

## Replay and scope

```sh
python -B research/workstreams/char2_field_spin/theory/regional_sparse_cover.py --check
python -B research/workstreams/char2_field_spin/theory/regional_sparse_cover.py --screen
python -B research/workstreams/char2_field_spin/theory/regional_sparse_cover.py --run-screen
```

The exact checks cover the diagonal rank step, actual packet covers, archived period emissions, and the repaired packet identity. The two numerical screens concern the earlier relaxed cover and remain explicitly uncertified. No native implementation, remote command, or benchmark is included.
