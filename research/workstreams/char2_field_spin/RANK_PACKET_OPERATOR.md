# A finite packet rank operator

This note supplies a positive five-state operator for the GRS[32,22] sparse-input proof. It is uniform over the common equation history and the regional orientations. It does not assume that rank determines the exact transition law. Uncertain transitions receive separate upper bounds, which can lose mass but preserve unconditional entrywise dominance.

The construction has sixteen inputs per packet and four state coordinates over GF256. Input diagonals are independent and uniform in \(\mathbb F_{256}^*\). One regional permutation orders the source rows in each outer-coordinate region. The outer support determines which input forms are fresh; the operator does not sample freshness independently.

The default operator concerns the unscaled dual-MDS repair from [DUAL_MDS_PHASE.md](DUAL_MDS_PHASE.md). A separately labeled option concerns the changed emission law multiplied by a fixed \(\eta\notin\mathbb F_{16}\). It uses the all-output-support invertibility theorem in [DUAL_MDS_OBSTRUCTIONS.md](DUAL_MDS_OBSTRUCTIONS.md); it is not valid for the unscaled construction.

## State and interface

Fix an exact nonempty outer support and its shortened message space \(V\), of dimension \(D\). Maintain the common span \(\mathcal W\subseteq V^*\) of all imposed zero equations. Every required input form must remain nonzero modulo \(\mathcal W\); otherwise discard the branch.

State 0 means that the current state vector has been constrained to zero. State \(r\in\{1,2,3,4\}\) means that its form matrix has rank \(r\) modulo \(\mathcal W\), and the witness state is required to be nonzero. Formal rank being positive is only a necessary condition for that requirement.

A zero-state checkpoint imposes all next-state equations. It adds their actual quotient rank to \(\dim\mathcal W\). A checkpoint is not a request that the unconstrained state map already have rank zero.

For a fixed support pair \(I,J\), split its sixteen labels into six classes:

| Class | Count |
|---|---:|
| Fresh input, zero output | \(f_0\) |
| Dependent input, zero output | \(d_0\) |
| Fresh input, allowed output | \(f_1\) |
| Dependent input, allowed output | \(d_1\) |
| No input, allowed output | \(o\) |
| No input or output | \(z\) |

Put \(a=f_1+d_1\), \(w=f_0+d_0+a\), and \(v=a+o\). An allowed output belongs to the specified exact output support. The analysis may omit its nonzero test, but does not condition on that test before revealing its diagonal.

For \(1\le\tau\le255\), an operator entry bounds the expected surviving potential \(\tau^{-\Delta\dim\mathcal W}\) that reaches its destination state. This statement holds conditionally on each fixed prior history. Consequently, any nonnegative continuation vector may follow the operator; no monotonicity assumption on that vector is needed.

The implementation exports:

```python
packet_matrices(tau, y, normalized=False, law='unscaled')
coefficient_table(t=16, s=4, law='unscaled')
```

The first function returns \(K[f,d,r,k]\), summing all fresh/dependent placements and output supports, with output mark \(y^v\). Here \(f=f_0+f_1\) and \(d=d_0+d_1\). Setting `normalized=True` divides by \(\binom{16}{f}\binom{16-f}{d}\), giving the uniform-placement average \(N_{f,d}\). The exact table has integer multiplicities indexed by \((f,d,v,r,k,c,j)\), representing \(y^v\tau^{-c}255^{-j}\).

## The zero-output phase

Suppose the incoming branch has positive rank \(r\). Four zero-input, zero-output labels force the state to zero by the emission MDS property. Thus the branch is dead when \(z\ge4\).

Otherwise let \(u\) be the state rank after the \(z\) deterministic root equations. Their equation rank is exactly \(r-u\), and

\[
 \max(1,r-z)\le u\le\min(r,4-z)
\]

on any surviving incoming-nonzero branch. Each designated-fresh zero-output input adds one equation without reducing the old-state rank. Each dependent zero-output equation reduces that rank by at most one.

If the rank after all zero-output equations is \(h\), then

\[
 \max(1,r-z-d_0)\le h\le\min(r,4-z).
\]

The alive diagonal lemma gives one factor \(\tau^{-1}\) for every input in \(I\setminus J\). For each possible \(h\), a valid entrywise charge is therefore

\[
 c(r,h):=f_0+d_0+r-\min(r,4-z,h+d_0).
\]

To justify the last term, the actual \(u\) is at most \(\min(r,4-z,h+d_0)\). Bounding each final-rank event by the full alive potential can count the same upper mass several times. This is an explicit relaxation, not an assumption about its distribution.

## Allowed outputs and checkpoints

After the zero-output phase, the exact quotient update is

\[
 S'=F_JS+\sum_{j\in I\cap J}c_jd_j\ell_j,
 \qquad F_J=I_4+\sum_{j\in J}c_ja_j.
\]

All diagonals in this sum remain independent and uniform. This remains true after conditioning on the zero-output phase, because that phase uses only other diagonals.

For \(v=0\), the next state equals the current state. A positive-state branch cannot return.

For \(v=1\), \(F_J\) is an invertible transvection fixing the sole injection column. If there is no input at the allowed output, the rank remains \(h\) and return is impossible. A designated-fresh input gives rank \(h\) or \(\min(h+1,4)\). An arbitrary dependent-designation input gives those same possibilities, plus rank \(h-1\) with probability at most \(1/255\). The designation “dependent” does not imply dependence on the current state: its form can instead be independent of the state span.

There is a stronger restriction when \(v=1\) and \(z\ge1\). A deterministic root row annihilates the post-zero state image but not the allowed input's injection column, since every distinct-label pairing is nonzero. Therefore the fresh input raises rank exactly to \(h+1\). A dependent-designation input gives \(h\) or \(h+1\), with no rank-loss exception. The code uses this restriction.

For unscaled packets with \(v\ge2\), the operator deliberately uses only

\[
 k_{\min}=\max\{0,\min(4,f_1),h-v-d_1\},
 \qquad k_{\max}=\min(4,h+a).
\]

Indeed, \(F_J\) loses at most \(v\) ranks; each dependent input update loses at most one. The designated-fresh forms are independent of all remaining current forms and prior forms, so their columns guarantee rank at least \(\min(4,f_1)\). Every positive integer in this interval receives a separate unit upper coefficient.

For a return, select \(m=\min(a,4)\) remaining injection columns. They are independent. Reveal all other remaining diagonals and left-invert the selected columns in the equation \(S'=0\). This produces \(m\) equations

\[
 d_i\ell_i+\text{offset}_i=0
\]

whose offsets contain no selected diagonal. Applying the alive diagonal lemma sequentially gives an extra factor \(\tau^{-m}\). Additional state equations can only reduce the potential.

There is also a pointwise next-state rank lower bound \(k_{\min}\), or \(h\) in the fresh one-output case. Imposing the checkpoint gives its corresponding factor. These two arguments use overlapping randomness and equations: the operator charges their **maximum**, not their sum.

A return with exactly one output requires all sixteen inputs to be nonzero. Indeed, the incoming state is \(c_jY_j\); every off-label cubic pairing is nonzero, and the on-label input equals \(Y_j\). The scaled law also requires at least thirteen inputs when there are two outputs. If their ratio is in GF16, both on-label inputs survive and the cubic has at most three off-label roots. If their ratio is outside GF16, no off-label input vanishes; invertibility prevents both on-label inputs from vanishing. Both filters are included.

At an incoming zero state, \(I=J\). A nonempty zero-to-zero packet has at least five active inputs and receives the four-pivot factor \(\tau^{-4}\). A zero-to-positive packet has rank between \(\max(1,\min(4,f_1))\) and \(\min(4,w)\). The empty packet stays at zero.

## Optional scaled emission law

Under the separate `eta-scaled` law, every \(F_J\) is invertible. Thus the initial allowed-output step preserves \(h\), even when \(v\ge2\). A return with \(I\cap J=\varnothing\) is impossible and is excluded.

Reveal the \(d_1\) remaining dependent-designation diagonals first. For one update from rank \(r>0\), the conservative transition is

\[
 r\longrightarrow r\text{ or }\min(r+1,4)
 \quad\text{with coefficient }1,
 \qquad
 r\longrightarrow r-1\quad\text{with coefficient }1/255.
\]

From rank zero, a required nonzero input form gives rank one. The code composes this positive transition \(d_1\) times. It then adds all \(f_1\) independent fresh forms. If the intermediate rank is \(h'\), the final rank lies between \(\max(h',\min(4,f_1))\) and \(\min(4,h'+f_1)\).

The return charge uses the deterministic lower bound \(\max(\min(4,f_1),h-d_1)\) together with the pivot charge, again taking their maximum. No general large-field arithmetic appears in this proof operator. The changed encoder's cost and all-\(J\) invertibility require their separate construction audit.

## A Bellman action without uncertain-rank duplication

The entrywise matrix can lose substantial mass by counting several possible ranks separately. The additional function `packet_action(tau,y,V,normalized=True)` avoids that loss for the scaled law. It accepts any nonnegative continuation vector \(V\), not only a monotone one. Its output is nonlinear in \(V\), so it must not be treated as a matrix and exponentiated.

For a nonnegative rank payoff \(B\), one dependent allowed-input update has the valid upper action

\[
 (\mathcal D B)_0=B_1,
 \qquad
 (\mathcal D B)_r=M_r+\frac{(B_{r-1}-M_r)_+}{255},
 \quad
 M_r=\max\{B_r,B_{r+1}\},
\]

where the nonexistent entry \(B_5\) is omitted. Outside the exceptional rank-loss event, only the two ordinary ranks are possible. The exceptional event has probability at most \(1/255\). Its contribution is therefore needed only when its payoff exceeds both ordinary payoffs. This proves the action for arbitrary \(B\).

For \(f\) independent fresh allowed inputs, use

\[
 (\mathcal F_f B)_r=
 \max_{\max(r,\min(4,f))\le k\le\min(4,r+f)}B_k.
\]

The updates are exposed in the opposite order to this backward action: dependent diagonals first, fresh inputs last. Thus \(\mathcal D^{d_1}\mathcal F_{f_1}B\) bounds their combined expected payoff.

To include both next-state branches, put

\[
 G_0=0,\quad G_k=V_k\ (k>0),\qquad
 B_k=G_k+\tau^{-k}V_0.
\]

The second term imposes the zero-state checkpoint with its actual rank charge. When \(a>0\), the allowed-output continuation is bounded by the minimum of

\[
 \mathcal D^{d_1}\mathcal F_{f_1}B
 \quad\text{and}\quad
 \mathcal D^{d_1}\mathcal F_{f_1}G+
 \tau^{-\min(a,4)}V_0.
\]

Each expression is independently valid; taking their minimum does not multiply overlapping charges. When \(a=0\), invertibility excludes return and the continuation is simply \(V_h\).

The zero-output phase also uses a maximum, not a sum. For each possible deterministic root rank \(u\), it takes the largest allowed-output continuation over \(\max(1,u-d_0)\le h\le u\), multiplied by \(\tau^{-(r-u+f_0+d_0)}\). It then maximizes over \(\max(1,r-z)\le u\le\min(r,4-z)\). This is valid because the total alive diagonal potential already bounds every one of these possible final ranks jointly.

The incoming-zero row retains \(I=J\), the minimum return weight five, and the four-pivot return charge. Its uncertain positive ranks are also maximized. Summing these actions over support placements and output supports gives \(\mathcal N_{f,d}(V)\).

For any positive vector \(V\) and positive marks \(u,v\), define

\[
 \rho=\max_{f,d,r}
 \frac{\mathcal N_{f,d}(V)_r}{u^fv^dV_r}.
\]

This is a common supersolution, uniform over all conditional histories. Across \(T\) packets, the terminal factor is at most \(\rho^T V_0/\min_rV_r\), times the category marks. No state-rank Markov assumption is needed. This action is intended to improve the evaluated matrix envelope, not to assume away its failure.

## Optional fresh orientation boundaries

The additional argument `fresh_orientation=True` changes the construction: the packet starts with a fresh independent uniform GL4 orientation. It must not be used at later packets under an already exposed orientation.

Condition on the entire preceding common equation history and the outer support. Before the new orientation is sampled, these determine the state-form map and all source forms. The orientation makes its rank-\(r\) image a uniform \(r\)-subspace of \(\mathbb F_q^4\), while preserving its row span.

For \(z\) distinct zero-input, zero-output labels, the root-row matrix has rank \(z\). If \(u\) is the state rank after imposing these root equations, its exact law is

\[
 p_u=
 \frac{{4-z\brack u}_q{z\brack r-u}_q
 q^{(4-z-u)(r-u)}}{{4\brack r}_q}.
\]

The brackets denote Gaussian binomial coefficients. The action discards \(u=0\), because the incoming witness state is required nonzero. For each other \(u\), it charges \(\tau^{-(r-u)}\), then processes the input-only zero equations as above. It sums these conditional bounds with weights \(p_u\), replacing the plain maximum over \(u\).

This is the only use of the fresh Grassmann law. Subsequent zero-output survival tests can bias the orientation. The remaining updates use bounds uniform over that conditioned history.

For a refresh every four packets, the backward continuation action is

\[
 \mathcal A_{\rm GL}\bigl(\mathcal A_{\rm plain}
 (\mathcal A_{\rm plain}(\mathcal A_{\rm plain}(V)))\bigr).
\]

A common positive supersolution bounds this whole four-packet action. The implementation does not replace it with a product of averaged rank matrices. The rational reference `packet_action_exact` also implements the fresh-orientation law.

## Whole-support composition

Each nonempty shortened outer row with support size \(w>10\) has \(w-10\) designated-fresh inputs followed by ten dependent inputs in stream order. The exact backward row automaton retains that ordering. A regional permutation samples the resulting three categories without replacement.

For positive marks \(u_j,v_j\), take the entrywise envelope

\[
 M_j=\max_{f,d}\frac{N_{f,d}}{u_j^fv_j^d}.
\]

For each fixed row-support pattern, successive packet category counts have hypergeometric laws. Their probabilities sum to one. Thus each region contributes at most its category mark times \(M_j^{L/16}\). No independent-freshness approximation or regional-profile union bound is needed.

Let \(H\) be the exact row-support polynomial evaluated with fresh marks \(\tau xu_j\) and dependent marks \(xv_j\). For the class \(W_{\min}\le\operatorname{wt}(X)\le W_{\max}\), the resulting unclipped bound is

\[
 \frac{H^L-1}{\tau x^{W_*}y^e}
 e_0\left(\prod_{j=1}^{32}M_j^{L/16}\right)\mathbf1,
 \qquad
 W_*=
 \begin{cases}W_{\max},&x\le1,\\W_{\min},&x>1.\end{cases}
\]

Here \(0<y\le1\) and \(e=\lceil N/10\rceil-1\). The subtraction removes the empty support before the matrix envelope is applied. This is a whole-class upper expression, not a union over individual field-valued words. It can still be loose because it drops clipping and duplicates uncertain rank transitions.

The independent implementation is [outer_freshness_transfer.py](theory/outer_freshness_transfer.py). A floating evaluation alone is a screen; a successful probability claim requires outward-rounded evaluation and sufficient coverage of all remaining weights.

## Executed screens

The first exact-freshness, entrywise-matrix composition was evaluated for the entire remaining K20 input interval \(1\le W\le289894\), with \(L=47680\), \(N=1525760\), and \(e=152575\). Its floating log upper envelope was approximately \(+0.71684N\). This proves no decay. It exposed substantial loss from the monomial profile envelope and uncertain-rank duplication.

The Bellman diagnostics retain the input-support count and the existing uniform-profile row bound, while conservatively treating all source forms as dependent. They remove uncertain-rank duplication. The following values use the independently optimized two-state refresh witnesses, then solve for a positive Bellman continuation vector:

| Fixed inverse-input density | Refresh every 16 symbols | Refresh every 64 symbols |
|---:|---:|---:|
| .01 | \(+.02108654N\) | not replayed at a tuned witness |
| .10 | \(-.00286359N\) | \(+.06162835N\) |
| .18 | \(-.12595539N\) | \(-.02668884N\) |

These are floating upper-envelope logarithms, with finite boundary factors omitted. Negative values at individual densities do not establish interval coverage. Positive values neither reject the construction nor prove that every parameter choice fails. The gain from the rank state over the two-state refresh screen is modest at these witnesses; the remaining regional-profile relaxation is now a concrete bottleneck.

The next proof step is a count-aware route bound that retains the canonical fresh/dependent row pattern. Another local transition identity alone will not supply that missing averaging theorem.

## Replay and scope

```sh
python -B research/workstreams/char2_field_spin/theory/rank_packet_operator_tests.py
python -B research/workstreams/char2_field_spin/theory/rank_packet_operator.py
python -B research/workstreams/char2_field_spin/theory/rank_packet_operator.py --diagnostics --refresh-packets 0 1 4
python -B research/workstreams/char2_field_spin/theory/rank_packet_operator.py --diagnostics --refresh-packets 1 4 --refresh-tuned
```

The exact tests check 1,536 GF4 support/quotient cases, enumerating every nonzero diagonal and testing \(\tau\in\{1,2,3\}\). Their four-symbol, two-state analogue has the exact same dual-MDS identities. A further 512 matrix/form cases check the exceptional rank-loss law over every GF4 diagonal, including nonmonotone continuation payoffs. A genuine scaled GF16 packet with GF4 labels checks every diagonal and output mask in three support classes against the exact rational continuation \((3,1,7)\). All 180 matrices in GL2(GF4) verify the Grassmann intersection law. Additional tests check homogeneity, floating/rational agreement, improvement over the matrix operator, and the scaled zero-overlap exclusion. These tests supplement the conditional proof; they do not establish the GF256 construction theorem by enumeration.

Both matrix tables are generated from exact integer multiplicities. The return-pivot proof, unscaled rank bounds, and fresh-orientation conditioning have independent mathematical audits. Neither the matrices nor the current screens establish the missing low-density distance certificate.
