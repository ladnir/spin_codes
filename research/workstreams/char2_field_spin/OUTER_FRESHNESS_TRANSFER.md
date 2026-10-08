# Route averaging with exact outer freshness

The fresh inputs supplied by GRS[32,22] cannot be sampled independently in each packet. This note keeps their rowwise dependence exactly. It gives a finite upper bound without enumerating regional count profiles. Numerical closure remains a separate question.

The construction has length \(N=bL\), with \(L\) outer rows and \(b\) coordinate regions. Each region has an independent uniform permutation of its \(L\) coordinates. Packets have width \(t\), where \(t\mid L\). The current parameters are \(b=32\), outer parity \(p=10\), \(t=16\), state dimension four, and \(q=256\).

The packet law must be specified separately. Both the unscaled repaired emission law and its eta-scaled variant have an interface in [rank_packet_operator.py](theory/rank_packet_operator.py). Changing that law changes the numerical operator, not the route-averaging argument below.

## The row polynomial

Fix an exact outer support. A nonempty row has support size \(m>p\) and shortened dimension \(m-p\). Designate its first \(m-p\) active coordinates in stream order as fresh. Its last \(p\) active coordinates are dependent. The fresh forms are independent of all earlier input forms, as established in [REGIONAL_SPARSE_COVER.md](REGIONAL_SPARSE_COVER.md).

For positive region marks \(A_j,B_j\), define

\[
 H_{b,p}(A,B)
 :=1+\sum_{\substack{S\subseteq\{0,\ldots,b-1\}\\|S|>p}}
       \prod_{j\in\operatorname{first}_{|S|-p}(S)}A_j
       \prod_{j\in\operatorname{last}_{p}(S)}B_j.
\]

The constant term represents the empty row. The exponent of the fresh marks is exactly the shortened dimension, not a probabilistic estimate.

A positive automaton evaluates \(H\) in \(O(bp)\) arithmetic operations. Scan regions backwards. States \(0,\ldots,p\) count the dependent coordinates already encountered. One additional state records that a fresh coordinate has occurred. Inactive coordinates preserve the state. An active coordinate advances a count below \(p\), with weight \(B_j\). From count \(p\), it enters the additional state with weight \(A_j\); later active coordinates remain there with weight \(A_j\). Accept only the empty state and the additional state. This gives twelve states when \(p=10\).

Across the \(L\) labelled outer rows, the support polynomial is \(H^L\). This factorization counts row choices. It does not assert independence of packet freshness flags.

For comparison, fix a routed ternary pattern in each region, with \(f_j\) fresh slots and \(d_j\) dependent slots. Its exact probability of inducing valid row patterns is

\[
 \frac{[\prod_j A_j^{f_j}B_j^{d_j}]\,H_{b,p}(A,B)^L}
      {\prod_j\binom{L}{f_j,d_j,L-f_j-d_j}}.
\]

Each regional permutation assigns labelled rows uniformly to the three slot classes. The numerator counts the assignments compatible with every row. The denominator counts all assignments. This identity is useful for small checks, but its direct coefficient evaluation is unnecessary for the bound below.

## Required packet interface

Fix an outer support and all regional orientations. Only routes and the independent nonzero input diagonals remain random. Maintain the common span of imposed zero equations. State zero means that the state-zero equations have actually been imposed. States one through four record a required-nonzero state map of that form-rank.

For \(1\le\tau\le q-1\) and \(0<y\le1\), let \(K_{f,d}(\tau,y)\) be a nonnegative matrix. It sums over every disjoint placement of \(f\) fresh inputs and \(d\) dependent inputs in one packet, and over every output support. An output support of size \(v\) contributes \(y^v\). Matrix entries upper-bound the conditional contraction of the alive rank potential \(\tau^{-\operatorname{rank}\mathcal W}\).

The required dominance is entrywise, conditional on every possible revealed history and compatible collection of input forms. It must include state-zero checkpoint equations and any exceptional-diagonal probabilities. It must not reuse a diagonal after conditioning on an event involving that diagonal. A statement only for monotone continuation vectors is insufficient for the matrix composition used here.

Define the uniform-placement average

\[
 N_{f,d}:=\frac{K_{f,d}}
                    {\binom{t}{f,d,t-f-d}}.
\]

Choose positive region marks \(v_j,w_j\), and take the maximum entrywise:

\[
 M_j:=\max_{f,d\ge0,\ f+d\le t}
              \frac{N_{f,d}}{v_j^f w_j^d}.
\]

If packet phases require different upper operators, each phase uses its own matrix in the product. The marks \(v_j,w_j\) must remain the same throughout that region. The current packet operator is uniform over the sixteen distinct labels in every phase, so one \(M_j\) suffices.

## Finite bound without regional profiles

Let \(1\le W_-\le W_+\le N\), and let \(0\le e\le N\). We bound the probability that a nonzero word has inner input weight in \([W_-,W_+]\) and output weight at most \(e\). Choose \(x>0\), and put

\[
 W_x:=\begin{cases}W_+,&x\le1,\\W_-,&x>1.\end{cases}
 \qquad
 \mathcal M:=\prod_{j=0}^{b-1}M_j^{L/t},
\]

where matrices occur in stream order. Let \(e_0\) be the row vector supported on state zero, and let \(\mathbf1\) be the all-ones column vector. Under the packet interface above, the desired failure probability is at most

\[
 \boxed{\min\left\{1,
   \frac{H_{b,p}(\tau x v,xw)^L-1}
        {\tau x^{W_x}y^e}\;
      e_0\mathcal M\mathbf1\right\}.}
\]

This is a support-and-rank union bound. It is not a union over field-valued source words. It drops conditional clipping and may therefore be loose.

### Proof

Fix one nonempty outer support, including its canonical freshness labels. Within a region, expose the route packet by packet. Conditional on the revealed route, the unplaced row identities remain uniformly assigned to the remaining positions. After also fixing the next packet's counts \((f,d)\), its ternary placements are uniform. The conditional packet interface therefore gives \(N_{f,d}\), even though the input forms and incoming state depend on the history.

Replace \(N_{f,d}\) by \(v_j^f w_j^dM_j\). For every possible count sequence through the region, the scalar factor becomes \(v_j^{f_j}w_j^{d_j}\). The factors depend only on the fixed regional totals. Summing the conditional hypergeometric probabilities over all count sequences gives one. No Bernoulli comparison or conditioning penalty is needed.

The same argument composes over regions because each regional permutation is independent. Entrywise dominance permits multiplication by every nonnegative future continuation. The initial state is zero. The final state is unrestricted because it is discarded, giving the terminal vector \(\mathbf1\).

For this fixed outer support, write \(D\) for its shortened dimension. A surviving nonzero witness leaves equation rank at most \(D-1\). Multiplying the expected rank potential by \(\tau^{D-1}\) bounds existence. Every fresh coordinate contributes one to \(D\), while every active coordinate contributes one to the input weight. Thus summing support weights gives

\[
 \sum_{\text{nonempty valid outer supports}}
   \tau^D x^{\operatorname{wt}(I)}
   \prod_j v_j^{f_j}w_j^{d_j}
 =H_{b,p}(\tau xv,xw)^L-1.
\]

For an admitted input weight \(W\), \(x^W\ge x^{W_x}\). For an admitted output weight \(E\), \(y^E\ge y^e\). Dividing by these two thresholds and by \(\tau\) proves the expression. A final minimum with one is always valid.

The subtraction of the empty source occurs in \(H^L-1\), before multiplying by the matrix envelope. Subtracting only one from the final expression would leave fictitious empty-source paths introduced by the envelope.

### A scalar continuation bound

The entrywise maximum can lose information by using a different pair \((f,d)\) for each successor rank. A fixed positive continuation vector avoids that particular loss. Choose \(V\in\mathbb R_{>0}^5\), and define

\[
 \rho_j:=\max_{f+d\le t}\max_r
       \frac{(N_{f,d}V)_r}{v_j^f w_j^d V_r}.
\]

Then the same theorem holds after replacing \(e_0\mathcal M\mathbf1\) by

\[
 \frac{V_0}{\min_r V_r}\prod_j\rho_j^{L/t}.
\]

Indeed, \(N_{f,d}V\le v_j^f w_j^d\rho_jV\) entrywise. Iterate this inequality through all packets, and use \(\mathbf1\le V/\min_rV_r\) at the terminal state. No ordering among the entries of \(V\) is assumed. The implementation exposes this alternative through `continuation=`.

## Implementation and remaining numerical gate

[outer_freshness_transfer.py](theory/outer_freshness_transfer.py) implements the positive row automaton, exact rational matrix envelope, and exact bound. Its adapter reads the integer terms of the separately owned packet operator. The `law` argument distinguishes `unscaled` from `eta-scaled`.

`floating_log_bound` uses scaled matrix powers and log-domain row evaluation. Its results are exploratory. Underflow of the nonempty row term returns an unusable bound, rather than claiming zero failure probability.

Seven exact small-case tests cover row enumeration, the coefficient/probability identity, packet normalization, and exhaustive route averaging. They also compare rational and floating matrix and scalar-continuation evaluation for both packet laws. These checks establish the implementation identities; the conditional packet theorem remains a separate obligation.

```sh
python -B research/workstreams/char2_field_spin/theory/outer_freshness_transfer_tests.py
python -B research/workstreams/char2_field_spin/theory/outer_freshness_transfer.py --screen
python -B research/workstreams/char2_field_spin/theory/outer_freshness_transfer.py --active-screen
```

### Executed sparse-class screen

The corrected eta-scaled packet operator was evaluated at \(K=2^{20}\), with \(L=47680\), \(N=1525760\), and \(e=152575\). The input interval was the entire unresolved class \(1\le W\le289894\), namely \(W<0.19N\).

A bounded search for the scalar continuation expression gave normalized base-two logarithm \(+0.7168375344N\) with constant region marks. It therefore yields only the trivial probability bound one. Allowing special marks at the two deterministic endpoints gave \(+0.7203545241N\) in that bounded search. No fresh coordinate occurs in the last ten regions, and the first region has no dependent coordinate. These values are not global optimization results or evidence that the construction fails.

A separate diagnostic isolates some packet loss. Treat every input as dependent, so that the packet polynomial is \(\sum_d x^dK_{0,d}\). Replace the two-state packet spectral radius in the root's return-only screen by this five-state radius, keeping its outer-row bound and parameters unchanged. The resulting normalized logarithms are:

| Input density | Return-only screen | Unscaled rank operator | Eta-scaled rank operator |
|---:|---:|---:|---:|
| .01 | +.090644 | +.091720 | +.088018 |
| .10 | +.074191 | +.082884 | +.151488 |
| .18 | +.002001 | +.015054 | +.076879 |

This diagnostic omits the root screen's finite boundary factors. It is a comparison at fixed witnesses, not a certificate or a comparison of globally optimized bounds. At densities .10 and .18, eta's conservative rank transitions add about \(.075N\) to that tested exponent. The executed whole-class monomial envelope loses substantially more. The current numerical result does not distinguish the best achievable monomial envelope from the bounded search's optimization loss.

The Bellman packet action now avoids duplicated successor-rank mass. At the constant-mark witness above, it changes the whole-class expression only from \(+.71683753N\) to \(+.71683682N\). Thus that witness's large loss is not caused by successor-rank duplication alone. A bounded retuning gave \(+.64733311N\), still unusable.

## Count-aware comparison using active outer rows

The monomial envelope ignores how many input positions a region can contain. An active-row count supplies a deterministic restriction. If exactly \(A\) outer rows are nonempty, every region has at most \(A\) active coordinates. This restriction is useful even when \(A\ge t\), where a local bound \(f+d\le\min(A,t)\) adds nothing.

Fix an interval \(1\le A_-\le A_+\le L\), and put \(\epsilon:=A_+/L\). Its exact source factor is

\[
 S_A:=\sum_{a=A_-}^{A_+}\binom La
       [H_{b,p}(\tau xv,xw)-1]^a.
\]

The subtraction now occurs separately in every active row. This preserves the condition that each selected row has more than \(p\) active coordinates.

### Conditioning on regional counts

Fix one actual row-support pattern and one region. Its routed fresh and dependent counts are \(f,d\), with \(f+d\le A_+\). Compare its uniform fixed-count placements to iid ternary categories with probabilities

\[
 \pi=(f/L,d/L,1-(f+d)/L).
\]

Conditioning this iid law on its exact counts gives the required uniform placement law. The conditioning cost is the inverse of

\[
 p_L(f,d):=\binom L{f,d,L-f-d}
       (f/L)^f(d/L)^d(1-(f+d)/L)^{L-f-d}.
\]

Let \(\nu(L,A_+)\) be the largest inverse cost over \(f+d\le A_+\). It has an exact closed form. Set

\[
 z_*:=\max\{L-A_+,\lfloor L/3\rfloor\},\quad
 f_*:=\lfloor(L-z_*)/2\rfloor,\quad d_*:=L-z_*-f_*.
\]

Then \(\nu=1/p_L(f_*,d_*)\). In particular, when \(A_+\le2L/3\), the minimizing counts are the balanced split of \(A_+\), followed by \(L-A_+\) zeros.

To verify the minimizer, write the mode probability as

\[
 p_L(f,d)=\frac{L!}{L^L}\prod_{n\in\{f,d,L-f-d\}}\frac{n^n}{n!}.
\]

The successive differences of \(\log(n^n/n!)\) are \(n\log(1+1/n)\), which increase with \(n\). Moving one unit from a larger count to a smaller count therefore decreases the product. Balance the counts subject to the required zero-count lower bound. Use \(0^0=1\).

### The restricted iid packet bound

Write \(\mathcal K_{f,d}(V)\) for the packet agent's positive Bellman action on a fixed vector \(V>0\). It is nonlinear in \(V\), so it must not be inserted into a matrix power. Its conditional upper-bound contract is the same as before: it is uniform over all compatible form histories.

For each region, bound

\[
 \rho_j\ge
 \max_{\substack{u,v\ge0\\u+v\le\epsilon}}\max_r
 \frac1{V_r}\sum_{f+d\le t}
  \left(\frac{u}{v_j}\right)^f
  \left(\frac{v}{w_j}\right)^d
  (1-u-v)^{t-f-d}\mathcal K_{f,d}(V)_r.
\]

The iid dependent probability is zero in the first region. The iid fresh probability is zero in the last \(p\) regions. These endpoint restrictions may be imposed in the displayed maximum.

There is a positive finite bound on this maximum. For a fixed incoming rank, its full-simplex Bernstein coefficients are

\[
 c_{f,d,r}:=
 \frac{\mathcal K_{f,d}(V)_r}
 {\binom t{f,d,t-f-d}v_j^f w_j^dV_r}.
\]

On the restricted triangle, the Bernstein coefficient at \((i,j)\) becomes

\[
 \widehat c_{i,j,r}
 =\mathbb E[c_{F,D,r}],\qquad
 F\sim\operatorname{Bin}(i,\epsilon),\quad
 D\sim\operatorname{Bin}(j,\epsilon),
\]

with independent binomials. Take the largest coefficient over \(i+j\le t\) and incoming rank. At a restricted edge, retain only the coefficients on that edge. The polynomial is a convex combination of these coefficients, so this is a proved upper bound, not a sampled maximization.

### Active-row bound and its scope

Intersect the input interval with

\[
 (p+1)A_-\le W\le bA_+.
\]

Define \(W_x\) using the resulting interval as before. If the intersection is empty, the class is empty. Otherwise the probability of the class is at most

\[
 \boxed{\min\left\{1,
 \frac{S_A\,\nu(L,A_+)^b}{\tau x^{W_x}y^e}
 \frac{V_0}{\min_rV_r}\prod_j\rho_j^{L/t}\right\}.}
\]

For the comparison, first fix an actual row-support pattern. Every routed history has the conditional packet bound. Extend this positive abstract Bellman game to all iid category histories. Histories incompatible with any actual outer pattern add only nonnegative terms. They are not asserted to be physical source words or independent fresh forms. Conditioning on the actual regional counts recovers the original route law, with cost at most \(\nu\) per region. The positive supersolution \(V\) bounds the extended iid game by \(\rho_j^{L/t}\). Summing the actual row-support patterns gives \(S_A\). This proves the expression without an iid-freshness assumption.

`active_row_log_bound` evaluates the expression with the eta-scaled Bellman action. The mode constant and Bernstein restriction have exact exhaustive tests. A further test compares the complete fixed-count/iid route inequality with all routes of a four-row example. There are ten tests in total.

An executed five-bin cover of all possible active-row counts at \(K=2^{20}\) gave:

| Active-row interval | Tested normalized log2 upper expression |
|---:|---:|
| 1--16 | +.072909 |
| 17--256 | +.090749 |
| 257--4096 | +.189134 |
| 4097--16384 | +.567055 |
| 16385--26354 | +.643528 |

These bounded searches used the entire input interval \(1\le W<.19N\). Every tested expression is positive. The cover is complete, but its evaluated upper bound is trivial. Wide input and active-row intervals also lose their correlation; fixed-weight slices provide a separate diagnostic. No numerical value here is directed or globally optimized.

At the fixed weight \(W=152576=.1N\), three further active-count slices gave \(+.131885N\), \(+.355067N\), and \(+.390493N\) for \(A=4768,9536,13870\). These remain weaker than the existing return-only middle-class screen. The `--active-screen` replay includes these witnesses and the additive components of every exponent.

The component breakdown identifies a concrete loss. The iid maximum includes the empty profile \(u=v=0\). Its state-zero self-loop forces \(\rho_j\ge1\) in every region. At the first four whole-class witnesses, all packet contributions are zero, meaning that this lower limit is attained. The separated bound can select the empty profile in each region even though an actual active row occupies at least eleven regions. The source polynomial still counts valid row patterns, but it no longer forces their state excursions inside the packet factor.

The next useful refinement should retain regional activity or chronology inside the packet composition, instead of maximizing each region separately. The exact row automaton, active-row restriction, and modal conditioning constants remain available for that refinement. They do not yet close the sparse-input class.
