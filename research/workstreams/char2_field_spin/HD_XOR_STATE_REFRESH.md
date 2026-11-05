# Diagonal-XOR-diagonal state refresh: exact support law

Replacing the fixed inner interpolation map by an XOR map could remove its
20 fixed GFNI products per packet. The two random diagonals would remain.
This note checks the resulting support law and identifies a possible proof
interface. It does not establish a distance bound. The separate
[native experiment](HD_XOR_STATE_NATIVE.md) now implements this refresh:
isolated inner time falls about 20%, but complete encoding improves only
about 1--5% in the matched runs.

## Candidate and remaining randomness

Work over $\mathbb F_q$ with $q=256$. Keep the current 64-symbol packet,
16-symbol state, and four RS[16,4] blocks in $A$, with $C=A^T$ and $CA=0$.
Replace only the refresh in

$$
Y=X+AQ,\qquad W=Q+CX,\qquad Q'=D_2HD_1W,
\qquad H=I+\mathbf1\mathbf1^T.
$$

Each packet has fresh independent diagonal entries from $\mathbb F_q^*$.
They are independent of the outer code, route, and other refresh keys.
Setup is sampled once and reused by the encoder. The proof averages over
that setup; the diagonals are not resampled separately for each message.

Since 16 is even and the characteristic is two, $H^2=I$. Thus every refresh
is invertible, and $Q'=0$ exactly when $W=0$. The existing causal inverse
still applies. Invertibility alone does not give the existing distance bound.

The static cost comparison retains eight variable GFNI products and removes
the fixed map's 20 GFNI products per packet. A straightforward XOR realization
on four ZMM registers uses nine vector XORs and two 128-bit-lane shuffles:
reduce corresponding lanes across the four registers, reduce the four lanes,
and add the replicated sum to each register.
At $K=2^{20}$ there are 32,768 packets, so the proposed
deletion is 655,360 fixed GFNI products. Random-diagonal key storage is unchanged.

## Exact conditional support kernel

Fix $W$ and let $S=\operatorname{supp}(W)$ and $s=|S|$. If $S$ is empty,
the output is zero. Otherwise $x=D_1W$ has independent uniform nonzero
entries on $S$. Write $\sigma=\sum_{i\in S}x_i$.

If $\sigma=0$, then $Hx=x$. The number of such inputs is

$$
Z_s=\frac{(q-1)^s+(q-1)(-1)^s}{q}.
$$

If $\sigma\ne0$, put $J=\{i\in S:x_i=\sigma\}$. Then
$\operatorname{supp}(Hx)=[16]\setminus J$. Fixing $J$, divide all entries
by $\sigma$. The remaining $n=s-|J|$ entries avoid both 0 and 1, and their
sum is $1+(|J|\bmod2)$ in characteristic two. For $n>0$, either required
sum, 0 or 1, has the same count

$$
L_n=\frac{(q-2)^n+(q/2-1)(-2)^n}{q}.
$$

For $n=0$, the count is one if $|J|$ is odd and zero otherwise; denote this
case-dependent count by $L_{0,|J|}$. For $n>0$, set $L_{n,|J|}=L_n$.
The exact support transition probability is

$$
P(S,T)=\frac{Z_s\,\mathbf1_{T=S}
 +(q-1)\sum_{J\subseteq S}
       L_{s-|J|,|J|}\,\mathbf1_{T=[16]\setminus J}}
 {(q-1)^s}.
$$

The zero-sum and nonzero-sum events are disjoint, but their output supports
can coincide. Their contributions must be added; in particular, this happens
when $S=T=[16]$. The formula is the existing
`hd_xor_column_probe.column_support_counts` kernel with group width 16.

Conditional on any output support $T$, the fresh $D_2$ makes the output
amplitudes independent uniform elements of $\mathbb F_q^*$. More precisely,
every vector $y$ with support $T$ has conditional probability
$(q-1)^{-|T|}$. This remains true for a measure weighted by the current and
earlier emissions: those weights are determined before the fresh $D_2$.

Two sparse cases illustrate the distinction from the rank-one refresh:

- A singleton input support becomes its 15-coordinate complement with probability one.
- A two-coordinate support stays unchanged with probability $1/(q-1)$ and becomes full support with probability $(q-2)/(q-1)$.

On the two-coordinate branch, the output amplitudes are independent, not
confined to the previous scalar line. For one fixed two-coordinate target
vector, the transition probability is $(q-1)^{-3}$. Hence the existing
dimension-16 diagonal-Cauchy atom cap cannot be reused.

## A possible 625-state interface

The exact support-flat cone has one component for each coordinate support
$S\subseteq[16]$: arbitrary nonnegative mass, uniform nonzero amplitudes on
that fixed support. Fresh $D_2$ restores this cone after every packet.
There is no additional arbitrary-amplitude component after refresh.

Partition the state coordinates into the four fixed RS blocks. Let
$\mathbf k=(k_0,k_1,k_2,k_3)$ count live state coordinates in each block.
There are $5^4=625$ such count vectors, including zero. The refresh alone
lumps exactly to these vectors. Its zero-sum contribution stays at
$\mathbf k$; for cancelled-coordinate counts $0\le j_i\le k_i$, the
nonzero-sum contribution goes to $(4-j_0,4-j_1,4-j_2,4-j_3)$ with probability

$$
\frac{(q-1)L_{s-|\mathbf j|,|\mathbf j|}}
     {(q-1)^s}\prod_{i=0}^3\binom{k_i}{j_i}.
$$

This exact lumping does **not** make the whole packet a 625-state Markov
chain. The columns of the RS emission and feedback maps distinguish state
coordinate positions. Incoming supports with the same block counts can have
different marked emission and return kernels.

A sound 625-state upper operator can instead retain arbitrary mixtures of
support-flat measures with each fixed count vector. Average fresh input and
route randomness for each exact incoming support, sum masses going to each
target count vector, then maximize over incoming supports in the source
count class. The resulting row dominates every such mixture. No assumption
of uniform support locations inside a block is justified by $D_2$ alone.

Computing useful packet entries remains the main proof obligation. Locally,
one needs bounds on

$$
\mathbb E\!\left[z^{\operatorname{wt}(X+A_SQ)}
                 \mathbf1_{\operatorname{supp}(Q+C_UX)=T}\right],
$$

with $Q$ uniform nonzero on an exact subset of four state coordinates.
Such a subset of RS coefficient columns need not itself be MDS or totally
nonsingular. One may dominate its nonzero-coefficient law by
$(q/(q-1))^{|S|}$ times the uniform-coefficient law, but must still use the
ranks of the actual selected columns.
The input law for $X$ must match the outer comparison: uniform-field values
for the Toeplitz comparison, or forced-nonzero values in the XOR support
analysis. These laws are not interchangeable. Exact support/zero constraints
can be counted by ranks and inclusion-exclusion in the four-state local
problem. Rank data are needed to calculate the kernel, not necessarily as
additional persistent proof states. A 16-state active-block mask could be
obtained by another maximum, but its loss has not been assessed.

## Obstructions, prior work, and next check

The sum-zero hyperplane is fixed by $H$, but fresh independent $D_1$ does
not preserve that condition. There is no nontrivial coordinate subspace
preserved by all refreshes: a singleton in that subspace produces nonzero
coordinates outside it. The two-coordinate branch is a genuine sparse
survival event of probability $1/255$, not a deterministic invariant.
These facts do not rule out other low-weight codewords.

The old uniform rank-one law preserves a scalar-line component. Its failed
three-state and block-mask bounds do not apply to this different law.
There was an older related diagnostic: `hd_systematic_mds_screen.py` uses
`feedback='lite'` for $D_2(I+c\mathbf1\mathbf1^T)D_1$, with $c\notin\{0,1\}$.
[HD_ENCODING_OPTIMIZATION.md](HD_ENCODING_OPTIMIZATION.md) records its failed
good-part/residual screen. That study used GF65536, global RS emission, and
separate carried-state/feedback maps; it did not analyze this common h4
refresh with the exact support-flat cone.

The five tests in `theory/hd_xor_state_refresh_tests.py` check every GF4
width-four input, every output diagonal, exact block-count lumping for the
refresh, overlapping full-support contributions, and width-16 sparse laws.
They are finite local checks, not a distance certificate.

The [exact local kernel](HD_XOR_STATE_LOCAL.md) now handles up to two
potential input positions within one four-state-coordinate block. Seven
additional tests check the actual basis, exact counts, and independent
enumeration. It uses uniform-field input values, not the XOR outer's
forced-nonzero comparison law. No regional or whole-code certificate follows.
The modest complete-encoder gain warrants prioritizing outer proof closure
and data-flow changes before a full 625-state investigation.
