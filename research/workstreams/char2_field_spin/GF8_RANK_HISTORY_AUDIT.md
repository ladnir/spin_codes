# Audit of the restricted-state rank-history transfer

The rank-cap transfer extends from eight state coordinates to any restriction
of the checked symmetric map with \(1\le s\le8\). Its joint return entry may
also retain the incoming cap. No additional independence assumption is
needed. This audit establishes the local composition contract; it does not
certify a distance bound for a new outer/state profile.

Fix one exact outer support, all routing choices, and the setup randomness
exposed before the current packet. Let \(V\) be the shortened source space,
and let \(W\subseteq V^*\) be the span of the zero equations imposed so far.
State and input coordinates are linear forms in \(V^*/W\). They may be
dependent. The current input and state diagonals remain mutually independent
and uniform on \(\operatorname{GF}(256)^*\).

The packet uses the restricted matrices \(A\in\operatorname{GF}(256)^{16\times s}\)
and \(C=A^T\). They are superregular and satisfy \(CA=0\). Its update is

\[
R:=DS,\qquad Y:=X+\eta AR,\qquad S':=R+CX.
\]

Here \(X\) includes the current input diagonals. Fix \(0<y\le1\) and
\(1\le\tau\le255\), the output mark and equation-rank tilt.

## State meaning and exposure order

Transfer state zero means the witness has zero state and its state-zero
equations have been imposed in the common span \(W\). A positive transfer
state \(c\) means the witness has nonzero state and the state-form rank is
at most \(c\). Positive states are caps, not mutually exclusive rank events.

On an incoming-nonzero branch, first select \(z\) zero outputs outside the
active input support. If the incoming state image has dimension \(r\le c\),
define

\[
H:=\dim\bigl(\operatorname{im}(S)\cap\ker(A_ZD)\bigr).
\]

The selected equations add rank \(r-H\). A branch with \(H=0\) cannot retain
an incoming-nonzero witness, so it is discarded. Current live inputs do
not turn that dead branch into a birth. Births belong to the separate
incoming-zero entry.

For \(H>0\), the \(i\) selected zeros at active input positions each expose
a fresh input diagonal. The alive diagonal-pivot bound supplies
\(\tau^{-i}\). Each selected equation expresses its input form as a linear
combination of the residual old-state forms. Consequently

\[
\operatorname{rank}(S'\bmod W')\le\min(s,H+w-i),
\]

where \(w\) is the active input count and \(W'\) includes all selected
output-zero equations. This argument permits dependent input forms and
singular allowed-output updates.

## Simultaneous tails and composition

The state diagonal satisfies

\[
0\le H\le\min(r,s-z),\qquad
\mathbb E_D[256^H-1]\le\frac{256^r-1}{255^z}.
\]

For \(1\le j\le s\), let \(T_j(c,z)\) be the largest value of
\(\mathbb E[\tau^{H-r}\mathbf1_{H\ge j}]\), over \(1\le r\le c\) and all
laws satisfying these constraints. Set \(T_{s+1}(c,z)=0\). These tail bounds
decrease with \(j\) and increase with \(c\). Thus

\[
\mu_h(c,z):=T_h(c,z)-T_{h+1}(c,z)\ge0
\]

defines one positive measure whose upper tails dominate every admissible
tilted law. Summation by parts gives domination for every nonnegative,
nondecreasing function of \(H\). The measure need not be an attainable
state-rank distribution.

Let the continuation vector satisfy \(v_j\ge0\) for \(0\le j\le s\) and
\(v_1\le\cdots\le v_s\). No order between \(v_0\) and \(v_1\) is required.
The successor map \(h\mapsto\min(s,h+w-i)\) is nondecreasing. Therefore the
measure can be applied to the continuation after the input pivots. Taking
unrelated maxima for individual successor masses would not justify this
step.

The return entries below increase with \(c\). Every packet matrix therefore
preserves this continuation cone. Positive matrix products and nonnegative
sums preserve the same cone. The all-ones terminal continuation belongs to
it, including when the final state is discarded. Regional occupancy sums
and density sums consequently preserve the local history contract.

## The cap-dependent joint return

Fix the exact incoming zero flat before exposing current diagonals. Let its
rank be \(k\), and let the actual live state support have size \(h\). A
nonzero incoming witness requires \(k<r\le c\). Return feasibility requires
\(w>s-h\), since otherwise the absent-state rows of \(C\) force every active
input to vanish. Inactive selected zeros require \(z<h\). These conditions
give

\[
k\le\min(w-1,c-1,s-1-z).
\]

There are at most \(\binom{s}{k}\) rank-\(k\) flats, by their distinct
canonical independent bases. Imposing the flat costs \(\tau^{-k}\) before
any current diagonal is used.

Put \(v:=16-i-z\), the number of unselected output coordinates. The numerical
matrix imposing both return and selected zeros has rank at least

\[
\rho_s(w,i,z):=\max\{s,\ i+z,\ w-16+i+2z+1\}.
\]

Return alone has rank \(s\) on every feasible branch. The selected output
rows have rank \(i+z\), using the active identity columns and inactive
Cauchy rows. Finally, the full packet map is invertible because \(CA=0\).
The joint kernel therefore embeds in the \(v\) unselected outputs, giving
rank at least \(h+w-v\ge w-16+i+2z+1\).

Let \(m_s(w)\) be a proved lower bound on nonempty return-output weight.
The conservative choice for \(s<8\) is \(\max(1,17-s-w)\). For \(s=8\), the
existing two-half argument gives the stronger recorded sequence. Define

\[
P_c(w):=y^{m_s(w)}\tau^{-s}
 \sum_{k=0}^{\min(w-1,c-1)}\binom{s}{k}\tau^{-k},
\]

\[
J_c(w):=
\sum_{\substack{0\le i\le w,\ 0\le z\le\min(s-1,16-w)\\
                 16-i-z\ge m_s(w)}}
\binom wi\binom{16-w}z
y^{16-i-z}(1-y)^{i+z}\tau^{-\rho_s(w,i,z)}
\sum_{k=0}^{\min(w-1,c-1,s-1-z)}\binom sk\tau^{-k}.
\]

Both are complete return covers. Their common successor is actual zero,
so the return coefficient may be \(\min(P_c(w),J_c(w))\). Both functions
increase with \(c\), and their minimum does too. The joint cover charges
the current diagonals once; it receives no additional marginal root factor.

Stronger deterministic return minima can replace \(m_s(w)\) in both covers.
On return, \(R=CX\), so \(Y=(I+\eta AC)X\). For input weight two, normalize
the first live coordinate to one and enumerate the second over GF256*.
Together with all input placements and packet phases, this exhausts the
possible output supports. A common nonzero scalar does not change support.
No probability or equation-rank factor changes when using such a minimum.

## Checks and remaining gate

The following existing suites passed together: 19 tests covering exact
GF4 tail bounds, all continuation-cone generators, regional composition,
and the GF16 joint-rank analogue.

```sh
python -B -m unittest state_rank_cap_transfer_tests state_rank_cap_composition_audit_tests state8_joint_return_tests
```

The nine tests in gf8_rank_history_tests also passed independently. They
check the generalized local transfer, cap-dependent joint returns, exact
cone preservation, and the batched aggregate against an exact local sum.
The aggregate's binomial identity agrees with the direct selected-zero
enumerator after substituting \(d=w-i\).

An additional independent check used seed 71063 to sample 500 source-form
configurations over GF16. Source dimensions ranged from one to four and
input weights from zero to four. It checked all sixteen selected-zero masks
for each configuration using exact field ranks. All 8,000 successor-rank
inequalities held; 3,805 were strict. This finite check supplements the
span argument above and is not an exhaustive claim over source forms.

## Removing impossible MDS row supports

The source comparison can omit nonzero row weights below the MDS distance.
Homogenizing the resulting polynomial requires a condition that the older
two-product comparison did not need.

Let the outer code have length \(b\), dimension \(b-p\), and \(L\) rows.
Here \(1\le p<b\), and \(L\) is even. Fix regional support counts
\(n_1,\ldots,n_b\), total weight \(W>0\), and a positive source tilt
\(\lambda\). Put \(\alpha_j:=n_j/L\), \(\alpha:=W/(bL)\), and
\(\beta:=\tau\lambda\). Under independent Bernoulli coordinates with
probabilities \(\alpha_j\), the exact feasible row factor is
\(\mathbb E[g(H)]\), where \(H\) is row occupancy and

\[
g(h):=
\begin{cases}
1,&h=0,\\
0,&1\le h\le p,\\
\tau^{-p}\beta^h,&p+1\le h\le b.
\end{cases}
\]

The zero row occurs once. A nonzero feasible row has shortened dimension
\(h-p\); its factor is therefore \(\tau^{h-p}\lambda^h\).

If \(\beta\ge2\), every second difference of \(g\) is nonnegative. When
\(p+2\le b\), the only potentially problematic boundary satisfies
\[
g(p)-2g(p+1)+g(p+2)
=\tau^{-p}\beta^{p+1}(\beta-2)\ge0.
\]
When \(p=b-1\), that boundary difference is absent. For \(p=1\), the initial second difference is
\(1+\tau^{-1}\beta^2>0\). Other boundary differences are nonnegative.

Fix all Bernoulli coordinates except two, and let \(K\) be their remaining
occupancy. Replacing probabilities \(a-d,a+d\) by \(a,a\) increases the row
expectation by \(d^2\mathbb E[g(K)-2g(K+1)+g(K+2)]\). Repeated pair averaging
therefore gives the homogeneous upper factor

\[
F(\alpha):=(1-\alpha)^b+
\tau^{-p}\sum_{h=p+1}^b\binom bh
(\beta\alpha)^h(1-\alpha)^{b-h}.
\]

Condition each regional Bernoulli column on its own count \(n_j\). That
conditional column is the required uniform subset. Each conditioning
penalty is at most \(M_L:=2^L/\binom L{L/2}\), the reciprocal central modal
probability. For shortened dimension \(D\), the resulting source bound is

\[
\mathbb E_{\mathrm{route}}
 [\mathbf1_{\mathrm{feasible}}\tau^D]
\le M_L^b\lambda^{-W}F(\alpha)^L.
\]

The all-empty support remains in the unconditional upper expectation.
Keeping that positive term is harmless for the fixed class \(W>0\).
The conditioning factor and source bound are positive scalars, so neither
changes the rank-cap continuation cone.

For \(\beta<2\), the homogeneous truncated formula needs another argument.
For example, take \(b=3,p=1,\tau=\lambda=1\) and probabilities \((1,1,0)\).
The exact feasible row expectation is one, while the homogeneous expression
at \(\alpha=2/3\) is \(7/9\). A safe implementation uses the original source
bound outside \(\beta\ge2\).

An independent finite check exhausted every binary support matrix for
\((L,b,p)=(2,3,1),(4,3,1),(2,4,2),(4,4,2)\). With \(\tau=3\) and
\(\beta\in\{2,5,15\}\), all 2,574 conditional profile comparisons satisfied
the displayed bound.

The next gate is numerical evidence that the refined transfer and source
factor improve the unresolved classes, followed by directed finite coverage
if they do. Existing certificates, fixed-map implementations, and runtime
claims are outside this audit.
