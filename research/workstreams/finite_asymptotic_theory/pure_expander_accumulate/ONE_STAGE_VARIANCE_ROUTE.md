# One-stage variance route

## Status

This document gives an exact reformulation of the one-stage shell variance.
It does not prove the required length-512 variance bound.

The reformulation replaces the sum over message-pair types by a covariance
operator indexed by subsets of the 512 sparse-map rows. The operator is
invariant under row permutations and has an explicit finite block
decomposition. The remaining work is a scalable outward evaluation that
retains the signed block cancellation.

## Objects and probability space

Fix \(k=256\), \(n=512\), and \(r=33\). Setup samples independent rows

\[
 R_1,\ldots,R_n
 \mathrel{\leftarrow}
 \{R\in\mathbb F_2^k:\operatorname{wt}(R)=r\}.
\]

The sparse map is

\[
 (Ex)_j=\langle R_j,x\rangle.
\]

Let \(A\) be the zero-initialized accumulator. For \(0\leq w\leq n\), define

\[
 N_w(E)=
 \bigl|\{x\in\mathbb F_2^k:
              \operatorname{wt}(AEx)=w\}\bigr|.
 \tag{1}
\]

If \(E\) is injective and \(w>0\), then \(N_w(E)\) is the ordinary weight-
\(w\) multiplicity of the realized code. The proposed setup samples a
bounded number of matrices and accepts the first injective matrix. An
unconditioned bad-event bound transfers to the accepted distribution by
division by a lower bound on the one-draw rank-success probability.

## Exact dual formula

Let \(K_w^{(n)}(h)\) be the binary Krawtchouk polynomial. For
\(t\in\mathbb F_2^n\), put

\[
 c_w(t)=K_w^{(n)}
        \bigl(\operatorname{wt}(A^{-\mathsf T}t)\bigr)
 \tag{2}
\]

and

\[
 I_t(E)=
 \mathbf 1\left\{
       \bigoplus_{j:t_j=1}R_j=0
       \right\}.
 \tag{3}
\]

Fourier inversion of the Hamming-shell indicator and character
orthogonality on \(\mathbb F_2^k\) give

\[
 N_w(E)=2^{k-n}\sum_{t\in\mathbb F_2^n}c_w(t)I_t(E).
 \tag{4}
\]

Let \(C\) be the covariance matrix of \((I_t)_t\). Then

\[
 \operatorname{Var}(N_w)
 =2^{2(k-n)}c_w^{\mathsf T}C c_w.
 \tag{5}
\]

Equations (4)--(5) are exact. Equation (5) retains the cancellation that is
lost by separately upper-bounding the second factorial moment.

## Radial rank-257 factorization

Define the row-character bias

\[
 \beta(x)=
 \mathbb E_R[(-1)^{\langle R,x\rangle}]
 =\frac{K_r^{(k)}(\operatorname{wt}(x))}{\binom kr}.
 \tag{6}
\]

For \(0\leq h\leq k\) and \(m\geq0\), let \(v_h(m)\) be the probability that
the XOR of \(m\) independent weight-\(r\) rows equals one fixed vector of
weight \(h\). A \((k+1)\)-state radial random walk computes every \(v_h(m)\).

For uniform independent \(X,Y\in\mathbb F_2^k\), Fourier convolution gives

\[
 \mathbb E[
   \beta(X)^a\beta(Y)^b\beta(X+Y)^d]
 =\sum_{h=0}^{k}\binom kh v_h(a)v_h(b)v_h(d).
 \tag{7}
\]

Thus every scalar covariance entry is determined by a
\(257\times513\) table. The implementation need not enumerate
\(\mathbb F_2^k\times\mathbb F_2^k\), nor revisit every message-pair type for
every row-subset triple.

The script verify_radial_triple_factorization_small.py checks all 165 triples
with \(a+b+d\leq8\) at \((k,r)=(4,3)\) in exact rational arithmetic.

## Exact covariance kernel

For a row-index subset \(s\), the return probability is

\[
 p_{|s|}:=\mathbb E[I_s]
 =2^{-k}\sum_{x\in\mathbb F_2^k}\beta(x)^{|s|}
 =v_0(|s|).
 \tag{8}
\]

For two subsets \(s,t\), put

\[
 a=|s\setminus t|,\qquad
 b=|t\setminus s|,\qquad
 d=|s\cap t|.
\]

Two applications of character orthogonality give

\[
 \mathbb E[I_sI_t]
 =2^{-2k}\sum_{x,y\in\mathbb F_2^k}
      \beta(x)^a\beta(y)^b\beta(x+y)^d.
 \tag{9}
\]

Hence

\[
 C_{s,t}=\mathbb E[I_sI_t]-p_{a+d}p_{b+d}.
 \tag{10}
\]

Equations (7), (9), and (10) give a rank-\((k+1)\) evaluation of the exact
kernel. The value depends only on \((|s|,|t|,|s\cap t|)\). Therefore \(C\)
commutes with the symmetric-group action on the \(n\) row indices.

## The parity mode

The degree \(r=33\) is odd. Thus

\[
 E\mathbf 1_k=\mathbf 1_n
 \quad\text{and}\quad
 AE\mathbf 1_k=(1,0,1,0,\ldots,1,0).
 \tag{11}
\]

Every realized code contains this deterministic word of weight \(n/2\).
Equivalently,

\[
 \beta(x+\mathbf 1_k)=-\beta(x).
 \tag{12}
\]

In particular, \(I_t=0\) whenever \(\operatorname{wt}(t)\) is odd. A proof
must retain this parity mode. Treating the row distribution as uniform on
all of \(\mathbb F_2^k\) is invalid.

The moderate-size primal calculation shows another consequence. Covariance
contributions indexed by message-difference weights \(h\) and \(k-h\) can be
large with opposite signs. At \((k,n,r)=(24,48,5)\), individual normalized
contributions near \(28\) cancel after complement pairing, while the complete
off-diagonal contribution is about \(1.806\). A termwise absolute-value
bound destroys this cancellation.

## Uniform-odd comparator

For orientation, replace the constant-weight row by a uniform vector in the
odd-weight affine hyperplane. This comparator preserves (11). Its only
nonzero row-character biases are the trivial character and the all-ones
character.

For every nonempty even subset \(t\), its row sum is uniform in the even-
weight subspace. Hence

\[
 \Pr[I_t=1]=2^{-(k-1)}.
\]

For distinct nonempty even subsets \(s,t\), the two row sums are independent.
The comparator covariance matrix is diagonal on the nonempty even subsets,
with diagonal entry

\[
 2^{-(k-1)}(1-2^{-(k-1)}).
 \tag{13}
\]

This is the reference operator. The certificate must bound the effect of
replacing a uniform odd row by a uniform weight-33 row.

## Why one global operator norm cannot close

For a two-element subset \(t\), the event \(I_t=1\) says that two sampled
rows are equal. Thus

\[
 \Pr[I_t=1]=\binom{256}{33}^{-1}
\]

and

\[
 \log_2\binom{256}{33}=138.1822087505\ldots.
\]

Consequently, the largest eigenvalue of \(C\) is at least approximately
\(2^{-138.18}\). A global use of
\(c_w^{\mathsf T}Cc_w\leq\lVert C\rVert\lVert c_w\rVert_2^2\) would need an
eigenvalue on the order of \(2^{-247}\) for a variance factor \(512\). This
route cannot close.

## Exact accumulator energy

Parseval gives

\[
 \sum_t c_w(t)^2=2^n\binom nw.
 \tag{14}
\]

The joint counts

\[
 \bigl|\{t:\operatorname{wt}(t)=u,
              \operatorname{wt}(A^{-\mathsf T}t)=v\}\bigr|
\]

have an elementary run decomposition. The script
analyze_dual_walk_and_accumulator_energy.py computes these counts, every
Krawtchouk coefficient, and (14) with exact integers.

For the central shell \(w=256\), the fraction of the nonconstant energy in
row-subset weights \(1\) through \(64\) is

\[
 2^{-7.4676\ldots}.
\]

The same fraction is at most \(2^{-11.5602}\) at \(w=192,320\), at most
\(2^{-29.7758}\) at \(w=128,384\), and at most \(2^{-104.8255}\) at
\(w=42,470\). These are binary64 logarithms of exact integer ratios.

Define the nontrivial row-walk spectral mass

\[
 S_m=\sum_{x\notin\{0,\mathbf 1_k\}}|\beta(x)|^m.
 \tag{15}
\]

The exact character values give:

| target for \(S_m\) | first \(m\) meeting the target |
|---:|---:|
| \(1\) | 22 |
| \(2^{-10}\) | 45 |
| \(2^{-20}\) | 68 |
| \(2^{-40}\) | 114 |
| \(2^{-80}\) | 207 |

These data show that low-step covariance modes are large but occupy little
accumulator Fourier energy, while high-step modes approach the uniform-odd
comparator. They do not by themselves bound (5).

## Exact row-permutation blocks

The standard decomposition of functions on all row subsets is

\[
 \mathbb R^{2^{[n]}}
 =\bigoplus_{j=0}^{\lfloor n/2\rfloor}
   S^{(n-j,j)}\otimes\mathbb R^{n-2j+1}.
 \tag{16}
\]

Because \(C\) is invariant, it acts in sector \(j\) as the identity on
\(S^{(n-j,j)}\) tensored with a matrix \(C^{(j)}\).

Write

\[
 F(p,q,\ell)=C_{s,t}
 \quad\text{when}\quad
 |s|=p,\ |t|=q,\ |s\cap t|=\ell.
\]

Put \(m=n-2j\). In the normalized lift basis of sector \(j\), the exact
block entry is

\[
\begin{split}
 C^{(j)}_{p,q}
 ={}&\sqrt{\frac{\binom m{p-j}}{\binom m{q-j}}}
 \sum_{u=0}^{j}(-1)^{j-u}\binom ju\\
 &\quad\cdot\sum_v
 \binom{p-j}{v}
 \binom{m-(p-j)}{q-j-v}
 F(p,q,u+v).
\end{split}
\tag{17}
\]

The inner sum runs over the integers for which both binomial coefficients
are defined.

For a proof of (17), choose \(j\) disjoint coordinate pairs and put
alternating signs on the \(2^j\) sets that choose one coordinate from every
pair. This is a \(j\)-harmonic vector. Its lift to level \(p\) is nonzero
exactly on sets that choose one coordinate from each distinguished pair, and
its square norm is \(2^j\binom m{p-j}\). Counting the number \(u\) of pairs
on which the level-\(p\) and level-\(q\) choices agree gives (17).

The script verify_subset_covariance_blocks_small.py checks the derivation in
two ways. At \((k,n,r)=(3,4,1)\) and \((4,5,3)\), the rational dual
quadratic form equals exhaustive sparse-map shell variance exactly. The
eigenvalue multiset of the blocks in (17), with the representation
multiplicities, matches the full covariance eigenvalue multiset to binary64
residuals below \(3\times10^{-16}\).

## Rejected weight-level relaxation

Define

\[
 e_{w,p}^2=\sum_{t:|t|=p}c_w(t)^2
 \tag{18}
\]

and

\[
 b_{p,q}=\max_{0\leq j\leq
      \min(p,q,n-p,n-q)}|C^{(j)}_{p,q}|.
 \tag{19}
\]

Cauchy--Schwarz at the level boundaries gives the valid bound

\[
 c_w^{\mathsf T}Cc_w
 \leq\sum_{p,q=0}^{n}e_{w,p}b_{p,q}e_{w,q}.
 \tag{20}
\]

The smallest tests initially made (20) look useful. At
\((k,n,r,w)=(4,8,3,6)\), the exact variance factor is \(1.2657\), (20)
gives \(8.2071\), and a global norm gives \(17.6220\). At
\((5,10,3,10)\), the corresponding values are \(1.0392\), \(1.5177\), and
\(3.0481\).

The next sizes reject (20) as the main certificate. At
\((24,48,5,31)\), the exact variance factor is \(2.7906\), while (20) gives
\(4611.74\). At \((32,64,7,40)\), the exact factor is \(4.6473\), while
(20) gives \(15854.46\). These are nondirected binary64 diagnostics, but the
gap is too large to attribute to rounding. Maximizing the sector
independently at every pair of levels destroys the required cancellation.

The script probe_radial_block_bound.py performs this calculation through the
radial factorization. The relaxation (20) is retained as a documented
failure, not as evidence for closure.

## Exact sector projection energies

The sector energy in (41) also has a finite exact interface. Fix a row-subset
level \(p\). Its possible intersection sizes are

\[
 \max(0,2p-n)\leq\ell\leq p,
\]

and its sectors are \(0\leq j\leq\min(p,n-p)\). The adjacency matrix that
selects pairs of level-\(p\) subsets with intersection \(\ell\) has
eigenvalue

\[
 \Lambda_{p,j}(\ell)
 =\sum_{u=0}^{j}(-1)^{j-u}\binom ju
   \binom{p-j}{\ell-u}
   \binom{n-p-j}{p-j-\ell+u}
 \tag{21}
\]

in sector \(j\), with undefined binomial coefficients interpreted as zero.
The square matrix \((\Lambda_{p,j}(\ell))_{j,\ell}\) is invertible.

Let \(\alpha_{p,j,\ell}\) be the solution of

\[
 \sum_\ell
 \alpha_{p,j,\ell}\Lambda_{p,j'}(\ell)
 =\mathbf 1\{j=j'\}.
 \tag{22}
\]

Define the accumulator correlation

\[
 Q_w(p,\ell)=
 \sum_{\substack{|s|=|t|=p\\|s\cap t|=\ell}}
 c_w(s)c_w(t).
 \tag{23}
\]

Then the exact level-\(p\), sector-\(j\) energy is

\[
 \lVert P_{j,p}c_w\rVert_2^2
 =\sum_\ell\alpha_{p,j,\ell}Q_w(p,\ell).
 \tag{24}
\]

Summing (24) over \(p\) gives \(\lVert P_jc_w\rVert_2^2\).

The script verify_sector_projection_energy_small.py computes (21)--(24) in
exact rational arithmetic. At \((k,n,r,w)=(4,8,3,6)\), the sector energies
reconstruct every row-subset level energy and the total Parseval energy
exactly. The sector-norm bound is \(6.3613\), compared with the exact
variance factor \(1.2657\) and the rejected weight-level bound \(8.2071\).
At \((5,10,3,10)\), the sector-norm bound is \(2.2252\), compared with the
exact factor \(1.0392\).

The remaining accumulator-side problem is therefore the scalable evaluation
of (23). It is a two-sequence, memory-one transfer problem: \(c_w(s)\) and
\(c_w(t)\) depend on the transition weights of the two subset indicators.
This is explicit, but no length-512 outward implementation is present.

## Walsh-conjugated low-shell route

There is a second exact route that avoids (23). Let \(W\) be the unnormalized
Walsh matrix on \(\mathbb F_2^n\), and define

\[
 f_w(z)=\mathbf 1\{\operatorname{wt}(Az)=w\}.
 \tag{25}
\]

Fourier inversion gives \(c_w=Wf_w\). Since \(W\) commutes with every row
permutation, the normalized conjugate

\[
 D=2^{-n}WCW
 \tag{26}
\]

has the same irreducible sectors as \(C\). If \(H^{(j)}\) is the Walsh
matrix in sector \(j\), its block is

\[
 D^{(j)}=2^{-n}H^{(j)}C^{(j)}H^{(j)\mathsf T}.
 \tag{27}
\]

The variance identity becomes

\[
 \operatorname{Var}(N_w)=2^{2k-n}f_w^{\mathsf T}Df_w.
 \tag{28}
\]

At rate one half, the scalar in (28) is one. Unlike \(c_w\), the vector
\(f_w\) has a short level profile in the low-distance shells. Put

\[
 g_{w,p}=|\{z:\operatorname{wt}(z)=p,
                 \operatorname{wt}(Az)=w\}|.
 \tag{29}
\]

The run count already used for (14) computes every \(g_{w,p}\) exactly. Its
support satisfies

\[
 p\le \min\{2w,\,2(n-w)+1\}.
 \tag{30}
\]

For

\[
 d_{p,q}=\max_j|D^{(j)}_{p,q}|,
\]

Cauchy--Schwarz at the level boundaries gives the valid exact-arithmetic
bound

\[
 f_w^{\mathsf T}Df_w
 \le \sum_{p,q}\sqrt{g_{w,p}}d_{p,q}\sqrt{g_{w,q}}.
 \tag{31}
\]

This relaxation still discards signed cross-level cancellation. Its advantage
is that (30) aligns the low-distance shell with a small part of the conjugated
operator. The script `probe_walsh_conjugated_level_bound.py` implements
(26)--(31). At \(n\le10\), a full Walsh-matrix check verifies
\(c_w=Wf_w\) and (28) directly.

The current nondirected diagnostics are:

| \(k\) | \(n\) | \(r\) | \(w\) | exact variance/mean | bound (31) |
|---:|---:|---:|---:|---:|---:|
| 4 | 8 | 3 | 6 | 1.26569 | 3.86162 |
| 24 | 48 | 5 | 31 | 2.79061 | 452.713 |
| 32 | 64 | 5 | 40 | 61.7570 | 25053.4 |
| 32 | 64 | 7 | 40 | 4.64727 | 1498.70 |
| 32 | 64 | 5 | 5 | 1.00217 | 1.68033 |
| 32 | 64 | 7 | 5 | 0.999994 | 1.07079 |
| 64 | 128 | 9 | 10 | not computed | 10.4285 |

The central-shell rows show that (31) is not a universal replacement for the
signed sector quadratic form. The low-shell rows are the relevant comparison:
they preserve almost all of the exact value at length 64, and the matched
length-128 diagnostic remains below the factor-512 target. These observations
do not extrapolate to length 512.

The early plan assumed that the exact constituent size could replace
shellwise central caps. The implemented finite transfer instead consumes a
central-shell envelope. The final proof therefore certifies the defect shells
first and later adds a cancellation-free central bound. Formula (30) limits
the primal levels for shell 42 to 84. The correlation table (23) remains a
valid fallback, but the completed certificate does not need it.

The conjugated blocks also have a direct primal factorization. For an ordered
message pair \((x,y)\), let \(P_{x,y}\) be its \(2\)-by-\(2\) one-row output
probability matrix. Independence of the sparse-map rows makes its full pair
kernel \(P_{x,y}^{\otimes n}\). Schur--Weyl decomposition gives the exact
sector block

\[
 \left(P_{x,y}^{\otimes n}\right)^{(j)}
 =\det(P_{x,y})^j
   \operatorname{Sym}^{n-2j}(P_{x,y}),
 \tag{32}
\]

where the symmetric power is written in the normalized Hamming-weight basis.
For the fixed-weight sparse row,

\[
 \det(P_{x,y})
 =\frac{\beta(x+y)-\beta(x)\beta(y)}4.
 \tag{33}
\]

Summing (32) over all ordered message pairs and subtracting the outer product
of the mean in sector zero gives \(D^{(j)}\), up to the general scalar
\(2^{n-2k}\). That scalar is one at rate one half. Formula (32) explains the
suppression of high sectors through powers of the determinant. It also
provides a route to the required low-level entries without first constructing
and conjugating every dense dual block.

The script `verify_primal_schur_power_factorization_small.py` aggregates all
256 ordered message pairs at \((k,n,r)=(4,8,3)\). Its blocks agree with the
Walsh-conjugated exact dual blocks to relative binary64 residual below
\(1.4\times10^{-15}\). The target implementation below aggregates (32) over
the \(\binom{k+3}{3}\) message-pair types without evaluating a dense
symmetric power separately for every type.

Only three sector diagonals are needed. For a nonnegative \(2\)-by-\(2\)
matrix \(P=(a,b;c,d)\), define

\[
 s_{m,u}(P)=
 \left(\operatorname{Sym}^{m}(P)\right)_{u,u}.
\]

Fix an output subset of size \(u+1\) among \(m+2\) coordinates. Sum over one
chosen coordinate inside the subset and one outside it. Restrict the input
subset to choose exactly one of those coordinates. Counting the same terms
inside \(s_{m+2,u+1}(P)\) gives

\[
 |\det P|s_{m,u}(P)
 \le (ad+bc)s_{m,u}(P)
 \le s_{m+2,u+1}(P).
 \tag{34}
\]

For fixed level \(p\), equation (34) compares adjacent sector diagonals in
(32). Every even sector has nonnegative summands. Hence each even sector
\(j\ge4\) is at most sector 2. For odd \(j\ge3\), positivity of the covariance
diagonal and the triangle inequality give

\[
 0\le D^{(j)}_{p,p}
 \le D^{(j-1)}_{p,p}
 \le D^{(2)}_{p,p}.
 \tag{35}
\]

Put

\[
 d_p=\max\{D^{(0)}_{p,p},D^{(1)}_{p,p},D^{(2)}_{p,p}\}.
\]

Each covariance block is positive semidefinite, so
\(|D^{(j)}_{p,q}|\le\sqrt{D^{(j)}_{p,p}D^{(j)}_{q,q}}\). Combining this fact,
(29), and (35) gives

\[
 f_w^{\mathsf T}Df_w
 \le\left(\sum_p\sqrt{g_{w,p}d_p}\right)^2.
 \tag{36}
\]

For sector one, the computation retains only positive determinant
contributions. Sector two already has nonnegative summands. Sector zero can
also be bounded by retaining positive locally centered pair contributions,
but that relaxation is useful only at the first defect shell.

The scripts `probe_primal_schur_diagonal_target.py` and
`evaluate_primal_schur_diagonal_shell_bound.py` implement this reduction.
They stream all 2,862,209 message-pair types. At
\((k,n,r,w)=(256,512,33,42)\), the monotone diagnostic from (36) is

\[
 \frac{\operatorname{Var}(N_{42})}{\mathbb E[N_{42}]}
 \le 31.0225.
\]

The reduction to (36) is exact. The displayed value uses nondirected
floating-point arithmetic and is not a certificate. Its factor-512 slack is
approximately \(16.5\). A direct signed diagnostic gives 28.2284, but global
sector-zero subtraction becomes numerically unstable at high primal levels.

The script `certify_primal_schur_diagonal_outward.py` now implements directed
binary64 interval arithmetic for the three one-sided formulas. At levels 80
and 84, the outward sector-zero bounds are 1.277282466 and 206.661452.
Substituting those two endpoints into the otherwise nondirected shell-42
calculation raises the variance factor to 33.8177. This intermediate mixed
calculation was only a diagnostic; later receipts cover levels 1 through 159.

The same sector-zero relaxation fails beyond the first shell. At levels 90
and 100, its outward bounds are approximately \(2.72\cdot10^6\) and
\(6.97\cdot10^{12}\). Stable signed diagnostics remain 1.00000023 and
1.00000048. Thus the growth is loss in the one-sided relaxation, not evidence
of large sector-zero variance.

Sector zero has a separate positive decomposition. Put

\[
 q_x=\frac{1-\beta(x)}2,
 \qquad
 \delta_{x,y}=\det(P_{x,y}),
\]

and define

\[
 h_{j,p}(q)=[z^p](1-z)^j(1-q+qz)^{n-j}.
 \tag{37}
\]

Expanding the product pair kernel around its independent-marginal kernel
gives the exact Hoeffding identity

\[
 \operatorname{Cov}\!\left(
   \mathbf1\{\operatorname{wt}(Ex)=p\},
   \mathbf1\{\operatorname{wt}(Ey)=p\}
 \right)
 =\sum_{j=1}^{n}\binom nj
   \delta_{x,y}^{\,j}h_{j,p}(q_x)h_{j,p}(q_y).
 \tag{38}
\]

The matrix \((\delta_{x,y})_{x,y}\) is a covariance matrix. The Schur product
theorem therefore makes each matrix
\((\delta_{x,y}^{\,j})_{x,y}\) positive semidefinite. After summing (38) over
all messages, every fixed-\(j\) term is nonnegative. Equation (38) is the new
sector-zero certificate target: it preserves the required cancellation
inside a positive quadratic form instead of discarding it pair by pair.
The script `verify_sector_zero_hoeffding_small.py` checks 2,304 pair-shell
instances of (38) in exact rational arithmetic at \((k,n,r)=(4,8,3)\). It
also checks all 72 fixed-order aggregates in that model.

The existing radial walk compresses each fixed-order quadratic form. For
\(0\leq t\leq k\), define

\[
 F_{j,\ell,p}(t)
 :=\sum_{a=0}^{k}K_a^{(k)}(t)
      h_{j,p}(q_a)\beta(a)^{j-\ell},
 \qquad q_a:=\frac{1-\beta(a)}2.
 \tag{39}
\]

If \(v_t(\ell)\) is the per-vector probability that \(\ell\) sampled sparse
rows XOR to one fixed vector of weight \(t\), then the order-\(j\) aggregate
in (38) equals

\[
 \binom nj4^{-j}
 \sum_{\ell=0}^{j}(-1)^{j-\ell}\binom j\ell
 \sum_{t=0}^{k}\binom kt v_t(\ell)
     F_{j,\ell,p}(t)^2.
 \tag{40}
\]

Thus (40) replaces \(2^{2k}\) ordered message pairs by \(k+1=257\) radial
states. The exact small verifier checks all 72 instances of (40) as well as
the pairwise identity. The alternating sum over \(\ell\) still requires Arb
intervals; fixed-order nonnegativity does not justify dropping its negative
terms.

There is also an exact polynomial form for the transform in (39). Substituting
\(q_a=(1-\beta(a))/2\) into (37) gives

\[
 h_{j,p}(q_a)
 =2^{-(n-j)}\sum_{s=0}^{n-j}
   \binom{n-j}{s}K_p^{(n)}(j+s)\beta(a)^s.
 \tag{40a}
\]

The Krawtchouk transform of \(\beta(a)^d\) is \(2^k v_t(d)\). Therefore

\[
 F_{j,\ell,p}(t)
 =2^{k-(n-j)}\sum_{s=0}^{n-j}
   \binom{n-j}{s}K_p^{(n)}(j+s)v_t(s+j-\ell).
 \tag{40b}
\]

The exact small verifier checks all 1,980 instances of (40b) at
\((k,n,r)=(4,8,3)\). A polynomial-convolution Arb implementation reproduces
the complete small-model diagonal exactly. At the target parameters, however,
it took 81.4 seconds through order 128 at level 100, compared with 44.2
seconds for the scalar Krawtchouk-transform implementation. Equation (40b)
is therefore a proved identity and a possible multi-level primitive, but the
current single-level polynomial implementation is not the production route.

The script `evaluate_sector_zero_hoeffding_radial_arb.py` evaluates (40)
with Arb intervals. At \((k,n,r,p)=(256,512,33,100)\), the complete sum over
orders 1 through 512 is

\[
 D^{(0)}_{100,100}
 \in
 1.0000035547087893253531579783751399012502409269071873\ldots
 \mathbin{\pm}2.14\cdot10^{-111}.
\]

The first run enclosed orders 1 through 511 but encountered a Python-FLINT
edge case at order 512: generic exponentiation returned `nan` when it squared
a zero-centered Arb interval. Direct interval multiplication computes the
same square and remains finite. A separate 2,048-bit order-512 receipt gives
a contribution below \(4.27\cdot10^{-46}\). The merger checks all 512 orders
and produces the displayed complete outward enclosure. This closes sector
zero at level 100, but not at the other levels required by the shell bound.

Odd degree supplies a useful diagnostic grouping. Averaging a pair
contribution over

\[
 (x,y),\ (x+\mathbf1,y),\ (x,y+\mathbf1),\
 (x+\mathbf1,y+\mathbf1)
\]

does not change the total covariance. At level 100, this four-term grouping
reduces the nondirected positive and absolute mass from about
\(1.3\cdot10^{12}\) to 1.0000032. The grouped terms are not always
nonnegative in small exact models. Orbit positivity is therefore false.
Nevertheless, an outward enclosure of each signed four-term orbit followed
by summation of its positive upper endpoint is a valid upper bound; it does
not assume orbit positivity. A directed binary64 prototype reproduces the
small \((4,8,3,p=2)\) positive-orbit bound as 0.8136835098268285, compared
with the nondirected value 0.8136835098266602 and exact diagonal
0.61776065826416015625. Canonicalizing under both message complements and
message swap reduces the target calculation from 2,862,209 types to 361,985
representatives. Constructing the determinant from its exact integer
numerator is essential; forming it by subtracting rounded biases gives a
spurious upper bound of \(7.51\cdot10^{18}\). With the exact determinant,
the directed target calculation proves

\[
 D^{(0)}_{100,100}\le 2.45522037042628.
 \tag{40c}
\]

The independent Arb calculation gives the sharper value
1.0000035547087893..., so (40c) is conservative but useful. The canonical
Python run takes about 480 seconds per level. The inequality is now viable;
the next implementation task is a compiled coefficient recurrence or a
multi-level calculation that amortizes the type traversal.

A full sector-zero trace bound was also tested and rejected. Trace is
preserved by the normalized Walsh conjugation, but the dual sector-zero trace
contains enormous low-level row-collision energy. Moreover, direct
long-double evaluation suffers catastrophic cancellation at the target
parameters. The trace cannot replace the required primal-band calculation.

## Remaining certificate target

A sharper valid bound is

\[
 c_w^{\mathsf T}Cc_w
 \leq
 \sum_j \lambda_{\max}(C^{(j)})
          \lVert P_jc_w\rVert_2^2,
 \tag{41}
\]

where \(P_j\) is the orthogonal projector onto sector \(j\). A signed
evaluation of each block quadratic form can be sharper still.

The first sufficient target is

\[
 \operatorname{Var}(N_w)\leq512\,\mathbb E[N_w]
 \tag{42}
\]

for every shell used by the existing 10.9-percent cap interface. The final
program may instead compare each certified variance directly with its actual
cap slack.

## Moderate-size evidence

The script probe_one_stage_variance_scaling.py evaluates the exact one-word
and pair formulas with nondirected floating-point arithmetic. Selected worst
positive-shell ratios are:

| \(k\) | \(n\) | \(r\) | maximum \(\operatorname{Var}(N_w)/\mathbb E[N_w]\) |
|---:|---:|---:|---:|
| 12 | 24 | 3 | 1.2349 |
| 16 | 32 | 3 | 2.0502 |
| 20 | 40 | 3 | 6.6102 |
| 24 | 48 | 3 | 38.6851 |
| 24 | 48 | 5 | 2.7906 |
| 32 | 64 | 5 | 61.7570 |
| 32 | 64 | 7 | 4.6473 |

These data show that fixed degree \(3\) or \(5\) is not an adequate proxy
for \(r=33\). Raising the degree sharply reduces the variance in the tested
models. No row in the table extrapolates to \(k=256\).

## Resolution and remaining extensions

The finite one-stage obligation is complete. Outward receipts cover sector
zero at levels 1 through 159, sector one at levels 1 through 159, and sector
two at levels 2 through 159. The diagonal relaxation built from those bounds
fails from the middle of the defect band; a lower-bound certificate first
crosses 512 at shell 65. The failed relaxation is retained as a proved dead
end.

The replacement argument retains one dominant row-character bias for each
ordered message pair. A pointwise likelihood ratio controls the two omitted
biases. Arb evaluates the reference laws. On central shells, Arb also
encloses each reference deviation from one before binary64 arithmetic. This
avoids subtracting nearly equal floating-point values.

The final receipts prove the factor-512 target on shells 42 through 79 and
433 through 470. They also bound every central shell from 80 through 432
well enough for Cantelli caps. Bounded rank-tested setup and the existing
RandomStepConv-M22 transfer yield 41.3621 bits of end-to-end margin at 10.9%
distance. See `../ONE_STAGE_SPARSE_EA_RANDOMSTEPCONV_CERTIFICATE.md`.

Two extensions remain outside this theorem:

1. transfer the realized-spectrum proof to the lower-XOR multistage outer;
2. replace RandomStepConv-M22 with RM2Sub or another practical inner.
