# Dense-occupation transfer for the RM2Sub-S19 inner

## Question

The fixed-occupation continuum controls every constant number of active
outer blocks. It does not control a sequence \(Q=Q_N\to\infty\). This note
gives a finite-\(L\) transfer that covers every such sequence.

The construction uses the random-outer RM2Sub ensemble from
`RANDOM_OUTER_RM2SUB_RAMP.md`. Set

\[
  n:=2^{19},\qquad M:=n-1,\qquad t:=128.
\]

Let \(A:\mathbb F_2^{19}\to\mathbb F_2^{128}\) be the frozen state-output
map, and let \(C=A^{\mathsf T}\) be the state-update map.

## Bernoulli reference law

Fix a support of \(Q\) active outer blocks and put

\[
  \alpha:=Q/L.
\]

For each active block, a fixed nonzero local message maps to a uniform
nonzero word in \(\mathbb F_2^B\). Replace that word by a uniform binary word.
This replacement costs

\[
  (1-2^{-B})^{-Q}.
  \tag{1}
\]

In each transposed region, the region permutation makes the \(Q\) active
block positions a uniform \(Q\)-subset of \([L]\). The values at those
positions are independent fair bits.

Fix \(p\in(0,1)\). Under the reference law, mark each region position
independently with probability \(p\). Give every marked position a fair bit,
and give every unmarked position the value zero. Conditioned on exactly
\(Q\) marks, this law is the required region law. The conditioning cost for
all \(B\) regions is

\[
  b_{L,Q}(p)^{-B},
  \qquad
  b_{L,Q}(p):=\Pr[\operatorname{Bin}(L,p)=Q].
  \tag{2}
\]

Before conditioning, every inner input bit is independent
Bernoulli-\(\beta\), where \(\beta=p/2\). The inner transfer is therefore
homogeneous across all \(N/t\) epochs.

## Exact constituent moments

Fix \(z\in(0,1)\). Define

\[
  u:=1-\beta+\beta z,
  \qquad
  v:=\beta+(1-\beta)z,
  \qquad
  r:=1-\beta-\beta z.
  \tag{3}
\]

The image of \(A\) has weight classes \(d_0,\ldots,d_6\), with

\[
  (d_i)_{i=0}^6=(0,48,56,64,72,80,128)
\]

and multiplicities

\[
  (c_i)_{i=0}^6
  =(1,5040,110848,292510,110848,5040,1).
  \tag{4}
\]

For \(a,b\in\mathbb F_2^{19}\), use the standard binary dot product. Define
the exact association matrix

\[
  \mathsf K_{ij}
  :=
  \sum_{\substack{a:\operatorname{wt}(A(a))=d_i}}
  \sum_{\substack{b:\operatorname{wt}(A(b))=d_j}}
  (-1)^{\langle a,b\rangle}.
  \tag{5}
\]

An integer Walsh transform gives

\[
\mathsf K=
\begin{pmatrix}
1&5040&110848&292510&110848&5040&1\\
5040&-6408&41280&-37596&-320&-3256&1260\\
110848&41280&-39680&-93216&-33024&-64&13856\\
292510&-37596&-93216&-76862&-80352&-4484&0\\
110848&-320&-33024&-80352&19712&-3008&-13856\\
5040&-3256&-64&-4484&-3008&7032&-1260\\
1&1260&13856&0&-13856&-1260&-1
\end{pmatrix}.
\tag{6}
\]

For class \(i\), set

\[
  D_i:=u^{128-d_i}v^{d_i},
  \qquad
  F_i:=u^{128-d_i}r^{d_i}.
  \tag{7}
\]

Define

\[
  k_{00}:=\frac1n\sum_i c_iF_i,
  \qquad
  k_{01}:=u^{128}-k_{00},
  \tag{8}
\]

\[
  \overline d
  :=
  \frac{n^{-1}F^{\mathsf T}\mathsf KD-k_{00}D_0}{k_{01}},
  \tag{9}
\]

and

\[
  m:=\frac{\sum_i c_iD_i-D_0}{M},
  \qquad
  \kappa:=\frac{M}{M-1}.
  \tag{10}
\]

Here \(k_{00}\) is the tilted mass of \(C(X)=0\) from state zero.
The quantity \(k_{01}\) is the complementary mass. Conditioned on that
transition, \(\overline d\) is the exact output moment in the next epoch.
The quantity \(m\) is the exact output moment for a uniform nonzero entering
state.

## Three-state envelope

Use the state classes \(Z,D,R\).

- \(Z\) is the zero state.
- \(D\) is the tilted nonzero state produced when an epoch leaves \(Z\).
- \(R\) is a uniform or punctured-uniform nonzero state.

Rows index the entering class and columns index the leaving class. The matrix

\[
  T_3(p,z):=
  \begin{pmatrix}
    k_{00}&k_{01}&0\\
    \overline d/M&0&\overline d\\
    \kappa m/M&0&\kappa m
  \end{pmatrix}
  \tag{11}
\]

is an entrywise envelope.

The first row is exact. From any nonzero entering state, the fresh multiplier
makes \(\alpha Q\) uniform on \(\mathbb F_{2^{19}}^*\). If \(C(X)\ne0\), the
next state is zero with probability \(1/M\); conditioned on survival, its law
omits at most one nonzero value. Removing one value increases a nonnegative
uniform-live average by at most \(\kappa\). The duplicated termination mass
in (11) only enlarges the transfer.

## Four-state small-density envelope

The factor \(\kappa\) in (11) is harmless at fixed positive density. It is a
fixed loss as \(p\to0\). A four-state envelope charges it only after a
nonzero syndrome.

Use \(Z,D,U,P\), where \(U\) is uniform nonzero and \(P\) is
punctured-uniform nonzero. Define

\[
  a(p)
  :=\Pr[C(X)\ne0]
  =1-\frac1n\sum_i c_i(1-p)^{d_i}.
  \tag{12}
\]

For any normalized entering-state law independent of the current input,

\[
  z^{128}a(p)
  \le
  \mathbb E[z^{\operatorname{wt}(Y)}\mathbf1_{C(X)\ne0}]
  \le a(p).
  \tag{13}
\]

For \(0<p\le8/(5\cdot10^4)\) and \(z=1-(8/5)\alpha\), the total live
moments exceed \(a(p)\). Indeed, each total moment is at least \(z^{128}\),
where
\(z^{128}\ge(1-1.6\times10^{-4})^{128}>0.97\), whereas
\(a(p)\le\Pr[X\ne0]\le64p\le0.01024\). Thus

\[
T_4(p,z):=
\begin{pmatrix}
k_{00}&k_{01}&0&0\\
a/M&0&\overline d-z^{128}a&a(M-1)/M\\
a/M&0&m-z^{128}a&a(M-1)/M\\
\kappa a/M&0&\kappa(m-z^{128}a)&\kappa a(M-1)/M
\end{pmatrix}
\tag{14}
\]

is an entrywise envelope. The fresh multiplier makes the termination event
independent of the current output after \(Q,X\) are fixed. This gives the
factors \(1/M\) and \((M-1)/M\) in (14).

## Finite first-moment bound

Let \(Z_{D,Q}\) count nonzero messages supported on exactly \(Q\) outer
blocks whose output weight is at most \(D\). For either valid transfer
\(T\in\{T_3,T_4\}\),

\[
\begin{split}
\mathbb E[Z_{D,Q}]
\le{}&
\binom LQ(2^{B/2}-1)^Q(1-2^{-B})^{-Q}
b_{L,Q}(p)^{-B}z^{-D}\\
&\cdot e_Z^{\mathsf T}T(p,z)^{N/128}\mathbf1.
\end{split}
\tag{15}
\]

The expectation is over the independent random outer injections, coordinate
permutations, region permutations, and nonzero RM2Sub multipliers. Equation
(15) follows from (1), (2), the Chernoff inequality, and transfer
multiplication. It is finite and requires no continuum approximation.

## Exact small-density certificate

Set

\[
  p:=\frac85\alpha,
  \qquad
  z:=1-\frac85\alpha,
  \qquad
  0<\alpha\le10^{-4}.
  \tag{16}
\]

Let

\[
  w:=
  \left(1,2^{-10},2^{-10},\kappa2^{-10}\right)^{\mathsf T}.
\]

Exact rational Bernstein coefficients prove

\[
  T_4(p,z)w\le(1-99\alpha)w
  \tag{17}
\]

throughout the interval in (16). The three nontrivial row inequalities have
degrees \(256,256,511\). After removing their zero factors, every Bernstein
coefficient is strictly negative.

Let \(D_{\rm KL}\) denote binary relative entropy in natural units. The
normalized exponent in (15), before the support and local-limit prefactors,
is at most

\[
  E(\alpha)
  :=\frac{\alpha\ln2}{2}
    +D_{\rm KL}(\alpha\Vert p)
    +\frac1{128}\ln(1-99\alpha)
    -\frac{11}{100}\ln z.
  \tag{18}
\]

The same exact receipt proves

\[
  E(\alpha)
  \le
  -\frac{780897}{6665600}\alpha
  <-0.11715\alpha.
  \tag{19}
\]

`certify_rm2sub_dense_small.py` reproduces (17)--(19).
`rm2sub_dense_small_exact_d11.json` records the exact certificate.

## Compact-density certificate

For fixed \(p,z\), define

\[
  \Psi(\alpha;p,z)
  :=\frac{\alpha\ln2}{2}
    +D_{\rm KL}(\alpha\Vert p)
    +\frac1{128}\ln\rho(T_3(p,z))
    -\frac{11}{100}\ln z.
  \tag{20}
\]

The map \(\alpha\mapsto\Psi(\alpha;p,z)\) is convex. Therefore a fixed
Collatz witness whose outward exponent bounds are negative at both endpoints
certifies the complete interval between them.

`certify_rm2sub_dense_compact_interval.py` partitions
\([10^{-4},1]\) into 31 intervals. Each interval has rational \(p,z\) and a
rational positive Collatz vector. All transfer operations and logarithms use
100-decimal-digit outward interval arithmetic. The largest accepted endpoint
upper bound is

\[
  -4.368933504109312\times10^{-7}
\]

natural units per output bit. The complete receipt is
`rm2sub_dense_compact_interval_d11.json`.

## Proof status

- **Proved:** the Bernoulli reduction, both entrywise transfers, the finite
  first-moment bound, the exact small-density interval, and the outward
  compact-density interval.
- **Imported exact facts:** \(C=A^{\mathsf T}\), the image spectrum in (4),
  and the fixed state dimension.
- **Diagnostic only:** optimized binary64 saddles in
  `rm2sub_dense_occupation_d11.json`. They propose witnesses but certify no
  interval.
