# Directed certificate for the middle input-weight class

The symmetric width-sixteen, state-eight candidate with GRS[256,176] has a directed bound below `2^-344` for every input weight in `[.03N,.24N)` at K20. This note certifies that class only. The separate [aggregate replay](FIELD256_K20_CERTIFICATE.md) now joins it with the matching sparse and dense bounds, giving 119 whole-code bits.

The replay is [state8_density_certificate.py](theory/state8_density_certificate.py). Its eight tests are in [state8_density_certificate_tests.py](theory/state8_density_certificate_tests.py). The ordinary replay uses frozen rational data and does not invoke an optimizer.

## Construction and event

Use the symmetric maps, four packet phases, and scalar `eta=2` defined in [STATE_CAPACITY_REPAIR.md](STATE_CAPACITY_REPAIR.md). In each packet, the emission matrix `A` is `16 x 8`, the injection matrix is `C=A^T`, and `CA=0`. Both maps are superregular. The input diagonal and eight state-diagonal coefficients are fresh independent elements of `GF(256)^*` in every packet. The route consists of independent regional permutations, applied in one global routing stage. The state starts at zero; the final state is discarded.

The fixed outer code has length `b=256`, dimension `176`, and redundancy `p=80`. At requested source length `2^20`, the replay uses

\[
L=6016,\qquad K_{\mathrm{pad}}=176L=1058816,
\qquad N=256L=1540096,
\qquad e=\lceil N/10\rceil-1=154009.
\]

The padded rate is `11/16`. Restricting the source to the requested length cannot introduce a new low-weight codeword. Here `W` denotes the support size of the vector entering the inner encoder. Routing and the input diagonal preserve this support size.

The bounded event is existence of a nonzero source vector whose inner input has `46203 <= W <= 369623` and whose output has weight at most `e`. Probability is over the ideal independent routing and diagonal choices just specified. The rank argument counts existence in a source kernel; it is not an ordinary word first moment.

## Positive packet cover

The two transfer states mean actual zero and actual nonzero witness state. Let `T_w(y,tau)` be the packet cover for a fixed input placement of size `w`, where `0<y<=1` and `1<=tau<=255`. Its zero-state entries are

\[
(T_w)_{00}=\mathbf1_{w=0}+\mathbf1_{w>8}y^w\tau^{-8},
\qquad (T_w)_{01}=\mathbf1_{w>0}y^w.
\]

For selected output zeros, write `i` for zeros at active input positions and `z` for zeros outside them. The nonzero-state continuation uses

\[
(T_w)_{11}=\sum_{i=0}^{w}\sum_{z=0}^{\min(7,16-w)}
\binom wi\binom{16-w}z
y^{16-i-z}(1-y)^{i+z}\tau^{-i}\Gamma_z(\tau).
\]

Here `Gamma_z(tau)` is the maximum of `rank_potential_bound(r,z,tau,256,8)` over `1<=r<=8`. That function evaluates the exact rational rank-moment polygon for diagonal state refresh. It requires superregular emission and makes no uniform-subspace assumption.

For `w>0`, the return entry is the minimum of two complete bounds. The pure bound is

\[
P_w=y^{m_w}\tau^{-8}\sum_{k=0}^{\min(w-1,7)}\binom8k\tau^{-k},
\]

where `(m_1,...,m_16)=(9,9,9,9,8,6,4,2,1,1,1,1,1,1,1,1)`.
The joint bound is

\[
J_w=\sum_{\substack{0\le i\le w,\ 0\le z\le\min(7,16-w)\\16-i-z\ge m_w}}
\binom wi\binom{16-w}z
y^{16-i-z}(1-y)^{i+z}
\tau^{-\max(8,i+z,w-16+i+2z+1)}
\sum_{k=0}^{\min(w-1,7-z)}\binom8k\tau^{-k}.
\]

Thus `(T_w)10=min(P_w,J_w)`; for `w=0` this entry is zero. The joint matrix receives one diagonal-pivot charge. No separate root factor is multiplied into this return entry. Its proof, including the weighted incoming state flats, is recorded with [state8_joint_return.py](theory/state8_joint_return.py) in `STATE_CAPACITY_REPAIR.md`.

Set

\[
T(x,y,\tau)=\sum_{w=0}^{16}\binom{16}w x^wT_w(y,\tau).
\]

Each packet placement is counted once. All entries are nonnegative, so entrywise upper bounds remain valid under products and sums.

## Finite source and route comparison

For a fixed exact input weight `W`, put `alpha=W/N`. The previously audited finite comparison gives

\[
B_W=\frac{(L+1)^{256}}{\tau(x\lambda)^W y^e}
\left[(1-\alpha)^{256}
+\tau^{-80}(1-\alpha+\tau\lambda\alpha)^{256}\right]^L
e_0T(x,y,\tau)^{N/16}\mathbf1.
\]

The factor `(L+1)^256` retains the cost of conditioning regional Bernoulli subsets on their modal counts. The factor `tau^-1` comes from the nonzero-kernel rank threshold. It is not an additional projective-word normalization. [STATE8_DENSITY_AUDIT.md](STATE8_DENSITY_AUDIT.md) records the conditioning and source-row comparison.

For an integer cell `[a,c]`, the replay fixes one rational quadruple `(x,y,tau,lambda)`. Put `beta=tau*lambda`. The row bracket is bounded by

\[
R_{a,c}=(1-a/N)^{256}
+\tau^{-80}\bigl(1+(\beta-1)\alpha_*\bigr)^{256},
\]

where `alpha_*=c/N` if `beta>=1`, and `alpha_*=a/N` otherwise. Also choose `W_*=a` if `x*lambda>=1`, and `W_*=c` otherwise. The cell contribution is at most

\[
(c-a+1)\frac{(L+1)^{256}R_{a,c}^L}
{\tau(x\lambda)^{W_*}y^e}
e_0T^{N/16}\mathbf1.
\]

This bounds the sum over every integer weight in the cell; it does not interpolate between point evaluations.

## Rational witness selection and replay

The source contains twelve frozen anchor densities and rational quadruples. For each cell, the replay considers its three nearest anchors. If an anchor has density `d` and tilt `lambda_d`, define `kappa=lambda_d*d/(1-d)`. At the cell midpoint `m=(a+c)/(2N)`, choose

\[
\lambda=\kappa(1-m)/m.
\]

This is an exact rational choice of a free Chernoff parameter. It does not change the finite comparison. The packet matrix and its prefactor do not depend on `lambda`, so they can be reused from the anchor. The replay chooses the smallest outward upper endpoint among the three valid cell bounds.

The initial cells contain at most `1541` consecutive weights. A cell is bisected until its directed upper endpoint is at most `2^-256`. Explicit contiguity checks require accepted cells to partition the entire interval, including when Python optimization is enabled. Every operation uses eighty-digit decimal intervals with downward and upward rounding; the enlarged decimal exponent range prevents premature overflow or underflow. Exact `Fraction` arithmetic supplies all witness values, rank factors, and state-flat sums.

Run:

```sh
python -B research/workstreams/char2_field_spin/theory/state8_density_certificate.py
python -B research/workstreams/char2_field_spin/theory/state8_density_certificate_tests.py
```

The replay accepted `216` cells, each containing between `770` and `1541` weights. It covered all `323421` weights and returned

\[
\Pr[\text{bad word in the middle class}]
\le 1.6089746439007743265160334370839809292483925603594688520394855486970406809360837\cdot10^{-104}
<2^{-344}.
\]

The partition SHA-256 is `01c230076663d1bff17808ee26577d6cc1304f2761a0e599a2e124e072abecff`. All eight tests passed. They include an independent exact-Fraction packet evaluation, both endpoint-monotonicity cases, exact small matrix powers, gap-free refinement, and the complete replay.

For the intended class union, every nonzero outer row has at least `81` active coordinates. Hence more than `600` active rows imply `W>=601*81=48681`, above this interval's lower endpoint. The sparse certificate can therefore overlap this middle class without a missing input-weight range. The independent dense class begins at `W=369624`. Matching and summing those two other certificates is a separate whole-code step. No implementation-speed claim follows from this replay.
