# Uniform GF(2^32) diagonals: valid rank extension, unresolved middle bounds

Fresh GF(2^32) diagonals substantially strengthen the one-row bound, but they
do not move the sampled rate/state frontier of the existing Cauchy cover.
GRS[64,48] still fails the middle-weight screen with four or eight state
symbols. GRS[128,87] still fails that screen after reducing its state to four
or six symbols. These are unsuccessful upper bounds, not proofs that the
corresponding codes fail distance.

The new [screen](theory/medium_field_screen.py) changes the diagonal alphabet
and its admissible rank tilt. Existing certificates and production code are
unchanged. All numerical results below are floating evaluations of sampled
classes; none is a whole-code certificate or runtime measurement.

## Construction and probability space

Fix an embedding of F0=GF(256) into F=GF(q), where q is either 256 or 2^32.
Fix a systematic GRS[b,k] generator over F0, and interpret its entries in F.
A message consists of L source rows in F^k. Each outer coordinate has one
independent uniform permutation of the L rows. Together these permutations
form exactly one regional routing stage. Every routed coordinate receives
an independent scalar from F without zero.

Use the symmetric width-sixteen maps from
[STATE_CAPACITY_REPAIR.md](STATE_CAPACITY_REPAIR.md). For state size s<=8,
retain the first s columns of its emission matrix A and the matching rows
of C=A^T. Retain the four packet phases and eta=2. The maps are superregular
and satisfy CA=0. Every packet independently samples a diagonal D with s
entries from F without zero. The packet recurrence is

\[
 R:=DS,\qquad Y:=X+\eta AR,\qquad S':=R+CX.
\]

The state starts at zero and the final state is discarded. All permutations,
input scalars, and state scalars are independent setup randomness. Setup is
sampled once and then fixed for every message. Failure means existence of a
nonzero message whose output weight is at most e=ceil(bL/10)-1.

For q=2^32, the diagonals range over all 2^32-1 nonzero field elements.
This distribution differs from sampling GF256 diagonals and extending the
resulting code. The new randomness can change which outer messages vanish
under a given set of output equations.

A 128-bit application record can hold four F coordinates. Apply the same
sampled scalar encoder to all four coordinates, sharing its route and
coefficients. Every nonzero record word has a nonzero component word.
Its record support contains that component's support. Thus a scalar distance
theorem would imply the same record distance without a factor-four failure
union. This representation needs no general GF(2^128) multiplication.

## Which proof steps extend

Every fixed matrix rank over F0 remains equal over F. A nonzero maximal
minor remains nonzero under the embedding, and larger minors remain zero.
Consequently the GRS shortening dimensions, Cauchy minors, CA=0 identity,
and previously verified support-rank laws remain valid. Restricting the
state takes a principal submatrix of CA, so it preserves CA=0.

The full packet map from (R,X) to (Y,S') is invertible. Given Y and S',
recover R=S'+CY and then X=Y+eta AR. The identity CA=0 justifies the first
equation. The nonzero state diagonal is also invertible, so the causal
encoder is injective for every setup.

The random steps require a separate argument because their distribution
has changed. Fix an exact outer support and a preceding common span W of
zero equations. All current forms are linear forms over F modulo W.
The following arguments condition only on previously exposed randomness.

* A fresh input equation has the form d*l+c=0, with l required nonzero.
  If c is independent of l, the equation adds rank one. Otherwise at most
  one nonzero value of d leaves the required form alive. Therefore the
  expected surviving rank potential contracts by tau^-1 whenever
  1<=tau<=q-1. The argument uses the actual alphabet size q-1.
* Let U be an r-dimensional incoming state image and A_Z contain z selected
  emission rows. For every nonzero u in U, superregularity gives
  Pr[A_Z D u=0] <= (q-1)^(-z). If u has at most z nonzero coordinates,
  this probability is zero. Otherwise fixing all but z active diagonal
  entries leaves at most one solution for the remaining entries.
* Write H=dim(U intersect ker(A_Z D)). Counting nonzero vectors gives
  E[q^H-1] <= (q^r-1)/(q-1)^z. The concave polygon through
  (0,0) and (q^h-1,tau^(h-r)) bounds
  E[tau^(H-r); H>0]. Its slopes decrease for 1<=tau<=q-1.
  This is exactly `rank_potential_bound(r,z,tau,q,s)`.
* Return covers fix the incoming coordinate-zero flat before exposing
  current diagonals. Their numerical joint equation matrix has the same
  rank bounds over F. Each fresh diagonal pivot then uses the first
  argument. No second root probability is multiplied into a joint return.

Define Gamma_z as the maximum polygon bound over 1<=r<=s. The maximum is
uniform over the common equation history. It makes no assertion that D
produces a uniform distribution on state subspaces. Positive transfers
therefore compose under the same history argument as the GF256 reference.

For completeness, let w be the packet input weight and 0<y<=1 an output
mark. The continuation entry is

\[
 (T_w)_{11}=
 \left(y+\frac{1-y}{\tau}\right)^w
 \sum_{z=0}^{\min(s-1,16-w)}
 \binom{16-w}{z}y^{16-w-z}(1-y)^z\Gamma_z.
\]

The zero-state entries are

\[
 (T_w)_{00}=\mathbf1_{w=0}+\mathbf1_{w>s}y^w\tau^{-s},
 \qquad (T_w)_{01}=\mathbf1_{w>0}y^w.
\]

For w>0, let m_w=max(1,17-s-w) when s<8. For s=8, retain the verified
sequence (9,9,9,9,8,6,4,2,1,1,1,1,1,1,1,1). The return entry is the
minimum of the following two complete bounds:

\[
 P_w=y^{m_w}\tau^{-s}
       \sum_{a=0}^{\min(w-1,s-1)}\binom{s}{a}\tau^{-a},
\]

\[
 J_w=\sum_{\substack{0\le i\le w,\ 0\le z\le\min(s-1,16-w)\\
                         16-i-z\ge m_w}}
 \binom wi\binom{16-w}z y^{16-i-z}(1-y)^{i+z}
 \tau^{-\max(s,i+z,w-16+i+2z+1)}
 \sum_{a=0}^{\min(w-1,s-1-z)}\binom{s}{a}\tau^{-a}.
\]

The entry is zero for w=0. Both covers retain one common equation span;
their minimum does not multiply probabilities from two uses of a diagonal.
Smaller states do not inherit the special eight-state returned-weight law.

The finite source and route factors also extend. An exact outer row support
of size h has shortened dimension max(0,h-(b-k)) over F. The input rank
argument uses this dimension, rather than counting q^k words. Independent
regional permutations have the same without-replacement law for both fields.
Thus the active-row and density formulas of
[STATE8_DENSITY_CERTIFICATE.md](STATE8_DENSITY_CERTIFICATE.md) apply after
substituting this packet transfer, b, k, and s. The script retains their
finite combinatorial factors and the rank-threshold factor tau^-1.
The density screen uses the older (L+1)^b conditioning factor, not the
sharper modal factor in the completed small-size certificates. Applying
that refinement would improve both fields equally; it cannot account for
the tens-of-thousands-of-bits deficits displayed below.

Only Gamma and the allowed interval for tau depend on q in these formulas.
In particular, Gamma_z>=tau^-z: take the full state image U=F^s, for which
the z root equations have deterministic rank z. Increasing q cannot remove
this contribution. A construction change that weakens a rank bound needs
a new transfer; increasing q does not automatically justify the old one.

## Paired finite screens

All screens use requested size K=2^20, L rounded to a multiple of 64, and
e=ceil(N/10)-1. GRS[64,48] has L=21888 and N=1400832.
GRS[128,87] has L=12096 and N=1548288. An active-row class includes every
support with exactly that many nonzero outer rows. A density entry evaluates
the single integer weight round(alpha*N).

The table reports margin -log2(B) for the evaluated bound B. Negative
margins do not establish probability decay. The final column compares q32
and q256 at the *same q256 witness*, isolating the field change from the
numerical optimizer. The q32 search also starts from that witness.

| Outer | State | Class | GF256 margin | GF(2^32) margin | Gain at same witness |
|---|---:|---|---:|---:|---:|
| [64,48] | 4 | 1 active row | 60.362 | 293.025 | 0.002 |
| [64,48] | 4 | 4 active rows | 103.753 | 103.759 | 0.006 |
| [64,48] | 4 | 16 active rows | -266.629 | -266.608 | 0.020 |
| [64,48] | 4 | alpha=.03 | -32656.277 | -32655.241 | 1.036 |
| [64,48] | 4 | alpha=.10 | -138421.734 | -138421.724 | 0.009 |
| [64,48] | 4 | alpha=.20 | -23649.661 | -23649.339 | 0.322 |
| [64,48] | 8 | 16 active rows | 76.499 | 76.518 | 0.020 |
| [64,48] | 8 | alpha=.03 | -20682.813 | -20681.280 | 1.529 |
| [64,48] | 8 | alpha=.10 | -64351.478 | -64349.156 | 2.322 |
| [64,48] | 8 | alpha=.20 | -23917.946 | -23917.656 | 0.290 |
| [128,87] | 4 | 16 active rows | 135.678 | 135.724 | 0.046 |
| [128,87] | 4 | 64 active rows | -551.722 | -551.559 | 0.163 |
| [128,87] | 4 | alpha=.03 | -9330.583 | -9329.414 | 1.169 |
| [128,87] | 4 | alpha=.10 | -68089.567 | -68089.556 | 0.011 |
| [128,87] | 4 | alpha=.20 | 169372.944 | 169373.475 | 0.531 |
| [128,87] | 6 | 16 active rows | 569.992 | 570.038 | 0.046 |
| [128,87] | 6 | 64 active rows | 203.511 | 203.671 | 0.160 |
| [128,87] | 6 | alpha=.03 | -209.837 | -208.293 | 1.543 |
| [128,87] | 6 | alpha=.10 | -26088.992 | -26087.329 | 1.663 |
| [128,87] | 6 | alpha=.20 | 168878.116 | 168878.646 | 0.531 |

The first row exploits tau approximately 2^32 in the q32 search. The
difficult classes use tau approximately 7 through 14 near densities .03
and .10, and below 74 at density .20. Those tilts already fit GF256.
Additional diagonal entropy therefore does little for the existing bound.
Bounded numerical search does not prove these are global optima.

The old geometric failures also remain. In
[RATE_DESIGN_REVIEW.md](RATE_DESIGN_REVIEW.md), the GRS[32,24], state-four
window event costs at most 36 exit equations but supplies 37 source
dimensions. That event is independent of all coefficient values. The
[GRS[32,22] realized witnesses](WINDOW_DIMENSION_SEARCH.md) have the same
field-independent character. Fresh GF32 diagonals cannot repair either
obstruction. No absence of such witnesses is claimed for the new geometries.

## Reproduction and next gate

```sh
python -B research/workstreams/char2_field_spin/theory/medium_field_screen_tests.py
python -B research/workstreams/char2_field_spin/theory/medium_field_screen.py --compare --outer 64 48 --state 4 --active 1 4 16 --density .03 .1 .2 --iterations 40
python -B research/workstreams/char2_field_spin/theory/medium_field_screen.py --compare --outer 64 48 --state 8 --active 16 --density .03 .1 .2 --iterations 40
python -B research/workstreams/char2_field_spin/theory/medium_field_screen.py --compare --outer 128 87 --state 4 --active 16 64 --density .03 .1 .2 --iterations 40
python -B research/workstreams/char2_field_spin/theory/medium_field_screen.py --compare --outer 128 87 --state 6 --active 16 64 --density .03 .1 .2 --iterations 40
```

Seven tests pass. They compare exact GF256 packet entries with the existing
cover, check exact large-field factors and tilt domains, independently
convolve a tiny regional transfer, reproduce the finite global formulas,
and verify CA=0 for every phase and tested state restriction.

The next useful gate is a changed packet/route cover that improves the
multirow and middle classes, accompanied by a costed GF32 arithmetic model.
A cheaper GF32 Cauchy circuit with the same rank laws can reuse these
screens, but the present results do not justify a certificate campaign.
A different state basis that weakens those rank laws first needs its own
placement-dependent transfer. Singleton return margins alone cannot select
that design.
