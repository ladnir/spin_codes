# Block Expand constituents for Structured SPIN

## Question and conclusion

The relevant use of Expand--Convolute is local. Replace the BA constituent by
one sampled sparse \([512,256]\) code, and reuse that code in every Structured
SPIN row. The global Expand--Convolute certificate remains useful evidence,
but it is not the proposed construction.

A direct local EC constituent has a promising first moment. More importantly,
its second-moment problem starts from the full message space. The input pair
enumerator is therefore explicit. This removes the unknown genus-two BCH
enumerator that blocked the BA proof.

An iterated Block Expand construction may give an even cleaner proof. Its
layers can be sparse invertible shears. Their one-word and two-word kernels
factor across small regions and require no convolution state.

## Direct block EC ensemble

Set \(K=256\) and \(B=512\). Partition the \(B\) intermediate coordinates
into fourteen consecutive regions. Eight regions have length 37, and six
regions have length 36.

For every input coordinate and every region, sample one neighbor uniformly
and independently. The resulting binary matrix

\[
 E\in\mathbb F_2^{K\times B}
\]

has left degree 14. Different edges can enter the same right coordinate and
cancel over \(\mathbb F_2\).

Independently sample a wrapped binary convolution \(C\) of memory 15. Define
the constituent generator

\[
 G_{\mathrm{BEC}}:=EC.
\]

Setup samples \((E,C)\) once. Structured SPIN reuses the realized
\(G_{\mathrm{BEC}}\) in all \(L=4096\) rows. This ensemble has the exact
power-of-two dimensions \([512,256]\); it needs no shortening or puncturing.

The fourteen unequal regions are deliberate. They preserve the exact output
length while retaining one independent edge per input coordinate in every
region. Two equal input rows collide in all regions with probability

\[
 37^{-8}36^{-6}.
\]

## First-moment diagnostic

For \(0<z<1\), define

\[
 M(z):=mathbb E_{E,C}\!\left[
   \sum_{x\in\mathbb F_2^{256}\setminus\{0\}}
   z^{\operatorname{wt}(xEC)}
 \right].
\]

The exact regional parity law and the exact reduced convolution transfer give
the following binary64 values.

| \(z\) | \(\log_2 M(z)\) | positive-output part | random \([512,256]\) |
|---:|---:|---:|---:|
| 0.03 | -57.6955 | -65.7742 | -234.1661 |
| 0.05 | -57.6914 | -64.9501 | -219.9607 |
| 0.10 | -57.6789 | -63.7288 | -185.5982 |
| 0.20 | -57.6411 | -62.2664 | -121.3264 |
| 0.50 | 43.5023 | 43.5023 | 43.5008 |

The expected number of nonzero kernel messages is below
\(2^{-57.7008}\) in binary64. The low-marker total is dominated by the
kernel event, especially by collisions of weight-two messages. After removing
zero outputs, the \(z=0.1\) moment is below \(2^{-63.7288}\).

For comparison, the two-sided \([510,255]\), degree-\(10/5\), memory-15
parent has only 44.9038 bits against the analogous row-collision event. Its
positive-output \(z=0.1\) moment is about \(2^{-51.1146}\). The degree-14
power-of-two constituent is therefore the better starting point.

These numbers are diagnostics. They do not prove a spectrum event or a SPIN
distance theorem. They do show two useful facts. The dense first moment is
almost random, and the sparse defect has more than the requested 40-bit
budget before the Structured SPIN transfer.

## Finite proof target

For a realized constituent \(G\), let

\[
 A_w(G):=\left|\left\{
 x\in\mathbb F_2^{256}\setminus\{0\}:
 \operatorname{wt}(xG)=w
 \right\}\right|.
\]

The existing all-occupation Structured SPIN transfer accepts deterministic
upper caps on every \(A_w(G)\). The direct block-EC proof should construct a
cap event \(\mathcal G\) such that

\[
 \Pr_{E,C}[\neg\mathcal G]\le 2^{-40-\eta}
\]

for enough reserve \(\eta>0\), and then reuse the existing conditional
transfer. The cap event is paid once because the same constituent is reused.

First moments alone do not provide useful central caps at 40-bit confidence.
The next object is the second factorial moment

\[
 \mathbb E[A_w(G)(A_v(G)-\mathbf 1_{w=v})]. \tag{1}
\]

Unlike the EBCH--BA case, the input to (1) is explicit. For two distinct
nonzero messages \(x,y\), their coordinate pairs have a four-symbol type

\[
 (n_{00},n_{01},n_{10},n_{11}),
 \qquad \sum_{a,b}n_{ab}=256.
\]

The number of ordered pairs of each type is its multinomial coefficient,
after excluding the three rank-one faces \(x=0\), \(y=0\), and \(x=y\).
No BCH pair enumerator is required.

Each expander region has an exact four-symbol parity-occupancy kernel. The
remaining technical task for direct EC is a pair transfer for the shared
wrapped convolution. A successful transfer yields (1), after which
Cantelli's inequality can reproduce the shell-cap interface used by the
random-code comparator.

## Block Expand--\(t\): expander plus accumulators

The most direct iterated construction keeps the expander as the rate-half
map and uses accumulators as rate-preserving mixers. Let \(E\) be the
degree-14 regional map above. For independent uniform permutations
\(P_1,\ldots,P_t\), define

\[
 G_{\mathrm{BE}\text{-}t}:=EP_1A\cdots P_tA, \tag{2}
\]

where \(A\) is the zero-state prefix accumulator. Setup samples \(E\) and
the \(t\) permutations once. Structured SPIN reuses the resulting
\([512,256]\) constituent in all rows.

Every accumulator is invertible. Hence the kernel of (2) equals the kernel
of \(E\), whose current first-moment diagnostic has 57.7008 bits of margin.
The exact one-word expander spectrum composes with the standard accumulator
input--output enumerator.

| Accumulators \(t\) | \(\log_2\sum_{w>0}\overline A_w0.1^w\) | Excess over random at \(z=0.5\) |
|---:|---:|---:|
| 0 | -38.5068 | 0.2075 bits |
| 1 | -53.0832 | 0.4842 bits |
| 2 | -64.8489 | 0.1607 bits |
| 3 | -72.4758 | 0.0189 bits |
| 4 | -78.4515 | 0.0004 bits |
| 5 | -85.8289 | 0.000003 bits |

Thus \(t=2\) slightly improves the direct memory-15 EC value
\(-63.7288\) at \(z=0.1\). The \(t=3\) candidate has more sparse reserve and
an almost random dense moment.

This construction retains the accumulator proof machinery developed for BA.
The accumulator pair-type kernel is known exactly. The missing initial pair
table now comes from \(E\), not BCH. For every four-symbol input type, the
regional edge choices give an explicit nonnegative generating function for
the output pair type. Composing that table with two or three accumulator
kernels gives the second factorial moments of all output shells.

This is the central advantage over BA. The EBCH pair table was an unknown
genus-two enumerator. The expander pair table is defined by the sampled
ensemble and can be computed from its regional law.

### Exact regional pair kernel

Fix an ordered pair of distinct nonzero input messages. Write

\[
 p=n_{10},\qquad q=n_{01},\qquad r=n_{11}.
\]

The zero coordinates (n_{00}=256-p-q-r) place no edges and do not affect
one region. In a region of length \(\ell\), let

\[
 m=(m_{00},m_{01},m_{10},m_{11}),\qquad \sum_{a,b}m_{ab}=\ell,
\]

be the output pair type. Introduce the local exponential generating
functions

\[
 F_{ab}(X,Y,Z)
 =\frac14\sum_{s,t\in\{-1,1\}}
   s^a t^b\exp(sX+tY+stZ). \tag{3}
\]

The number of labeled edge assignments that produce type (m) is exactly

\[
 p!q!r!\binom{\ell}{m_{00},m_{01},m_{10},m_{11}}
 [X^pY^qZ^r]\prod_{a,b\in\mathbb F_2}F_{ab}(X,Y,Z)^{m_{ab}}. \tag{4}
\]

Dividing (4) by \(\ell^{p+q+r}\) gives the regional transition
probability. Formula (3) is the parity projector: a type-10 edge contributes
the sign (s), a type-01 edge contributes (t), and a type-11 edge
contributes (st). Formula (4) then assigns the labeled edges to the
unlabeled occupancy of each bin and finally chooses which bins have the four
output symbols.

The fourteen regions are conditionally independent. Their regional kernels,
with eight copies at \(\ell=37\) and six at \(\ell=36\), convolve to the
exact pair-type law of (E). This is a proved finite identity, not a mixing
hypothesis. `verify_block_expand_pair_region_small.py` checks (4) against
direct enumeration in 176 small cases and 3,036 output-type entries using
exact integers and rationals.

### Exact second-moment pipeline

For an input pair type
\(n=(n_{00},n_{01},n_{10},n_{11})\), the number of ordered message pairs is

\[
 \binom{256}{n_{00},n_{01},n_{10},n_{11}}. \tag{5}
\]

The sum includes precisely the types satisfying

\[
 n_{10}+n_{11}>0,\quad
 n_{01}+n_{11}>0,\quad
 n_{10}+n_{01}>0, \tag{6}
\]

which say (x\ne0), (y\ne0), and (x\ne y). Compose (5), the fourteen
regional kernels (4), and (t) copies of the exact common-permutation
accumulator pair kernel. Summing final pair types with first and second
weights (w) and (v) gives

\[
 \mathbb E[A_wA_v-\mathbf1_{w=v}A_w] \tag{7}
\]

exactly. In particular, only the diagonal (w=v) is needed for the
per-shell Cantelli caps. Equations (3)--(7) completely specify the finite
second-moment calculation; no constituent enumeration and no unknown BCH
pair spectrum occurs.

The open problem is computational rather than definitional. A length-512
four-symbol type vector has

\[
 \binom{512+3}{3}=22{,}632{,}705
\]

entries, and a materialized pair transition is much larger. The certificate
implementation must therefore evaluate (7) by nonnegative coefficient
bounds, sparse/on-demand convolution, or an equivalent shell-specific
compression. Binary64 evaluation alone will remain diagnostic until an
outward implementation covers every shell.

### Finite cap gate

The one-shot random-constituent proof sets the shell caps to zero outside
weights 42 through 470. Keeping that support gives a quick necessary test:
the expected positive spectrum outside the interval, plus the expected
number of nonzero kernel messages, must already fit the 40-bit setup budget.
The binary64 results are

| Accumulators (t) | positive mass outside 42--470 | kernel term |
|---:|---:|---:|
| 2 | (2^{-32.4286}) | (2^{-57.7008}) |
| 3 | (2^{-51.6918}) | (2^{-57.7008}) |
| 4 | (2^{-52.4060}) | (2^{-57.7008}) |
| 5 | (2^{-52.6725}) | (2^{-57.7008}) |

Thus Block Expand--2 fails this gate regardless of its variance. It could
only be revived by admitting lower-weight shells and redoing the conditional
SPIN transfer. Block Expand--3 is the first candidate that fits the existing
support.

The positive random-code caps also cannot be copied literally. In dense
shells they sit only about (2^{25.5}sqrt{mu_w}) above the random mean, so
even a small multiplicative change in the mean eventually crosses them. For
a hypothesized variance bound

\[
 \operatorname{Var}(A_w)\le F\mathbb E[A_w], \tag{8}
\]

the diagnostic instead recenters each positive cap on the Block Expand mean
and gives it per-shell Cantelli failure (2^{-51}). The resulting event
margin is about 42.25 bits for every tested (F\in\{1,2,4,16,512\}); the
tail and kernel terms are included. The price is a change in the three band
majorants used by the conditional SPIN proof.

At the uniform per-shell allocation (2^{-51}), for (F=1), the
low/central/high band-majorant inflations are respectively
0.3122/0.3134/0.3086 bits at (t=3), 0.0711/0.0699/0.0695 bits at (t=4),
and 0.0041/0.0039/0.0039 bits at (t=5). For (F=2), they are
0.3653/0.3378/0.3618 bits at (t=3), 0.1285/0.0964/0.1270 bits at (t=4),
and 0.0627/0.0309/0.0625 bits at (t=5). These figures do not prove that the
conditional all-occupation transfer closes; that transfer must be rerun with
the recentered caps.

The uniform allocation is unnecessarily conservative. The certified random
proof spends only about (2^{-42.26}) on its outer event, while the combined
target permits almost (2^{-40}). For Block Expand--5 and (F=2), allocate
only half of the available combined 40-bit budget to the 429 positive shells.
The resulting per-shell allowance is (2^{-49.7456}). With this allocation,
the low, central, and high band majorants are respectively 0.0082, 0.0018,
and 0.0084 bits *smaller* than the frozen random-code majorants.

This domination closes the complete conditional transfer. The exact-shell
outward rerun gives 67.9441 bits for (Q=1) and 105.3578 bits for (Q=2).
The existing outward bounds for (Q=3,ldots,159) and
(Q=160,ldots,4096) apply unchanged because they use only the three larger
frozen band majorants. The resulting conditional distance-failure bound is

\[
 2^{-51.6439589890}. \tag{9}
\]

The outward mean calculation is now complete. It uses exact integer regional
occupancy polynomials, five exact binomial accumulator transitions, and
256-bit Arb enclosures. For the frozen integer caps, (8) with (F=2) implies
an outer-event bound of (2^{-42.6282605259}). Combining it with (9) gives

\[
 \Pr[d_{min}<228{,}590]<2^{-42.6254759461}. \tag{10}
\]

Statement (10) is conditional only on the length-512 variance inequality
(8). The shell means, kernel term, zero-cap tail, cap arithmetic, and complete
SPIN transfer are outwardly certified.

Occupation one is not the bottleneck. At (t=3), the exact-shell binary64
RandomStepConv-M22 transfer has 75.17 bits of margin under (F=1), 74.47
bits under (F=2), 72.93 bits under (F=16), and 70.39 bits even under
(F=512). In every case the dominant admitted shell is weight 42.

Exact rational small models give evidence about (8). With four degree-four
regions, the maximum positive-shell variance-to-mean ratio after five stages
is 1.3505 at ((K,B)=(4,8)), 1.2397 at ((6,12)), and 1.1659 at
((8,16)). The latter two models decrease across the last accumulator
stages. These calculations use the exact regional kernel and exact
accumulator pair kernel. They support the target (F=2), but they are not a
length-512 bound and supply no monotonicity theorem in (B) or (d).

## Pure invertible expansion layers

The wrapped-convolution pair transfer may still be unnecessarily difficult.
An iterated construction can retain exact dimensions and remove the recursive
state.

Start from the systematic rate-half embedding

\[
 S(u):=(u,0)\in\mathbb F_2^{512}.
\]

One Block Expand layer samples a permutation of the 512 coordinates and
splits the permuted word as \((a,b)\in\mathbb F_2^{256}\times
\mathbb F_2^{256}\). It independently samples a sparse regional map

\[
 H:\mathbb F_2^{256}\longrightarrow\mathbb F_2^{256}
\]

and applies the shear

\[
 X_H(a,b):=(a,b+Ha). \tag{11}
\]

The map in (3) is invertible for every \(H\). Its transpose is another sparse
shear. For independent layers \((P_i,H_i)\), define

\[
 G_{\mathrm{BX}\text{-}t}
 :=S P_1X_{H_1}P_2X_{H_2}\cdots P_tX_{H_t}. \tag{12}
\]

Equation (12) defines a binary \([512,256]\) constituent for every realization.
It is sampled once and repeated in every outer row. Ordinary and transposed
encoding require \(O(tdB)\) bit operations when each \(H_i\) has left degree
\(d\). This pure-shear variant is secondary to (2).

The one-word weight kernel of a layer is exact and low dimensional. Condition
on input weight \(h\). The permutation gives a hypergeometric split
\((\operatorname{wt}(a),\operatorname{wt}(b))\). Regional parity occupancy
gives \(\operatorname{wt}(Ha)\). A second hypergeometric variable gives the
overlap of \(b\) and \(Ha\). These three variables determine the output
weight.

The two-word kernel has the same structure over four-symbol types. It factors
across regions and has no memory-15 state.

## Scalable theorem target

Let \(B_N\) be an admissible even block length and let
\(K_N=B_N/2\). Define (2) at length \(B_N\), with fixed \((d,t)\), and reuse
one sample in all \(L_N=N/B_N\) rows.

The asymptotic target is a certified mean-spectrum envelope
\(\overline A_{B_N}(w)\) and an event

\[
 A_w(G_{\mathrm{BE}\text{-}t})
 \le B_N^2\overline A_{B_N}(w)
 \quad\text{for every }w,
\]

whose failure probability is \(o(1)\). If the resulting envelope satisfies
the existing Structured SPIN variational inequalities, then
\(B_N=c\log N+O(1)\) gives linear distance and linear encoding work.

For the finite \(k=2^{20}\) theorem, the polynomial cap is insufficient. The
finite proof should instead use the exact second factorial moments to allocate
at most \(2^{-40-\eta}\) across all shells.

## Recommended order

1. Freeze Block Expand--5 and the target (F=2) in (8). Its complete
   conditional transfer is already closed.
2. Evaluate the exact second-moment pipeline (3)--(7) by a non-materialized
   length-512 method and certify (8).
3. Convert the current conditional receipt into an unconditional outward
   certificate. All numerical components already close with 42.6254 bits;
   only the variance premise must be discharged.
4. Benchmark ordinary and transposed encoding after the finite theorem
   closes. Compare (t=3,4,5) if the fifth accumulator is material.

Block Expand--3 is the cheapest candidate that passes the present zero-cap
gate. Block Expand--5 is now the proof target because it is the first tested
candidate whose factor-two, budget-optimized caps are dominated by the
already certified three-band interface. Block Expand--4 remains the likely
performance fallback if a sharper cap allocation or transfer proof can absorb
its small central inflation. Direct block EC remains a comparison point, but
its shared memory-15 pair state is less attractive.

## Diagnostic artifacts

- `evaluate_block_ec_constituent_mgf.py` and
  `block_ec_510_255_d10_m15_mgf_diagnostic.json` record the two-sided parent.
- `block_ec_518_259_d14_m15_mgf_diagnostic.json` records a second two-sided
  comparison.
- `evaluate_block_ec_balanced_left_regular_mgf.py` and
  `block_ec_512_256_d14_m15_mgf_diagnostic.json` record the exact-size
  constituent proposed above.
- `evaluate_block_expand_accumulate_spectrum.py` and
  `block_expand_accumulate_512_256_d14_t0_5_spectrum.json` record the
  Block Expand--\(t\) expected-spectrum sweep.
- `evaluate_block_expand_cap_budget.py` and
  `block_expand_cap_budget_diagnostic.json` record the zero-cap gate and the
  recentered-cap cost under five variance-factor hypotheses.
- `evaluate_block_expand_q1_cap_transfer.py` and
  `block_expand3_q1_cap_transfer_diagnostic.json` record the exact-shell
  occupation-one transfer for Block Expand--3.
- `verify_block_expand_pair_region_small.py` and
  `block_expand_pair_region_small_exact.json` verify the regional coefficient
  identity.
- `analyze_block_expand_pair_moments_small.py` and
  `block_expand_pair_moments_K4_B8_exact.json`,
  `block_expand_pair_moments_K6_B12_exact.json`, and
  `block_expand_pair_moments_K8_B16_exact.json` give exact small-model shell
  variances.
- `certify_block_expand5_F2_conditional_transfer.py` and
  `block_expand5_F2_conditional_transfer_outward.json` close every occupation
  for the frozen Block Expand caps.
- `certify_block_expand5_mean_cap_outward.py` and
  `block_expand5_mean_cap_outward.json` certify the complete mean spectrum and
  reduce the conditional theorem to the factor-two variance inequality.
- `build_block_expand5_conditional_manifest.py` and
  `FINITE_K20_BLOCK_EXPAND5_F2_CONDITIONAL_MANIFEST.json` hash-bind the
  conditional theorem artifacts and name the remaining hypothesis.
