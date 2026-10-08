# Higher-rate Field-SPIN screens

The utility target is rate at least 3/4 with 10% relative symbol distance and competitive encoding cost.
The alternative target is substantially more than 10% distance at rate 1/2.
The half-rate, 10% rectangular-outer work stopped when this target changed.
Its unfinished files remain separate.

The current higher-rate screens do not establish a useful profile.
Every requested rectangular Toeplitz profile has a sampled middle class with a negative margin.
The inexpensive fixed GRS[16,12] outer also has weak sparse margins under the stated inner.
These are failures of upper-bound screens, not distance counterexamples.
No native timings or new directed certificates are reported here.

## Rectangular Toeplitz source law

Let F=GF(16), fix an outer dimension k and parity length p, and write b=k+p.
Each source row u in F^k is encoded as (u,Du).
For each row, independently sample the k+p-1 diagonals of a p-by-k Toeplitz matrix D:

    D_(i,j) = theta_(i-j),     -(k-1) <= i-j <= p-1.

Fix any nonzero u before sampling D.
Let j_* be the largest index with u_(j_*) nonzero.
Condition on every diagonal except theta_(r-j_*) for 0<=r<p.
The coefficient of theta_(r-j_*) in output coordinate i is u_(i-r+j_*), when that index exists.
It vanishes when r<i, because i-r+j_*>j_*.
The diagonal coefficient is u_(j_*), which is nonzero.
Thus the remaining p-by-p map is upper triangular and invertible.
All selected diagonal indices lie in the allowed range, including when j_*>=p.
Consequently Du is uniform in F^p for every nonzero u.

For any nonnegative function f on F^k times F^p, linearity of expectation gives

    E_D sum_(u != 0) f(u,Du)
      = 16^(-p) sum_(u != 0, v in F^p) f(u,v).

Independent matrices make this identity multiply across active source rows.
One may fix the route and inner randomness before choosing f to mark the output weight.
The resulting value-counting measure makes nonzero input values independent within each fixed support.
This construction therefore uses no input diagonal.

Each outer coordinate has one independently permuted region of L symbols.
The width-16 monomial inner carries its state across all b regions.
It uses independent uniform GL(s,16) refresh matrices, zero initial state, and discarded final state.
These are the same inner maps as in the preserved Toeplitz64 construction.

For nominal source size 2^ell, the screen uses

    K = 16k ceil(2^ell/(16k)),     L=K/k,     N=bL.

The largest bad weight is ceil(N/10)-1.
Thus the rate is exactly k/b at the rounded size; the screen does not silently identify K with 2^ell.

### Sparse source sum without zero systematic rows

Let R_h(z) be the two-state regional transfer for h uniformly placed nonzero inputs, as defined in the existing regional argument.
For a active source rows, the unrestricted region matrix is

    M_a(z) = sum_(h=0)^a binom(a,h) 15^h R_h(z).

The new sparse calculation retains the requirement that every active source row have a nonzero systematic coordinate.
It tracks t, the number of source rows already seen in the first k regions.
From t seen rows, choose i previously seen rows and j unseen rows at the next coordinate.
The transition to t+j has matrix weight

    binom(a-t,j) sum_(i=0)^t binom(t,i) 15^(i+j) R_(i+j)(z).

After the k systematic regions, retain only t=a, then apply M_a(z)^p for parity.
Multiply by 16^(-pa), binom(L,a)/15, and z^(-ceil(N/10)+1).
The division by 15 removes scalar multiples once per complete source word.
This positive calculation sums exact systematic support restrictions under the averaged value measure.
Its local matrices still use the existing inner envelope.

The middle-class bound uses the unrestricted positive source envelope:

    H(x,z) = sum_(w=0)^16 binom(16,w) (15x)^w (1+x)^(16-w) T_w(z),

    B_a <= binom(L,a)^(-(b-1)) 16^(-pa) / 15
           * z^(-ceil(N/10)+1) x^(-ba)
           * [1,0] H(x,z)^(bL/16) [1,1]^T.

All searches below evaluate these bounds in floating arithmetic.

## Rectangular screening results

The following values are -log2 of the numerically optimized class bound at nominal K20.
The sparse column retains nonempty systematic support for every active row.
The middle witness column identifies one sampled class where the bound fails.
No table row claims complete finite coverage.

| Outer | State s | Sparse a=1 | Sparse a=16 | Middle a | Middle margin |
|---|---:|---:|---:|---:|---:|
| [32,24] | 3 | 2.33 | -345.11 | 8739 | -249518.89 |
| [32,24] | 4 | 6.15 | -259.87 | 8739 | -209084.03 |
| [32,24] | 5 | 9.34 | -196.38 | 8739 | -165371.28 |
| [32,24] | 6 | 12.23 | -142.63 | 8739 | -123554.37 |
| [64,48] | 3 | 18.21 | -470.18 | 4371 | -234043.36 |
| [64,48] | 4 | 22.82 | -297.54 | 4371 | -193592.47 |
| [64,48] | 5 | 26.30 | -168.86 | 4371 | -149863.26 |
| [64,48] | 6 | 29.35 | -59.38 | 4371 | -108033.04 |
| [128,96] | 3 | 50.85 | -734.27 | 2185 | -226519.74 |
| [128,96] | 4 | 57.52 | -386.74 | 2185 | -186055.53 |
| [128,96] | 5 | 61.40 | -128.55 | 2185 | -142321.04 |
| [128,96] | 6 | 64.65 | 92.34 | 2185 | -100511.07 |

The rounded (K,N,L) triples are (1048704,1398272,43696), (1049088,1398784,21856), and (1049088,1398784,10928).
The middle scan also evaluated 16, 64, 256, 1024, and approximately 10%, 20%, 33%, 50%, 80%, and 100% of L.
The displayed negative classes suffice to stop directed closure with this envelope.

The requested rate-2/3 comparison also fails the sampled middle screen through s=6.
For [48,32], s=6 gives 37.67 bits at a=1 but -22377.77 bits at a=6553.
For [96,64], s=6 gives 79.62 bits at a=1 but -10855.49 bits at a=3276.
Rate 2/3 would not itself meet the revised rate target.

## Fixed GRS[16,12] target

The native cost study identified a cheap parity circuit for systematic GRS[16,12].
The fixed code interpolates evaluations at field labels 0 through 11 to labels 12 through 15 in GF(16).
Its minimum distance is five.
For each prescribed support of size w, its number of codewords is

    A_w = sum_(j=0)^w (-1)^j binom(w,j) 16^max(0,w-j-4).

Unlike the random Toeplitz construction, this fixed outer requires independent uniform nonzero GF16 input multipliers for the local value law used here.
The new screen sums the exact MDS support spectrum through four active rows.
For middle classes it optimizes rho>0 in

    eta(rho) = max_(5<=w<=16) A_w / (16^(w-4) rho^w).

The source envelope is eta(rho)^a 16^(-4a), with coordinate factor 16rho.
The middle bound consequently substitutes 16rho x for 15x in H and multiplies by eta(rho)^a.

At nominal K20, K=1048704, N=1398272, and L=87392.
The preserved local envelope gives:

| State s | Sparse a=1 | Sparse a=2 | Sparse a=4 | Middle a=29130 |
|---:|---:|---:|---:|---:|
| 3 | -1.10 | -14.82 | -47.64 | -312412.96 |
| 4 | 4.85 | -6.77 | -32.27 | -260646.68 |
| 5 | 10.44 | -0.18 | -22.16 | -223630.86 |
| 6 | 15.83 | 5.89 | -14.46 | -172448.48 |

For one active row, each regional input has weight zero or one.
The preserved local transfer is exact for both weights.
Thus a more accurate return term at larger input weights cannot change the one-row transfer.
The displayed one-row value remains a Chernoff upper bound, not an exact failure probability.

### A diagnostic for local slack

The existing nonzero-state row sums to g+J_upper-J_lower.
Here g is the exact output moment, whereas the return event is bounded separately in the two state transitions.
The excess can be measured exactly when the input X is uniform in F^16.

Let Z(z)=E_X z^wt(X), and let Z_0(z)=E_X[1_(CX=0) z^wt(X)].
Then

    Z(z) = ((1+15z)/16)^16,
    Z_0(z) = 16^-16 sum_w binom(16,w) birth_zero_count(w,s,16) z^w.

Because CA=0, Y=X+AR is uniform and Q'=R+CY.
For a nonzero previous state, R is uniform over nonzero state vectors.
The exact aggregated transfer is therefore

    [[Z_0, Z-Z_0], [(Z-Z_0)/(16^s-1), Z-(Z-Z_0)/(16^s-1)]].

Both row sums equal Z.
For s=3, the preserved envelope's nonzero-state row sum is 2.81485 times Z at z=0.1 and 1.36551 times Z at z=0.2.
This calculation diagnoses a possible improvement; it does not supply a general input-fugacity bound.

The screen accepts an aggregate callback F(v,z,s), where

    F(v,z,s) = sum_w binom(16,w) 15^w v^w T_w(z).

For the GRS source envelope, H=(1+x)^16 F(16rho x/[15(1+x)],z,s).
For the Toeplitz source envelope, H=(1+x)^16 F(x/(1+x),z,s).
This normalization is tested against the preserved local calculation.

### Executed exact-return middle screens

The separate `exact_inner_character.py` engine computes the joint return term from an exact integer character histogram.
Its floating aggregate was used in new s=3 searches, with the same source envelopes and size rounding.
The improvement does not close either candidate:

| Outer | Active rows | Preserved local envelope | Exact local aggregate, floating |
|---|---:|---:|---:|
| Toeplitz[128,96] | 2185 | -226519.74 | -58992.19 |
| Toeplitz[128,96] | 3642 | -233555.24 | -3430.02 |
| Toeplitz[128,96] | 5464 | -158344.32 | 114798.11 |
| GRS[16,12] | 17478 | -284318.27 | -118589.82 |
| GRS[16,12] | 29130 | -312412.96 | -83629.15 |
| GRS[16,12] | 43696 | -249190.89 | 26106.71 |

For Toeplitz at a=2185, the optimized witness is approximately x=0.135651613267 and z=0.285520658494.
For GRS at a=17478, it is approximately x=0.136311693947, z=0.290217277884, and rho=0.915754398398.
Exact Fraction evaluation of the local matrices at these decimal witnesses gives positive entries.
The floating entries agree within 2.6e-15 relative error.
A dense Toeplitz witness also agrees within 2.1e-15.
These checks address cancellation in the local character sum; they do not certify the global floating bounds.

The exact-aggregate middle deficit also persists at the smaller rounded sizes.
Toeplitz[128,96] gives -4271.10 bits at K=66048, a=137, and -15285.37 bits at K=262656, a=547.
GRS[16,12] gives -7515.64 bits at K=65664, a=1094, and -29738.86 bits at K=262272, a=4371.
All four checks use s=3 and approximately 20% active source rows.

## Reproduction and next step

From the repository root:

```sh
python -B research/workstreams/char2_field_spin/theory/field_rate_tradeoff.py --length 128 --dimension 96 --state 3 4 5 6 --active 1 2 4 16
python -B research/workstreams/char2_field_spin/theory/field_rate_tradeoff.py --length 128 --dimension 96 --state 3 4 5 6 --active 2185 --middle
python -B research/workstreams/char2_field_spin/theory/field_rate_tradeoff.py --mds --state 3 4 5 6 --active 1 2 4
python -B research/workstreams/char2_field_spin/theory/field_rate_tradeoff.py --mds --state 3 4 5 6 --active 29130 --middle
python -B research/workstreams/char2_field_spin/theory/field_rate_tradeoff.py --length 128 --dimension 96 --state 3 --active 2185 3642 --middle --exact-local
python -B research/workstreams/char2_field_spin/theory/field_rate_tradeoff.py --mds --state 3 --active 17478 29130 --middle --exact-local
python -B -m unittest discover -s research/workstreams/char2_field_spin/theory -p field_rate_tradeoff_tests.py -v
```

Ten tests pass.
They exhaust the rectangular law over small GF4 shapes, verify the complete value identity, and compare the positive source DP with exhaustive support enumeration.
They also check size rounding, preserved formula agreement, MDS counts, input-uniform row sums, and aggregate normalization.

The exact return term removes substantial middle-class slack, but the sampled deficits remain too large for directed closure.
The next useful step is a stronger source-composition argument or a changed construction, while separately addressing the already-exact one-row transfer.
Neither a larger outer nor a larger inner state should advance to directed closure without a competitive encoding-cost estimate.
