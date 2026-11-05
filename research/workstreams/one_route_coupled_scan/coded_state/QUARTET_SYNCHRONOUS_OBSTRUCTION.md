# Exclude all coordinate-synchronous dimension obstructions

2026-10-31. Full-state bursts can dominate a scalar comparison even though
there are few such histories. This note bounds their entire geometric
obstruction class directly through the random regional routing.

Use `n=512`, `T=4096`, `M=8192`, `N=nT=256M`, and `W=943718`.
Each of the 256 regions routes its `M` slots by an independent uniform
permutation of the outer rows. Quartet partitions and field maps may be
fixed arbitrarily throughout this proof.

A coordinate-synchronous support has `S_t` equal to either the empty set
or all `n` coordinates at every update. Define its exact live sets,
footprint, state volume, deaths, and footprint rank as usual:

```
U_t = Gamma_t(S_(t-1)),          S_0=empty,
V_t = U_t union S_t,
F = sum_t |V_t|,                L = sum_t |S_t|,
D = sum_t |S_(t-1) \ S_t|,
R = sum_rows min(128, covered degree).
```

There is no death after `T`. For these supports `U_t` is empty or full
according to `S_(t-1)`, independently of the quartet partition.

**Theorem.** With probability less than `2^-119` over the regional
routing, some coordinate-synchronous support satisfies

```
L>R  and  R+D+1<=W.                                  (1)
```

The theorem excludes this whole sufficient dimension-obstruction class.
It does not exclude every bad codeword with a synchronous support, nor
the remaining nonsynchronous obstruction classes. It requires no union
over eligible row sets, ranks, footprints, field maps, or state amplitudes.

## Count the histories and partial regions

Let `r` be the number of closed full-state runs. Then

```
D = nr,             F = L+D.                         (2)
```

The footprint is obtained by extending each full-state run one update
to the right, except at the free terminal. Extensions may merge runs.
Thus the footprint has at most `r+1` interval components.

Call a region partial when its footprint volume is strictly between
zero and `M`. A partial region contains an endpoint of a footprint
component. Therefore its number `m` satisfies

```
m <= min(256,2r+2).                                  (3)
```

There are exactly

```
binom(T,2r) + binom(T,2r+1) = binom(T+1,2r+1)         (4)
```

binary histories with `r` deaths. The two terms distinguish an empty
or full terminal state. Changes alternate from the initial empty state,
so choosing their update positions proves (4).

When `r=0`, a nonempty history is a single full suffix. Its footprint
contains some whole suffix regions and at most one partial region.
Consequently `R=min(F,128M)`. Since `L=F` and `W<128M`, (1) is impossible.
The empty history also fails (1).

## A route bound for one fixed history

Fix a synchronous history. Let `v_j` be its footprint volume in region
`j`, and put `p_j=v_j/M`. In that region, the covered rows form a uniform
`v_j`-subset of the `M` rows. These subsets are independent across regions.

For each row `i`, let `Z_i` be its covered degree. Every `Z_i` has the
same marginal distribution as

```
Z = sum_(j=1)^256 Bernoulli(p_j),       mu=E Z=F/M.    (5)
```

The summands in (5) are independent. The different row degrees need not
be independent.

We use three standard negative-association facts: indicators of a uniform
subset are negatively associated; independent unions preserve that
property; and increasing functions of disjoint variable blocks preserve
the corresponding product inequality. Here each row is one disjoint
block of the regional indicators. For nonnegative increasing row
functions, repeated application gives a product bound on their joint
expectation. These facts are from
[Joag-Dev and Proschan, *Negative Association of Random Variables, with
Applications*, Annals of Statistics 11(1), 286–295 (1983)](https://doi.org/10.1214/aos/1176346079).

Define the total excess above the RS cap by

```
A = sum_i (Z_i-128)_+ = F-R.
```

For any `a>1`, negative association therefore gives

```
E a^A <= [E a^((Z-128)_+)]^M.                        (6)
```

Set `a=3/2`. Only the `m` partial regions contribute randomness to `Z`.
For a centered Bernoulli variable, the logarithm of its moment generating
function has second derivative at most `1/4`. Integrating twice from zero
gives the bound `exp(lambda^2/8)`. Independence in (5) yields

```
E a^(Z-mu) <= exp(m*(ln a)^2/8) <= 2^(m/32).          (7)
```

For the last inequality, `ln(3/2)<13/32` and `ln 2>=2/3` imply
`(ln a)^2/(8 ln 2)<507/16384<1/32`. The first logarithm bound follows
from the seven-term alternating series for `ln(1+1/2)`. The second
follows by integrating the convex function `1/x` on `[1,2]`.

The pointwise inequality
`a^((Z-128)_+) <= 1+a^(Z-128)` now gives

```
E a^((Z-128)_+) <= 1+2^(m/32)*a^(mu-128).             (8)
```

If (1) holds, (2) and integer ranks imply

```
A >= nr+1+(F-W)_+.
```

Apply Markov's inequality to (6), and write `omega=W/M`. The resulting
bound is maximized over `F` at `F=W`: below that value the bracket
increases with `F`; above it the external factor cancels the growing
term and leaves a decreasing term. Explicitly, with
`u=(mu-omega)_+`, the bracket after cancellation is

```
a^(-u) + 2^(m/32)*a^(mu-128-u)
  <= 1+2^(m/32)*a^(omega-128).
```

Thus every fixed history with `r>=1` satisfies the uniform bound

```
Pr[(1)] <= a^(-nr-1)
             * [1+2^(m/32)*a^(omega-128)]^M.          (9)
```

No footprint or rank union is hidden in (9). The maximization is an
inequality for each fixed history. We will use the simple bound
`omega<116`.

## Sum all histories with one fixed rational marker

First take `1<=r<=57`. By (3), the nonconstant term in the bracket of
(9) is bounded by

```
nu_r = 2^((r+1)/16)*(2/3)^12 <= r*nu_1,
nu_1 = 2^(1/8)*(2/3)^12 < 1/112.                     (10)
```

To justify the first inequality, `ln r-(r-1)ln(2)/16` is concave on
`[1,57]` and is nonnegative at both endpoints. At the second endpoint
use `57^2>2^7`. For the second inequality in (10), use

```
2^(1/8)<9/8,
(9/8)*(2/3)^12 = 512/59049 < 1/112.
```

Bernoulli's inequality and `ln 2>=2/3` now give

```
(1+nu_r)^M <= (1+nu_1)^(Mr) < 2^(128r).              (11)
```

Indeed, `M log2(1+nu_1)<8192*3/224=768/7<128`.
Also `(3/2)^12>2^7`, so `log2 a>7/12`. The number in (4) is less than
`8192^(2r+1)=2^(26r+13)`. Summing (9) using these inequalities gives

```
Pr[some history with 1<=r<=57 satisfies (1)]
  < sum_(r=1)^57 2^(13-(434/3)r)
  < 2^(-130).                                       (12)
```

For `r>=58`, use `m<=256` and

```
2^(m/32)*a^(omega-128) < 256*(2/3)^12 < 2.
```

The bracket in (9) is less than three. There are at most `2^4096`
binary histories altogether. Since `3^5<2^8` and `log2 a>7/12`,

```
Pr[some history with r>=58 satisfies (1)]
  < 2^4096 * 3^8192 * (3/2)^(-(512*58+1))
  < 2^(4096+8192*8/5-(512*58+1)*7/12)
  = 2^(-2401/20)
  < 2^(-120).                                       (13)
```

The impossible `r=0` case, (12), and (13) prove the theorem:
`2^-130+2^-120<2^-119`.

## A neighborhood with sixteen perturbations per closed run

The same argument covers partially occupied states near these histories.
Measure distance by the number of coordinate-time positions at which two
support masks differ. For a synchronous base history with `r` closed
runs, allow distance at most `16r`. When `r=0`, the allowed distance is
zero.

**Neighborhood theorem.** With probability less than `2^-73` over the
same routing, some support within this neighborhood satisfies (1).
The probability bound is uniform over all quartet partitions. There is
no union factor for the locations or values of the perturbations.

To prove this, first change one coordinate-time support bit. The state
volume changes by at most one. Its two incident temporal edges show
that the death count changes by at most one. The footprint can change
only at its current slot and at the four slots of its next quartet.
Consequently the footprint rank changes by at most five. At the free
terminal there is no next quartet, which only improves the bound.

Let `S` be an obstructing mask at distance `k` from a synchronous base
`B`. Use subscripts to distinguish their statistics. Repeating these
one-bit comparisons gives

```
|L_S-L_B|<=k,       |D_S-D_B|<=k,       |R_S-R_B|<=5k.
```

Since ranks and volumes are integers, (1) implies

```
R_B <= L_B-1+6k,
R_B <= W-D_B-1+6k.
```

For the base, `D_B=512r` and `F_B=L_B+D_B`. Hence its rank excess must
satisfy

```
F_B-R_B >= 512r+1+(F_B-W)_+-6k
         >= 416r+1+(F_B-W)_+,                        (14)
```

where the last inequality uses `k<=16r`. The event that any perturbation
works is contained in this single event for its base. This containment
explains why perturbation locations require no union bound.

The proof of (9) now gives the same bracket with prefactor
`a^(-416r-1)`. For `1<=r<=75`, inequality (10) still holds:
`ln r-(r-1)ln(2)/16` is nonnegative at both endpoints of `[1,75]`.
At the second endpoint, the exact comparison `75^8>2^37` suffices.
Equations (4) and (11) therefore give

```
Pr[some allowed perturbation of a base with 1<=r<=75
   satisfies (1)]
  < sum_(r=1)^75 2^(13-(266/3)r)
  < 2^(-74).                                        (15)
```

For `r>=76`, use the bracket bound by three and all `2^4096` bases:

```
Pr[some allowed perturbation of a base with r>=76
   satisfies (1)]
  < 2^4096 * 3^8192 * (3/2)^(-(416*76+1))
  < 2^(4096+8192*8/5-(416*76+1)*7/12)
  = 2^(-24801/20)
  < 2^(-1200).                                      (16)
```

The `r=0` bases are unchanged and remain impossible. Adding (15) and
(16) proves `2^-74+2^-1200<2^-73`.

The distance is a total count of changed coordinate-time positions.
It does not mean sixteen arbitrary coordinate changes at every update.
The theorem covers a union of neighborhoods of all synchronous bases;
it imposes no restriction on which perturbation coordinates are chosen.

## What this resolves

There is no deterministic exclusion valid for every routing. For example,
fill the first eight updates of every region and leave its last eight
empty. This support is saturated and has

```
L=1048576,     D=f=131072,     F=1179648.
```

If the same 4608 rows occupy the nine covered updates in every region,
then `R=589824` and `R+D+1=720897<W`. This is a genuine dimension mask
for those routes, independently of the quartet partitions. The probability
of that particular alignment, allowing any common row set, is exactly
`binom(8192,4608)^(-255)`. Its existence does not imply typical failure.

The theorem bounds every synchronous obstruction, including less extreme
alignments. Full-state phase histories therefore cannot be treated as a
typical obstruction merely because they dominate a relaxed scalar bound.
The next step is to retain this small-history-count advantage for broader
families containing partial states, rather than to retune the same scalar
comparison.

## Verification status

The derivation uses exact inequalities and one fixed marker `a=3/2`.
Negative association, the partial-region accounting, and the complete
synchronous proof have independent derivation review. The standalone
[constant and structural checker](quartet_synchronous_obstruction_checks.py)
also covers the neighborhood constants and one-bit comparisons. The
parent ran it serially in 0.152 seconds. The exact replay passed:

* 26 constant comparisons and 75 small-range inequalities;
* 4,432 synchronous histories and 12,960 partial-region comparisons;
* 19 history-count identities and 81 suffix-rank checks;
* 57,384 one-bit comparisons and 4,218 obstruction-threshold comparisons;
* six negative-association controls over 1,929 dependent route outcomes,
  together with six centered-moment controls.

The written synchronous proof, radius-`16r` extension, and checker source
all passed independent final audit. The review confirmed the threshold
constants, physical one-bit comparisons, dependent route controls, and
the absence of a perturbation union. No numerical job was run by the
note's author. These results certify the stated geometric classes; the
full-code certificate remains open.
