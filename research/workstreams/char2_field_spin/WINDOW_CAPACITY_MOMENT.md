# A moment bound for a specified window-capacity event

The negative-association and equal-length steps in
`theory/window_dimension_mgf.py` have the correct inequality direction.
They give an upper bound on the probability that an aligned window family
has shortened outer dimension greater than the charged exit capacity.
This is a bound on a structural obstruction event, not on all low-weight
codewords.

The event must specify which exits are charged. Charging one exit for
every used region gives the formula audited below. A window ending at the
global stream end has no required exit. Windows touching across a regional
boundary can also share a single connected output interval. Excluding only
the first event, even with a final-suffix correction, does not exclude all
such merged-window obstructions.

## Routing experiment and shortened dimension

Fix an outer GRS[b,k] code, put p=b-k, and repeat it independently across
L source rows. There are b coordinate regions, each containing one
coordinate of every source row. In each region, sample an independent
uniform permutation of the L row labels.

Fix a set H of h used regions. In region j in H, fix one nonempty window
of length ell_j at a specified start. The starts and lengths in this
paragraph are fixed before sampling the route. Define

\[
 X_{ji}=\mathbf1\{\text{row }i\text{ is routed into window }j\},
 \qquad c_i=\sum_{j\in H}X_{ji},
 \qquad D=\sum_{i=1}^{L}(c_i-p)_+.
\]

For fixed windows, `(X_ji)_{i=1}^L` is the indicator vector of a uniform
ell_j-subset of the row labels. These subset vectors are independent
between regions. The rowwise shortened GRS space has dimension
`(c_i-p)_+`, so D is the dimension of the full outer space supported
inside these routed windows. Independent nonzero coordinate scalings do
not change D.

Let the causal inner have s state coordinates and zero initial state.
Impose zero state after every selected window. These homogeneous equations
have total rank at most sh. If D>sh, their common kernel contains a
nonzero outer word. Its inner output vanishes outside the selected windows,
because its input and state both vanish there. Thus total window length at
most e makes `D>sh` a sufficient bad-code event, independently of the
chosen linear inner coefficients. It is not a necessary condition for
all bad codewords, or even for all words supported in the same windows.

## Negative association of the regional subsets

For completeness, the needed fixed-size sampling fact has an elementary
proof. Let S be a uniform ell-subset of a finite set, and let I and J be
disjoint coordinate sets. Let f and g be increasing functions of the
indicators in I and J, respectively. Put `T=|S intersect I|`.

Conditioned on T=t, the chosen subsets in I and its complement are
independent and uniform, with respective sizes t and ell-t. Consequently,

\[
 \mathbb E[fg\mid T=t]=F(t)G(t),
 \quad F(t)=\mathbb E[f\mid T=t],
 \quad G(t)=\mathbb E[g\mid T=t].
\]

Coupling uniform subsets by a random ordering shows that F is increasing
and G is decreasing. For an independent copy T' of T,

\[
 2\operatorname{Cov}(F(T),G(T))
 =\mathbb E[(F(T)-F(T'))(G(T)-G(T'))]\le0.
\]

Hence `E[fg]<=E[f]E[g]`: fixed-size subset indicators are negatively
associated. This is a special case of the sampling-without-replacement
results of [Joag-Dev and Proschan (1983)](https://doi.org/10.1214/aos/1176346079).

Two closure steps apply. First, the union of independent negatively
associated vectors is negatively associated. To see this for two vectors,
condition on one vector, apply the inequality within the other, and then
apply it to the two resulting increasing conditional expectations.
Second, increasing functions of disjoint groups of coordinates preserve
negative association. Both properties are stated in Proposition 7 of
[Dubhashi and Ranjan, *Balls and Bins: A Study in Negative Dependence*](https://www.brics.dk/RS/96/25/BRICS-RS-96-25.pdf).

Apply these facts to the independent regional subset vectors, regrouped
by row. For z>=1, each row factor `z^((c_i-p)_+)` is nonnegative and
increasing in its indicator coordinates. Products of such factors are
also increasing, so induction gives

\[
 \mathbb E[z^D]
 \le\prod_{i=1}^{L}\mathbb E[z^{(c_i-p)_+}]
 =\left(\mathbb E[z^{(B_1+\cdots+B_h-p)_+}]\right)^L,
 \tag{1}
\]

where the B_j are independent Bernoulli variables with probabilities
`q_j=ell_j/L`. The equality uses identical row marginals, not independence
of the actual row hit counts. Replacing actual rows by independent rows
is an upper-moment comparison only.

## Why equalizing the lengths gives an upper bound

Define `f(c)=z^((c-p)_+)`. For integer c>=0, its second differences are
nonnegative: they are zero before the threshold, z-1 at c=p-1 when that
index exists, and `z^(c-p)(z-1)^2` for c>=p.

Fix the other h-2 Bernoulli variables and write their sum as T. For two
remaining probabilities a and d,

\[
 \mathbb E[f(T+B_a+B_d)\mid T]
 =f(T)+(a+d)\Delta f(T)+ad\Delta^2f(T).
\]

Holding a+d fixed and replacing a,d by their average increases ad.
It therefore cannot decrease the expectation. Repeated pair averaging,
followed by continuity, yields

\[
 \mathbb E[f(B_1+\cdots+B_h)]
 \le \mathbb E[f(B)],
 \qquad B\sim\operatorname{Binomial}(h,\bar q),
 \quad \bar q=\frac{\sum_j\ell_j}{hL}.
 \tag{2}
\]

This upper bound does not require equal actual window lengths. It also
does not approximate their sampling law by independent regional slots.
Only the one-row moment is enlarged after applying negative association.

Let a be the alignment, assume a divides L, and put `F=floor(e/a)`.
For total window length at most e,

\[
 \bar q\le q_h:=\min\left(1,\frac{aF}{hL}\right).
\]

The binomial expectation is increasing in its success probability.
Equations (1) and (2) thus give `E[z^D]<=M_h(z)^L`, where

\[
 M_h(z)
 =1+\sum_{c=p+1}^{h}\binom hc q_h^c(1-q_h)^{h-c}
             (z^{c-p}-1).
 \tag{3}
\]

For h<=p, D=0 deterministically; those families cannot satisfy either
positive dimension threshold used below.

## Counting and bounding the closed-window event

Choose the h used regions in C(b,h) ways. Write each positive window
length as a times a positive integer. The number of ordered positive
length vectors with sum at most F is C(F,h). This follows by summing
`C(t-1,h-1)` over t=h,...,F. Imposing the individual length cap L can
only reduce this number.

For a length a*m, the number of aligned starts is `L/a-m+1<=L/a`.
Consequently, the number of fixed window families is at most

\[
 N_h=\binom bh\binom Fh(L/a)^h.
\]

D is integer, so `D>sh` means `D>=sh+1`. Markov's inequality and a union
bound give, for any independently chosen real marks z_h>1,

\[
 \Pr[\exists\text{ such family with }D>sh]
 \le
 \sum_{h=p+1}^{\min(b,F)}
 N_h\min\{1,z_h^{-(sh+1)}M_h(z_h)^L\}.
 \tag{4}
\]

The final sum may also be capped at one. Each h may use its own mark.
Taking an infimum over marks improves the bound; using a finite rational
mark is already valid. Fixed starts and lengths were chosen before the
route only to establish each summand. The union bound covers a family
chosen after observing the route, without an independence assumption
between different families.

## A global-final suffix and connected windows

An exit can be omitted when the last selected window is in region b and
ends at the global stream end. Its start is fixed by its length. Such
families have count at most

\[
 N_h^{\rm tail}
 =\binom{b-1}{h-1}\binom Fh(L/a)^{h-1}.
\]

Charging the other h-1 exits gives the additional event bound

\[
 \sum_{h=p+1}^{\min(b,F)}
 N_h^{\rm tail}
 \min\{1,z_h^{-[s(h-1)+1]}M_h(z_h)^L\}.
 \tag{5}
\]

The moment proof is unchanged. Adding (4) and (5) is conservative even
though the events overlap. The free-exit condition refers to the global
stream end, not merely to the last selected window.

There is a separate connection issue. Suppose a selected window ends at
the end of region j, and the next selected window starts at the beginning
of region j+1. No output coordinate separates them. The inner state need
not vanish between those two windows. More generally, if the selected
windows form r connected intervals in the global stream, only r exits
need be charged; this becomes r-1 when the final interval reaches the
stream end. The relevant sufficient condition is then

\[
 D>s(r-\nu),\qquad
 \nu=\mathbf1\{\text{last connected interval reaches the stream end}\}.
\]

Equations (4) and (5) do not alone cover all these smaller thresholds.
A broader capacity bound must enumerate or dominate the touching-boundary
patterns and use their corresponding exit count. This does not invalidate
(4); it limits the structural event that (4) excludes.

## Code audit and scope

At the inspected revision, `window_dimension_mgf.py` computes (3), then
the h-sum in (4). Its `sh+1` crossing, positive-composition count, start
bound, and equal-probability direction match this proof. The loop assumes
the budget permits every h being visited; a generic small-geometry caller
should omit h>F rather than take the logarithm of C(F,h)=0.

This audit did not rerun the numerical optimizer. Floating minimization
does not certify a bit margin. A directed replay must upper-bound the
moment power and the union sum at specified marks. Exact or directed
evaluation of (4) or (5) certifies only the stated obstruction event.
It does not imply minimum distance, bound arbitrary low-output codewords,
or complete a proof for any particular inner coefficient family.

The next proof step is therefore to settle the connected-window counting
and the exact terminal convention before attaching a numerical capacity
margin to the broader obstruction family. The remaining whole-code
distance argument is a separate requirement.

## Connected-window extension and directed replay

`theory/connected_window_capacity.py` now includes boundary joins. Fix h
used regions whose windows form c connected intervals. Let v of these
intervals occupy a single region, and put G=L/a. For the moment, charge
an exit at every interval end, including the global end if it is reached.

There are

\[
 \binom cv\binom{h-c-1}{c-v-1}
\]

ordered component-length patterns when c-v>0. If every component is a
singleton, there is one pattern precisely when h=c=v. Placing the ordered
components into b regions, with zero or more unused regions between
components, gives at most `C(b-h+c,c)` placements. Adjacent components
are allowed in this count; constraints preventing their windows from
touching are dropped conservatively.

A component spanning at least two regions has two variable endpoint
lengths and full interior regions. Thus there are `f=h-2c+v` forced full
regions and `d=2c-v` variable lengths. A singleton component has at most G
starts. Each other endpoint is anchored to its regional boundary. The
total family count is consequently at most

\[
 \binom cv\binom{h-c-1}{c-v-1}
 \binom{b-h+c}{c}\binom{F-fG}{d}G^v. \tag{6}
\]

Terms with insufficient remaining length are omitted. The length count
allows variable lengths greater than G, so it is only an upper bound.
The singleton-only convention above replaces the first two factors when
their displayed expression is undefined.

For a component ending at the global stream end, remove its exit and
fix the placement's final endpoint. This replaces the placement factor
by `C(b-h+c-1,c-1)`. There are two cases:

- If the last component is a singleton, replace `C(c,v)` by
  `C(c-1,v-1)` and `G^v` by `G^(v-1)`.
- Otherwise replace `C(c,v)` by `C(c-1,v)`. Its last endpoint region is
  full, so increase f by one and decrease d by one.

Every row hits all f full regions. In (3), use d trials, parity p-f,
and probability `min(1,(F-fG)/(Gd))`. The integer crossing is `sc+1`
without a free final exit and `s(c-1)+1` with one. Sum the resulting
Chernoff bounds over all feasible h,c,v and the three endpoint cases.
Every family with at most one aligned window per region is included,
using its actual connected components. Overlap between the counted cases
only enlarges the bound.

For GRS[64,44], L=23872, s=8, a=16, and e=152780, this procedure checks
17,160 grouped terms. Rational-mark replay with outward 80-digit decimal
intervals gives

\[
 \Pr[\exists\text{ counted family with }D>s(c-\nu)]
 <5.375\cdot10^{-48}<2^{-157}.
\]

The floating optimizer suggests the marks; it does not certify the
inequality. The replay evaluates every moment, power, count, and final
sum with outward intervals. Five independent tests verify component
patterns, count domination against all small window families, tail
capacity, and interval containment of exact rational expressions.
A separate read-only audit checked the counting formulas.

This excludes the specified dimension-forced mechanism with 157 bits of
margin. It does not exclude every codeword supported in these windows,
multiple windows in one region, or other low-weight outputs. It is not
a whole-code distance certificate and does not justify a native kernel.

```sh
python -B -m unittest discover -s research/workstreams/char2_field_spin/theory -p connected_window_capacity_tests.py
python -B research/workstreams/char2_field_spin/theory/connected_window_capacity.py --rows 23872 --outer 64 44 --state 8 --replay
```
