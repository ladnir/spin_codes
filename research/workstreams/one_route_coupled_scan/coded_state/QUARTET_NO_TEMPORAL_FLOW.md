# No-temporal flow: a valid field bridge and a Hall obstruction

2026-11-01. Removing temporal edges from the proof network gives a
strong all-death field census and an exact 129-state adaptive operator.
However, disjoint Hall witnesses prove that this network cannot meet
the required half-row deficit tail. This is a limitation of the proof
network, not an encoder change or a counterexample to the code.

Use `n=512`, `T=4096`, `N=2097152`, `M=8192`, `W=943718`, and
`q=2^32`. Fix a row set `H` of size 4096. Each region has exactly
4096 `H` inputs and 4096 non-`H` inputs, uniformly positioned by its
independent permutation. Write `J_t=H_t^c`.

## The valid all-death field calculation

Take `(a,b,c)=(1,0,2)` in
[QUARTET_GENERAL_FLOW_COVER.md](QUARTET_GENERAL_FLOW_COVER.md).
Sources have capacity one, non-`H` inputs capacity two, and each
source can use its own input or its next quartet. There are no temporal
edges. The integer deficit is

```
delta_H=max_S [L-2|V(S) outside H|].                  (1)
```

For an actual support whose eligible set is `H`, the rank envelope is

```
R >= L/2+524288-delta_H/2.
```

If every half-row set satisfies `delta_H/2<=52429-s`, the fixed-support
field bound at `x=q^-1/2` becomes

```
gamma/(q-1) * q^-s * A^(L-D-a_T) * z^D,
A=1+2^-16,       z=2^-16.                            (2)
```

Here `a_T` is the terminal state population. Equation (2) retains all
deaths and output cancellations. Summing over every binary state
history, starting from zero and with a free endpoint, gives the exact
per-coordinate factor

```
Z_T=e_0 [[1,A],[z/A,A]]^T (1,1/A)^transpose.
```

Thus the conditional probability of any bad output with eligible-row
count 4096 is at most `gamma/(q-1)*q^-s*Z_T^512`, with no death or
volume restriction and no additional row-set union.

Put `t=2^-8`, so `A=1+t^2`. Conjugating by
`diag(1,t/A)` gives the symmetric transfer

```
[[1,t],[t,1+t^2]],
```

with endpoint `(1,256)`. Its largest row sum is strictly below
`256/255`: after subtracting one, this is
`t+t^2=257/65536<1/255`, since `257*255=65535<65536`.
Equivalently, `65793*255=16777215<16777216`. Consequently

```
Z_T^512 < 2^4096*(256/255)^N < 2^15996.              (3)
```

For the last inequality use `ln(1+1/255)<=1/255`,
`ln2>56/81`, and
`N*81=169869312<169932000=11900*56*255`.
The parent verified these integer inequalities, the accepted-field
factor `gamma/(q-1)<2^-31`, and the following consequences exactly:

```
s=512  gives  15996-31-32s=-419,
delta cap=2*(52429-512)=103834.                      (4)
```

Therefore a sufficiently small probability of any half-row set having
deficit at least 103834 would yield a conditional field contribution
below `2^-419` for the entire half-row class. Using threshold 103834
is conservative: strict violation of the displayed cap begins at
103835. The next section proves that even the weaker threshold 104860
has substantial probability for this proof network.

## Disjoint Hall witnesses rule out the required tail

For each `t<T` and each quartet `Q` in the partition of update `t+1`,
let

```
r_(t,Q)=|H_t intersect Q|,
j_(t,Q)=|J_(t+1) intersect Q|.
```

Consider the `r_(t,Q)` sources at time `t` whose own input is in `H`.
Their only positive-capacity neighbors are the `j_(t,Q)` inputs in
that next quartet. Their Hall deficit is therefore
`(r_(t,Q)-2j_(t,Q))_+`. Retain this source set only when its deficit
is positive. These retained groups have disjoint sources and disjoint
positive-capacity input neighborhoods across both quartets and times.
Terminal sources with own input in `H_T` have no neighbor at all.
Taking the union of these witnesses proves, deterministically,

```
delta_H >= |H_T|+
 sum_(t=1)^(T-1) sum_Q (r_(t,Q)-2j_(t,Q))_+.        (5)
```

Under independent fair labels, `r,j` are independent `Binomial(4,1/2)`.
Only `j=0,1` contribute, and

```
E(r-2j)_+
 = (1/16)*2+(4/16)*(6/16)=7/32.                    (6)
```

The exact regional quotas preserve a sufficient lower bound. Condition
on any quartet partition. Each term in (5) depends on eight distinct
input positions: four at time `t` and four at time `t+1`. For any
prescribed pattern on these positions, each successive prescribed
label has conditional probability at least `4089/8192`. This holds
whether the positions lie in one region or cross a region boundary.
Thus every pattern has probability at least

```
2^-8*(1-7/4096)^8 >= 2^-8*(505/512).                 (7)
```

Multiplying the independent-label expectation in (6) by this uniform
pattern ratio is valid because the summand is nonnegative. Also
`E|H_T|=256` exactly. Equations (5)-(7) give

```
E delta_H >= 256+(505/512)*(7/32)*128*4095
            =14508593/128
            =113348.3828125.                        (8)
```

Since `0<=delta_H<=N`, for every integer `d` one has

```
Pr(delta_H>=d) >= (14508593/128-(d-1))/(N-(d-1))
```

whenever the numerator is positive. In particular,

```
Pr(delta_H>=103834)>1/256,
Pr(delta_H>=104860)>1/256.                           (9)
```

The parent verified the expectation and both rational tail bounds by
exact arithmetic. The respective lower bounds are
`1217969/255144832` and `1086641/255013504`, each greater than `1/256`.

These are lower bounds for a single fixed `H`, averaged over the
actual regional routing and fresh quartet partitions. No independence
between different Hall summands is required. They rule out the small
uniform-half-row deficit tail needed by (4), even before the row-set
union. They do not imply a low-weight codeword: the Hall witnesses need
not have this `H` as their eligible set or satisfy the encoder's field
equations.

Two independent mathematical reviews found the same disjoint-star
family and passed (5)-(9), including exact quota conditioning and the
integer event threshold. No numerical jobs or no-temporal operator
implementation were needed.

An independent full-note audit also passed the field transfer and its
endpoint, all scalar constants, and the population kernel below.

## The exact 129-state operator, for reference

The operator itself remains mathematically valid. In reverse time let
`A=V_(t+1)` be the next hole set; only its population `ell=|A|` is
needed. It is divisible by four. For a current quartet put

```
u=|A intersect Q|, j=|J_t intersect Q|, r=|J_t intersect A intersect Q|,
B=2j-r,       A_loss=max(0,4-B).
```

Use the adaptive probabilities from the general-flow note with these
loss/gain bounds. Define nonnegative quartic label polynomials

```
T_u(x)=sum_(j,r) binom(u,r)binom(4-u,j-r)
                    *x^j*y^(u-r)*g(j,r),
S_u(x)=sum_(j,r) binom(u,r)binom(4-u,j-r)
                    *x^j*y^(u-r)*p_select(j,r),
T(x,z)=sum_(u=0)^4 binom(4,u)T_u(x)z^u,
S(x,z)=sum_(u=0)^4 binom(4,u)S_u(x)z^u.
```

With `m=n/4`, the exact fresh-partition average is

```
K_x[ell,4k]=binom(m,k)/binom(n,ell)
             *[z^ell] T(x,z)^k S(x,z)^(m-k).        (10)
```

The denominator counts old marked positions; the factor `binom(m,k)`
chooses the hole quartets. The physical start is `ell=n`, and the
first-layer endpoint is

```
q_x(ell)=(y+x)^ell*(1+x)^(n-ell).
```

The exact quota operators are
`B=[x^4096]K_x^16/binom(M,4096)` and
`b=[x^4096]K_x^15 q_x/binom(M,4096)`, with one row-set union outside
`e_n B^255 b`. There is no extra terminal factor `y^n`. At fixed
positive `x`, (10) needs only positive one-variable polynomial products
and a 129-by-129 matrix; no population bins or large Fourier transform
are needed. Nonetheless any valid upper majorant of this deficit
moment must respect (9), so evaluating this operator cannot prove the
proposed half-row tail.

Next work should retain enough temporal capacity to repair these Hall
witnesses while accounting for the resulting death weight, or impose
additional eligible-support constraints. Further tuning of the
no-temporal uniform-row deficit bound cannot close the target.
