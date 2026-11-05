# A three-state obstruction to the regional Cauchy comparison

2026-11-01. The adaptive population majorant with `y=2` and quota
evaluations `x=1,2,4` cannot produce the proposed threshold-60518
geometry bound through one scalar supersolution. This remains true
when its physical population action is evaluated exactly, with no
bins or exponential projection. A three-state subkernel proves a
regional growth floor of `2^518` and a complete comparison floor of
`2^79750`. An additional return path proves the same obstruction
even after splitting destinations, at this fixed evaluation grid.
A subsequent exact certificate gives a `2^400` growth floor for the
rowwise Cauchy comparison over all positive evaluations. That stronger
scope does not include destination-separated quota comparisons.

This is an obstruction to the specified Cauchy comparison, not to the
exact quota operator, the actual flow-tail probability, or the code.
The encoder and setup distributions remain unchanged.

## The comparison and its three retained states

Use the original `(1,1,2)` proof network and the exact adaptive
population operator `K_x` from
[QUARTET_EXACT_POPULATION_HORNER.md](QUARTET_EXACT_POPULATION_HORNER.md).
Put `m=128`, `C0=binom(8192,4096)`, and define the positive homogeneous
map

```
F(phi)(s)=min_(x in {1,2,4}) x^-4096*(K_x^16 phi)(s)/C0.
```

The exact quota operator `B=[x^4096]K_x^16/C0` satisfies
`B phi<=F(phi)`. No reverse inequality is assumed.

Retain the physical states

```
O=(0,0),       P=(512,0),       H=(512,512).
```

Their old quartet types `(u,w)` are respectively `(0,0)`, `(4,0)`,
and `(4,4)`. Selecting every quartet gives destination `O`. Choosing
every quartet as a hole gives destination `(512,a)`, where `a` is
the old hole population. Thus the exact submatrix on these states is

```
L_x = [[S00(x)^128, T00(x)^128, 0],
       [S40(x)^128, 0,           T40(x)^128],
       [S44(x)^128, 0,           T44(x)^128]].       (1)
```

Here `S` denotes selection and `T` a hole. The adaptive probabilities
give the following exact rational polynomials:

```
T00(x)=1+4x/5,
S00(x)=(1+x)^4-T00(x),
T40(x)=[256(1+x/4)^4-(1+x)^4]/255,
S40(x)=(1+x)^4-T40(x),
T44(x)=16*[256(1+x/4)^4-(1+x/2)^4]/255,
S44(x)=(2+x)^4-T44(x).                              (2)
```

All coefficients are nonnegative. Pure old states have only one
quartet type, so (1) has no residual partition average or population
multinomial divisor. Other population paths are discarded positively:
`(K_x^16)|_{O,P,H}>=L_x^16` entrywise.

## An exact lower bound on all nine row/evaluation pairs

Let `v` be one on `O,P,H` and zero on every other physical state.
The following paths retain sixteen transitions. A repeated endpoint
means that the path stays there for the remaining transitions.

| Evaluation | Start `O` | Start `P` | Start `H` |
| --- | --- | --- | --- |
| `x=1` | `O -> P -> H`, then 14 `H` loops | `P -> H`, then 15 `H` loops | 16 `H` loops |
| `x=2` | 16 `O` loops | `P -> O`, then 15 `O` loops | 16 `H` loops |
| `x=4` | 16 `O` loops | `P -> O`, then 15 `O` loops | `H -> O`, then 15 `O` loops |

At `x=1`,

```
T00=9/5>1,       T40=203/85>2,       T44=9919/255>32.
```

At `x=4`,

```
S00=3104/5>512,  S40=51968/85>512,  S44=53248/51>1024.
```

At `x=2`, all three factors used by the paths exceed `2^(25/4)`:

```
S00=392/5,       S40=1296/17,       T44=4096/51.
```

Their fourth-power comparisons reduce to small integer inequalities:

```
49^4=5764801 >5120000=625*2^13,
3^16=43046721>42762752=512*17^4,
51^4=6765201 <8388608=2^23.                         (3)
```

We use `C0<2^8186`. One elementary proof is the induction

```
binom(2n,n)^2 <=16^n/(3n+1).
```

For its successive-step comparison, the difference between
`4(n+1)^2(3n+1)` and `(2n+1)^2(3n+4)` is `n`. At `n=4096`,
`3n+1>2^12`, which gives the stated bound on `C0`.

After including `x^-4096/C0`, the selected paths give these strict
base-two exponent lower bounds:

| Evaluation | Row `O` | Row `P` | Row `H` |
| --- | ---: | ---: | ---: |
| `x=1` | 902 | 1542 | 2054 |
| `x=2` | 518 | 518 | 518 |
| `x=4` | 2054 | 2054 | 2182 |

For example, every `x=2` path has numerator greater than
`(2^(25/4))^(128*16)=2^12800`; subtracting 4096 and 8186 gives 518.
Consequently

```
F(v)>=2^518*v.                                     (4)
```

Outside the retained states the right side is zero, so (4) is a
valid statement on the entire physical state space.

## Every positive scalar potential inherits the obstruction

Suppose a strictly positive physical potential `phi` satisfies
`F(phi)<=rho*phi`. Since the state space is finite, `v<=c*phi` for
some finite positive `c`. Monotonicity, homogeneity, and iteration give

```
2^(518k)*v <= F^k(v) <= c*rho^k*phi.
```

At any of the three retained states, taking `k` to infinity proves
`rho>=2^518`. This argument does not require the zero-extended `v`
to be strictly positive, and does not require a stationary distribution.

The endpoint cannot cancel this growth floor. Let

```
b=[x^4096]K_x^15 q_x/C0,
q_x(a,w)=(2+x)^w(1+x)^(512-w).
```

For fixed labels and quartet partitions, the adaptive hole/selection
probabilities sum to one, including forced choices. Each complete
coin history has cover cost at least zero. Its normalized expectation
of `2^cost` is therefore at least one. Averaging over the exact balanced
regional label quota and partitions gives `b(s)>=1` for every old
state. This argument uses the conditional normalized expectation,
not just a lower bound on the endpoint polynomial.

For any valid endpoint upper bound `b_upper` and
`nu=max_s b_upper(s)/phi(s)`, it follows that
`phi(H)*nu>=b(H)>=1`. The full scalar comparison is consequently at
least

```
C0*2^-60518*phi(H)*rho^255*nu
 >=C0*2^(-60518+255*518)
 >2^79750.                                         (5)
```

The final strict inequality uses
`C0>=2^8192/8193>2^8178`, obtained from the largest coefficient of
`(1+x)^8192`. Thus neither further potential iteration nor a different
positive potential can make this particular comparison close.

## The fixed grid also defeats destination splitting

A destination-separated quota bound changes the order of the minimum
and destination sum. Define the three-by-three matrices

```
A_x=x^-4096*L_x^16/C0,
M_ij=min_(x in {1,2,4}) (A_x)_ij.                   (6)
```

The all-ones argument alone would not establish a floor for `M*v`.
However, the `H,H` entry itself has the same lower floor. At `x=1`
and `x=2`, the sixteen `H` loops already give exponents 2054 and 518.
At `x=4`, retain the path

```
H -> O, then13 O loops, then O -> P -> H.
```

This has exactly sixteen transitions. In addition to the preceding
selection bounds, use `T00(4)=21/5>4` and
`T40(4)=1157/85>8`. Its normalized exponent is strictly greater than

```
128*(10+13*9+2+3)-8192-8186=518.
```

It follows that `M_HH>2^518`. The same lower bound holds for the
`H,H` entry of the full physical regional matrix after an entrywise
minimum over this grid. Every positive scalar supersolution of that
matrix therefore has growth at least `2^518`.

For a full physical implementation, disjoint destination groups
`J_1,...,J_k` permit

```
B phi <=sum_j min_x x^-4096*K_x^16(1_(J_j)*phi)/C0. (7)
```

Each group is propagated through sixteen fixed-evaluation updates
before its quota minimum. The minima are then summed. Equation (7)
does not authorize a per-update quota choice or an overlapping
destination decomposition without accounting for its multiplicity.
The group containing `H` retains the diagonal path above, so no such
split removes the obstruction at `x=1,2,4`.

This does not exclude other positive quota evaluations. The next
bounded diagnostic is the single scalar function

```
x -> x^-4096*(L_x^16)_HH/C0,          x>0.
```

A rigorous lower bound for its infimum would test whether the retained
return family obstructs arbitrary positive Cauchy evaluations. A small
value would only remove this particular obstruction, not prove a useful
full-population bound. The diagonal diagnostic alone makes no claim
about all positive evaluations or the exact quota coefficient.

## Bounded continuous-evaluation diagnostic

The parent completed a serial floating minimization of the three-state
`H,H` expression. It found `x=2.756067424` with logarithmic diagonal
value `177.10420`. Direct evaluations at `x=2.75` and `x=2.8` gave
approximately `178.97` and `209.03`, respectively. The available bulk
growth budget is about `205.2268` bits before endpoint costs.

Thus a refined evaluation escapes the diagonal obstruction proved for
`x=1,2,4`. These floating values are not directed certificates or a
bound on the full population operator. They do not show that a useful
positive potential exists.

The next bounded diagnostic compares two three-state objects:

```
G(v)_i=inf_(x>0) sum_j (A_x)_ij*v_j,
M_ij=inf_(x>0) (A_x)_ij,
```

where `A_x=x^-4096*L_x^16/C0` as above. A scalar growth lower bound
for `G` tests whether finer quota evaluations alone can remove the
retained three-state obstruction. The Perron growth of `M` tests the
stronger entrywise destination-separated comparison. Always
`M*v<=G(v)` for `v>=0`; changing this order can therefore matter.

Both are lower diagnostics for their corresponding full physical
comparisons. A value above the available budget identifies an
obstruction only after rigorous validation. A value below the budget
merely leaves that comparison viable.

The parent subsequently ran these tiny three-state controls serially.
After twelve half-damped steps, the rowwise diagnostic reported growth
logs between `412.813` and `413.129`, with potential proportional to
approximately `(2^-676.5,2^-363.276,1)`. The entrywise-continuous matrix
had floating Perron growth log `177.1042`. The following certificate
turns the former observation into an all-evaluation theorem. The
entrywise result remains a floating diagnostic and supplies no upper
bound on the full physical operator.

## Exact obstruction over all positive rowwise evaluations

Extend the positive three-state vector

```
v=(1,2^314,2^677)
```

by zero outside `O,P,H`. For each retained row `i`, define

```
P_i(x)=(L_x^16*v)_i,
f_i(t)=log2(P_i(2^t)/(C0*v_i))-4096*t.
```

The polynomial `P_i` has nonnegative coefficients. Consequently
`f_i` is convex, with exact slope at a rational evaluation `x`

```
s_i(x)=x*P_i'(x)/P_i(x)-4096.
```

Choose rational brackets `x_l,x_u` for which
`s_i(x_l)<0<s_i(x_u)`. Every global minimum then lies between their
logarithms. The decreasing supporting line at the left endpoint gives

```
inf_t f_i(t)
 >=f_i(log2 x_l)+s_i(x_l)*log2(x_u/x_l).             (8)
```

The root-owned exact checker
[quartet_cauchy_phase_obstruction.py](quartet_cauchy_phase_obstruction.py)
evaluates `P_i` and its logarithmic derivative by positive integer
three-by-three matrix powers. It bounds logarithms with directed
rational intervals from a dyadic mantissa and a positive logarithm
series. The displayed decimal values summarize rational lower bounds;
the actual comparisons with 400 use exact fractions.

| Row | Rational bracket | Lower bound from (8), decimal display |
| --- | --- | ---: |
| `O` | `1819/1000,1822/1000` | 410.755936 |
| `P` | `2112/1000,2116/1000` | 410.358917 |
| `H` | `2316/1000,2320/1000` | 411.295089 |

The parent exact replay passed all three opposite-slope tests and all
three strict lower-bound comparisons. Therefore the continuum rowwise
map `G` satisfies `G(v)>=2^400*v`. Its full physical counterpart

```
F_cont(phi)=inf_(x>0) x^-4096*K_x^16(phi)/C0
```

has the same zero-extended subeigenvector inequality. The comparability
proof above implies that every strictly positive scalar supersolution
has growth `rho>=2^400`. The endpoint lemma then gives the complete
comparison floor

```
C0*2^-60518*phi(H)*rho^255*nu >2^49660.              (9)
```

Thus finer quota evaluations or further positive-potential iteration
cannot close the unsplit rowwise comparison, even when all positive
evaluations are available. This theorem still does not refute the
exact quota operator. It also does not refute the destination-separated
comparison in (7): placing the minimum before the destination sum
changes the positive operator. The next justified full-size test must
retain that distinction, and should first use the small entrywise
diagnostic to choose its quota evaluations.

The proof, its endpoint/comparability argument, the all-ones constants,
and the diagonal return-path extension passed independent mathematical
audit. The parent completed an exact serial replay of all nine exponent
bounds, both binomial comparisons, and the complete exponent 79750;
all passed. The parent also passed the exact continuum certificate
giving (9). Independent source review passed its integer matrix and
derivative powers, logarithm intervals, slope signs, and final exact
comparisons. No numerical job or helper implementation was performed
by the note's author. The full-code certificate remains open.
