# A quantitative flow bound for all low-death output amplitudes

2026-11-01. Excluding a dimension obstruction alone does not prove the
distance target. A stronger flow bound, with positive rank margin, does
give a complete field-probability bound for every actual state history
with a bounded number of deaths. This note states that conditional
theorem, including output cancellations and the free terminal state.
No required flow-tail estimate is assumed to have been proved.

Use `q=2^32`, `n=512`, `T=4096`, `N=nT=2097152`, `K=1048576`, and
`W=943718`. The accepted-field density loss is `gamma<=1024/989`.
All setup and rank conventions are those in
[QUARTET_GLOBAL_RANK_RESET.md](QUARTET_GLOBAL_RANK_RESET.md) and
[QUARTET_GLOBAL_DEATH_COVER.md](QUARTET_GLOBAL_DEATH_COVER.md).

## Conditional geometry and field theorem

Fix all routing and quartet partitions. For an actual state support
`S`, let `L` be its volume, `D` its deaths before the free endpoint, and
`a_T` its terminal population. Its footprint has parity rank `R` and
eligible row set `H`, with `h=|H|` and `p=h/8192`.

Suppose the geometry satisfies, simultaneously for every row set of
this size,

```
tau_H <= p(K-W)-s,                                    (1)
```

where `tau_H` is the physical deficit defined in
[QUARTET_FIXED_ROW_FLOW.md](QUARTET_FIXED_ROW_FLOW.md). The cut identity
then gives every actual support with this eligible-row count

```
R >= (1-p)L-pD+pK-tau_H
  >= (1-p)L-pD+pW+s.                                 (2)
```

For `0<x<=1`, the fixed-support amplitude bound from the death-cover
note is

```
gamma/(q-1) * q^-R * x^-W * (1+qx)^L
 * [qx^2/(1+qx)]^D * [qx/(1+qx)]^a_T.                (3)
```

It bounds all accepted-field bad outputs at this support, after applying
the exact-live field rank and compatible birth equations. Required
nonzero states were relaxed only after that rank step. Thus it already
includes outputs whose weight is smaller than `L+D` through cancellations.

Set `x=q^-p` and `A=1+q^(p-1)`. Substituting (2) into (3) gives

```
gamma/(q-1) * q^-s * A^(L-D-a_T)
 <= gamma/(q-1) * q^-s * A^N.                        (4)
```

No condition `L+D<=W` is used. The terminal factor in (3) accounts for
the free endpoint and makes the exponent exactly `L-D-a_T`.

Each state support is uniquely encoded by the positions where its
binary state changes, starting from the empty state. It has `D+a_T`
birth transitions and `D` death transitions, hence `2D+a_T` changes.
For an integer `d>=0`, put

```
Census(d)=sum_(i=0)^min(N,2d+512) binom(N,i).
```

There are at most `Census(d)` supports with `D<=d`. Summing (4) proves
the conditional bound

```
Pr_field[bad output, eligible count h, D<=d | fixed good geometry]
 <= gamma/(q-1) * q^-s * (1+q^(p-1))^N * Census(d).   (5)
```

This covers all amplitudes and all actual state volumes in that class.
There is no additional row-set union in (5): condition (1) already holds
for every row set of size `h`. Different eligible-row counts can be
summed once, together with the probability that their geometry bounds
fail. The existing large-footprint and other field classes retain their
separate failure allowances.

The class statement is over `GF(2^32)`. Once every bad-output class is
covered over that field, a nonzero basis component gives the usual
transfer to the 128-bit extension. A restricted low-death class does
not by itself transfer with the same death or eligible-row label.

## Explicit sufficient margins for the half-row case

For every `p<=1/2`,

```
(1+q^(p-1))^N < 2^47,
gamma/(q-1) < 2^-31.
```

For the first inequality, `ln(1+u)<=u` gives a natural-log exponent at
most 32, and `ln2>56/81` gives `32/ln2<47`. Therefore any integer
`C_d` satisfying `Census(d)<=2^C_d` gives the simple bound

```
Pr_field[class in (5)] < 2^(C_d+16-32s).             (6)
```

The parent verified the following constants by exact integer arithmetic.
They suffice for a per-`h` bound strictly below `2^-55`:

| Death bound d | Census exponent C_d | Required margin s | Physical tau cap at h=4096 | Integer deficit cap |
| ---: | ---: | ---: | ---: | ---: |
| 65536 | 709338 | 22170 | 30259 | 60518 |
| 98304 | 943002 | 29472 | 22957 | 45914 |
| 131072 | 1141364 | 35670 | 16759 | 33518 |

Here `p(K-W)=52429`, the physical cap is `52429-s`, and the integer
network doubles that cap. Each row satisfies `32s>=C_d+71`; the actual
exponents in (6) are respectively `-86,-86,-60`. These are sufficient
flow targets, not measured or proved deficit tails.

For `k=2d+512<N/2`, consecutive binomial coefficients give the exact
geometric-series upper bound

```
Census(d) <= binom(N,k)*(N-k+1)/(N-2k+1).
```

The parent used this standalone integer recipe, without a floating
logarithm:

```
from math import comb
N = 1 << 21
for d in (65536, 98304, 131072):
    k = 2*d + 512
    upper_floor = comb(N, k)*(N-k+1)//(N-2*k+1)
    C_d = upper_floor.bit_length()
    s = (C_d+71+31)//32
    tau_cap = 52429-s
    print(d, C_d, s, tau_cap, 2*tau_cap, C_d+16-32*s)
assert 11**32 < 2**111
assert 1024*2**32 < 1978*(2**32-1)
```

The last assertions also verify the two scalar factors used in (6):
`e<11/4` proves `A^N<2^47`, and the second inequality proves
`gamma/(q-1)<2^-31`. The exact check ran serially in 2.41 seconds on
the parent-controlled host. No flow probability was computed by it.

## Exact run and death transfers for sharper comparisons

If a sharper census is useful, the per-coordinate state-support
transfer with volume marker `u` and death marker `z` is

```
M(u,z)=[[1,u],[z,u]],
support polynomial=e_0 M(u,z)^T (1,1)^transpose.
```

At `u=1`, its death polynomial is exactly

```
A_T(z)=sum_d binom(T+1,2d+1)z^d
      =[(1+sqrt(z))^(T+1)-(1-sqrt(z))^(T+1)]/(2sqrt(z)).
```

The full support census with `D<=d` is bounded by
`inf_(0<z<=1) z^-d A_T(z)^512`. This keeps the alternating changes on
each coordinate, whereas (5) only uses their total number.

The field and output markers can also remain in the transfer. Put

```
a=(1+qx)/q^(1-p),
b=q^(p+1)x^2/(1+qx),
t=qx/(1+qx).
```

Under (1), the low-death field bound is, for `0<z<=1`,

```
gamma/(q-1) * q^(-pW-s) * x^-W * z^-d
 * [e_0 [[1,a],[bz,a]]^T (1,t)^transpose]^512.        (7)
```

It retains all output cancellations and all deaths up to `d` without
enumerating supports. At `z=1`, `x=r q^-p`, a diagonal similarity gives
the symmetric matrix

```
[[1,r],[r,r+q^(p-1)]]
```

with free endpoint `(1,1)`, and prefactor
`gamma/(q-1) q^-s r^-W`.

## Why the death restriction matters

Even an ideal half-row bound `tau_H=0` does not make the unrestricted
support census close. Consider formal supports where every coordinate
has volume 1152, exactly 688 closed runs, and zero terminal state.
Their number is

```
[binom(1151,687)binom(2944,688)]^512.
```

They have `L=589824`, `D=352256`, and `L+D=942080<W`. The flow rank
envelope at `tau=0` gives surplus `(K-L-D)/2=53248`. The highest-output
term in the relaxed enumerator (3), after multiplication by `x^-W`,
contributes at least `q^-53248` per support for every `0<x<=1`, before
the common factor `gamma/(q-1)`.

The parent verified the support count directly with exact integers:

```
c = comb(1151, 687)*comb(2944, 688)
assert c.bit_length()-1 == 3417
assert 512*3417 - 32*53248 - 32 == 45536
```

Thus `c>=2^3417`. The rank penalty is `32*53248` bits, and
`gamma/(q-1)>q^-1` costs fewer than another 32 bits. The unrestricted
relaxed count is therefore larger than `2^45536`, even at `tau=0`.

These formal supports are not asserted to have eligible count 4096 or
to realize compatible codewords. Many may already lie in a previously
closed footprint class. The example identifies a loss in dropping
those constraints; it does not disprove a distance theorem. The new
bridge closes low-death classes only after their quantitative flow
tails are proved. Larger-death histories still require an explicit
coverage argument using the established rank, footprint, or amplitude
bounds.
