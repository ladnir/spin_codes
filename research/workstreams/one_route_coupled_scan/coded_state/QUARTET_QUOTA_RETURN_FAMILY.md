# A partial-population return family is too small to obstruct exact quotas

2026-11-01. A proposed lower test retained three-update returns through
partial populations in the original adaptive Hall operator at `y=2`.
The family preserves the exact regional quota, but an analytic upper
bound makes its coefficient computation unnecessary. Even all renewal
layouts of these returns contribute less than `2^103` per region.
Their complete restricted contribution at threshold 60518 is less than
`2^-21971`, including the row-set union and the unrestricted endpoint.

This bounds only the specified return family. It does not upper-bound
the full adaptive operator or certify the code. The construction and
setup distribution remain fixed.

## The exact return family

Use the population state `(a,w)` and the original `(1,1,2)` proof
network from [the Horner note](QUARTET_EXACT_POPULATION_HORNER.md).
Set `n=512`, `m=128`, and `C0=binom(8192,4096)`.
The polynomial `T_uw(x)` sums the local hole branch for an old quartet
with `u` positions in `A` and `w` positions in `W`. The polynomial
`S_uw(x)` sums its selected branch. Both include current label
multiplicities and cover weight `2^(w-r)`, where `r` counts non-H labels
in old `W`. Neither includes the old-color multinomial.

For `0<=k<m`, define

```
H=(n,n),       D_k=(4k,4k),       P_k=(n,4k).
```

Choosing exactly `k` hole quartets from `H` gives `D_k`. An all-hole
update from `D_k` gives `P_k`; another all-hole update gives `H`.
The three transition polynomials are

```
A_k(x)=binom(m,k) T_44(x)^k S_44(x)^(m-k),

B_k(x)=[z^(4k)] [sum_(u=0)^4 binom(4,u) T_uu(x) z^u]^m
        /binom(n,4k),

C_k(x)=[z^(4k)] [sum_(w=0)^4 binom(4,w) T_4w(x) z^w]^m
        /binom(n,4k).                                      (1)
```

Each denominator averages a fresh uniform partition. Both all-hole
transitions have deterministic destination populations, so no further
destination normalization occurs. Define the cycle and loop polynomials

```
P(x)=sum_(k=0)^(m-1) A_k(x) B_k(x) C_k(x),
h(x)=T_44(x)^m.                                            (2)
```

The intermediate states of every cycle are non-H. Thus a specified
sequence of one-update loops and three-update cycles has no duplicate
representation. With `r` cycles and `16-3r` loops, its exact regional
contribution is

```
[x^4096] h(x)^(16-3r) P(x)^r / C0,       0<=r<=5.           (3)
```

The initially proposed balanced extraction retained 256 labels in one
loop and 768 in each of five cycles. Positivity makes that a valid
lower bound on (3) for `r=5`. The upper bound below applies directly
to (3), including every allocation of the regional quota.

## A uniform evaluation at x=4

Let `j` be the number of non-H labels in a quartet, and let `r` count
those labels in old `W`. When both adaptive choices are allowed, the
hole probability satisfies

```
g(u,j,r) <= 2^(-(2j-r)) / (1-2^(-(4+u))).                   (4)
```

The inequality also holds for the forced choices. At `x=4`, multiplying
(4) by the label marker and cover cost cancels its dependence on `j,r`:

```
4^j 2^(w-r) g(u,j,r) <= 2^w/(1-2^(-(4+u))).
```

Summing the sixteen label masks gives

```
T_uu(4) <= 16*2^u*(16/15),
T_4w(4) <= 16*2^w*(256/255).                               (5)
```

Apply (5) inside the two coefficient formulas in (1). The binomial
normalizers cancel exactly, yielding, for `c=4096/3825`,

```
B_k(4) C_k(4) <= 2^1024 * 2^(8k) * c^128.                 (6)
```

The full-state coefficients obey the exact identities

```
T_44(4)=64240/255 <256,
S_44(4)+T_44(4)=6^4,
S_44(4)+256*T_44(4)=65536.                                 (7)
```

Extend the sum in (2) positively to `k=m`. Equations (6)-(7) imply

```
P(4) <= 2^1024*c^128 * [S_44(4)+256*T_44(4)]^128
     = 2^3072*c^128,
h(4) <2^1024.                                             (8)
```

For any layout with `r` cycles, its numerator in (3) is at most
`4^-4096 h(4)^(16-3r) P(4)^r`. Since `C0>2^8178`, (8) gives

```
[x^4096]h^(16-3r)P^r/C0 <2^14*c^(128r) <=2^14*c^640.
```

Two small integer comparisons give a convenient uniform ceiling:

```
4096*14=57344 <57375=3825*15,
15^8=2562890625 <2951578112=2*14^8.
```

Therefore `c<2^(1/8)`, and each layout in (3) contributes less than
`2^94`. There are

```
sum_(r=0)^5 binom(16-2r,r)=1+14+66+120+70+6=277<2^9
```

layouts. The complete regional renewal family consequently contributes
less than `2^103`.

## Global scope and the endpoint

The global comparison starts at `H`, uses 255 bulk regions, and ends
with the exact free-endpoint region. In that final region, the cover
cost is at most its 4096 H labels. The adaptive coins form a probability
distribution, so the conditional endpoint value is at most `2^4096`.
There is no additional `2^512` terminal factor.

Retain only the renewal family in each bulk region, but allow the
entire endpoint region. Its contribution to the H-unioned comparison
at threshold 60518 is strictly less than

```
2^8186 * 2^-60518 * (2^103)^255 * 2^4096 = 2^-21971.        (9)
```

Thus this family cannot supply an obstruction to the exact-quota
adaptive majorant. Other population paths remain outside the argument.
In particular, a stable band of partial populations need not return
through two all-hole updates.

The return identities passed two independent mathematical reviews.
The analytic upper bound was derived before any coefficient job ran.
The isolated, untested coefficient helper was then removed. No
production source changed and no numerical result is asserted here.
