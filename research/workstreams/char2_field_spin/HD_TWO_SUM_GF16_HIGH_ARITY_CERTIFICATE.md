# Directed high-total-arity certificates for three g16 source weights

The three classes below have a combined failure-probability bound of
$2^{-62.24877612169}$. Fresh 256-bit and 384-bit replays agree on all 555
ordered matching-profile ceilings, all three class ceilings, and the aggregate.
A separate 512-bit sum of the rounded class ceilings gives the same aggregate.
This is a partial certificate: the low-total-arity complements and all other
source weights are excluded.

## Construction and event

Work over $\mathbb F_{256}$ with AES polynomial `0x11b`. Each of the 8192
source groups contains sixteen RS$[16,8]$ constituent rows. Its column mixer
is $H=I+\mathbf1\mathbf1^T+uu^T$, where the ordered labels are

```
06 07 0a 0b 56 57 5a 5b b6 b7 ba bb e6 e7 ea eb
```

These are the coset $06+\mathbb F_{16}$ in the established byte basis. The
ideal outer setup uses independent nonzero input and output diagonal entries.
The gauge-normalized implementation is interpreted through its established
lift to that ideal distribution. Routing uses the alternating adjacent-pair
matchings, independent uniform pair permutations in all 256 regions, and
independent fair within-pair swaps. The h4 inner construction is unchanged.
The ensemble identifier is
`two-sum-gf16-coset-g16-rs16x8-h4-swapped-pairs-k20-rate-half`.

For a scalar codeword, let $G$ be its active source-group set, let $r=|G|$,
and let $A$ count its nonzero constituent RS rows across $G$. The certified
event is existence of a nonzero scalar codeword of output weight at most
$L=734003$ in one of the three classes $(r,A\ge A_0)$ listed below. Probability
is over the ideal randomized encoder setup. Every indexed set $G$ of the
specified size, every even/odd source balance, and every constituent-row-set
assignment satisfying the threshold are included.

This restricted scalar claim is not a claim about a record codeword's own
active-group or arity class: a scalar component can have smaller support.
The usual scalar-extension equality of minimum distance becomes sufficient
for a whole record-code certificate only after all scalar classes are covered.

## Uniform comparison without a row-set union factor

Fix $G$ and a target outer array $y$ on its $256r$ potential output coordinates.
Condition on all outer output diagonals $D_{\mathrm{out}}$. In each column,
$H^{-1}D_{\mathrm{out}}^{-1}y$ determines one pre-input-diagonal array $v$.
Because input diagonal entries are nonzero, this array determines one
nonzero constituent-row set and hence one arity $A(v)$. The row sets are not
resampled across columns.

For a prescribed nonzero RS$[16,8]$ row of support size $w$, its shortened
code has dimension $w-8$ when $w\ge9$. An information set therefore bounds
the number of full-support words by $255^{w-8}$. Independent nonzero input
diagonal entries make the prescribed scaled row uniform on its $255^w$
amplitude assignments. Its expected membership is at most $255^{-8}$;
an impossible support contributes zero and a zero row contributes one.
Independence across rows gives, for the expected outer counting measure,

$$
\mu_{G,A\ge A_0}(y\mid D_{\mathrm{out}})
\le \mathbf1_{A(v)\ge A_0}\,255^{-8A(v)}
\le255^{-8A_0}.
$$

Averaging the conditioned output diagonals preserves this bound. There is
no row-set, cohort, or deficit-binomial multiplier: at each fixed target and
conditioned output diagonal, only one row-set assignment can contribute.
Relative to uniform field amplitudes on all $256r$ potential coordinates,
the density factor is

$$\beta_{r,A_0}=256^{256r}/255^{8A_0}.$$

This comparison allows zero amplitudes. It therefore uses the established
**uniform-amplitude** swapped-pair operator, not the forced-nonzero operator.

## Reused operator and new scalar factor

Let $R_{r,f}(z,\alpha)$ be the directed uniform-amplitude regional matrix
returned by `regional_directed` in `hd_shifted_swap_sparse.py`. Its local
packet caps, good-good composition before the entrywise fractional power,
fair-swap averaging, positive coefficient extraction, and unpowered
matching-histogram normalization are unchanged. The adapter rebuilds this
matrix at the requested precision; it reads no saved numerical receipt.

Let $N_r(a,b)$ denote the exact-integer upper bound `profile_count_upper`
on indexed source supports with $a$ full pairs in the even matching and
$b$ in the odd matching. With $e_0^T=(1,0)$ and $\mathbf1=(1,1)^T$, the bound is

$$
\left(\frac{\beta_{r,A_0}z^{-L}}{255}\right)^\alpha
\sum_{a,b}N_r(a,b)\,
e_0^T\bigl(R_{r,a}R_{r,b}\bigr)^{128}\mathbf1.
$$

The initial state and regional order are retained. The source-support union
is outside $\alpha$ and is already counted by $N_r$; no additional
$\binom{8192}{r}$ factor is multiplied. Projective normalization occurs once.
The new logarithmic scalar factor is evaluated directly as

$$
\alpha\bigl(2048r-8A_0\log_2 255-\log_2 255-L\log_2z\bigr).
$$

As an algebraic control, this is the legacy inherited-diagonal factor plus
$8\alpha(16r-A_0)\log_2 255$. The sharper invertible-Toeplitz comparison is
not used or transferred.

## Fresh directed results

All three witnesses use $\alpha=1/5$ and upward-rounded 60-bit packet caps.

| $r$ | Included $A$ | $z$ | Ordered profiles | Exact logarithmic upper ceiling |
|---:|---:|---:|---:|---:|
| 18 | 247 through 288 | $691/700$ | 99 | $-15985512996797/250000000000$ |
| 24 | 317 through 384 | $172/175$ | 168 | $-64512660498557/1000000000000$ |
| 32 | 463 through 512 | $171/175$ | 288 | $-31649989891019/500000000000$ |

The raw aggregate and the 512-bit rounded-class aggregate both have ceiling
$-6224877612169/100000000000$. The complete replay arithmetic took 2.543
seconds at 256 bits and 2.670 seconds at 384 bits. Each anchor had a
cooperative 60-second computation guard; no guard was approached.

The eight focused tests pass in 0.295 seconds. They include exact GF4
enumeration of the unique-row-set comparison, the legacy-prefix identity,
ordered noncommuting matrix/profile controls, fresh-provider precision checks,
scope failures, and independently rounded dual-receipt verification.

Files and SHA-256 identifiers:

- `theory/hd_two_sum_gf16_high_arity_certificate.py`:
  `9c9c2e4d1bfc95f2952d6bc743b9047489128d5cc26a43f5186ea721b1da4027`
- `theory/hd_two_sum_gf16_high_arity_certificate_tests.py`:
  `b7885e7d68f24c9725e4c768ce1fad1784fd2618e34fb19bfc0b568f74c65eef`
- `theory/hd_two_sum_gf16_high_arity_receipt256.json`:
  `195c0bfe9637d44d8f57cf35011d5897bbab8d4c7c3d7b27d709e3c571256605`
- `theory/hd_two_sum_gf16_high_arity_receipt384.json`:
  `f8ff61fb61253faaf9f1eff1378e6185db8d608a046254a48e15295b71b26248`
- `theory/hd_two_sum_gf16_high_arity_dual_verified.json`:
  `4bee3257cafb59ff44d6ffbdb76e7ef7d890e2d449a41e78a9263347e6d932ed`

The receipts pin the reused dependencies. Hashes preserve provenance; the
numerical evidence is the observed fresh directed arithmetic and its
independent rounded union, not the hashes alone.

## Remaining proof task

The complementary classes are $A\le246$, $A\le316$, and $A\le462$ at the
three respective source weights. The separately recorded bounded low-arity
FLOAT control did not close them. The next useful step is to strengthen that
complementary fixed-type operator or its support-floor constraints while
preserving these high-side receipts. These three anchors are not a continuous
weight band, a complete weight-class cover, or a whole-code certificate.
