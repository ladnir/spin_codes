# Binary two-sided regular expanders

This note records the proof route, the initial rate-one-half diagnostics, and
the resulting interval certificates.  Exploratory calculations use
`scripts/binary_biregular_diagnostic.py`; the final claims use the separate
outward-rounded verifiers named below.

## Exact regional law

Let a region contain `ell` right coordinates of degree `d_R`, so
`k = d_R ell`.  For a fixed message support of size `r`, a random regional
permutation makes its slots a uniform `r`-subset of the `k` slots.  Define

```text
E_d(X) = ((1+X)^d + (1-X)^d)/2,
O_d(X) = ((1+X)^d - (1-X)^d)/2.
```

If `Q_0(Y)` and `Q_1(Y)` are the recursive-map transfers for binary inputs
zero and one, respectively, then the exact regional transfer is

```text
R_{ell,d,r}(Y)
  = binom(k,r)^(-1) [X^r]
    (E_d(X) Q_0(Y) + O_d(X) Q_1(Y))^ell.
```

Independent regional permutations compose by taking the `d_L`-th matrix
power.  This formula applies to the accumulator and the wrapped convolution.

For small `r`, the diagnostic first computes the exact regional parity-shell
law

```text
Pr[U=u]
  = binom(ell,u) [X^r] E_d(X)^(ell-u) O_d(X)^u / binom(k,r).
```

Conditional on `U=u`, the regional input is uniform among weight-`u` binary
vectors.  The existing uniform-slice transfer then gives the exact output
moment.

## Parity restriction

For the all-one message, every right coordinate receives exactly `d_R`
ones.  Therefore

```text
1 B = (d_R mod 2) 1.
```

An even right degree puts the all-one message in the kernel.  At rate one
half, `d_L=2d_R`, so a valid binary profile requires `d_L = 2 mod 4`.

## EA result

Choose

```text
k = 1,048,575
d_L/d_R = 62/31
ell = 33,825
n = 2,097,150
L = 230,729
```

The exact weight-two contribution is `2^-18.244216`.  For comparison, the
left-regular calculation at the same degree gives `2^-18.243922`.  Fixing the
right degree gains only `0.000294` bits.  Thus two-sided regularity does not
reduce the EA degree from 64 to 62 at this cutoff and failure target.

## EC results

The first conservative candidate was `d_L/d_R=26/13` with memory 9.  Its
exact contributions for supports one through eight range from `2^-106.14`
to `2^-318.38`.  Sparse supports therefore have much more margin than the
existing left-regular degree-28 certificate suggests.

The more useful candidate is

```text
k = 1,048,575
d_L/d_R = 10/5
ell = 209,715
n = 2,097,150
L = 230,729
memory = 15
```

Every exact contribution through support 64 is below `2^-50`; support one is
the largest of those terms.  Sampled coefficient-saddle terms outside the
central band are also strongly negative.

Degree `6/3` with memory 60 does not close by the same exact method.  Its
support-64 term is approximately `2^1.94`.  Memory 79 repairs this term.  An
outward-rounded verifier covers every nonzero support and bounds the total by
`3.558444e-7`, or more than `21.4222` failure bits.  It uses exact transfer
blocks through support 383, positive-coefficient blocks on the outer ranges,
and degree-three local point-mass blocks in the central range.  Because some
exact transfer entries are smaller than the ordinary binary64 range, the
verifier stores entries in separate power-of-two exponent layers.

Two intermediate profiles fill the degree--memory curve.  Both use
`k=1,048,635`, `n=2,097,270`, and cutoff `230,742`:

```text
degrees   memory   total bound       failure bits   d_L + memory
14/7      6        1.402384e-8       26.0875        20
18/9      4        2.987318e-7       21.6746        22
```

For any odd right degree, each parity-conditioned group enumerator has only
negative real roots.  It is therefore the probability-generating polynomial
of a sum of independent Bernoulli variables.  Summing endpoint lower bounds
on their variances extends the central local-limit proof to these profiles.

The next smaller memories do not close under the same sparse-support bound.
For degrees `14/7` at memory `5`, support two is bounded only by `2^-17.47`.
For degrees `18/9` at memory `3`, support two is bounded only by `2^-9.06`.
These diagnostics identify the current frontier but do not prove that another
output-tail argument cannot improve it.

## Dense-band slack

The generic positive coefficient saddle drops the fixed-support condition in
each region.  Near support `k/2`, this costs approximately one square-root
factor per region.  It reports `2^44.87` for the degree-10 candidate even
though the central parity distribution is almost uniform.

There is an exact check at the two central supports.  Write `k=2R+1` and fix
support `r=R`.  For a Fourier character involving `a` of the `k` slots, its
normalized coefficient has magnitude

```text
|K_R(a;2R+1)| / binom(2R+1,R)
  = binom(R,floor(a/2)) / binom(2R+1,a).
```

Only characters whose slot sets are unions of right groups occur.  Summing
their absolute Fourier coefficients gives `2` to numerical precision for
the candidate parameters.  Hence one regional point mass is at most
`2^(1-ell)`.  Invertibility of the convolution then gives a first-moment term
of approximately `2^-58.35` at each central support for degree 10.  The
generic saddle's positive value is therefore proof slack.

The paper now proves a uniform replacement for this central identity when
`d_R=5`.  Tilt the uniform subset by selecting every slot independently with
probability `p=r/k`, then condition on the total count.  Conditional on one
group's parity, half its selected-slot count is a sum of two independent
Bernoulli variables.  A Fourier integral bounds its largest point mass.

For the degree-10 candidate, the resulting bound over
`k/5 <= r <= 4k/5` has largest individual support term `2^-51.8319`.
A monotone endpoint partition gives a sharper bound on the sum.

The outward-rounded verifier now covers every nonzero support.  Its total is
`1.6474999656e-10`, or `32.4990025` failure bits.  The contributions are:

```text
exact supports 1..16       9.4364693e-16
111 outer saddle blocks    1.6336666e-10
47 central local blocks    1.3823976e-12
full support               < 1e-341856
```

The optimizer-selected markers and all block endpoints are frozen in
`results/binary_biregular_ec_rate_half_d10_m15_gv.json`.  The standalone
standard-library checker validates the schema, parameter identities, marker
domains, and complete support partition before the Arb verifier evaluates the
stored choices.

## Recommended next step

Replace the simple `d_L + memory` proxy with measured encoding costs.  Also
determine whether a sharper exact-support partition can close memory 78 for
the degree-6 profile.
