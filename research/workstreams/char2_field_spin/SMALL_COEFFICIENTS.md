# A certified small-coefficient starting point

We can close a first characteristic-two family without full-width field
products. The starting construction uses a genuine GRS outer, two uniform
symbol permutations, and two input-weighted accumulators. It is a baseline
for the new family, not yet the intended one-permutation batched SPIN design.

At rate 1/2, the inexpensive GRS[8,4] outer already certifies 10% distance:

| Message symbols K | Output symbols N | Minimum distance at least | Setup-failure margin |
|---:|---:|---:|---:|
| 65,536 | 131,072 | 13,108 | 41.290077942075 bits |
| 262,144 | 524,288 | 52,429 | 45.298335145044 bits |
| 1,048,576 | 2,097,152 | 209,716 | 49.300334745625 bits |

These are complete finite all-message bounds, not sampled distances.
They concern the ideal independent setup distribution defined below.
Matched implementation measurements are recorded separately in
[RESULTS.md](RESULTS.md).

## Field and coefficients

Fix F = GF(2^m), with m in {8,32,64,128}. Fix an irreducible polynomial
P_m of degree m and represent F as GF(2)[theta]/P_m(theta). The coefficient
alphabet is

    S = {the polynomial bit patterns 1,...,128}.

These are 128 distinct nonzero elements because their degrees are below m.
At m > 8, S is not contained in the GF(256) subfield: it contains theta, which generates
F over GF(2). Multiplication by an element of S is multiplication by a
polynomial of degree at most seven, followed by field reduction. This
permits a small-coefficient kernel; it does not require a general
full-width-by-full-width field multiplier.

The proof uses only the size and nonzero property of S. Other sets of 128
distinct nonzero elements have the same guarantee. The arithmetic costs,
however, depend on the selected representation and coefficient set.

Fix the field before sampling setup. The guarantee holds separately for
each requested m. Reusing one bit-pattern setup across all four native
polynomial fields is not scalar extension: their defining polynomials and
the meanings of theta differ. A union bound would lose at most two margin
bits if one requires that setup to succeed simultaneously in all four.

## A small GRS outer

Partition the message into blocks of four field symbols. For a block
(u_0,u_1,u_2,u_3), define

    S_1(X) = X^2 + X,
    f(X) = u_0 + u_1 X + u_2 S_1(X) + u_3 X S_1(X).

Evaluate f at the eight distinct field elements whose polynomial bit
patterns are 0,...,7. The basis polynomials have respective degrees
0,1,2,3 and leading coefficient one. Thus they form a basis of all
polynomials of degree at most three. A nonzero f has at most three roots,
so this evaluation code is [8,4,5] MDS.

There is no interpolation or division in this encoder. To evaluate a pair
of points v and v+1, use S_1(v+1)=S_1(v) and compute

    a = u_0 + S_1(v) u_2,
    b = u_1 + S_1(v) u_3,
    f(v) = a + v b,
    f(v+1) = f(v) + b.

For v in {0,theta,theta^2,theta^2+theta}, every displayed constant has
polynomial degree at most four. This gives a concrete implementation
without dense full-field constants. Further common-subexpression or
additive-FFT optimization must preserve this field map or explicitly
declare the corresponding change of message basis.

## Two scans with multiplication outside the recurrence

Setup samples two independent uniform permutations pi_1 and pi_2 of the
N symbol coordinates. Independently, sample c_(j,i) from S for every
round j in {1,2} and position i. Setup is public and fixed for all messages.

Starting from the concatenated outer output, round j computes

    x_i = v_(pi_j(i)),
    y_(-1) = 0,
    y_i = c_(j,i) x_i + y_(i-1).

The final output is the complete second-round vector. Both scans are
invertible because every c_(j,i) is nonzero. Therefore the composed
encoder has dimension K for every setup and rate 1/2.

Writing A for the ordinary prefix-sum matrix and D_c for the diagonal
coefficient matrix, one scan is A D_c. Its field transpose is D_c A^T:
compute suffix sums, then multiply each result by its local small
coefficient. Permutations are transposed in the usual reverse order.
The small products can be batched outside the prefix/suffix dependency
chain. No long product of sampled coefficients is needed.

## Why the companion proof transfers

The companion's generic bound for feedback-weighted scans also applies
to input-weighted scans. Fix a source support and quotient the message
space by all output-zero equations imposed so far. At a return to zero,
let ell be the required nonzero input form and r the required nonzero
state form in this quotient. The new zero equation is

    c ell + r = 0.

At most one c in S can make this equation identically zero in the quotient.
Indeed, two distinct such coefficients would imply ell=0. Every other
coefficient either adds an independent zero equation or kills the branch.
For a rank weight tau >= 1, each return therefore multiplies expected
branch potential by at most

    g = 1/128 + (127/128)/tau.

The argument conditions only on already exposed coefficients and labels.
The next coefficient remains independent. It uses no ordering of the
field, no odd-characteristic opposite-pair argument, and no assumption
that separate zero equations are independent.

The inverse scan satisfies x_i = c_i^(-1)(y_i-y_(i-1)), so an output of
weight h has input weight at most 2h. Hence final weight at most e restricts
intermediate weight to at most 2e and source weight to at most 4e.

For a prescribed support of j outputs in a local [8,4] MDS block, the
shortened message space has dimension max(0,j-4). Its dimension-weighted
support polynomial is exactly

    Q(z,tau) = 1 + 56 z^5 tau + 28 z^6 tau^2 + 8 z^7 tau^3 + z^8 tau^4.

These are the same source dimensions and return factors used by the
companion's four-input certificate. Consequently its whole-run support
count, sparse return moment, and rational bulk endpoint bound transfer
without changing their numerical inequalities. The integer Newton
implementation and its addition counts do not transfer.

For completeness, let A_g(s,h) be the coefficient of x^s y^h in

    [1,0] [[1,xy],[g*x,(1+x)*y]]^N [1,1]^T,

and H_g(h,e) = sum_(j=1)^e A_g(h,j). The resulting all-message bound is

    Pr[d_min <= e]
      <= (1/tau) sum_(s=1)^(4e) ([z^s] Q(z,tau)^(N/8) / binom(N,s))
           * sum_(h=1)^(2e) A_g(s,h) H_g(h,e) / binom(N,h).

The first permutation gives the binom(N,s) normalization. The second
gives binom(N,h). Rank potential is retained across both scans. Different
source ranges may use different tau values. These normalizations require
the two independent uniform permutations; they do not justify replacing
either permutation with SPIN region routing without a new argument.

## Exact finite evaluation

[theory/return_certificate.py](theory/return_certificate.py) adapts the
companion's verified return bound. It does not modify that checkout.
The sparse sum covers every source weight through 128. The bulk intervals
cover every integer source weight from 129 through 4e. Each bulk interval
is bounded by 2^-80 using rational endpoint inequalities and outward
logarithm bounds. Floating optimization proposes witnesses only.

All three rows above use tau=512 for the sparse sum and e=floor(N/10).
Their bulk covers contain 78, 96, and 114 intervals, respectively.
Generation ends by replaying the complete certificate. A later independent
invocation of --verify does not optimize or trust the recorded bounds.

Example:

```sh
python -B research/workstreams/char2_field_spin/theory/return_certificate.py \
  --generate --log-k 16 --alphabet 128 --distance 0.10 --tau 512 --margin 40 \
  --output tmp/char2-field-spin/k16-grs8-q128-d10.json

python -B research/workstreams/char2_field_spin/theory/return_certificate.py \
  --verify tmp/char2-field-spin/k16-grs8-q128-d10.json
```

The initial witnesses are in the ignored tmp directory. Their SHA-256 values
are recorded here to identify this run:

| K | Witness filename under tmp/char2-field-spin | SHA-256 |
|---:|---|---|
| 2^16 | k16-grs8-q128-d10.json | 0333bd06af9961c82d86a09933e0cc9c6e15440ed29d0de8d2f3f431e7b9ad00 |
| 2^18 | k18-grs8-q128-d10.json | a5bc5aa4e547b4b950270addc80df7cfd2c22e82cdf3b41b3be83e8cd0ec26a6 |
| 2^20 | k20-grs8-q128-d10.json | 9d4e224ca41e611374b4c78824c276ee023fb6c879b2ea73256a48f7fe7b9376 |

The report also records hashes of the companion verifier dependencies.
The companion path can be supplied explicitly with --companion. These
research scripts currently depend on that adjacent checkout; they are not
yet a self-contained released proof artifact.

## A stronger outer already has a transferred certificate

The companion's generic direct certificate for a [16,8] MDS outer was also
replayed on 2026-10-03. With K=2^20, the same two input-weighted scans, and
the same 128-element alphabet, it gives

    d_min >= 639632,
    Pr[failure] <= 2^(-58.865325164206).

The relative distance exceeds 30.5%. This transfer uses only the MDS
shortening polynomial and the generic rank continuation operator. It
does not use the companion's later odd-characteristic pair refinements.
There are 159 sparse intervals and 158 dense leaves.

Replay it with [theory/direct_certificate.py](theory/direct_certificate.py):

```sh
python -B research/workstreams/char2_field_spin/theory/direct_certificate.py \
  --verify ../gen-BAA/field-universality/expander-lift/results/convolution-width-one-direct-r8-q128-e639631.json
```

The source witness hash is
`4ba36d670c74fc9a7c3100e057a3f7bf7266417c51b278d33dca52e2fe6f4404`.
Use a genuine GRS[16,8] map in characteristic two, not the companion's
integer Newton circuit. The implementation report measures both outer choices.

## Next steps

The GRS[8,4] and GRS[16,8] implementations now pass scalar-reference,
field-adjoint, SIMD, aliasing, and sanitizer checks. Complete precomputed
timings expose the two global permutations as the main large-size cost.
Keep these proved controls while developing the separate
[batched region-route construction](BATCHED_PROOF.md). Higher-rate variants
still require a new source polynomial and finite certificate.
