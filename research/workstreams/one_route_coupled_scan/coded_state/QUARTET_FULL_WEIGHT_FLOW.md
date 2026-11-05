# From the certified band to all state weights

2026-10-27. The construction is unchanged. The goal remains an end-to-end
45% distance bound with at least 40 bits of margin. This round develops
positive transfer bounds for all state weights and all turnover counts.
The numerical results below are floating-point screens, not certificates.

## What retaining weight flow resolves

The previous scalar band bound could maximize over a different incoming
weight at every update. `quartet_weight_flow.py` instead retains the zero
state and each weight from 216 through 232. Transitions between live states
still have 24--40 deaths; zero runs have arbitrary lengths.

Both comparisons use the same coefficient tables and proof-p grid:
0.42, 0.54, 0.60, and 112627/167924. At the proposed tau=0.092430:

| Inactive rows r | Retained weight flow | Scalar envelope |
| --- | ---: | ---: |
| 1850 | -1083.72 | +2409.49 |
| 1974 | +1067.31 | +4559.33 |
| 2100 | +7226.80 | +12264.14 |

Entries are base-two logarithms of probability upper bounds. Positive
entries provide no nontrivial probability guarantee. Retaining weight flow
helps, but does not close even this restricted family at every tested r.
These three row counts do not certify the intervening interval.

A marked-coefficient calculation also measures active input coverage under
the artificial one-update coefficient distribution. At p=112627/167924,
the conditional means for (incoming weight, outgoing weight, deaths) equal
321.30 for (216,216,31), 330.56 for (224,224,32), and 345.14 for (232,232,40).
These are not means under the exact regional quota or the actual encoder.
They suggest that insufficient aggregate coverage is not the main loss
inside this dense band. They say nothing comparable about sparse states.

## A full one-update kernel

Use the construction and conditional message count from
[the band proof](QUARTET_HISTORY_BAND.md). Fix an inactive-row set of size r.
For one quartet let U and V be its incoming and outgoing supports, and let
T mark positions assigned to inactive rows. Write a=|U|, v=|V|,
d=|U\V|, and t=|T|. The existing feasibility checks give a selected-zero
polynomial P(ell). Let c count covered inactive positions.

Introduce a coverage marker gamma>=1. The number of covered active input
positions in a feasible quartet is

    c_active = 4-t if a>0, and v if a=0.

When a=0, feasibility forces V and T to be disjoint. The local coefficient
for incoming count a and outgoing count v is the sum, over V and T, of

    binom(4,a) p^t (1-p)^(4-t) [P(ell)/q^c]^theta
      * X^(v-a) ell^(-theta*a) (q*ell)^(-theta*d)
      * gamma^c_active.

Here q=2^32, theta=7/40, and q*ell>=1. The potential is
f(w)=X^w*ell^(theta*w). Multiplying 128 copies of this polynomial and
dividing incoming row u by binom(512,u) gives a 513-by-513 kernel.
It sums every outgoing weight and every death count, without a band cutoff.
Uniform quartet partitions justify the binomial denominator. The row sum
depends only on incoming weight, so the kernels compose over shared states.

The helper `quartet_full_weight_kernel.py` computes this kernel in the
logarithmic domain. A build takes approximately 6--7 seconds here. Boundary
rows agree with independent activation and reset formulas. Its 289 entries
on weights 216--232 dominate the separately computed turnover-restricted
kernel, as required.

## A necessary outer constraint

An active RS[256,128] row has at least 129 nonzero coordinates. Every such
coordinate must be covered by the input support allowed by the history.
Thus a realizable history satisfies

    C_active >= 129*(8192-r).

After the conditional message bound and its theta moment, multiply by
the feasibility indicator. For gamma>=1, that indicator is at most
gamma^[C_active-129*(8192-r)]. This gives the local marker above and a
global factor gamma^[-129*(8192-r)]. The marker is outside the theta power.
It is a necessary aggregate condition, not the stronger row-by-row test.

For example, repeated transitions 0 -> weight one -> 0 cover at most five
input coordinates per two updates. Over 4096 updates they cover at most
10240 coordinates. At r=2100, active rows require at least 785868 covered
coordinates. Such histories cannot represent messages in this row class.
This example identifies a genuine overcount, not a low-distance codeword.

## Regional composition and proof parameters

Grade the 0-to-0 edge as quiet and every other edge as nonquiet. After
16 updates, apply the existing regional cap kappa_h(p,r) to degree h.
Let Gamma_p denote this positive regional operator. It bounds the same
fixed-quota route for every artificial Bernoulli parameter p.

For any positive vector h, define F(h) coordinatewise by min_p Gamma_p(h).
Each alternative bounds the same exact regional continuation operator.
Consequently, if F(h)<=lambda*h, monotonicity gives a lambda^256 bound,
with the initial coordinate h(0) and the maximum terminal ratio included.
The terminal vector is 1/f(w), and the common global logarithmic factor is

    log2 binom(8192,r) + theta
      + 32*theta*(W-K+128*r) + theta*W*log2(ell)
      - 129*(8192-r)*log2(gamma).

All alternatives use the same theta, ell, X, and gamma. Choosing p in the
proof does not change or resample the route. Six potential iterations
produce candidate upper bounds; they are not a convergence certificate.

The helper also groups weights into bins and retains both regional endpoint
bins and h before minimizing over p. Each grouped transition maximizes over
every incoming weight after summing the outgoing bin. It never substitutes
a representative weight. This remains an upper envelope, but loses
within-bin flow. Choosing X=ell^(-theta) makes f constant and removes
artificial variation from the diagonal potential inside a bin.

## Full-weight screens and what they imply

At r=2100, the all-weight positive-potential bound remains far from closing:

| ell | gamma | Proof-p choices | Best log2 upper bound in six iterations |
| --- | ---: | --- | ---: |
| 13592/664659 | 1.01 | 112627/167924, 2100/8192 | +1029971.07 |
| 0.00377 | 1.5 | same | +1224530.31 |
| 0.00377 | 1 | same | +685155.11 |
| 0.00377 | 1 | same, plus 0.94 | +685155.11 |

The coverage-marker comparison is not an optimization over gamma.
In particular, the attempted gamma=1.5 does not help at these other
parameters. Sparse impossible paths are one source of loss, not an
established explanation for the largest eigenvalue or final bound.

A separate constant-weight screen isolates another issue. With ell=0.00377
and gamma=1, minimizing the complete fixed-weight bound over the first two
p values gives its largest value at weight 489: +483685.50. Adding p=0.94
moves the largest value to weight 103: +148438.89. Weight 224 is already
negative at -9791.04. All these histories allow unrestricted turnover.

Nearly full states cover almost every inactive-row coordinate. The proof
parameter inherited from the 224/32 band handles this regime poorly.
The extra p value improves that restricted screen, but not the full
positive-potential bound. The unresolved loss includes how different
weight regimes are summed within a region; constant-weight checks alone
do not close the full transfer.

Using bins of 16 with the original X gives +827219.13. This is also open.
The constant-potential choice gives +1042598.82, so it is not an improvement.
No positive screen here establishes a structural obstruction in the code.

## Reproduction and remaining work

Run the source helpers serially, with Python -B and one BLAS/OMP thread.
They have a 60-second watchdog, require at least 20 GiB free, and write no
generated files. An eight-p weight-flow collection exceeded the watchdog;
the recorded matched comparison uses four p values instead.

Useful commands, relative to this directory:

    python -B quartet_weight_flow.py --optimize-tau
    python -B quartet_full_weight_kernel.py --tests-only
    python -B quartet_full_weight_kernel.py --compare-band
    python -B quartet_full_weight_kernel.py --regional --adaptive-p --extra-p 0.94 --ell 0.00377
    python -B quartet_full_weight_kernel.py --coarse-regional 16 --adaptive-p --extra-p 0.94 --ell 0.00377 --neutral-potential

Independent checks cover nine exact small-power cells, 117 explicitly
enumerated short paths, 81 potential-conjugacy cells, and 81 coverage
monotonicity cells. These check the implementation identities, not directed
rounding of the large bound.

Next retain regional coverage and inactive occupancy jointly with weight
flow. In particular, avoid charging all 512 positions of a sparse update
as though their inactive labels were constrained. Then reintroduce the
rowwise RS feasibility test. Other inactive-row counts and directed
rounding remain necessary before whole-code closure.
