# Exact regional routing for the first sparse classes

The unchanged randomized-GF65536 outer and shared beta64 inner now have a
directed bound for **every word with one through five active outer rows**.
This includes quiet and nonquiet histories, all row choices, and all messages.
The failure probability for their complete union is below

    2^-56.51820720499855.

The first unclosed class in the requested one-through-sixteen ladder is
six active rows. Every occupancy from one through sixteen was screened.
The calculations below do not certify six through sixteen or the whole code.
No encoder parameter or setup distribution changed.

## Exact coefficients preserve the regional route

Use the setup and the nonnegative three-component operators P_j(z) from
[the state-density proof](HD_RANDOMIZED_OUTER_STATE_PROFILE.md). The components
are zero-state mass, a density cap relative to the uniform nonzero state,
and arbitrary residual nonzero mass. The current mixer is the actual
invertible rank-one law over q=65536, with four state coordinates.

Fix a set of r active rows. In each region its potential positions form
a uniform r-subset of 8192 positions. Uniform outer comparison averages
the amplitudes and within-packet subsets before any fractional power.
The 128 packets of a region retain their chronological state evolution.

At fractional exponent one, the regional upper operator is exactly

    R_r(z) = [x^r] (sum_(j=0)^64 binom(64,j) P_j(z) x^j)^128
               / binom(8192,r).

This is an exact expectation of the chosen upper operators under regional
placement. It does not claim that the upper operators equal the true process.
In particular, the denominator is the regional subset count; no Bernoulli
marker or conditioning penalty is introduced.

For a fractional exponent 0<alpha<=1, first group b chronological packets.
For every count tuple (j_1,...,j_b), multiply the operators before taking
their entrywise alpha power. Define the block polynomial

    H_(b,alpha)(x) = sum_(j_1,...,j_b)
        product_i binom(64,j_i) * x^(sum_i j_i)
        * (P_(j_1)(z) ... P_(j_b)(z))^(entrywise alpha).

Only tuples whose total is at most the desired r contribute to [x^r].
The corresponding regional operator is

    R_(r,b,alpha)(z) = [x^r] H_(b,alpha)(x)^(128/b) / binom(8192,r).

The implementation uses b in {1,2,4,8,16}. All divide 128. The binomial
placement weights remain outside alpha. Positive matrix-path expansion
and subadditivity justify the entrywise powers, even though the measure
components need not be disjoint probability states.

Put beta=q^256/(q-1)^128 and D=838860. Conditional Markov, the route average,
and the union over active-row sets give the complete r-row failure bound

    B_r = binom(8192,r) * beta^(alpha*r) * z^(-alpha*D)/(q-1)^alpha
              * e_0 R_(r,b,alpha)(z)^256 1.

The row-set union stays outside alpha. Division by q-1 counts projective
nonzero witnesses and stays inside alpha. The power 256 preserves state
between regions. Summing all three terminal components leaves the final
state unconstrained.

## Averaging a sparse return before taking a state supremum

The earlier residual-to-zero bound charges the entire nonuniform mixer
mass whenever its incoming state might lie in the feedback image. Averaging
the fresh potential-slot subset first gives a smaller valid bound.

Let C have 64 columns, with every four independent. Fix a nonzero state h.
For a uniform j-subset U of columns, define

    p_j(h) = Pr_U[h belongs to span(C_U)].

For j=1,2,3,

    p_j(h) <= j/64.

To prove this, expose an ordered sample without replacement. After j-1
columns, suppose h is not in their span. At most one remaining column
can extend that span to contain h. Two such columns would make their
j+1<=4-column union dependent. If h is already in the old span, every
extension succeeds. Consequently

    p_j(h) <= p_(j-1)(h) + (1-p_(j-1)(h))/(65-j),       p_0(h)=0.

Induction gives j/64. The bound is attained when h is a column of C.
No claim of this form is used at j=4: every four columns span the state space.

On nonzero states, write the actual rank-one law as

    M = epsilon * (a mixture of I and scalar-line averaging)
        + c * (uniform nonzero refresh),
    epsilon=1/q,       c=1-1/q,

where the equality denotes equality of transition kernels. The first part
preserves the scalar line of h. The uniform part hits span(C_U) with
probability (q^j-1)/(q^4-1). A uniform j-coordinate input then matches that
state with probability q^-j. Thus the return probability, averaged over U,
is at most

    q^-j [epsilon*j/64 + c*(q^j-1)/(q^4-1)],       1<=j<=3.

This statement remains uniform over every incoming nonzero h. The fresh
within-packet subset is independent of the earlier weighted state measure
after conditioning on packet counts, so the bound applies to the residual
component of the state measure.

At most j FFT8 tiles contain potential inputs. Their complement emits A h
without input cancellation and has at least max(0,61-8j) nonzero symbols.
For the output mark z, the joint emission-and-return bound is therefore

    kappa_j <= z^max(0,61-8j) q^-j
                 [epsilon*j/64 + c*(q^j-1)/(q^4-1)].

The emission factor is pointwise, so this multiplication requires no
independence between return and output weight. The helper takes the smaller
of this bound and the previously proved complete bounds for kappa_j.

For every j>0 it also uses a weaker general refinement. Average
z^max(0,61-8u) over potential subsets, where u is their occupied-tile count,
and multiply by the old unweighted return bound. That return bound is
uniform for each fixed subset. The j=0 return remains exactly zero.

## Directed endpoints and their union

The following witnesses use the membership refinement. The weight mark is
the exact positive real z=2^ell, where ell is the displayed rational number.
Arb evaluates z and all later operations with enclosures.

| Active rows r | ell=log2 z | alpha | Packets per block b | log2 failure upper bound |
| ---: | ---: | ---: | ---: | ---: |
| 1 | -1/512 | 1 | 1 | -391.9847734422 |
| 2 | -1/256 | 1 | 1 | -777.5398361082 |
| 3 | -3/512 | 11/16 | 4 | -535.3231446326 |
| 4 | -1/128 | 5/16 | 4 | -213.8444812453 |
| 5 | -5/512 | 1/4 | 4 | -56.5182072049 |

Each row includes every choice of its active-row set. Summing the five
directed endpoints gives

    log2(sum_(r=1)^5 B_r)
      <= -56.518207204998555118899149411951304647522802240880485.

The acceptance test compares the enclosing Arb value against -40.
Rounded table entries are not used for acceptance. As a separate regression,
the same coefficient engine with the unmodified local operator certifies
one through five rows with a weaker union below 2^-51.5378269956.

## What fails after five rows

The first screen used alpha=1 and exact regional coefficients for all
integer occupancies one through sixteen. One and two rows closed, but
three rows gave only about -2.75 bits. At four rows the expression exceeded
one by over a thousand bits. These are failed upper bounds, not construction
counterexamples.

Fractional four-packet blocks move the frontier to five rows. The following
floating diagnostics use the membership refinement and a bounded grid;
only the successful first five were replayed into the certificate above.

| Active rows | Best four-packet expression in the grid |
| ---: | ---: |
| 6 | +84.17 |
| 7 | +191.85 |
| 8 | +314.60 |
| 9 | +439.50 |
| 10 | +544.44 |
| 11 | +623.83 |
| 12 | +701.55 |
| 13 | +778.59 |
| 14 | +856.09 |
| 15 | +933.73 |
| 16 | +1011.31 |

A separate bounded refinement at r=6 used sixteen-packet blocks. Its best
floating expression was +47.701781, at ell=-47/4096 and alpha=13/64.
The new flat-state and coupled-column moments changed this value negligibly
near that weak tilt. Larger blocks improve the bound, but these tests do not
identify a closed witness for r=6.

The observed failure is specific: after charging the r-row outer count and
the row-set union, the current route-conditioned weighted-history envelope
does not have forty bits of slack. Fractional composition substantially
reduces the first-moment loss. The experiments do not prove which remaining
inequality is dominant or whether the actual code has a low-distance word.

The next useful step is to retain more information about successive births
and returns across sparse input packets, then test that improvement at r=6
before extending a larger numerical grid. Increasing the grid alone has
not closed the requested one-through-sixteen range.

## Implementation and verification

The new files are
[theory/hd_randomized_outer_sparse_ladder.py](theory/hd_randomized_outer_sparse_ladder.py)
and [its tests](theory/hd_randomized_outer_sparse_ladder_tests.py).
The floating path uses positive matrix-polynomial convolutions to select
witnesses. The directed path reconstructs the local operators independently
and uses positive Arb products throughout. It truncates only degrees above
the target occupancy, which cannot affect a lower coefficient.

Six tests pass. They compare every local occupancy through sixteen with the
directed reconstruction, check noncommuting coefficient products, recover
the existing one-row endpoint, and verify regrouping at alpha=1. They also
check every small-field span in an eight-column GF8 MDS example, all potential
subset normalizations, and both complete one-through-five unions.

    python -B -m unittest hd_randomized_outer_sparse_ladder_tests -v
    python -B hd_randomized_outer_sparse_ladder.py --receipt
    python -B hd_randomized_outer_sparse_ladder.py --provider base
    python -B hd_randomized_outer_sparse_ladder.py --steps 4 --fractional

No previous proof file, encoder implementation, benchmark, remote host,
or commit was changed for this calculation.
