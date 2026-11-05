# Joint regional occupancy and the remaining state bound

2026-10-27. Follow-up to [long quiet intervals](QUARTET_QUIET_REGIONS.md).
The whole-code 45%/40-bit goal remains active. Construction and old directed
certificates are unchanged. Numerical results here are FLOAT diagnostics.

## Separate bounds before summing history classes

The previous regional calculation used one interpolation parameter tau for
all histories. Now retain the total number h of nonquiet updates through
all 256 regions. For each h, choose the smallest complete bound among the
tested tau values, then sum h. This does not select proof weights after
observing setup: it is a partition of the candidate history class.

For the unchanged live-weight band 216--232, deaths 24--40 between live
states, and unrestricted zero-state runs, the three screened row counts give:

| Exact inactive-row count | Best FLOAT log2 bound after summing h | Hardest h |
|---|---:|---:|
| 1850 | -17566.92 | 3668 |
| 1974 | -2866.61 | 3626 |
| 2100 | +14865.59 | 3626 |

These are not directed certificates or coverage of intermediate row counts.
Extending the regional p candidates from at most 0.67070 to 0.76 did not
change these three combined results. The tested p range was therefore not
the source of this remaining gap.

The interpolation can also be simplified. Put

    ell' = ell (q ell)^(-tau).

The interpolated transfer's state/death factors, prefactor, and potential
are exactly those obtained by replacing ell with ell'. Its local selected-zero
polynomials still use ell. Because ell'<=ell and coefficients are nonnegative,
using P(ell') instead of P(ell) is at least as strong. Future proof searches
should vary ell directly rather than treat tau as essential new machinery.

## Retain inactive occupancy through each entire region

The new helper `quartet_joint_quota.py` keeps a polynomial in the number
of inactive-row labels. For fixed incoming weight u, it computes the joint
coefficients in u and target count t of 128 local quartet factors, divided
by binom(512,u). This division averages the uniform quartet partition.
All coefficients are positive; the diagnostic uses log-domain products,
not an FFT or a cutoff of small coefficients.

The live-to-live envelope takes the coefficientwise maximum over incoming
weights 216--232. Reset coefficients are treated similarly. Activation and
quiet coefficients are explicit: quiet coefficients are the Bin(512,p)
probabilities, and activation coefficients additionally sum
binom(512-t,w) X^w over w=216,...,232.

These four polynomials form a two-state matrix. Multiply it 16 times,
extract coefficient r, and divide by the regional binomial probability
pi_r(p). This retains inactive occupancy jointly with the state-class paths
within the region. Compose 256 such regional matrices, with initial state
zero and an unrestricted terminal class.

The artificial Bernoulli parameter p cancels from this exact coefficient
expression. This provides an independent numerical check. No optimization
of p is needed in this version. Repeating the complete calculation at p=1/2
and p=112627/167924 changes the three final logarithms by at most
4.66e-10 bits, consistent with floating-point rounding.

The cheap live envelope initially sums all outgoing weights and turnover.
To test that relaxation, additional coefficient bounds insert markers
z^(v-u) y^d, where v is next weight and d is deaths. For v in [216,232],
divide by the appropriate endpoint power of z; for d in [24,40], do the
same for y. Take the smallest bound for each incoming-weight/occupancy
coefficient before maximizing the incoming weight. These bounds retain
the two interval constraints without a large three-dimensional array.

## What the comparison says

At tau=0.1 with the preceding frozen markers:

| Joint regional envelope | r=1850 | r=1974 | r=2100 |
|---|---:|---:|---:|
| Outgoing weight/deaths relaxed | +16744.67 | +19539.35 | +23068.31 |
| Outgoing constraints bounded by eight extra marker choices | +15226.41 | +18055.39 | +21347.91 |

These rows use a single tau and sum histories directly. They are not directly
comparable to the preceding best-per-h table. More importantly, exact regional
occupancy did not by itself close the scalar state envelope. The outgoing
constraints help, but the tested bounds still fail.

Two checks of global proof weights also remained inconclusive. With tau=0,
ell=0.00328, theta=0.2, and X=0.22, the three screens were approximately
-5953, +7271, +21132. At theta=0.15 and X=0.32 they worsened. These are
bounded diagnostic choices, not an optimization result or a distance ceiling.

Separately, [mixed minors](QUARTET_MIXED_MINORS.md) give a proved setup-event
bound for the actual construction. Their numerical gain in this screen is
less than one bit, smaller than their union cost. They are not the next lever.

## Next discriminating calculation

Keep the joint occupancy count, but retain transitions between distinct
state weights or small weight intervals. The current coefficientwise maximum
collapses the entire incoming band, so successive updates can use different
worst incoming weights without accounting for the state flow that connects
them. The outgoing marker test limits each endpoint but does not retain that
flow. This is a remaining relaxation to test, not a demonstrated obstruction
in the code.

Do not claim that closing this band would finish the goal. Other live weights,
turnover ranges, row counts, and their transitions still require a global
coverage argument and a directed union. The reusable target is a composable
weight-state operator, not indefinitely adding isolated profile certificates.

## Resource and validation record

The joint helper uses arrays below 16 MiB, a 60-second watchdog, and the
20-GiB disk guard. Baseline joint runs took about 8--10 seconds; the eight
extra outgoing markers took about 45 seconds. Runs were serial and file-free.
Tests compare 52 polynomial coefficients, 16 short matrix products, and
10 normalized occupancy coefficients with independent exact calculations.
No benchmark, native/library/paper edit, setup change, commit, or push.
