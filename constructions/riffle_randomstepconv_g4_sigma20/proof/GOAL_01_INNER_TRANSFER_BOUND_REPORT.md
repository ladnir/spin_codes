# Goal 01 report: exact inner transfer bound

## Result

Goal 01 is complete.  The RandomStepConv inner admits an exact two-state
transfer enumerator.  A rank-one reduction converts its fixed-support moment
into a positive sum indexed by the number of terminated live episodes.  The
implementation bounds each summand with an independent coefficient tilt.

The resulting target-size calculation does not expand a polynomial in the
packet count (N=524352).  Its work for one pair ((h,D)) is linear in the
input packet support (h).

## Exact identities

Fix the nonzero values of the (h) active input packets.  The probability is
over the global packet permutation and the independent setup matrices.

For a live step, a uniform random matrix maps the nonzero input-state vector
to a uniform element of \(\mathbb F_2^{g+\sigma}\).  Therefore, the output
packet is uniform in \(\mathbb F_2^g\), and the next state is independently
uniform in \(\mathbb F_2^\sigma\).  This fact proves the two-state transfer
matrices in the Goal 01 statement.

The active transition matrix has rank one.  For

\[
b(z)=2^{-g}(1+z)^g,
\qquad
d(z)=(1-2^{-\sigma})b(z),
\qquad
\alpha(z)=\frac{2^{-\sigma}b(z)}{1-d(z)},
\]

one post-activation gap of length \(\ell\) contributes

\[
f_\ell(z)=\alpha(z)+(1-\alpha(z))d(z)^{\ell+1}.
\]

Expanding the product of the (h) gap factors gives the positive
termination-count decomposition recorded in the Goal 01 statement.  No
probabilistic approximation enters these identities.

## Verification

The checker exhaustively enumerated every matrix sequence in three small
instances:

| (g) | \(\sigma\) | (N) | active value | matrix sequences | result |
|---:|---:|---:|---:|---:|---|
| 1 | 1 | 3 | 1 | 4,096 | exact match |
| 2 | 1 | 2 | 1 | 262,144 | exact match |
| 2 | 1 | 2 | 3 | 262,144 | exact match |

Each instance matched both the transfer coefficients and the rank-one gap
identity.  The two (g=2) instances also check invariance under the nonzero
packet value.

## Target-size bounds

The following table reports outward-rounded upper bounds.  The binary output
length is 2,097,408.  The terminal state is discarded.

| Relative threshold | (h=18) | (h=24) | (h=92) | (h=512) | (h=2048) |
|---:|---:|---:|---:|---:|---:|
| 5% | \(2^{-31.5902}\) | \(2^{-50.2693}\) | \(2^{-267.2045}\) | \(2^{-1639.4995}\) | \(2^{-6733.5861}\) |
| 9% | \(2^{-16.7937}\) | \(2^{-30.4413}\) | \(2^{-190.0895}\) | \(2^{-1206.2087}\) | \(2^{-4981.6084}\) |
| 12% | \(2^{-9.6899}\) | \(2^{-20.8920}\) | \(2^{-152.6298}\) | \(2^{-994.9681}\) | \(2^{-4128.7712}\) |

The outer BCH construction activates at least 18 packets.  Bounds below this
support are not required for the inherited outer code.

## Numerical status

Floating-point arithmetic selects the output tilt and the component radii.
The verifier converts the selected values to exact decimal rationals.  It then
evaluates every positive summand with 90-digit decimal arithmetic and directed
rounding.  The receipt stores both the floating diagnostic and the
outward-rounded result.

The calculation proves only the conditional inner bounds in the table.  It
does not yet sum the outer spectrum or prove an end-to-end minimum distance.
