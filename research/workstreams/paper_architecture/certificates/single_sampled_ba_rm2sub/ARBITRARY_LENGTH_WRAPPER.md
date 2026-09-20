# Arbitrary-length wrapper

## Objective

Each SPIN construction is easiest to define at lengths that satisfy its local
block, inner-step, and routing constraints.  This document gives one wrapper
for all other finite lengths.  The wrapper is deterministic and is applied
after proving or certifying an admissible-length code.

The wrapper does not infer a distance certificate at a new length.  It
transports an existing certificate with an explicit loss.

## Admissible lengths

For a construction parameter record

\[
  \vartheta=(B,K_B,t,s,\mathsf{route}),
\]

let \(\mathcal A_\vartheta\subseteq\mathbb Z_{>0}\) be the set of admissible
intermediate and output lengths.  A length \(N\) belongs to
\(\mathcal A_\vartheta\) only if all of the following applicable conditions
hold:

1. **Outer blocks:** \(B\mid N\), and the outer input length is
   \(K=(K_B/B)N\in\mathbb Z\).
2. **Inner steps:** \(t\mid N\).
3. **State termination:** any required flush or tail length is included in
   \(N\), or the manifest defines the final-state rule explicitly.
4. **Factored interleaver:** every transpose, region, lane, and bucket count is
   integral.
5. **Permutation domain:** the sampled permutation family is defined on the
   exact domain used at length \(N\).
6. **Implementation-only constraints:** any power-of-two mask, tile size, or
   vectorized batch constraint is recorded separately from mathematical
   admissibility.

For a plain uniform interleaver with no tail, conditions 1 and 2 reduce to

\[
  \operatorname{lcm}(B,t)\mid N.
\]

The frozen Structured SPIN implementation has

\[
  (B,K_B,t,s)=(256,128,128,19)
\]

and the single frozen length \(N=2^{21}\) in 128-bit blocks.  Its source also
uses 256 regions, 8192 outer blocks, 2048 outer blocks per tile, and
power-of-two addressing.  These facts authenticate that instance.  They do
not yet define an arbitrary admissible-length family.

Every construction theorem should publish \(\mathcal A_\vartheta\), not leave
the divisibility rules implicit in code.

## Exact-input padding wrapper

Suppose an admissible construction gives a linear encoder

\[
  E^+:\mathbb F_2^{K^+}\to\mathbb F_2^{N^+}
\]

with minimum distance at least \(D^+\).  For any requested input length
\(K\le K^+\), define

\[
  E_{\mathrm{pad}}(x):=E^+(x\mathbin\|0^{K^+-K}).
\]

**Lemma 1 (input padding).** The padded encoder is injective and its image is
a binary \([N^+,K,d_{\min}\ge D^+]\) subcode.

*Proof.* The padding map is injective.  Every nonzero padded word is a nonzero
input to \(E^+\), so its output weight is at least \(D^+\). \(\square\)

This wrapper supports every message length.  Its output length is the next
admissible output length, so it may have a small rate loss.

## Exact-output puncturing wrapper

Let \(C^+\) be a binary \([N^+,K^+,D^+]\) linear code.  Fix a requested
output length \(N\le N^+\), and set \(q:=N^+-N\).  Puncture any fixed set of
\(q\) output coordinates.

**Lemma 2 (puncturing).** The punctured code has minimum distance at least

\[
  D^+-q.
\]

If \(q<D^+\), puncturing is injective on \(C^+\), so the punctured code has
dimension \(K^+\).

*Proof.* Puncturing deletes at most \(q\) nonzero coordinates from each
codeword.  If its kernel contained a nonzero codeword, that word would have
weight at most \(q<D^+\), a contradiction. \(\square\)

If the punctured code has dimension at least the requested \(K\), choose any
\(K\)-dimensional subcode.  Taking a subcode cannot decrease minimum distance.

The output coordinates to puncture must be fixed independently of the sampled
setup unless the probability experiment explicitly includes that choice.
Selecting coordinates after inspecting the code requires a separate proof.

## One wrapper for exact input and output lengths

Let \((K^+,N^+,D^+)\) be a certified admissible instance.  Given requested
\((K,N)\) with \(K\le K^+\) and \(N\le N^+\), perform these steps:

1. Restrict the admissible encoder to a \(K\)-dimensional input subspace.
2. Puncture a fixed set of \(q=N^+-N\) output coordinates.

**Theorem 3 (arbitrary finite length).** If \(q<D^+\), the wrapped code is a
binary

\[
  [N,K,d_{\min}\ge D^+-q]
\]

linear code.  If \(q\ge D^+\), the valid transported lower bound is only
\(\max\{D^+-q,0\}\); injectivity then requires a separate check.

For an ensemble with certified failure probability \(\varepsilon^+\), the
same failure bound applies to the deterministic wrapper.  The wrapper adds no
randomness and cannot invalidate a good admissible setup when \(q<D^+\).

This theorem separates two choices:

- The admissible-length certificate determines \(D^+\) and
  \(\varepsilon^+\).
- The wrapper determines the explicit rate and distance loss.

## Selecting the neighboring admissible instance

For each requested \(N\), define

\[
  N^+(N):=\min\{M\in\mathcal A_\vartheta:M\ge N\},
  \qquad q(N):=N^+(N)-N.
\]

Use the certified code at \(N^+(N)\), then apply Theorem 3.  If admissible
lengths are multiples of \(L\), then \(0\le q(N)<L\).

Suppose admissible instances have rate at least \(R-o(1)\), relative distance
at least \(\delta-o(1)\), and \(q(N)=o(N)\).  Then the wrapped family has

\[
  \frac{K}{N}\ge R-o(1),
  \qquad
  \frac{d_{\min}}{N}\ge\delta-o(1),
\]

provided the requested dimension obeys \(K\le K^+(N)\).  For constant
\(L\), the distance loss is at most \(L-1\).  For
\(L=\Theta(\log N)\), the relative loss is
\(O(\log N/N)\).

## Adjacent block sizes

When \(B=B(N)=\Theta(\log N)\), rounding every local block to one size can
leave \(\Theta(B)\) coordinates.  A two-size outer layer can remove this
remainder without puncturing.

Fix adjacent allowed output block sizes \(B\) and \(B+1\).  If

\[
  N=uB+v(B+1)
\]

for nonnegative integers \(u,v\), use \(u\) blocks of the first constituent
and \(v\) blocks of the second.  This method is valid only after proving:

- compatible local rates and spectrum envelopes for both constituents;
- a common interleaver and inner-transfer theorem for the mixed block list;
- explicit ordering and setup independence for the two block types; and
- the inner-step and routing constraints at the resulting total length.

The Frobenius bound implies that every sufficiently large integer has such a
representation because \(\gcd(B,B+1)=1\).  This arithmetic fact does not prove
the spectrum or transfer conditions.  Until those conditions are discharged,
puncturing the next admissible instance is the theorem-safe default.

## Shortening

Shortening a binary \([N^+,K^+,D^+]\) code on \(q\) coordinates means taking
codewords that are zero on those coordinates and deleting the coordinates.
The shortened code has length \(N^+-q\), dimension at least \(K^+-q\), and
minimum distance at least \(D^+\).

Shortening is attractive when the desired dimension tolerates its possible
dimension loss.  It is not interchangeable with input padding: a chosen set
of output coordinates need not correspond to padded input coordinates in a
nonsystematic SPIN encoder.

## Rules for theorem statements

An arbitrary-length theorem should state:

- the requested \((K,N)\);
- the selected admissible \((K^+,N^+)\);
- the certificate bound \((D^+,\varepsilon^+)\);
- the operation used: padding, puncturing, shortening, or two block sizes;
- the resulting dimension and distance bounds; and
- whether the wrapper changes the encoder interface or only its parameters.

No construction proof should duplicate this boundary logic.  It should prove
its claim on \(\mathcal A_\vartheta\) and invoke Theorem 3.

For `RANDOM_OUTER_RM2SUB_CERTIFICATE.md`, take \(L_m=128m\) and let \(B_m\)
be the least positive even integer with
\(B_m\ge9\log_2(L_mB_m)\). The resulting lengths satisfy

\[
  \frac{N_{m+1}-N_m}{N_m}=o(1).
\]

Thus next-admissible puncturing preserves its asymptotic rate and relative
distance. The public interface must still choose whether exact output length
or exact message dimension has priority.

## Open obligations

1. Define mathematical admissibility for each remaining scalable SPIN
   variant.
2. Separate those rules from frozen implementation tiling constraints.
3. Prove a bounded-gap property for each admissible-length set.
4. Decide whether exact output length or exact message length is the public
   interface.
5. If adjacent block sizes are used, prove both local spectrum envelopes and
   the mixed-block transfer theorem.
