# RM(4,9) Outer Route

Status: active finite-\(n\) investigation.  The first target is a direct
exact-spectrum transfer to the already defined uniform-routing,
RandomStepConv-M22 proof model.

The occupation-one gate has now been certified.  It gives only

\[
  Q_1 < 2^{-26.5608944032},
\]

and is dominated by local weight 32.  Consequently the current
occupation-one transfer does not by itself meet a 40-bit target at 10.9%.
This is a limitation of the present bound, not a lower bound on the true bad
probability and not a refutation of the RM construction.

## Outer code

The local constituent is the binary Reed--Muller code

\[
  \operatorname{RM}(4,9)=[512,256,32].
\]

Its 256 input coordinates index the monomials in nine variables of degree at
most four.  Its 512 output coordinates are the evaluations at every point of
\(\mathbb F_2^9\), in a fixed order.  The global outer is the direct sum of
4096 copies of this one fixed constituent.  The code is not resampled between
blocks and has no outer setup-failure event.

The exact spectrum is `scripts/rm512_256_spectrum.csv`, with SHA-256
`995aab561da18f22074b5c6f5413882f492084510aa1cd19b848355f1fcd4ed7`.
It has total mass \(2^{256}\), minimum nonzero weight 32, complement symmetry,
and support on 107 weights.  `scripts/verify_rm_gleason_formula.py` checks the
table against the Type-II Gleason expansion.  This workstream treats that
authenticated integer table as exact input.

## Current transfer target

The target parameters are

\[
  k=2^{20},\qquad N=2^{21},\qquad D=228{,}590,
  \qquad D/N=10.9000205994\%.
\]

For each of the 4096 outer rows, routing samples one independent uniform
permutation of its 512 coordinates.  For each of the 512 transposed regions,
routing then samples one independent uniform permutation of the 4096 row
positions.  RandomStepConv-M22 independently samples one uniform binary
\(23\)-by-\(23\) linear map at every output position.  The probability in the
desired theorem is over these routing and inner maps; the RM constituent is
fixed.

The first proof gate was occupation \(Q=1\), beginning with the 52,955,952
weight-32 codewords.  `certify_rm49_q1_outward.py` sums every nonzero shell
with outward arithmetic and records the result in
`rm49_q1_randomstepconv_M22_d109_outward.json`.  The certificate obtains
26.5608944032 bits, not 40 bits.  Binary64 diagnostics at memories 23, 27,
and 40 obtain 29.3433, 33.9054, and 34.5109 bits, respectively.  The last
value shows that increasing memory alone has little remaining leverage in
this transfer.

Occupation two could use the exact product spectrum.  Occupations at least
three require a fresh band or shell-sensitive transfer.  The sparse-EA caps
cannot simply be reused: they vanish below weight 42, whereas RM has nonzero
weight 32, and RM's tail multiplicities are much larger than random-code
means.

## Prior proved checkpoint

`scripts/certificates/rm512_256_sig32_delta009.json` records a proved 9%
checkpoint with certified total log probability at most \(-40.520728\) in a
different inner ledger.  `scripts/verify_rm_block_checkpoint.py` rechecks its
prefix, tail, and exact spectrum.  This establishes that the spectrum was
used successfully before.  It does not prove the present 10.9%,
RandomStepConv-M22 target.

## Implementation hypothesis

A length-512 Boolean Möbius transform evaluates all nine-variable Boolean
functions using at most \(9\cdot 2^8=2304\) XORs.  RM encoding restricts the
256 input coefficients to degree at most four.  The transposed map is another
butterfly followed by selection of those 256 coefficients.  This gives a
credible path below the 5,868-XOR two-BCH comparison, but it is presently a
cost bound, not a benchmark.

## Proof obligations

1. Decide whether to accept a margin below 40 bits, lower the distance target,
   or strengthen the routing/inner argument for the weight-32 shell.
2. If the revised occupation-one gate closes, certify \(Q=2\) from the exact
   repeated-constituent product spectrum.
3. Replace the sparse-EA high-occupation cap domination with an RM-compatible
   exact-spectrum transfer.
4. Combine all occupations and state the probability only over routing and
   RandomStepConv-M22.
5. Only after proof viability is established, implement and benchmark the
   forward and transposed Möbius kernels sequentially against the existing
   BCH outer.
