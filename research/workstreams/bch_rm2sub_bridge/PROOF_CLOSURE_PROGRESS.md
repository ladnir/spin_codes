# Closure attempt after the implementation study

Latest update: the selected t128_s19 map is fully certified at K=2^18,
with a 52.3463883689-bit distance/setup margin and coverage of Q=1,...,2048.
The preceding K=2^16 certificate has a 53.9443672720-bit margin.
See [T128_S19_M18_CLOSURE.md](T128_S19_M18_CLOSURE.md) for the exact union,
512-bit replay, and scalable all-one-split search. The initial search results
below describe an earlier checkpoint.

The selected t128_s19 map was the first target for a new full BCH-256
certificate. The implementation study established its speed and correctness.
The proof attempt started at message dimension K=2^16,
with 512 outer rows and bad-output cutoff floor(2K/10)=13107.

## Current t128_s19 coverage

The serial search has certified every occupancy Q=1,...,111, with 512-bit
numerical replays and an exact rational union ledger. Their combined upper
bound corresponds to 54.0067322321 bits. Occupancies 112,...,512 remain
uncovered; the reported 54 bits are not a full-code margin.

The accepted ledger is
`generated/closure_t128_s19_m16_v1/report_0001.json`.
The first 12-job run covered Q=1,...,6. A resumed run with a 150-second
budget extended coverage through Q=111. No bounds were transferred from
a different message size or inner map.

At Q=512, seven fixed-weight trials used tilt indices -4,-2,0,2,4,6,8,
where lambda=exp(index/10). All produced weak upper bounds. The best
diagnostic margin was -8826.4 bits, at index 2. These weak trials were
not accepted or replayed as certificates.

A separate scalar-density screen tested Q=32,64,128,170,256,384,512 with
three row-witness banks and 41 tilts. It also remained weak, including
-13209.1 bits at Q=512. Those values use nearest binary64 arithmetic and
sampled densities; they are search diagnostics, not certificates or
evidence that the code fails. The receipt is
`generated/closure_t128_s19_m16_dense_screen_v1.json`.

## A first-moment obstruction for t256_s14

The separate t256_s14 diagnostic changes its proof priority. At K=2^20,
the expected number of bad messages exceeds 2^150700. Therefore, tightening
an unconditional first-moment upper bound cannot prove a 40-bit margin
for this parameter choice. This conclusion does not establish a lower
bound on the probability that a setup has any bad message.

The expectation is over the existing setup model: independent uniform
row and region permutations, with the selected inner map and zero initial
state. The multiplier randomness does not affect the event used below.
The fixed outer is the project's BCH [256,128,d>=38] subcode.

Let T contain its words of weights 38 through 80. The retained exact BCH
constraints prove

\[
|T|\ge a_0:=443405397513809550622719915749.
\]

Choose Q=2620 occupied positions among L=8192 rows, with a word from T
in each position. Every resulting message has input weight at most
80Q=209600, below the bad-output cutoff 209715.

We count messages whose input to every inner epoch lies in ker B. The
state then remains zero, and the output equals the permuted input. Thus
each counted message is bad, independently of the multipliers.

The even-row convolution and range-restriction argument in
[ZERO_STATE_FIRST_MOMENT_OBSTRUCTION.md](ZERO_STATE_FIRST_MOMENT_OBSTRUCTION.md)
applies unchanged. When the words from T are sampled for this counting
argument, the probability that every region has an even input weight
between 64 and 1536 exceeds 2^-256. This auxiliary sampling averages over
the counted messages; it does not change the setup distribution.

Let H(z) be the exact weight enumerator of ker B for the selected
t256_s14 map. Each region has L/t=32 epochs. Given its input weight j,
a uniform region permutation places every epoch in ker B with probability

\[
\beta_j=\frac{[z^j]H(z)^{32}}{\binom{8192}{j}}.
\]

Outward arithmetic verifies beta_j>2^-447 for every even j from 64 to
1536. Conditional on the row outcomes, the region permutations are
independent. Consequently, the average probability of the zero-state
event is greater than

\[
2^{-256}(2^{-447})^{256}=2^{-114688}.
\]

If Z counts all bad messages, linearity of expectation gives

\[
\mathbb E_{\rm setup}[Z]
>\binom{8192}{2620}a_0^{2620}2^{-114688}
>2^{150700}.
\]

The final comparison uses exact integers. Its diagnostic logarithm is
150741.8982893 bits. The proof concerns this fixed map at K=2^20;
it does not transfer an obstruction to other state sizes or message sizes.

`certify_t256_s14_obstruction.py` rechecks the BCH tail's exact primal and
dual solution, computes the region probabilities, and verifies the integer
comparison. It stores the compact formula rather than its expanded integer.
The producer uses 256-bit Arb arithmetic; `--verify` recomputes at 512 bits.
Existing producers and hash-bound receipts remain unchanged.

The first replay attempt found the tail receipt absent from this workspace.
The four files `lower.json`, `h_38.lp`, `h_38.sol`, and `warm.bas` were
restored byte-for-byte from the original ba80 worktree's
`generated/joint_tail_lower_80` directory. Their total size is about 8.9 MB;
they remain ignored local replay inputs. The new obstruction receipt records
the dependency hashes without changing the frozen migration manifest.

From this directory:

```text
python -B certify_t256_s14_obstruction.py
python -B certify_t256_s14_obstruction.py --verify
python -B test_zero_state_parity.py
```

The 256-bit producer and 512-bit replay both passed, including fresh
exact BCH-tail checks. The parity and certificate-search regression suites
passed all 13 tests. The separate report verifier accepted the t128_s19
partial-coverage ledger; it did not rerun those already-recorded numerical
replays a second time.

## Next proof priority

Continue t128_s19 at K=2^16. Close the sparse occupancies with the existing
outward engine, then refine mixed row-weight types in the dense range.
The entrywise maximum over bands and the scalar-density relaxation both
lose too much information there. Weak upper bounds are not evidence of
an obstruction for t128_s19.

Keep t64_s16 as the implementation-backed fallback and t64_s20 as the
already-certified reference at K=2^20. Keep the t256_s14 performance result,
but do not spend further effort trying to close its unconditional first
moment at K=2^20. Returning to that configuration requires a different
probability argument, rather than a tighter implementation of this bound.
