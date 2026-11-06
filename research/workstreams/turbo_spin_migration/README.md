# Turbo SPIN results before manuscript revision

This workstream asks whether Turbo's packet-route/s22 architecture can
carry the finite, asymptotic, encoding, and application results now
attributed to Structured SPIN. Its finite outer is RS8; its asymptotic
outer must grow, as the paper's current Golay--BAA outer does. Work
here before making further manuscript edits.
The provisional paper draft is frozen at `0f3f0c4a` on
`codex/turbo-spin-paper`; this research branch starts from `f0656bf4`.
Neither the draft's new claims nor a timing ratio substitutes for a matched
certificate and implementation result.

**First result.** The [setup audit](SETUP_AUDIT.md) found no algebraic
map or sampling-interface mismatch in the RS8 profile, while preserving
the ideal-setup versus seeded-implementation distinction. The
[one-group obstruction](ONE_GROUP_OBSTRUCTION.md) proves that the fixed
RS8/64-region family has asymptotic minimum distance at most 5/64 with
probability tending to one. The direct 10%/11% asymptotic transfer gate
is closed negatively; finite-length migration remains open. A scalable
Turbo-like family needs a new local packet-support floor or a growing
outer. No manuscript claim is changed by this result.
The [selected shared-outer candidate](SHARED_OUTER_CANDIDATE.md) uses
the exact RS8 expected enumerator to show that fixed local profiles
with packet minimum distance at least 19 exist. It has no global
distance certificate yet; sharing a profile changes the outer proof.
The [pair-return obstruction](PAIR_RETURN_OBSTRUCTION.md) shows that
even a selected good fixed outer cannot support an asymptotic
`1-o(1)` positive-distance theorem with the current fixed-width route
and uniform-GL inner: a fixed one-group message has a small but
length-independent probability of pairwise state returns and output
weight at most `N/16`. The shared-profile candidate remains relevant
only to finite-size or redesigned-inner questions.
The [family bridge](TURBO_FAMILY_BRIDGE.md) fixes the replacement target:
finite RS8 and growing Golay--BAA outers with Turbo's packet route and
s22 inner. The [growing RS8 base](SCALABLE_RS8_BRIDGE.md) is a later
outer-performance option, not a prerequisite for the asymptotic proof.
The [one-active-block packet transfer](Q1_PACKET_TRANSFER.md) derives an
exact BAA packet-support enumerator and an s22 finite-step upper transfer.
A preliminary screen using a selected 16% packet-support floor leaves
room for the old `39/4` outer-growth schedule in the one-block class.
The [packet-floor analysis](PACKET_FLOOR_TAIL.md) closes the sparse-source
part; its dense-source exponent still needs a uniform certificate.
The other occupation classes still need proofs;
the numerical screen does not establish a scalable distance theorem.

**Next decision.** Keep the existing Structured theorem as the control.
Certify the dense-source packet-tail exponent for the shared BAA outer
and audit the finite-step one-block transfer.
Then test the other occupation classes. The finite K=2^22/2^24 certificate sweep is
independent of that asymptotic proof.

## Baseline and gaps

| Result | Retained evidence | Migration gap |
|---|---|---|
| Half-rate finite distance | RS8 packet code has >10% distance with 61.884, 58.246, and 51.822 bits of setup-failure margin at K=2^16, 2^18, and 2^20 | K=2^22 and 2^24, plus any intermediate lengths we claim |
| Transposed 128-bit encoding | Public RS8 profile: 0.099526, 0.461416, and 2.861814 ms at those three lengths | A matched BCH/IMT and BAA campaign; scaling at larger lengths |
| Scalable distance | Growing Golay--BAA outer with bit route and IMT: rate 1/2, >11% distance with probability 1-o(1), O(N) ordinary and transposed work, Lean proof | Keep the growing outer as control; certify it with Turbo's packet route and s22 inner |
| Quarter-rate distance | Older BCH/IMT profile at K=2^20: >16.5% with 40 bits and >19% with 30 bits | A binary RS-packet quarter-rate construction, certificate, and timing |
| Ordinary and wide encoding | Older profile has tuned forward and wide kernels; RS8 forward is a scalar fallback | Optimized Turbo forward kernels and matched measurements |
| Applications | Existing Silent OT, PCS, and Flock results use BCH/IMT | Re-run with the exact certified Turbo profile and report complete workload boundaries |

The [RS8 profile](../../../spin/RS8.md) identifies the current public code,
setup distribution, certificates, and benchmark protocol. The
[scalable theorem](../../paper/structured_proof.tex),
[BCH finite theorem](../../paper/finite_certificates.tex), and
[application sections](../../paper/applications.tex) define the results
that could be displaced. The older asymptotic construction has a growing
outer; its theorem does not transfer merely because both encoders use one
structured permutation.

## Work order and decision gates

1. **Freeze the comparison object.** Pin the exact RS8 generator, column
   randomizers, packet route, inner maps, setup independence, and forward
   adjoints. Verify that the public encoder and each proof use the same
   distribution. Keep the existing three finite certificates as controls.
2. **Test the common Turbo architecture.** Use the paper's growing
   Golay--BAA outer first, with `b=Theta(log N)` and `b/4` packet regions.
   Keep Turbo's packet route and s22 inner. Start with fixed occupation,
   including the binary-weight-to-packet-support and label interface.
   Then cover growing sparse and positive occupation. The gate is a
   complete >11% / O(N) theorem or a precise obstruction for this
   route/inner. Do not require an RS8-based growing outer for this gate.
3. **Extend the finite ladder.** For K=2^22 and 2^24, first evaluate the
   existing RS8/s22 distribution at the 10% cutoff. If its margin falls
   below 40 bits, change one parameter or component at a time and record
   the effect on both the certificate and online cost. Verify every active
   group count, exact union, and map/route identity. Do not infer a
   certificate for intermediate lengths from the three known endpoints.
4. **Design the quarter-rate member.** Specify its outer dimension,
   packet route, and inner map before screening. Seek the older 16.5%/40-bit
   and 19%/30-bit targets at K=2^20, allowing a different RS constituent.
   Check its actual distribution and implementation separately.
5. **Build and time the required kernels.** Optimize 128-bit forward
   encoding first, then wide/generic paths needed by consumers. Preserve
   the fast transposed path. Use correctness and adjoint checks for each
   implementation change. Benchmark one process at a time; never run two
   benchmarks concurrently. Compare candidates in both orders with fixed
   CPU, compiler, memory policy, setup exclusion, and input-reset policy.
6. **Migrate application measurements.** Integrate the certified profile
   into the precomputed Silent OT mode, then the PCS row encoder and Flock.
   Report encoder time and complete workload time separately. Recompute
   each application's conditional security calculation with the selected
   distance and setup-failure bound; do not carry over BCH/IMT numbers.

## Manuscript decision

If Turbo SPIN closes the scalable theorem and the finite, forward, and
application gates, make it the main construction and retain only older
proof ingredients that explain or support it. If the scalable theorem
remains unique to the older family, present one SPIN architecture with a
scalable theorem family and a faster finite RS-packet family. In either
case, update the abstract, comparison table, and applications together
only after their evidence is matched. Preserve complete older proofs in
the appendix or research record when removing them from the main story.

Each gate should leave a compact source pin, exact claim, proof/benchmark
protocol, and result summary here. Keep raw timing logs, large receipts,
build products, and generated plots outside version control.
