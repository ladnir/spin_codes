# BCH spectrum continuation status

## Original ideal-M22 goal closed, 2026-09-05

Read [RANDOM_INNER_M22_CLOSURE.md](RANDOM_INNER_M22_CLOSURE.md).
The fixed BCH-derived [256,128] outer at k=20 now has a deterministic-spectrum
computer-assisted distance theorem under ideal RandomStepConv-M22:
failure <= U < 0.61 * 2^-40, with diagnostic margin 40.72497637925699 bits.
The event is a nonzero message with output weight <= 209716, at output
length 2^21. No statistical shell assumptions or variance estimate are used.

The final improvement jointly bounds the paired weights 38, 40, and 42
under the certified OA29/hull constraints. Exact dual verification checks
all 1163 rows and the general objective. The full bound includes every shell
and all 8192 occupations. It does not require A38 <= 10^13 or a full spectrum.

Frozen receipts: generated/shift_rank_oa29_joint/audit.json and
generated/bch256_m22_unconditional_closure_audit.json. Final audit and replay
passed: 78 JSON evidence files, 92 independently recomputed M22 transfer
coefficients, exact higher-occupation reaggregation, and reconstructed LP.
The separate cell-DP replay matched every array bit-for-bit. Base rank,
refined rank, OA21, signed-weight, and Fourier toy audits also replayed.

All computations are terminal. No benchmark ran concurrently with another.
The statement does not transfer to RM2Sub, a PRG, equation (19) of the draft,
or 11% distance. Next: independent review and paper integration; keep the
practical-inner transfer as a separate task. Earlier entries below are history.

## Refined rank cover gives OA29; M22 now39.72180bits, 2026-09-05

An exhaustive split of Q Fourier case15 on whetherF23iszero gives two
rank29witnesses. All other cases already reach29orhigher, so independent
auditing proves extendedQdual>=30 and Q/cosetsOA29. No literature input.
See BCH_SHIFT_RANK_CERTIFICATES.md. New exact1163-row model gives
A38(C)<=24866368872377, fullM22 margin39.7217972337603bits. Goal ACTIVE.

Replay certify_bch_shift_rank_split.py --verify; audit_bch_extended_hull_cap.py
prepare_bch_shift_rank_oa29_probe shift_rank_hull_probe --verify.
Joint optimization of weights38/40/42is now inprogress; the current separate
caps miss the target byabout21%. Existingfrozen evidence retained.

## Independent dual-rank certificates; M22 now38.64383bits, 2026-09-05

Read BCH_SHIFT_RANK_CERTIFICATES.md. Direct Fourier independent-set and
root-run witnesses cover all17Q and16P nonzero Fourier cases. Independent
auditing reconstructs generators, root sets, parity, every witness step,
and affine invariance. This proves extendedQdual>=24 andPdual>=26 without
the inaccessible SchaubPlus table. The formerly provisional equations now
have independent proofs.

New certified model generated/shift_rank_hull_probe reuses the exact old
solution only after exported-LP identity checks. Rational audit gives
A38(C)<=90005732934389 and fullM22 margin38.64383138557605bits. Goal ACTIVE.
Replay certify_bch_shift_rank.py --verify and audit_bch_extended_hull_cap.py
prepare_bch_shift_rank_probe oa21_hull_dual18_probe --verify.

Further complete-cover auditing proves Qdual>=28 (not30): receipt
generated/bch256_shift_rank_q28.json, witnesses shift_rank_q30_complete.
Onlycase15stalls at28; the other cases reach29orhigher. Further branching
and stronger moment models are in progress. No higher M22 cap yet claimed.

## Canonical orbit method and provisional stronger moments, 2026-09-05

Read BCH_AFFINE_CANONICAL_AND_SCHAUBPLUS_LEAD.md. Ordered-pair normalization
now provides exact affine canonical keys and stabilizer sizes without full
orbit expansion. Independent Python verification matches all97 previously
completed weight18 orbits; nontrivial-stabilizer regression cases pass.

At weight20,2040verified orbits contain133171200Dwords, none inHQdual.
This proves only A20(HQdual)<=832972800, not zero. Two candidate methods
give the SAME orbit union. No new cap from this partial result is claimed.

A newly found indexed Ponchio-Sala2003 table reports stronger dual-distance
lower bounds24(Q0) and26(P0). Primary PDF retrieval failed; these inputs
remain unaudited. A clearly PROVISIONAL1157-row exact LP, adding Kq22=Kh22=0
and K(q+255h)24=0, gives conditional cap90005732934389 and margin38.64383.
The certified margin remains37.57062; goal ACTIVE. Do not promote this
conditional result or use its folder as an audited envelope.

Next: audit or independently certify the new dual distances, then try a
joint weighted objective and full intersection constraints under that input.
Frozen receipts retained; all substantial runs sequential and terminal.

## Exact weight18 shell; full intersection obstruction, 2026-09-05

Read BCH_HULL_DUAL_WEIGHT18.md. A complete 97-orbit enumeration exhausts
all6332160 weight18 words of the known supercode D=hull(L71)^perp. None
lies in HQdual, so A18(HQdual)=0 exactly. The candidate construction is
not assumed exhaustive: matching the published exact count proves that.

The new196-variable1154-row exact certificate gives
A38(C)<=226659579388193 and fullM22 margin37.570616954093566bits.
Goal ACTIVE, not closed. Exact feasibility witnesses also establish that
the full intersection constraints plus A14/A16, and plus A14/A16/A18,
remain insufficient even with exact inner tails. These hypothetical
spectra do not lower-bound actual BCH failure. See the new note for replay.

Next: new code-specific shell/structural constraints, not more cuts of the
exhausted relaxation. Preserve frozen evidence; all jobs sequential.

## Exact dual-hull shells; M22 now37.56826 bits, 2026-09-05

Read BCH_HULL_DUAL_LOW_SHELLS.md. All65280 weight14 words of the known
supercode D=hull(L71)^perp form one explicitly verified affine orbit; none
lies in HQdual. Complete weight16 enumeration combines3212592 affine
four-flats with195840 skew unions of affine three-flats. Exactly16592 of
the full3408432-word shell belong to HQdual. Hence d(HQdual)=16 exactly.
These are deterministic code-specific counts, with the published spectrum
anchor supplying the completeness cardinalities.

Adding both exact shell equations gives a196-variable1153-row rational
certificate: A38(C)<=227085449947170; fullM22 margin37.56825962339644bits.
The40bit goal remains ACTIVE. The smaller model's true-tail primal
obstruction still holds. All experiments were sequential and terminal.

Earlier this turn, Farkas rounds2/3 passed exact checks and mapped warm
starts, but gains were negligible. A separate rational feasible point for
the full OLD intersection model has h38>=5e12 and true-tail lower first
moment>2^-40. It proves that OLD relaxation insufficient, not BCH unsafe.
It does not establish insufficiency after the new weight14/16 equations.

    python -B code/bch_hull_intersection_feasibility.py verify
    python -B code/certify_bch_hull_dual_weight14.py --verify
    python -B code/certify_bch_hull_dual_weight16.py --verify
    python -B code/audit_bch_extended_hull_cap.py prepare_bch_hull_dual16_count_probe oa21_hull_dual16_probe --verify

Next: combine full intersection constraints with the new exact dual shells;
constructively enumerate the known supercode's weight18 shell. Preserve all
hashed scripts and receipts. Do not iterate old cuts indefinitely.

## Intersection constraints converted to an exact cut, 2026-09-05

Read BCH_HULL_INTERSECTION_CUTS.md. J=HP intersect hull(L71) has dimension61
and a primitive quotient of size256 in the known dimension69 anchor. New
primal/dual inclusion constraints give a229-variable1329-row model; exact
fixed-coordinate substitution gives213variables1297rows. Bounded90s/120s
solver attempts ended without certificates. They are terminal, not proofs
of infeasibility. No old frozen model or receipt was edited.

The reciprocal containment HQ/HP into Qdual/Pdual was verified and added
to the196-variable model. Exact1148-row optimization proves no change to
its optimum. Fixing that exact primal leaves a17-variable J-extension LP.
An explicit exact Farkas witness proves this extension infeasible. Its16
nonzero multipliers eliminate t and yield one valid integer-coefficient cut.

The1149-row model with that cut solved and verified exactly. NewA38 cap:
228870861644645 (only8959 below the preceding cap); fullM22 margin:
37.55841850570869. The40bit goal remains ACTIVE. The numerical improvement
is negligible, but exact cut generation now avoids the unstable full model.
Next: iterate extension tests and cuts; obtain a feasible extension to prove
when this relaxation has been exhausted. All jobs sequential and terminal.

    python -B code/audit_bch_hull_intersection_structure.py --verify
    python -B code/audit_bch_extended_hull_cap.py prepare_bch_hull_benders_cut oa21_hull_reciprocal_probe --verify

## Coupled hull spectra and a derived exact anchor, 2026-09-05

Read BCH_HULL_SPECTRUM_CONSTRAINTS.md. Exact binary elimination constructs
nested doubly-even hulls HP[256,85] and HQ[256,93], both with minimum >=40.
Their primitive cyclic quotient gives common nonzero coset spectra. Coupled
MacWilliams and signed Fourier inequalities yield a 196-variable model.

The first exact 1017-row certificate improves the M22 margin to 36.62779.
The published dimension-71 code has polar rank two and negative Gauss sum;
its dimension-69 hull spectrum is exactly its weights divisible by four.
That hull lies in HQ and supplies an additional exact anchor. The augmented
1115-row model proves A38(C)<=228870861653604, yielding 37.55842 M22 bits.
Both rational primal/dual certificates and 60 exhaustive Fourier toy cases
passed. No statistical caps or orbit-lattice rounding are used.

    python -B code/audit_bch_hulls.py --verify
    python -B code/check_bch_hull_fourier.py --verify
    python -B code/audit_bch_hull_cap.py --verify
    python -B code/audit_bch_hull_anchor_cap.py --verify

The original M22 goal remains ACTIVE and not closed. The augmented exact
primal still exceeds the target when paired with the true-tail lower bound.
Next: couple the dimension-61 intersection of the new anchor with HP.
All substantial jobs ran sequentially and are terminal. Preserve receipts.

## Exact mod-four refinement and local-slack limitation, 2026-09-05

Read BCH_LOCAL_SLACK_AND_MOD4.md. Exact quadratic-form elimination gives
S(P)=S(Q)=-2^108 and S(H)=0. The polar radicals have dimensions 85 and 93.
Adding the signed identities to OA21 yields an exactly solved 404-row LP:
A38<=674685258232216, versus 678661807513927 before. Full M22 margin is now
36.12681 bits. The target A38<=10^13 remains unproved.

A smaller independent-block correction of the old dual witness was also
tested. All 20 local moment LPs have exact zero-cost primal witnesses, so
that fixed-dual/OA6 correction cannot improve the cap. This does not rule out
the fully coupled split LP. New larger solver attempts also timed out; all
jobs are terminal, and no rational split solution was produced.

    python -B code/verify_bch_local_slack.py
    python -B code/audit_bch_quadratic_sums.py --verify
    python -B code/audit_bch_mod4_cap.py --verify

Receipts: generated/bch256_local_slack_certificate.json,
generated/bch256_quadratic_weight_sums.json, generated/oa21_mod4_probe/audit.json.
The exact true-tail primal obstruction still holds for OA21 plus mod-four.
Next: exploit the radical/hull structure or coupled split constraints. Goal
ACTIVE; original M22 parameters unchanged. Preserve all frozen evidence.

## Code-specific split constraints verified; no new A38 cap, 2026-09-04

Read BCH_SPLIT38_CONSTRAINTS.md. For the explicit Wambach weight-38 support,
Q's projection has exact rank 38. Every inside pattern has 2^85 extensions.
The shortened [218,85] code has dual minimum at least 7, proved by exhaustive
integer checks of 218 columns, 23,653 pairs, and 1,703,016 triples. Therefore
every inside pattern's outside fiber is OA6. These are new exact local
moment equalities beyond the generic OA21 constraints. Replay and a separate
exhaustive small-code identity check passed.

    python -B code/audit_bch_wambach_shortening.py --verify
    python -B code/check_bch_split38_identities.py

The split LP has 1696 variables/720 rows. The code-specific reduced LP has
1656 variables/633 rows. Three sequential exact-solver attempts ended at time
limits (50,240,180 seconds); none produced a solution or an infeasibility
certificate. All jobs are terminal. See generated/bch_split38_attempts.json.
The prepared solution checker has not run, because no solution exists.

Current M22 bound remains 36.11873 bits, or 40.14372 conditional on A38<=10^13.
Next: stabilize the code-specific LP or extract smaller inequalities from the
certified local moments. The goal is ACTIVE; the A38 obligation is unchanged.

## M22 reduced to one unproved shell cap, 2026-09-04

Read M22_ONE_SHELL_REDUCTION.md. Cellwise brackets on active time replace the
one-row Chernoff bound at weights 38,40,42. The new full deterministic M22
bound has about 36.11873 bits, improved from OA21's 32.90151 bits.
An explicit bit-level toy enumeration and full directed replay passed.

Only A38<=10^13 is now needed to certify the full target, with about 40.14372
bits. Weights 40 and 42 use their deterministic OA21 caps; neither requires
statistical acceptance. The largest sufficient integer A38 cap for these
saved coefficients is 14,584,006,710,777. This is NOT a proved BCH cap.

The true-tail lower bound and an independent integer recurrence both show
that the feasible OA21 primal spectrum already exceeds the first-moment
budget at weight 38 alone. Exact inner tails cannot make the current OA21
relaxation suffice. This limits the relaxation, not the actual code.

    python -B code/certify_bch_q1_activity_cells.py --verify
    python -B code/audit_bch_q1_cell_budget.py --verify

New receipts: generated/bch256_q1_activity_cells_outward.json and
generated/bch256_q1_cell_budget.json. Preserve their hashed sources and NPZ.
Next: deterministic A38<=10^13, with the slightly looser exact threshold if
needed. No new sampling, no model substitution, and no goal completion.

## OA21 refinement: original M22 goal remains active, 2026-09-04

Read BCH_OA21_REFINEMENT.md. A published dual-distance bound, together with an
explicit translation/puncturing argument, raises the usable OA strength from
15 to 21. Three exact LP solutions passed all 402 primal rows, all 130 dual
columns, sign checks, and exact objective equality. No empirical lower bounds
or statistical shell caps are used. The deterministic full M22 bound improves
to approximately 32.90151 bits, still short of 40.

An exact feasible pseudoenumerator also proves that these constraints with
the current transfer coefficients cannot certify 40 bits, even under joint
objective optimization. This limits the relaxation, not the actual BCH code.
The published Schaub rank certificate itself has not been independently rerun.

    python -B code/certify_bch_oa21_probe.py --verify

New receipt: generated/oa21_closure_probe/audit.json. Preserve its hashed
sources and exact solutions. Next: quantify one-row Chernoff slack and seek
BCH-specific constraints beyond OA21. No user approval to substitute M25 has
arrived; the original goal remains active. All experiments remain sequential.

## Full-closure goal: M25 certified; original M22 remains open, 2026-09-04

The user set a goal to fully close the result under the random-inner model.
That goal remains ACTIVE pending the parameter choice: the established M22
target is still conditional, but a nearby M25 variant now has a complete
deterministic-spectrum certificate. No three-shell statistical assumptions
are used by the M25 result. Its full bound has approximately 41.90289 bits
and is exactly below 2^-41 at the same cutoff 209716 and 2^20 message bits.

Read RANDOM_INNER_M25_CLOSURE.md for the statement, proof, and interface audit.
The new checker re-derives all 396 BCH LP rows, checks the two published anchor
tables and generator containments, and removes every orbit-search lower bound.
The resulting caps are slightly weaker; a proved factor 128/127 safely accounts
for their effect on every reused higher-occupation bound. A marginal coupling
shows that M22 transfer bounds remain valid when memory increases to M25.

All 92 M25 Q1 shells were recomputed independently at 256-bit Arb precision.
The final read-only checker passed: 17 JSON dependencies, exact low-occupation
reaggregation, full tail coverage, and exact full sum below 2^-41.

    python -B code/audit_bch_closure_envelope.py --verify
    python -B code/certify_bch_m25_closure.py --verify
    python -B code/verify_bch_m25_closure.py

New receipts: generated/bch256_closure_deterministic_envelope.json,
generated/bch256_m25_unconditional_distance.json, and
generated/bch256_m25_independent_transfer_check.json. Their sources are frozen.
No original source or receipt was overwritten, and no new sampling was needed.

The current paper's general Theorem 3.1 and Section 6.4 route match our event
and route. Its equation (19) random convolution does NOT match RandomStepConv:
it retains trailing outputs and is invertible for every setup. Also, the
current finite certificate targets 10% distance, not the paper's 11% target.
No transfer to that convolution, RM2Sub, or a PRG implementation is claimed.

Binary64 diagnostics at 10% with deterministic caps give M22 32.716 bits,
M23 36.795, M24 39.801, and M25 41.903. These are bounds, not estimates of the
true bad-setup probability. M25 has a separate rigorous certificate above.

Next: ask whether to adopt M25 as the random-model endpoint or require the
original M22 parameters. Do not mark the user's full-closure goal complete by
silently substituting M25. The asynchronous question also asks whether the
user instead intends the draft's distinct convolution and 11% target.

## All-occupation random-inner goal completed: 2026-09-04

The q>=4 tail is now outward-certified below 2^-83, versus the required 2^-45.
All 8,189 remaining occupations are covered. A new deterministic cumulative
envelope uses squares of degree-at-most-seven Krawtchouk polynomials and the
existing OA15 identities. The reference argument uses monotonicity in the
averaged RandomStepConv law, not pointwise monotonicity of a fixed encoder.

The published k=71 anchor was checked directly against Table 7, its generator
containment in Q was verified, and its exact dual minimum is 16. Independent
checks passed for 68 polynomial squares, 5,397 prefix comparisons, 83 full-Arb
small-occupation recurrences, and 8,106 direct full-length matrix powers.

Combining all occupations gives an exact bound below 2^-40, with displayed
margin approximately 40.04769 bits, CONDITIONAL on the three tested shell caps.
The new statement explicitly defines the ideal RandomStepConv-M22 setup.
No deterministic proof of the three caps, RM2Sub transfer, or arbitrary PRG
replacement is claimed. No new shell sampling or adaptive retry was performed.

Read RANDOM_MODEL_FULL_DISTANCE_BOUND.md first. Final receipt:
generated/bch256_full_random_inner_conditional.json. Read-only verification:
python -B code/verify_bch_full_random_inner.py (passed, 30 JSON dependencies).
Preserve all frozen sources and receipts. Earlier status fields saying q>=4
is open or that only Q1-Q3 are covered are historical and superseded here.

Next recommended goal: audit the explicit conditional theorem against the
intended SPIN construction and package the proof for the paper. Deterministic
shell-cap proofs and transfer to RM2Sub remain separate research objectives.

## Larger goal completed: 2026-09-04

The full fixed-budget tests under code/run_higher_endpoint_certificate.py
completed successfully. Weight 40 had 436,628,033 zero-hit trials and 723 Python
audit checks. Weight 42 had 139,431,904 zero-hit trials and 427 Python checks.
Both full independent replays matched. The receipt, frozen manifest, fresh
tapes, and audit records are in generated/higher_endpoint_certificate_20260904.
Do not modify frozen sources, executables, or retained evidence.

Independent Arb arithmetic recomputed all 91 interior Q1 transfer coefficients
and included the previously omitted all-one shell, bounded below 2^-642.
The Q2 positive recurrence is now outward-certified at 73.4876764 bits using
only the existing deterministic shell caps. A direct capped-weight Q3
calculation also passed at 114.5938572 bits, without statistical caps.
The evidence-aware Q1 audit and exact Q1+Q2+Q3 aggregate passed. The combined
bound has 40.0476888685 bits (display approximation), leaving 3.2515045185%
of the 2^-40 budget. Occupations q>=4 remain open. The next sufficient target
is sum_(q=4)^8192 F_q<=2^-45; the exact residual exceeds that allowance.
See RANDOM_MODEL_NARROW_CERTIFICATE.md for the mathematical scope, and
generated/bch256_q1_full_arb_transfer.json and
generated/bch256_q2_positive_outward.json and
generated/bch256_q3_capped_outward.json for completed arithmetic receipts.
The joint claim that all three shell caps hold, accepted only if all three
predeclared tests pass, has false-accept error at most 2^-40 under ideal IID
sampling by an intersection-union argument. The older 2^-39 familywise bound
is conservative for this joint decision. See RANDOM_MODEL_NARROW_CERTIFICATE.md.

Final receipt: generated/bch256_q123_joint_evidence.json. The read-only command
python -B code/verify_bch_q123_evidence.py passed: all three contributions
were independently reaggregated and 23 JSON dependencies checked. Exact LP
dual witnesses and the Johnson certificates were also freshly rechecked in
generated/bch256_deterministic_caps_reaudit.json. No full-SPIN theorem or
RM2Sub transfer is claimed. The earlier entries below are historical milestones;
their prospective status and remaining-Q2/Q3 language are superseded here.

## Higher-shell samplers validated: 2026-09-04

The medium implementation milestone is complete. The new fast weight-40/42
C++ samplers and independent replay passed 4,096 Python fixtures, including
2,048 known endpoints and nine singular cases. Sequential 65,536-trial
preflights at both weights had zero endpoints and respectively 314 and 227
singular systems. Full replay matched all histograms, byte counts, and audit
records. Python full Gaussian elimination verified all 576 retained records.
No shell cap is established by these small runs.

See HIGHER_SHELL_PREFLIGHT.md and
generated/higher_endpoint_preflight_20260904/preflight.json. The new preflight
verifier and the old weight-38 certificate verifier both passed. Their hashed
source inputs and executables must now remain unchanged for reproducibility.

Next: add a separate full statistical driver, freeze its acceptance rules,
and run the prospective weight-40/42 experiments sequentially with fresh tapes
and independent replay. The full experiments have not started. Q2 outward
arithmetic and Q>=3 remain separate proof obligations.

## Higher-shell reduction and Q2 diagnostic: 2026-09-04

The weights 40/42 have a new exact endpoint-incidence reduction. Their locators
have form g(z^2)+z^37 h(z^2), with deg h=1 or 2. Fixing h(0) in the code's
32-element syndrome subspace leaves 20 or 22 unknown coefficients. Every
20- or 22-subset of a genuine endpoint's roots gives a nonsingular system;
the proof uses coprimality and an auxiliary polynomial with too many roots.
This yields an exact hit-probability formula, including the zero-syndrome
subcode, rather than assuming a random spectrum.

Exhaustive GF(16) analogues and 2,048 BCH positive-control reconstructions
passed. The proposed A40<=5e13 and A42<=5e15 zero-hit tests require respectively
436,628,033 and 139,431,904 IID trials at false-accept error 2^-41 each. These
full experiments have not run; neither cap is established. See
HIGHER_SHELL_ENDPOINT_REDUCTION.md and
generated/higher_shell_endpoint_reduction_checks.json. Combining those error
allocations with the existing A38 test would give error <=2^-39, not <=2^-40.
This confidence error is not the SPIN setup failure probability.

A separate fixed-code RandomStepConv Q2 diagnostic includes all ordered
weight pairs and the all-one shell. Independent row permutations make the
ordinary products A_a*A_b sufficient; no pair enumerator is needed here.
Nearest-binary64 upper-bound calculations give 73.7324 bits using the existing
LP-only spectrum envelope, 82.0810 with the accepted A38 target, and 87.5977
with the still-unproved A40/A42 targets. Thus Q2 is not an apparent numerical
obstruction. This is not an outward-rounded certificate. Q>=3 remains open.
See generated/bch256_q2_envelope_diagnostic.json and its retained source snapshot.
Its coarser tilt grid explains why its Q1 figures differ from the finer directed
application calculation; the newer diagnostic does not replace that calculation.

Next: implement and validate the fast weight-40/42 samplers with separate replay
before drawing their fixed-budget tapes; then certify Q2 with outward arithmetic
and attack Q>=3. Keep the statistical and deterministic routes distinct.

## Correction to the application reduction: 2026-09-04

The historical claim that weight 38 alone closes the 40-bit Q1 screen was
incorrect. The directed A_38 target 3,827,351,840,403 is an allocation from a
joint inflated low-shell model. It is not sufficient when combined with the
current other-shell LP caps. An exact rational aggregation of outward transfer
coefficients gives 36.8688 bits with this cap, or 36.9955 bits even if A_38
were zero. Weight 40 is the next dominant upper-bound obstruction.

This correction supersedes the statements below that call weight 38 the only
remaining application obligation. The standalone locator identities, exact LP
caps, and endpoint statistical test remain valid. Reproduce the budget audit
with `python code/audit_endpoint_application_budget.py` (python-flint required).
Its receipt is `generated/endpoint_application_budget_audit.json`.

## Direct endpoint test: 2026-09-04

The exact identity Pr[N_B=37] = L_37*C(37,18)/C(255,18) permits a direct
statistical test of the weight-38 cap, without estimating the sixth moment.
If the current application cap fails, the endpoint probability is at least
1.84395410537e-7. The predeclared 150,361,035-sample experiment completed with
zero endpoint hits and a full matching replay. All 407 retained audit vectors
also passed the independent Python checker. The exact false-accept bound is
at most 2^-40 under ideal IID bytes. Actual execution uses Windows CNG
pseudorandom bytes, retained in a 2.812 GB archive. The result is conditional
statistical acceptance, not an unconditional deterministic theorem.
The complete receipt is `generated/endpoint_certificate_20260904/certificate.json`.
The same sample gives a diagnostic sixth-moment ratio of 0.9753118802.

The existing five-million fixed-seed MT19937 sample is insufficient even
under the ideal IID interpretation: a just-violating shell would yield zero
hits with probability about 0.398. The RNG premise must also be audited before
claiming a statistical certificate. See ENDPOINT_SAMPLING_CERTIFICATE.md and
generated/endpoint_sampling_analysis.json. This route bounds A_38 directly;
it neither proves R_6<=6 nor covers higher SPIN occupations.

Next: attack weights 40 and 42. The conditional targets A_40<=5e13 and
A_42<=5e15, together with the accepted weight-38 cap and current higher-shell
bounds, yield 40.0476889 bits for the audited interior-shell sum. These two
new shell caps remain unproved. See ENDPOINT_APPLICATION_BUDGET_CORRECTION.md.

## Reproduced exactly

- The bundled anchor, MacWilliams, Wambach-word, and affine-stabilizer checks pass.
- The BCH generators have dimensions 131 and 123 and differ by
  `x^8 + x^5 + x^3 + x^2 + 1`.
- Multiplication by `x` modulo that factor has one orbit on all 255 nonzero
  quotient elements.
- The exact JSON model now enforces the Wambach lower bounds and exports the
  affine 1- and 2-design lattice conditions explicitly.

## Solver work

- An exactly equivalent, column- and row-scaled rational export has 130
  structural variables and 396 rows after replacing the OA moments by the
  independent even Krawtchouk basis.
- Project-local QSopt_ex and SoPlex binaries are unpacked under `.tools/`.
- The support/OA/Wambach core is feasible in rational SoPlex mode.
- Full-model floating and iterative-refinement statuses are not stable enough to
  certify either feasibility or infeasibility.  In particular, the apparent
  SoPlex infeasibility result has not been independently verified and must not be
  used as a mathematical claim.
- A full exact Z3 rational feasibility attempt remained unresolved after five
  minutes and was stopped.  Its status is `unknown`.
- The unscaled exact QSopt_ex objective run timed out after 120 seconds.  The
  scaled formulation solved in seconds.  Independent exact checkers verify
  primal feasibility, dual feasibility, and strong duality for all seven
  low-shell objectives.

## Calibrated low-shell search

- A fixed-width C++ scanner now enumerates up to five systematic generator
  rows over every inequivalent primitive-element representation.
- At radius four, affine closure exactly recovers the published minimum shells
  of the rate-half extended BCH codes of lengths 8, 32, and 128.
- At radius five it proves `A_38(P) >= 7507200`, hence
  `A_38(C) >= 912640`.
- Independent radius-five closures prove `A_40(Q) >= 12990720` and, after
  adding the discovered `P \ Q` orbits, `A_40(C) >= 18045952`.
- The same method proves `A_42(C) >= 42289664` and
  `A_44(C) >= 267281408`.  The latter directly probes the shell used by the
  current LP's pathological worst-case solution.
- Radius-five searches prove `A_46(C) >= 1247821056` and
  `A_48(C) >= 5525363968`.  A radius-four search proves
  `A_50(C) >= 1574477312`.
- These are rigorous lower bounds, not completeness or upper-bound claims.
  Every calibrated heuristic through weight 50 remains millions of times below
  its sufficient application cap.
- Exact rational Johnson-scheme Delsarte LPs prove
  `A_52(C) <= 2979928035058718446462766169` and
  `A_54(C) <= 47955805775385501045909339217`.  These are 4.35 and 4.44 bits
  tighter than the elementary constant-weight packing bounds.
- Combining those exact caps, shellwise packing through weight 128, and all
  seven calibrated low-shell estimates gives 45.799588 RandomStepConv bits.
  The low estimates can grow jointly by a factor of 4.47 million before the
  margin reaches 40 bits.
- The same exact Johnson LP at weight 38 gives only
  `A_38(C) <= 971390652389758748 < 2^59.753`, about 18 bits above the
  application cap.  Distance-only Delsarte bounds therefore do not settle the
  low shells; the next upper bound must use BCH linear or algebraic structure.
- A residual 2-design Johnson program improves this to
  `A_38(C) <= 128630323189940096 < 2^56.837`.
- The full exact BCH-sandwich programs prove
  `A_38(C) <= 773397229452928 < 2^49.459`; the remaining exact caps for weights
  40 through 50 have log2 values 51.903, 55.087, 59.989, 62.637, 65.606, and
  69.129.  The last two are within 0.83 and 0.21 bits of their simple
  application allocations.
- Substituting all seven exact caps leaves only 32.706 RandomStepConv bits.
  Weight 38 alone accounts for essentially the entire gap, so the remaining
  target is a 7.66-bit improvement in that shell rather than simultaneous
  progress on all seven shells.
- A native PB/XOR ApproxMCPB encoding exactly recovers fixed-pair counts for
  lengths 8 and 32, but stalls on the known length-128 count.  It is not yet a
  credible length-256 estimator.

## Next recommended reduction

Attack the normalized weight-37 locator-polynomial list implied by the BCH
Newton identities.  It is a structured Reed--Solomon coset-agreement problem;
if `L_37` denotes its list size, then `L_37=19*h_38/128`.  Bounding `L_37`
below `2^34.093` closes the current weight-38 target.
A directed-rounding transfer certificate now proves that a total incidence
factor below 6.098404530135757 is sufficient. After charging the
derivative-zero family, the corresponding smooth-locus budget is 6.069.
The clean standalone target remains an incidence factor at most 6. The exact
random-rank baseline is C(255,24)/256^6.
Every candidate has at most 37 agreements, since the agreement equation lifts
to a nonzero degree-37 locator polynomial.  Importance sampling five million
uniform 18-point interpolants estimates the 24-point incidence ratio at 1.032,
with a normal 95-percent upper endpoint of 1.136.  This directly supports the
needed constant-factor conjecture but is not a point-count proof.
An exact/sampled size ladder distinguishes this narrow claim from a broader
random-splitting claim. The exact \(q=32\) incidence ratios are at most 1.994,
and the sampled sixth-order ratio at \(q=128\) is 0.784. However, the exact
\(q=128\) fully split endpoint is enriched by a factor of 96.94. The proposed
heuristic therefore assumes only \(R_6\leq2\); a proof of \(R_6\leq6\) is
sufficient for the application.
The application threshold is no longer a floating-point proof obligation.
The exact integer sufficient cap is A_38(C) <= 3827351840403, and the lemma
R_6 <= 6 has 0.02347 bits of certified margin.

The ordinary Walsh ridge has also been isolated. On affine perturbations its
signed aggregate is -64, and exhaustive degree-at-most 1, 2, and 3 slices
have at most 4, 7, and 10 roots. Those slices therefore contribute no
18-base or 24-point incidences. However, the full character sum has a
six-dimensional Vandermonde kernel whose derivative denominators couple all
24 points, so the ordinary one-variable Walsh spectrum does not factor the
remaining problem.

The universal Reed--Solomon coset moments are now exhausted exactly. The
first 18 factorial moments and the 37-agreement endpoint permit an exact LP
extremizer with \(R_6=3.3072\cdot10^7\), supported on agreement counts
0 through 17 and 37. Generic MDS moment identities therefore cannot close the
proof.

The monomial-specific 24-point equations have been rewritten as
\(e_{24}h_{122}=1\) and \(h_{123}=\cdots=h_{127}=0\). The identity
\(H=EH^2\) halves all six indices through Frobenius squaring. This symmetric
form is the current highest-value exact reduction.

The binary Frobenius-fixed exceptional family is exhausted. All
\(2^{18}\) binary-coefficient interpolants contribute only
\(3.33\cdot10^{-9}\) to \(R_6\), even though the family contains ten endpoint
locators. The corresponding invariant-set stratum is highly enriched but
contributes only \(1.05\cdot10^{-17}\) globally. Four of its five affine
locator orbits are new. Their exact closures improve the rigorous floors to
\(A_{38}(P)\geq7{,}768{,}320\) and \(A_{38}(C)\geq944{,}384\).

The full \(\mathbb F_4\)-coefficient stratum is now exhausted as well. Among
267,516,561 invariant 24-sets, exactly 67,900 are admissible. They recover
40,712 interpolants with at least 24 agreements and exactly 34 endpoint
locators. The full stratum contributes \(1.134\cdot10^{-8}\) to \(R_6\), with
99.786 percent coming from the endpoints. The endpoints occupy 15 affine
orbits, including ten new orbits of size 65,280. Their exact closures improve
the rigorous floors to \(A_{38}(P)\geq8{,}421{,}120\) and
\(A_{38}(C)\geq1{,}023{,}744\).

The next nested stratum has been sampled with an unbiased invariant-set
estimator. Ten million invariant 18-point interpolations over
\(\mathbb F_{16}\) give ratio 1.00624 with standard error 0.00412. The normal
95-percent upper endpoint is 1.01432. The run also yields ten exact endpoint
polynomials in nine new affine orbits. Their closures add 587,520 words to
the rigorous \(P\)-shell floor. Thus
\(A_{38}(P)\geq9{,}008{,}640\) and
\(A_{38}(C)\geq1{,}095{,}168\).

The next recommended proof step is to exploit the nested
\(\mathbb F_2\subset\mathbb F_4\subset\mathbb F_{16}\subset\mathbb F_{256}\)
results in the Vandermonde-kernel character sum. All three structured strata
are harmless on their incidence scale. Endpoint concentration remains the
visible exception. The next useful split is therefore the complement of the
\(\mathbb F_{16}\)-coefficient locus, together with a separate endpoint term.
