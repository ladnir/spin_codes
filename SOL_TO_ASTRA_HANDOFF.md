# Sol to Astra handoff: BCH-256 spectrum evidence for SPIN

Updated: 2026-09-05

## Newest: original ideal-M22 goal CLOSED at 40.724976 bits

Read bundle/RANDOM_INNER_M22_CLOSURE.md, where bundle is
bch_spectrum_work/bch_spectrum_codex_bundle in this workspace.
The fixed BCH-derived [256,128] outer now meets k=20, rate 1/2, relative
distance > 0.1, with exact setup-failure upper bound U < 0.61 * 2^-40.
This is a computer-assisted theorem for ideal RandomStepConv-M22, not a
statistical shell conclusion. No variance estimate or A38 <= 10^13 premise.

Breakthrough: independently proved Q-dual distance >= 30 gives OA29, then
the exact 1163-row hull model bounds the paired weights 38/218, 40/216,
and 42/214 JOINTLY. The separate-cap result was only 39.721797 bits.
The joint witness plus all other shells and occupations gives 40.72497637925699.
Frozen joint folder: generated/shift_rank_oa29_joint. This is NOT an h38-cap
LP; never use the unit-h38 checker on it. The objective is rounded upward
to dyadics with 192 fractional bits and exactly normalized.

Final receipt: generated/bch256_m22_unconditional_closure_audit.json.
Run python -B code/verify_bch_m22_closure.py. Creation and read-only replay
both passed: 78 JSON evidence files, all 92 original M22 transfer coefficients
independently recomputed at 256-bit Arb precision, exact joint dual witness,
and exact all-occupation reaggregation. The cell-DP replay reproduced every
saved array; base/refined rank, OA21, signed-weight, and Fourier toy audits
also passed. Higher-occupation independent arithmetic receipts are retained
and hash-checked, not claimed rerun by the final checker.

The Q-rank proof covers 18 leaves after splitting pivot 15 on F(23)=0 or
nonzero. Both children have rank 29. Parity and affine extension yield
Q-dual >= 30; P-dual inherits it. No external dual-distance table is needed,
but the published exact endpoint spectra remain explicit mathematical inputs.

All computations are terminal; no live solver or verifier to resume.
All hashed sources/receipts are frozen. Preserve old failed/provisional work.
The P32 search did NOT prove 32; OA27 builder was never run. Neither is needed.

Scope: ideal refreshed-state M22, not RM2Sub, not a PRG, not equation (19)
of the draft, not 11% distance, and not a full protocol-security theorem.
Next: independent review and paper integration; practical-inner transfer is
separate. Older handoff entries below are historical and may say goal ACTIVE.

## Newest: independent Fourier certificates, CERTIFIED38.64383bits

Read bundle/BCH_SHIFT_RANK_CERTIFICATES.md. CERTIFIED cap90005732934389,
fullM22 margin38.64383138557605bits, goalACTIVE. No inaccessible table input
needed anymore: independent complete rankcases prove Qdual>=24/Pdual>=26.
code/certify_bch_shift_rank.py --verify rebuilds dualgenerators, zeros,
casecovers17Q/16P, parity, affineextensioninvariance. Completefolders
shift_rank_q24_complete andshift_rank_p26_complete; receipt
bch256_shift_rank_dual_distances.json. NEW model shift_rank_hull_probe
(196vars1157rows) byteidenticalexporttooldprovisional, sooldrationalwitness
reusedthenindependentlyaudited: audit_bch_extended_hull_cap.py
prepare_bch_shift_rank_probe oa21_hull_dual18_probe --verify.

Further proof: Qdual>=28, receiptbch256_shift_rank_q28.json. Replay
certify_bch_shift_rank_threshold.py39 shift_rank_q30_complete
bch256_shift_rank_q28 --verify (space after .py). Full17case coverminimum28:
pivot7rank29,pivot15rank28,pivot23rank29,othersrootrun>=30. Target30NOTproved.
P32probecomplete minimum28also, not yet independentlyaudited. The only
Qcasepreventing30is15; a binarysplit23is inprogress (current context handle).
Do notrunnewcomputationbefore checkinglivehandle. Allrunssequential.

Searchsourceprobe_bch_shift_rank.py: rankindependentsetrowsv_e=(alpha^ej)on
unknownsupport. StartpivotbwithF(b)!=0. Frobeniuse->2^ke andcommonshift
preserveindependence; ifalltransformedexponentsareknownFourierzeros, adding
b increasesrankby1. Saved(k,s)relativeatb means actualshift s+(1-2^k)b.
Verifiercheckseverycontainment. Binarycasecover firstnonzerocyclotomiccoset;
earliercosetszero. Forlatercases useexplicitBCHrootrunsthantoowidebeam.
Evenpunctureddualparityrounds23to24,25to26. Translationofzerocoordinateto
extensionparitypositiontransfersboundtoextendeddual. Proofinnewnote.

run_bch_shift_rank_cases.py isearlierinefficientallbeamrunner; partiallate
Qcaseshitlimits. complete_bch_shift_rank_cases.py addsrootruns andretained
witnesses BUThasa frozenbugwhenoldretainedBCHwitnessisbelownewtarget
(KeyErrorsteps), so DON'Tuseitforhigherthresholds. Correctednewsource
extend_bch_shift_rank_cover.py DIST SOURCEFOLDER NEWFOLDER TARGETRANK
[--width500 --seconds35]. Sourcecan'-'. Itretainssuccesses,coversallcases,
andreportsminimumranktruthfullyregardlesssearchtarget. Frozenq30_probe
has2casefilesonlyandnocover(summaryabsent); q30_completeisfull audited28.
Thresholdindependentauditor acceptsDIST FOLDER RECEIPTNAME [--verify].
probe_bch_shift_split.py SOURCEFOLDER PIVOT SPLIT TARGET NEWFOLDER tests
splitcoset nonzero/zero; eachleafindset independentlycheckablebutneedsfull
treecoverauditbeforeclaimingnewglobaldistance. CandidateLPbuilder
prepare_bch_shift_rank_oa27_probe.py existsbutnotrunyet; addsq/h24/26.

## Newest lead: provisional SchaubPlus moments; exact canonical orbits

Read bundle/BCH_AFFINE_CANONICAL_AND_SCHAUBPLUS_LEAD.md. CERTIFIED margin
still37.570616954093566/cap226659579388193. Goal ACTIVE, not closed.
New PROVISIONAL stronger moments give conditional38.64383138557605bits,
cap90005732934389. Do not promote until source/algorithm audited.

Web index of Ponchio-Sala, A lower bound on distance of cyclic codes,
2003-02-25, Table1,p16 reports Q0(dim123,d39)dual>=24, P0(dim131,d37)dual>=26
under SchaubPlus column. Primary PDF retrieval/visual table audit FAILED.
Source URL in new note; archive returned429; CiteSeer reset connection;
RICAM2006slides403; no local paper files created. 2015Kaida-Zheng reference11
identifies original as BCRI-07.ps. BCRI site unavailable. No independent
rank certificate reproduced. Provisional equations Kq22=Kh22=0 and
K(q+255h)24=0 via same affine extension argument. DO NOT impose Kh24=Kq24=0.
prepare_bch_schaubplus_provisional.py builds196vars1157rows and exact warm
solve in generated/schaubplus_provisional. audit_bch_schaubplus_provisional.py
[--verify] checks arithmetic and writes ONLY provisional_audit.json with
conditional keys, forces targetclosedfalse. This exact primal no longer
has the old weight38-only true-tail obstruction. Next prioritize auditing
new dual bounds, then jointweighted LP/fullJ and sharpenother40/42caps/tails.

New frozen canonical kernel bch_affine_canonical.cpp templates14/16/18/20/22.
For supportS, normalize all orderedpairs x!=y by(z+x)/(y+x), chooseleast
integerword, count minimizers=stabilizer, orbit65280/count. Independent
certify_bch_canonical_shell.py rebuilds all pairimages of everyoutput rep,
checks canonical/stabilizer/membership, totals against knownD shellcount.
check_bch_affine_canonical.py --verify crosschecksfullclosure/nontrivial
stabilizers and all97old18orbits. run_bch_affine_canonical.py W SEEDFOLDER
NEWFOLDER compiles and runs once. Outputcanonical.bin 36byte records:
32littlebyteword +4littlebytestabilizer, sortedintegerword.

Weight20 searches: translated3flattriples221336seeds and additions of
all3212592affine4flats to98known14/18orbitseeds10992seeds yield IDENTICAL
2040canonicalrecords/133171200Dwords, noneHQdual. Exact partial upper
832972800 onHQdual20, NOTzero. Not addedtoLP. Receipts
bch256_hull_dual_weight20_partial.json and _weight20_flat_partial.json;
canonicalfolders hull_weight20_canonical and hull_weight20_flat_canonical.
Oldcomplete18alsoindependentlycanonverified in hull_weight18_canonical,
receipt bch256_hull_dual_weight18_canonical.json. All frozen, jobssequential.

## Latest: exact weight18 shell; M22 margin37.57062bits

Read bundle/BCH_HULL_DUAL_WEIGHT18.md (bundle remains bch_spectrum_work/
bch_spectrum_codex_bundle). Best cap226659579388193, fullM22 margin
37.570616954093566bits. Goal ACTIVE, originalM22 unchanged.

New exact fact A18(HQdual)=0. Candidate search of triples of translated
3flats gives25492seeds; independent complete affine closure finds97
disjoint orbits and6332160Dwords, matching exact knownD18count. None inHQdual.
code/certify_bch_hull_affine_shell.py 18 hull_weight18_affine_triples
bch256_hull_dual_weight18_complete --verify replays. Older weight18 receipt
without _complete is partial49orbits only. NewLP builder
prepare_bch_hull_dual18_probe.py (196vars1154rows) replays all proofs.
Audit via audit_bch_extended_hull_cap.py prepare_bch_hull_dual18_probe
oa21_hull_dual16_count_probe --verify. All hashed files frozen.

New configurable bch_full_hull_feasibility_probe.py run/audit/verify
MODULE FOLDER proves FULLJ+14/16 and FULLJ+14/16/18 models insufficient.
Folders hull_full_exact16_feasibility and hull_full_exact18_feasibility;
modules prepare_bch_hull_dual16_count_probe and prepare_bch_hull_dual18_probe.
Exact rational witnesses with AUXILIARY h38>=5e12 satisfy1475original
intersection rows and all augmented rows, true-tail lower first moment
exceeds2^-40. This is NOT actual BCH spectrum/failure lower bound. Do not
aggregate auxiliary models as valid upper-bound models. New information
is required; old-cut iteration cannot close the goal. All jobs sequential.

## Latest: exact dual-hull shells14/16; M22 margin37.56826bits

Read BCH_HULL_DUAL_LOW_SHELLS.md. Best current cap227085449947170, fullM22
37.56825962339644bits. Goal ACTIVE, originalM22 parameters unchanged, no
statistical shell assumptions, no live processes. All runs sequential.

New exact code facts: HQdual has minimum distance EXACTLY16,
A14(HQdual)=0 and A16(HQdual)=16592. Let A=hull(L71), dimension69, with
known spectrum0mod4 part of L71. D=Adual is known [256,187,d14]; HQdual
subsetD. Its exact known counts areA14=65280 andA16=3408432.

Weight14 audit enumerates97155 linear3spaces of F2^8 via RREF. Indicator
syndromes againstA produce255 equal-syndrome pairs, each intersecting at0.
One pair's symmetric difference hasweight14; its affine orbit has65280
distinct words, exhaustingD's known shell. None is inHQdual. Affine
invariance ofA,HP,HQ underx->x+1 andx->2x is checked on full bases.

Weight16 audit enumerates200787 linear4spaces; all indicators inD. Their
16cosets give3212592 distinct affine4flat words. 1037 linear4spaces lie
inHQdual, giving16592 words. For the255 matching3space pairs, enumerate
all32x32 coset pairs: skew unions give195840 further distinctweight16
words, none inHQdual. Their affine span is7, so disjoint from4flatclass.
Totals match3408432 knownD shell, certifying completeness without assuming
a low-weight classification theorem. Both audits replayed successfully.

Frozen sources/receipts:
code/certify_bch_hull_dual_weight14.py --verify;
code/certify_bch_hull_dual_weight16.py --verify;
generated/bch256_hull_dual_weight14.json and bch256_hull_dual_weight16.json.
NewLPs: generated/oa21_hull_dual16_probe (A14zero,1152rows,196vars,
cap228602379051008,37.55989408948bits), and
generated/oa21_hull_dual16_count_probe (+A16exact,1153rows,bestcapabove).
Replay latest via audit_bch_extended_hull_cap.py
prepare_bch_hull_dual16_count_probe oa21_hull_dual16_probe --verify.
Prepare modules rebuild shell audits and carry source hashes; do NOT edit.
The latest SMALL model's exact primal still gives true-tail lower first
moment>2^-40. This says nothing yet about its fullJ intersection extension.

Earlier this turn: new generic bch_hull_cut_iteration.py runs exact
extension/Farkas/master iterations and mapsBAS row names by unchanged
normalized-row identities. Rounds2/3 completed, replaypassed, tinycapgains:
round2 cap228870861604159; round3 cap228870861596719. Folders:
generated/hull_cut_iterations/round_002 and round_003, each extension/,
farkas/, master/. `python -B code/bch_hull_cut_iteration.py verify 2` etc.
The source is now hashed/frozen. Warm master solves took29s and8s.

More importantly, bch_hull_intersection_feasibility.py found and checked
an exact feasible point for the FULL OLD intersection model, plus reciprocal
and cuts throughround2, with auxiliaryh38>=5e12. Its weight38 true-tail
lower contribution>2^-40 proves that OLD relaxation cannot close. All1475
raw original rows and1479 augmented raw rows checked (normalized LP row
counts differ). Folder generated/hull_intersection_feasibility; verify mode
replays without solver. Its auxiliarylowerbound is NOT a BCH assumption;
never aggregate it as a claimed upper bound. It predates newA14/A16 equations
and does NOT prove their full intersection model insufficient.

Next priorities:
1. Combine exactHQdual shell14/16 equations with fullJ model. Create a new
generic feasibility probe based on the frozen old helper, selecting newest
master as previous; same mapped-warm approach worked in25s before. Do not
alter old sourceglobals or edit hashed sources. Explicitly classify any
auxiliary h38 lower constraint as relaxation-only, not actualBCH.
2. Attack D weight18. Its known cardinality is6332160=97*65280. A proposed
search (UNTESTED) is syndrome-zero triples of8-point indicators. There are
97155 linear3spaces, exactly381 scalar orbits (stabilizer must divide7 and
255, hence1). Fix firstspace to381 representatives and scan all97155 second
spaces; syndromeXOR lookup for third is~37million pairs. Keepweight18
candidate words, then close affineorbits, comparing their total with6332160.
This candidate construction is NOT known exhaustive; matching total proves
completeness if achieved. Use efficient compiled/batched lookup if needed.
Represent orbit union by processedcandidate removal and canonical/digest
checks to avoid storing6.3million256-bit integers unnecessarily. Other
weight18 shapes may require nonconcurrent affine3flats. No weight18code yet.

A read-only multiplier-root search (probe_bch_hull_multiplier_bounds.py)
found no betterBCHdistance than40forHP/HQ and8for theirduals acrossall
invertible cyclicmultipliers. This probe is not the source ofd16; explicit
shell enumeration is. Avoid repeating that root-run search.

## Latest: intersection eliminated by exact Farkas cut; negligible numeric gain

Read BCH_HULL_INTERSECTION_CUTS.md. New best cap228870861644645 (only8959
improvement), fullM22 margin37.55841850570869. Goal ACTIVE, original40bit
target still open, same true-tail primal obstruction. No live processes.

J=HP intersect hull(L71), dim61, primitive anchor/J quotient factor0x169
verified exactly. New model229vars1329rows; exact fixedJ-coordinate reduction
213vars1297rows. QSopt90s/P256 and120s/P512 both timed out with no certificates.
Do not repeat those large solves blindly. All retained under
generated/oa21_hull_intersection_probe and oa21_hull_intersection_reduced.

Missing reciprocal inequality was added: s <=2^-131 K(q-h), since
HQ intersect Pdual=HP. New196var1148row model solved EXACTLY with unchanged
optimum. Receipt generated/oa21_hull_reciprocal_probe/audit.json.

Fixing its rational primal gives17var150row J-extension LP. Solver said
infeasible but saved only status, so a separate151var18row alternative LP
obtained explicit Farkas multipliers. Exact y>=0, y^T B<=0, y^T d>0 checks
passed. Replaying16 nonzero multipliers on the original un-fixed rows
eliminates unknownJ coordinates, giving one valid integer-coefficient cut.
New196var1149row model solved exactly and its certificate passed.

New folders: oa21_hull_fixed_extension (status-only infeasible),
oa21_hull_extension_farkas (explicit exact witness),
oa21_hull_benders_cut (new cap/certificate). All in generated/.
Key frozen scripts: prepare_bch_hull_intersection_probe.py,
prepare_bch_hull_intersection_reduced.py, prepare_bch_hull_reciprocal_probe.py,
prepare_bch_hull_fixed_extension.py, prepare_bch_hull_extension_farkas.py,
prepare_bch_hull_benders_cut.py. Each has guarded main and read-only build().
Do NOT edit frozen sources. New generic run_bch_exact_probe_once.py runs one
named folder once with --seconds / --precision, retains logs, uses h_38.lp
and h_38.sol filenames even when auxiliary objective is different.

Checks: audit_bch_hull_intersection_structure.py --verify; and
audit_bch_extended_hull_cap.py prepare_bch_hull_benders_cut oa21_hull_reciprocal_probe --verify.
Generic cap auditor reconstructs model from named module.build(), requires
module.FOLDER, checks exact primal/dual, and aggregates from named previous
folder/audit.json. It expects original h_38 objective, not auxiliary LPs.
All scripts now hashed; preserve. Research-table detour: hullBCHdelta61
dimension61 is notJ (intersection53); no equivalence under any invertible
cyclic multiplier, but arbitrary coordinate equivalence not ruled out.

Next: iterate fixed-primal extension tests and Farkas cuts at the new optimum.
The first cut only shaved~4e-11 relatively, so do not oversell. Seek a feasible
extension to certify exhaustion of the intersection relaxation. Reuse a
correctly mapped old simplex basis to speed196var re-solves: export_lp appends
16 OA rows, so a new cut shifts their c-index names; remap BAS row names using
normalized_rows(row['name']) identity before warm-starting. No warm start
has been attempted yet. Current fixed-extension and Farkas scripts hardcode
previous folders; create a new generalized engine rather than editing them.

## Latest: coupled hull spectra improve M22 to 37.55842 bits

Read BCH_HULL_SPECTRUM_CONSTRAINTS.md in the active bundle. Exact nullspaces
give doubly-even HP[256,85] subset HQ[256,93] subset Q; P intersect HQdual=Q.
Punctured generator quotient HP/HQ is primitive factor0x169. Use spectra
r=HP and s=common nonzero HQ/HP coset; HQ=r+255s, supported at multiples4.
Signed Fourier inequalities plus containments solved exactly (1017rows,
196variables), cap466844536111807, full36.62779bits.

New exact anchor: L71 has polar rank2, hull dimension69, signed sum -2^70.
Its hull is exactly its weight0mod4 words, so its spectrum is known from
published L71_HALF. This hull lies in HQ. Adding primal/dual containment
constraints gives1115rows and proves A38<=228870861653604, full M22 bound
37.55841850566503bits. Original40bit target remains open; feasible rational
primal still exceeds2^-40 using the true-tail LOWER coefficient. This is a
limitation of the relaxation, not actual code failure probability.

All checks passed and replayed sequentially:
code/audit_bch_hulls.py --verify; code/check_bch_hull_fourier.py --verify;
code/audit_bch_hull_cap.py --verify; code/audit_bch_hull_anchor_cap.py --verify.
New receipts: generated/bch256_hull_structure.json,
generated/bch256_hull_fourier_toy.json (60 exhaustive small codes),
generated/oa21_hull_probe/audit.json, generated/oa21_hull_anchor_probe/audit.json.
Models, exact solutions and logs are adjacent. Both solver runs terminal,
about36seconds each, no parallel benchmarks. Sources are now frozen by hashes.

Next: J=R_L intersect HP has dimension61, and J subset R_L subset HQ,
J subset HP subset HQ, R_L+HP=HQ. Introduce J spectrum t and exploit known
R_L=t+255u with cyclic quotient symmetry (must verify quotient first),
t<=r and u<=s, plus Jdual containment/difference constraints. The existing
structural receipt records L71 intersect HP dimension61 and L71 intersect
HQ dimension69, which equals R_L. No J model yet. Do not assume U187dual
lies in HP/HQ: exact read-only check found intersection dimension53, not69.
Goal ACTIVE; original M22 parameters unchanged. No M25 substitution approved.

## Latest, 2026-09-05: small exact improvement; remaining cap still open

Read BCH_LOCAL_SLACK_AND_MOD4.md in the bundle. Quadratic elimination proves
signed weight sums S(P)=S(Q)=-2^108 and S(H)=0, with polar radical dimensions
85 and 93. OA21 plus these two identities solved exactly: new A38 upper cap
674685258232216; full M22 bound 36.12681 bits. Still need A38<=10^13 for the
40.14372-bit conditional theorem. No goal closure.

All 20 independent local slack LPs have exact zero-cost primal witnesses.
Thus that correction of the fixed old dual cannot help; this does not refute
the fully coupled split formulation. Verifier uses no LP optimizer.

New split attempts: balanced QSopt_ex (120s), balanced SoPlex exact (120s),
cleaned SoPlex 100-digit exact (120s) all ended at time limits, no rational
solutions. All handles terminal. SoPlex exit0 is NOT success; logs say timeout.
Tools extracted locally under generated/split38_soplex_tools, not installed.
Models/logs under split38_balanced_probe and split38_clean_probe.

Checks: code/verify_bch_local_slack.py; code/audit_bch_quadratic_sums.py --verify;
code/audit_bch_mod4_cap.py --verify. Receipts and scripts are frozen by hashes.
Next: radical/hull structure or genuinely coupled split constraints. Goal
ACTIVE, M22 unchanged, experiments sequential. No model substitution approved.

## Latest exact facts: code-specific shortening; split LP inconclusive

Read BCH_SPLIT38_CONSTRAINTS.md in the active bundle. The Wambach support S
has size38 and Q projects onto S with rank38. Its shortening on the remaining
218 positions has dimension85 and dual minimum>=7. This was certified by
exhaustive integer XOR/sort checks through 1,703,016 triples, not sampling.
Thus each fixed inside pattern has 2^85 extensions and its outside fiber
has six exact OA moments. Receipt: generated/bch256_wambach_shortening.json.
Replay and the independent RM(1,3) split-identity toy check passed.

New models: generated/split38_probe (1696 variables/720 rows) and
generated/split38_local_probe (1656 variables/633 rows). All three solver
attempts are terminal TIME_LIMIT_REACHED, with no solution files and no new
A38 cap. Read generated/bch_split38_attempts.json before starting more work.
The exact solution checker exists but has not run; it is not evidence of a
bound. At split center (19,109), the variable denotes half the actual count.

Goal ACTIVE. M22 still has 36.11873 deterministic bits and closes conditionally
at 40.14372 if A38<=10^13. Next: stabilize this code-specific relaxation or
derive smaller exact inequalities from the new local moments. Preserve all
hashed sources. Experiments sequential; no model substitution authorized.

## Newest: M22 needs only A38<=10^13

Read M22_ONE_SHELL_REDUCTION.md in the active bundle. Directed cell-count
brackets on active time give true one-row tail bounds at weights 38,40,42.
The fully deterministic M22 bound improves to 36.11873 bits. With only the
additional A38<=10^13 assumption, the full bound has 40.14372 bits. No
statistical assumptions at weights 40 or 42 are needed. The exact largest
sufficient A38 cap is 14,584,006,710,777; neither cap has been proved.

An exact feasible OA21 primal spectrum and a certified true-tail lower bound
show that even exact inner tails cannot make this relaxation suffice.
An independent integer-truncated recurrence confirms the obstruction.
Toy checks enumerate 6,400 full bit-level paths. Directed full replay passed.

Receipts: generated/bch256_q1_activity_cells_outward.json and
generated/bch256_q1_cell_budget.json. Check code/certify_bch_q1_activity_cells.py
--verify and code/audit_bch_q1_cell_budget.py --verify sequentially.
Goal ACTIVE. Next: focus deterministic BCH work on A38<=10^13, not on all
three shells. No user approval to replace M22 by M25 has arrived.

## Latest refinement: OA21 exact bounds; M22 still open

Read BCH_OA21_REFINEMENT.md and generated/oa21_closure_probe/audit.json in the
active bundle. Published dual minimum distance 22 for BCH(255,123,39), plus
the extension argument, supplies OA21. New exact LPs give caps A38 <=
678661807513927, A40 <= 3350512214809133, A42 <= 27886467650792183.
No statistical caps, orbit lower bounds, or lattice rounding are used.
Full deterministic M22 margin is about 32.90151 bits. Exact primal witnesses
also show these constraints and fixed transfer coefficients cannot prove 40
bits by joint LP optimization. This is not an obstruction for the actual code.
Verify with python -B code/certify_bch_oa21_probe.py --verify.

Goal remains ACTIVE, with no approval to replace M22 by M25. Next useful work:
quantify one-row Chernoff slack and investigate BCH-specific constraints beyond
OA21. Do not run two substantial experiments simultaneously. The primary
1996 author PDF is retained and hashed; its rank certificate was not rerun.

## Latest: deterministic M25 alternative; M22 full-closure goal still active

Read RANDOM_INNER_M25_CLOSURE.md in the active bundle first. The user asked to
set a goal to fully close the random-inner result. A deterministic-spectrum
certificate now proves the full distance bound at memory 25 instead of 22,
with 2^20 message bits, 2^21 output bits, cutoff 209716, and an exact bound
below 2^-41 (about 41.90289 bits). No statistical shell assumptions remain in
this M25 variant. The original M22 target is still conditional.

The goal remains ACTIVE: ask whether M25 is an acceptable random-model
endpoint rather than silently changing the parameters and marking completion.
An asynchronous question also asks whether the user intended the draft's
different convolution and 11% target instead of our M22/10% model.

New mathematical work: independently reconstructed 396 BCH LP rows, verified
both published anchor spectra and generator containments, removed every search
orbit lower bound and weight-38 lattice rounding, and derived marginal
monotonicity with increasing inner memory. The weaker spectrum bounds cost
at most 128/127 on the reused complete higher-occupation sum. All 92 M25 Q1
coefficients were independently checked using normalized Arb recurrences.

New final receipt: generated/bch256_m25_unconditional_distance.json.
Read-only checker: python -B code/verify_bch_m25_closure.py (passed, 17 JSON
dependencies, exact sum below 2^-41). Recompute independent Q1 arithmetic:
python -B code/certify_bch_m25_closure.py --verify. Re-derive BCH constraints:
python -B code/audit_bch_closure_envelope.py --verify. Preserve these sources.

Application audit: the general Theorem 3.1 event and Section 6.4 route match.
Equation (19) random convolution does NOT match RandomStepConv; its state is
trailing outputs and it is invertible for every setup. Our ideal full-matrix
inner is different. Our finite target is 10%, not the draft's 11% theorem.
No paper edits, RM2Sub transfer, PRG argument, or performance claims were made.

Next recommendation: adopt the proved M25 variant if that parameter change is
acceptable, then treat M22 tightening as an optimization problem; otherwise
continue the open deterministic low-shell/transfer work at M22.

## Previous: the M22 all-occupation bound is complete, conditionally

Read RANDOM_MODEL_FULL_DISTANCE_BOUND.md in the active bundle first.
The q>=4 tail is outward-certified below 2^-83, covering all 8,189 occupations.
It uses deterministic cumulative moment envelopes, the published OA15 anchor,
and optimized coefficient bounds. Independent full-Arb recurrences and direct
matrix powers checked every occupation; all prefix and polynomial checks passed.

The complete first-moment bound is below 2^-40 (approximately 40.04769 bits),
CONDITIONAL on A38<=3,827,351,840,403, A40<=5e13, and A42<=5e15.
The three statistical tests remain the evidence for those hypotheses; their
joint ideal-IID false-accept bound is 2^-40, with the CNG qualification.
The theorem concerns the explicitly defined ideal RandomStepConv-M22 model.
No deterministic spectrum theorem, RM2Sub transfer, or PRG replacement is claimed.

Receipt: generated/bch256_full_random_inner_conditional.json in the active bundle.
Read-only checker: python -B code/verify_bch_full_random_inner.py (passed).
Next: audit this conditional theorem against the intended SPIN interface and
package the proof for the paper. Preserve frozen evidence. All earlier claims
that q>=4 remains open are historical and superseded by this result.

## Previous completed larger goal: shell tests and Q1-Q3

Read RANDOM_MODEL_NARROW_CERTIFICATE.md in the active bundle first.
The full A40<=5e13 and A42<=5e15 tests passed: 436,628,033 and 139,431,904
zero-hit trials, matching full independent replay, and 723 plus 427 Python
audit records. Together with the prior A38 test, all three caps have conditional
statistical acceptance. The joint all-tests-must-pass claim has false-accept
error at most 2^-40 under ideal IID bytes; actual CNG requires the stated
randomness qualification. This is not a deterministic spectrum theorem.

Independent outward computations now cover Q1, Q2, and Q3, including the
all-one shell. Their exact combined upper bound has 40.0476888685 bits as a
display approximation. It leaves 3.2515045185% of the 2^-40 target budget.
The next sufficient proof target is sum_(q=4)^8192 F_q<=2^-45. Q>=4 remains
unproved, so no full-SPIN 40-bit theorem is claimed. RM2Sub transfer remains
separate; this task keeps RandomStepConv-M22 as requested.

Final receipt: generated/bch256_q123_joint_evidence.json in the active bundle.
Read-only verification: python -B code/verify_bch_q123_evidence.py (passed).
The endpoint tests and the older A38 verifier were also fully rechecked.
Preserve frozen inputs and tapes. The entries below describe earlier milestones;
their statements that A40/A42 have not run or Q2/Q3 remain open are superseded.

## Earlier milestones (historical)

Latest implementation milestone: HIGHER_SHELL_PREFLIGHT.md records validated
fast weight-40/42 samplers and independent replay. Both paths passed 4,096
Python fixtures; two sequential 65,536-trial preflights passed replay and all
576 Python audit checks. These runs establish no shell cap. Their frozen
receipt is generated/higher_endpoint_preflight_20260904/preflight.json in the
active bundle. Next add a new full statistical driver without modifying frozen
sources, then execute the full fixed-budget experiments with fresh tapes.

Latest continuation: HIGHER_SHELL_ENDPOINT_REDUCTION.md in the active bundle
derives exact direct sampling identities for A40 and A42, with a nonsingularity
proof and exhaustive small-field checks. Prospective zero-hit budgets are
436,628,033 and 139,431,904 trials at error 2^-41 each. Neither experiment has
run and neither cap is established. The separate Q2 binary64 diagnostic gives
73.7324 bits even with existing LP-only caps; outward rounding and Q>=3 remain
open. See the new top section of PROGRESS.md for receipts and next steps.

## Astra correction and continuation, 2026-09-04

The historical application claim in this handoff is too strong. The cap
A_38(C)<=3,827,351,840,403 is an allocation from a joint low-shell envelope,
not a sufficient standalone condition for the full 40-bit Q1 screen. Combining
it with the current other-shell LP caps gives only 36.8688 directed bits.
Even setting A_38 to zero leaves 36.9955 bits. Weight 40 is the next dominant
upper-bound obstruction. This correction supersedes every later assertion
that weight 38 alone closes the application.

A direct endpoint sampling test has since been implemented. Its exact identity
is Pr[N_B=37]=L_37*C(37,18)/C(255,18). The fixed 150,361,035-sample experiment
and full independent replay completed with zero endpoint hits. All 407 retained
audit vectors passed an independent Python checker. Results are under the active
bundle's generated/endpoint_certificate_20260904 directory; certificate.json
records conditional statistical acceptance. Its 2^-40 false-accept bound assumes ideal independent bytes;
the execution uses retained Windows CNG pseudorandom bytes.

Read ENDPOINT_SAMPLING_CERTIFICATE.md, PROGRESS.md, and
generated/endpoint_application_budget_audit.json in the active bundle first.
The deterministic BCH upper-bound problem and the full application remain open.

Repository root:

    C:\Users\peter\.codex\worktrees\ba80\permute_conv

Active BCH work folder:

    bch_spectrum_work\bch_spectrum_codex_bundle

This document transfers the current working state. It is not a chronological
summary. Numerical claims below were reconciled against the repository and
the generated receipts on 2026-09-04.

## 1. Objective

### Ultimate goal

The application is a finite SPIN construction that needs a fast rate-half
outer constituent. The candidate studied here is the extended BCH-derived
binary code

\[
 C=[256,128,38].
\]

The engineering conjecture is that this BCH constituent preserves the strong
distance margin predicted by random-like spectrum models. The complete
spectrum of \(C\) is unknown.

The user does not require the complete spectrum if a narrower certificate
controls the application functional. The current proof-backed alternative is
Expanded Accumulant. The user reports that Expanded Accumulant is about four
milliseconds, or a factor of 1.4, slower than the BCH path.

### Current subproblem

Under the temporary RandomStepConv-M22 transfer, weight 38 is the only
material unresolved shell. The exact application-side reduction shows that
the following integer cap is sufficient:

\[
 A_{38}(C)\leq 3{,}827{,}351{,}840{,}403.
\]

The BCH algebra converts this shell question into a point-count problem for
degree-18 polynomials over \(\mathbb F_{256}\). The clean proof target is the
sixth incidence ratio

\[
 R_6\leq 6.
\]

The current work pursues two tracks:

1. collect tightly matched evidence for the exact functional \(R_6\);
2. develop a narrow proof using BCH and monomial-specific algebra.

### Success criteria

For the temporary RandomStepConv model, any of the following closes the
weight-38 application obligation:

- prove \(R_6\leq 6\);
- prove \(A_{38}(C)\leq3{,}827{,}351{,}840{,}403\);
- prove the normalized locator-list bound \(L_{37}<2^{34.0932149}\).

The final SPIN application is expected to use RM2Sub rather than
RandomStepConv. The application threshold must therefore be recomputed after
the final RM2Sub transfer artifact arrives. The user reported that RM2Sub is
at least as good as random convolution near \(m\approx20\), but that result is
not present in this folder.

## 2. Current state

### Repository state

- Git is at detached commit
  \(0a743a59f344db8f7426213bc1a562076c68e930\),
  whose subject is “Promote EBCH128 parity theorem.”
- The complete folder bch_spectrum_work is untracked.
- The root tmp folder is also untracked and contains unrelated or shared
  scratch state. Do not delete or normalize it.
- None of the BCH continuation work described here is committed.
- The aggregate BCH verifier passed on 2026-09-04:

      cd bch_spectrum_work\bch_spectrum_codex_bundle
      python code\check_all.py

The verifier checks the compact receipts. It does not rerun the 267-million
candidate GF(4) enumeration or the ten-million-sample GF(16) experiment.

### BCH chain

The structural codes are

\[
 Q=\operatorname{EBCH}(256,123,40),\qquad
 P=\operatorname{EBCH}(256,131,38),\qquad
 Q\subset C\subset P,\quad\dim C=128.
\]

The quotient \(P/Q\) has dimension eight. Cyclic shift is transitive on its
255 nonzero elements. Every nonzero coset of \(Q\) in \(P\) therefore has the
same weight distribution. Let \(H_w\) denote one such coset multiplicity.
Then

\[
 A_w(P)=A_w(Q)+255H_w,\qquad
 A_w(C)=A_w(Q)+31H_w.
\]

The spectrum is independent of which five-dimensional subspace of \(P/Q\)
defines \(C\).

### Application-side threshold

The directed certificate in
generated/random_inner_threshold_outward.json gives:

\[
\begin{aligned}
A_{38}(C)&\leq3{,}827{,}351{,}840{,}403,\\
R_6&<6.098404530135757,\\
R_6\leq6&\quad\text{suffices with }0.02347\text{ bits of margin}.
\end{aligned}
\]

After charging the derivative-zero family in full, the smooth-locus budget is
6.069270771430856. Use this directed value. An older binary64-derived receipt
contains a slightly different final decimal.

### Best exact lower inventory

The radius-five Wambach search originally certified 115 affine weight-38
orbits and

\[
 A_{38}(P)\geq7{,}507{,}200.
\]

Nested Frobenius searches subsequently found 23 new affine orbits:

- four from the exact GF(2) stratum;
- ten from the exact GF(4) stratum;
- nine from the sampled GF(16) stratum.

Every new orbit has trivial affine stabilizer and size 65,280. The current
exact lower floors are

\[
 A_{38}(P)\geq9{,}008{,}640,\qquad
 H_{38}\geq35{,}328,\qquad
 A_{38}(C)\geq1{,}095{,}168.
\]

These are lower bounds. They do not imply completeness and do not address the
needed upper bound.

The 138 certified affine orbits generate 5,244 normalized locator
polynomials. Their degrees are:

- degree 16: 7;
- degree 17: 31;
- degree 18: 5,206.

Every one has 37 agreements. All 20,976 sampled admissible 24-subsets used in
the current Jacobian receipt have rank six.

### Present confidence

The heuristic evidence for the narrow incidence claim is strong. The direct
five-million-sample estimator gives \(R_6=1.03172\), with a normal 95-percent
upper endpoint of 1.13608. Exact GF(2) and GF(4) strata are negligible in the
global sum. The GF(16)-invariant estimator is also near one.

The result is not yet a proof. The main unresolved risk is an exceptional
component outside the GF(16)-coefficient locus, especially at the
37-agreement endpoint.

## 3. Key conclusions established so far

### Proved or exactly verified

1. The BCH quotient relations above are exact. The choice of intermediate
   five-dimensional quotient subspace does not change the spectrum of \(C\).

2. Every normalized candidate has at most 37 agreements. If
   \(g(u)=u^{146}\) at \(N\) nonzero points, those points lift to roots of the
   nonzero degree-37 polynomial

   \[
   z^{37}+g(z^2).
   \]

3. For a uniform 18-set \(B\subset\mathbb F_{256}^{*}\), define the unique
   polynomial \(g_B\) by

   \[
   \deg g_B\leq18,\qquad g_B(0)=1,\qquad
   g_B(u)=u^{146}\quad(u\in B).
   \]

   Set

   \[
   N_B:=\#\{u\neq0:g_B(u)=u^{146}\},\qquad E_B:=N_B-18.
   \]

   Then the exact double-counting identity is

   \[
   R_6:=
   \frac{\mathbb E_B\binom{E_B}{6}}
        {\binom{237}{6}/256^6}
   =
   \frac{I_{24}}{\binom{255}{24}/256^6}.
   \]

4. The normalized endpoint list obeys

   \[
   L_{37}=\frac{A_{37}(P_{\rm punctured})}{255}
         =\frac{19H_{38}}{128},
   \qquad
   A_{38}(C)=\frac{3968L_{37}}{19}.
   \]

5. If \(g'\neq0\), the six incidence equations have Jacobian rank six.
   After scaling, six suitable columns form a Vandermonde matrix.

6. If \(g'=0\), then \(g=h^2\) with \(\deg h\leq9\). A nine-point count gives

   \[
   L_{37}^{\rm deg}\leq87{,}550{,}900<2^{26.384}.
   \]

   This family is more than seven bits below the locator target.

7. For a 24-set \(T\), let \(e_i(T)\) and \(h_i(T)\) denote elementary and
   complete homogeneous symmetric functions. Admissibility is equivalent to

   \[
   e_{24}(T)h_{122}(T)=1,\qquad
   h_{123}(T)=\cdots=h_{127}(T)=0.
   \]

   In characteristic two,

   \[
   h_n=\sum_{\substack{0\leq i\leq24\\i\equiv n\pmod2}}
       e_i h_{(n-i)/2}^{\,2}.
   \]

8. The Walsh spectrum of \(x^{146}\) is exact. Its values are
   \(\{-32,-16,0,16,32,64\}\). The value 64 occurs exactly at
   \((a,b)=(a,a^{131})\), for \(a\neq0\). All other absolute values are at
   most 32. The differential uniformity is six.

9. Perturbations of degree at most one, two, and three have at most 4, 7, and
   10 roots. They contribute nothing to the 18-base or 24-point incidence
   counts. The exceptional affine Walsh ridge has signed aggregate \(-64\).

10. The complete GF(2)-coefficient family is exact:

    - \(2^{18}=262{,}144\) polynomials;
    - 82 have at least 24 agreements;
    - 10 have 37 agreements;
    - no polynomial has 29 through 36 agreements;
    - total global contribution \(3.3276472\cdot10^{-9}\).

11. The complete GF(4)-coefficient family is exact:

    - 267,516,561 invariant 24-sets;
    - 67,900 admissible invariant sets;
    - 40,712 distinct polynomials with at least 24 agreements;
    - invariant-set hit ratio 1.0396306 against six independent GF(4)
      equations;
    - full GF(4)-coefficient contribution
      \(1.13382774\cdot10^{-8}\);
    - exactly 34 endpoint polynomials;
    - endpoints supply 99.7856 percent of the GF(4) contribution;
    - 15 affine endpoint orbits, including ten new ones.

12. The nine affine orbits obtained from the GF(16) sample are exact
    codeword witnesses even though the containing sample is statistical.
    Their affine closures add 587,520 exact words to \(A_{38}(P)\).

13. The exact rational BCH-sandwich LP gives

    \[
    A_{38}(C)\leq773{,}397{,}229{,}452{,}928.
    \]

    This is still 7.6587 bits above the application cap.

### Experimentally supported

1. Five million independent uniform 18-point bases over
   \(\mathbb F_{256}^{*}\) give

   \[
   R_6=1.0317155.
   \]

   The sample-variance normal 95-percent upper endpoint is 1.1360727.

2. Degree-18 interpolants contribute 1.02708 to the global sampled ratio.
   Degree-17 interpolants contribute only 0.004632.

3. Conditioning on the affine Walsh-ridge image changes the sampled
   conditional ratio from 1.02383 outside to 1.06338 inside. The normal
   95-percent upper endpoint inside is 1.29456.

4. Ten million invariant 18-point bases for GF(16) give an unbiased
   invariant-set estimate

   \[
   R_{6,\mathrm{GF}(16)\text{-inv}}=1.0062402
   \]

   with standard error 0.0041180 and normal 95-percent upper endpoint
   1.0143115.

5. The GF(16) run sampled ten exact endpoint polynomials in nine new affine
   orbits. Their closures contain 64 new GF(16)-coefficient locators. The
   certified GF(16) inventory is now 190 normalized locators in 30 affine
   orbits. This inventory remains a lower bound.

6. The size ladder supports the narrow incidence model:

   - \(q=8\): exact ratio 1.82857 at order two;
   - \(q=32\): exact ratio 1.99381 at order four;
   - \(q=128\): sampled order-six ratio 0.784;
   - \(q=256\): sampled order-six ratio 1.032.

   The exact \(q=128\) fully split endpoint is enriched by 96.94. Endpoint
   behavior is therefore not random-like even when the lower incidence moment
   is close to its reference.

### Plausible hypotheses

- The working heuristic is \(R_6\leq2\).
- The minimum required proof statement is \(R_6\leq6\).
- The smooth component appears random-rank-like after derivative-zero and
  proper-subfield components are isolated.
- Larger BCH constituents likely continue the favorable smaller-code trend.
  The present evidence supports this engineering choice but does not certify
  it.

### Speculation

- A signed character-sum argument may control the complement of the
  GF(16)-coefficient locus.
- Endpoint concentration may admit a separate orbit or locator bound sharper
  than a generic smooth-variety estimate.
- More GF(16) endpoint orbits almost certainly remain undiscovered. The
  ten-million sample was not designed as a completeness proof.

## 4. Current approach

### Application reduction

The downstream transfer consumes low-shell multiplicities, but its current
failure is effectively one-dimensional. Exact caps for weights 38 through 50
leave 32.706 RandomStepConv bits. Weight 38 accounts for essentially the
entire deficit. This justifies attacking \(A_{38}(C)\) alone.

The directed application certificate converts the weight-38 cap into the
incidence target \(R_6\leq6\). This reduction avoids computing the whole
spectrum.

### Locator reduction

The BCH Newton identities identify normalized weight-37 punctured locators
with polynomials \(g\) of degree at most 18 satisfying \(g(0)=1\) and agreeing
with \(u^{146}\) on 37 nonzero points.

The proof target is a 24-point incidence count. The sixth factorial moment of
the extra agreements is exactly the application functional.

### Smooth/exceptional split

The derivative-zero family is bounded separately and charged in full. The
remaining nondegenerate locus has the expected codimension six.

Proper subfields form explicit exceptional strata:

\[
\mathbb F_2\subset\mathbb F_4\subset\mathbb F_{16}
\subset\mathbb F_{256}.
\]

The GF(2) and GF(4) strata are exhausted exactly. The GF(16) stratum is too
large for exhaustion, so it uses invariant-base importance sampling.

### Character-sum form

Let \(\Phi(X)\in\mathbb F_{256}^6\) contain the constant and five
high-coefficient defects for an ordered 24-tuple \(X\) of distinct nonzero
points. Character orthogonality gives

\[
24!I_{24}=256^{-6}\sum_{a\in\mathbb F_{256}^6}S_a,
\qquad S_0=(255)_{24}.
\]

The target follows from

\[
\sum_{a\neq0}S_a\leq5(255)_{24}.
\]

An equivalent coefficient-first expansion uses the six-dimensional kernel

\[
K_X=\left\{\lambda:
\sum_i\lambda_i x_i^j=0,\ 1\leq j\leq18\right\}.
\]

Every \(\lambda\in K_X\) has the unique form

\[
\lambda_i=\frac{h(x_i)}{x_iP_X'(x_i)},\qquad\deg h\leq5.
\]

This representation exposes the key obstacle: \(P_X'(x_i)^{-1}\) couples all
24 points. Ordinary one-variable Walsh coefficients do not factor the sum.

### Implementation structure

The exact GF(4) kernel stores coefficients in GF(4), computes only the
complete-homogeneous indices needed by the Frobenius-halved equations, and
processes 32 candidates per AVX2 batch. It avoids generic abstractions in the
hot loop.

The GF(16) sampler samples invariant 18-point bases. Frobenius invariance and
interpolation uniqueness force the resulting polynomial coefficients into
GF(16). Each agreeing invariant six-point extension \(T\setminus B\) receives
weight \(1/m(T)\), where \(m(T)\) counts invariant 18-subsets of \(T\). This
weight makes the invariant 24-set estimator unbiased.

## 5. Important files and code locations

Paths in this section are relative to
bch_spectrum_work\bch_spectrum_codex_bundle.

### Read first

| File | Role | Status |
|---|---|---|
| PROGRESS.md | Compact current BCH status and next reduction | Current |
| LOCATOR_INCIDENCE_DUAL_TRACK.md | Central definitions, exact reductions, heuristic track, proof track | Authoritative research ledger |
| SPIN_HEURISTIC_EVIDENCE.md | Application interpretation, curve evidence, stress tests, RandomStepConv/RM2Sub boundary | Current but transfer remains provisional |
| README.md | Bundle map and first commands | Current |
| handoff\bch_spectrum_codex_handoff.md | Original BCH-chain and endpoint-spectrum handoff | Authoritative for initial algebra; predates recent incidence work |

### Core verification

| File | Important entry point | Role |
|---|---|---|
| code\check_all.py | main, line 31 | Aggregate receipt checker; passed |
| code\bch_quotient.py | verify_quotient_algebra | BCH dimensions and quotient action |
| code\spectrum_exact.py | MacWilliams utilities | Exact anchor transformations |
| code\verify_scaled_rational_solution.py | verify | Exact primal/dual LP checker |
| code\certify_random_inner_threshold_outward.py | main, line 229 | Directed application threshold |
| code\locator_incidence_bounds.py | main, line 17 | Locator relations, Jacobian theorem, derivative-zero bound |
| code\verify_symmetric_incidence.py | complete_homogeneous, monomial_remainder | Symmetric incidence identities |
| code\analyze_locator_incidence.py | incidence_jacobian_rank, normalized_locators | Expands certified affine orbits and checks sampled ranks |

### Exact exceptional-family work

| File | Role |
|---|---|
| code\analyze_binary_146_family.py | Exhausts all GF(2)-coefficient polynomials |
| code\analyze_frobenius_invariant_incidence.py | Exhausts \(x\mapsto x^2\)-invariant 24-sets |
| code\sample_frobenius4_incidence.cpp | AVX2 sampler and exact enumerator for \(x\mapsto x^4\)-invariant sets |
| code\analyze_frobenius4_invariant_incidence.py | Checks raw GF(4) enumeration and certifies endpoint orbits |
| code\sample_frobenius16_locator_extensions.cpp | GF(16) invariant-base importance sampler |
| code\analyze_frobenius16_locator_sample.py | Checks sample statistics and exactly certifies sampled endpoints |

### Other proof/evidence tools

| File | Role |
|---|---|
| code\sample_locator_extensions.cpp | Main uniform 18-base \(R_6\) sampler |
| code\analyze_locator_size_ladder.py | Exact/sampled field-size ladder |
| code\analyze_monomial_146.py | Exact differential and Walsh spectra |
| code\analyze_low_degree_146_slices.py | Exhaustive degree-one through degree-three perturbations |
| code\analyze_rs_coset_moment_lp.py | Exact universal moment LP and dual |
| code\analyze_wambach_search.py | Calibrated low-shell orbit-search evidence |
| code\enumerate_wambach_rows.cpp | Performance-sensitive fixed-width row scanner |

### Authoritative generated receipts

| Receipt | Key content |
|---|---|
| generated\random_inner_threshold_outward.json | Directed \(A_{38}\) and \(R_6\) thresholds |
| generated\exact_bch_lp_caps.json | Exact shell caps; transfer margins remain binary64 |
| generated\locator_incidence_bounds.json | Locator cap, Jacobian theorem, derivative-zero bound |
| generated\locator_extension_monte_carlo.json | Five-million main \(R_6\) sample |
| generated\locator_size_ladder.json | Size ladder |
| generated\monomial_146_spectrum.json | Exact Walsh and differential spectra |
| generated\low_degree_146_slices.json | Exact low-degree exclusions |
| generated\rs_coset_moment_lp.json | Exact primal/dual universal moment obstruction |
| generated\symmetric_incidence_identity.json | Exact symmetric equations |
| generated\binary_coefficient_146_family.json | Exact GF(2) family |
| generated\frobenius_invariant_incidence.json | Exact \(x\mapsto x^2\) set stratum |
| generated\frobenius4_invariant_incidence_exact_raw.json | Raw 267-million-set GF(4) enumeration |
| generated\frobenius4_invariant_incidence.json | Compact checked GF(4) receipt |
| generated\frobenius16_locator_extension_monte_carlo_raw.json | Raw ten-million GF(16) result |
| generated\frobenius16_locator_extension_monte_carlo.json | Checked GF(16) statistics and exact orbit records |
| generated\locator_incidence_jacobian.json | Current 138-orbit locator inventory |
| generated\wambach_search_evidence.json | Radius-search results; its weight-38 floor predates new Frobenius orbits |

Selected SHA-256 hashes:

- frobenius4_invariant_incidence.json:
  C5CD08C7579A0C70B2A2C940EC65DD214C8A4D691D74A2E71CB5D0317D918E92
- frobenius16_locator_extension_monte_carlo.json:
  E9C533B0DAB4D7C64EA07819B7E4842B1550CD3EB6EE64DE699B182934161F46
- random_inner_threshold_outward.json:
  4C3EFCE399480A50C9267A2C1807FB2C98E79C6D49734274E0E2E5D3E1F74B05

## 6. Experiments and evidence

### Main uniform incidence sample

Setup:

- field: GF(256), modulus 0x14d;
- target monomial: \(u^{146}\);
- base: uniform 18-subset of the 255 nonzero points;
- samples: 5,000,000;
- statistic: \(\binom{E_B}{6}\).

Result:

\[
R_6=1.0317155313,\qquad
\text{normal 95-percent upper}=1.1360727194.
\]

Receipt: generated\locator_extension_monte_carlo.json.

Limitation: the normal endpoint uses the empirical variance. It is not a
rigorous upper confidence bound against unseen heavy algebraic components.

### Exact GF(4) exhaustion

Build in a Visual Studio x64 developer shell:

    cl /nologo /O2 /std:c++20 /EHsc /arch:AVX2 ^
      code\sample_frobenius4_incidence.cpp ^
      /Fe:generated\sample_frobenius4_incidence.exe

Run:

    generated\sample_frobenius4_incidence.exe all 0 ^
      generated\frobenius4_invariant_incidence_exact_raw.json
    python code\analyze_frobenius4_invariant_incidence.py

Observed runtime: 159.878 seconds on the current machine.

The enumerator checked 267,516,561 sets, found 67,900 admissible sets, and
recovered 40,712 distinct polynomials. The compact analyzer independently
reconstructs and evaluates endpoint remainders and verifies BCH membership.

### GF(16) invariant-base sample

Build:

    cl /nologo /O2 /std:c++20 /EHsc ^
      code\sample_frobenius16_locator_extensions.cpp ^
      /Fe:generated\sample_frobenius16_locator_extensions.exe

Run:

    generated\sample_frobenius16_locator_extensions.exe 10000000 0xf160146 ^
      generated\frobenius16_locator_extension_monte_carlo_raw.json
    python code\analyze_frobenius16_locator_sample.py

Observed runtime: 145.492 seconds while the machine was sharing CPU with
unrelated jobs.

Result:

\[
R_{6,\mathrm{GF}(16)\text{-inv}}=1.0062402040,\qquad
\mathrm{SE}=0.0041179791,
\]

with normal 95-percent upper endpoint 1.0143114431.

The nine new affine orbits are exact. Only the incidence-ratio estimate is
statistical.

### Exact low-degree and Walsh work

The Walsh calculation exhausts all 65,280 nonzero-frequency pairs. The
degree-three calculation exhausts \(256^3=16,777,216\) perturbations. These
results rule out low-degree concentration but do not control general
degree-18 perturbations.

### Universal moment LP

The exact moment identities for orders zero through 18 permit an extremizer
with

\[
R_6=33{,}072{,}431.23375.
\]

The extremizer uses agreement counts \(0,1,\ldots,17,37\). Exact primal and
dual objectives agree. The calculation proves that generic MDS/Reed--Solomon
moments cannot close the target.

### Calibrated orbit search

Radius-four affine closure exactly recovers the published minimum shells at
lengths 8, 32, and 128. Radius five finds 115 weight-38 affine orbits in \(P\).
The Chao1 estimate is diagnostic only because discoveries are not independent.

The radius search also produced exact lower bounds at weights 40 through 50.
These shells are not the current proof bottleneck.

### Exact shell LPs

Seven scaled rational LPs have independently checked primal feasibility, dual
feasibility, and strong duality. The caps are:

| Weight | Exact cap for \(A_w(C)\) |
|---:|---:|
| 38 | 773,397,229,452,928 |
| 40 | 4,210,950,726,984,874 |
| 42 | 38,281,539,335,776,556 |
| 44 | 1,144,487,702,470,434,098 |
| 46 | 7,169,795,968,863,348,791 |
| 48 | 56,162,197,953,802,356,662 |
| 50 | 645,473,488,233,155,997,799 |

Only weight 38 needs a large additional improvement for the current transfer.

## 7. Failed approaches and dead ends

### Generic distance-only bounds

Johnson/Delsarte and residual-design bounds are exact but far too loose:

\[
A_{38}(C)\leq9.7139\cdot10^{17}
\]

from the direct Johnson LP, and

\[
A_{38}(C)\leq1.2863\cdot10^{17}
\]

from the residual two-design program. Do not expect another
parameter-only constant-weight LP to close a 7.66-bit BCH-specific gap.

### Generic BCH-sandwich LP alone

The full exact sandwich LP improves the cap to
\(7.73397\cdot10^{14}\), but it remains 7.66 bits too large. It needs a new
BCH-specific locator inequality.

### Universal Reed--Solomon moments

The first 18 exact factorial moments admit an endpoint-heavy extremizer with
\(R_6\approx3.31\cdot10^7\). This route is a genuine dead end unless a
monomial-specific constraint is added.

### Ordinary Walsh factorization

The ordinary Walsh spectrum is favorable, and the affine ridge cancels in
aggregate. However, the six-dimensional kernel contains
\(P_X'(x_i)^{-1}\), which couples all points. A product of one-variable Walsh
bounds is invalid. Revisit Walsh methods only after a reduction removes or
averages this coupling.

### Generic smooth-variety estimates

The nondegenerate locus is smooth of codimension six, but generic
degree-dependent Weil/Lang--Weil constants are too large at \(q=256\).
Smoothness is structural input, not a sufficient bound.

### ApproxMCPB

The native PB/XOR encoding exactly recovers fixed-pair counts at lengths 8
and 32. It stalls at the known length-128 instance. It is not credible for the
length-256 target without a major algorithmic change.

### Direct exact GF(16) enumeration

There are 348,770,148,711,057,905 \(x\mapsto x^{16}\)-invariant 24-sets.
Brute-force exhaustion is infeasible. The invariant 18-base estimator is the
current substitute.

### Unscaled or generic exact solvers

- Full exact Z3 feasibility remained unknown after five minutes and was
  stopped.
- Unscaled QSopt_ex timed out after 120 seconds.
- Full-model floating and iterative-refinement SoPlex statuses were unstable.

The scaled rational formulation is the reliable LP path.

### BA-3 as the proof fallback

The user rejected BA-3 for the present proof because BA-3 is randomized and
the required variance or concentration calculation is unavailable. Do not
silently restore BA-3 as a proved constituent. The current extracted SPIN
draft still says “sampled spectrum-certified BA-3”; that prose is stale or
belongs to another branch of the project.

### Deferred, not dead

- Expanded Accumulant has a known spectrum and computable variance. It is a
  valid proof route, but the user reports a roughly four-millisecond,
  1.4-times performance cost.
- Additional GF(16) sampling can discover exact new endpoint orbits, but
  sampling cannot certify completeness.
- A larger Wambach radius may find more words, but lower-bound discovery does
  not close the required upper bound.

## 8. Things that are easy to misunderstand

1. RandomStepConv is temporary. The final application is expected to use
   RM2Sub. Do not quote the RandomStepConv threshold as an RM2Sub theorem.

2. The current exact application threshold is
   generated\random_inner_threshold_outward.json. The older
   generated\locator_incidence_bounds.json uses a binary64-derived threshold
   and differs in the eighth decimal place for the smooth budget.

3. generated\wambach_search_evidence.json records the old
   \(A_{38}(C)\geq912{,}640\) radius-search floor. The current floor after
   Frobenius discoveries is \(1{,}095{,}168\).

4. The exact lower floor and the needed upper cap differ by more than six
   orders of magnitude. Finding more words can calibrate structure but cannot
   by itself prove the application.

5. The GF(4) number 1.03963 is a conditional invariant-set hit ratio. The
   number \(1.1338\cdot10^{-8}\) is the GF(4)-coefficient family's
   contribution to the global \(R_6\). They are different quantities.

6. The GF(16) ratio 1.00624 concerns the invariant 24-set stratum only. It is
   not a direct estimate of the entire GF(16)-coefficient contribution.

7. Endpoint polynomials found inside a Monte Carlo run are exact witnesses
   after finite-field reconstruction and BCH membership checks. Their
   existence and affine closures are not probabilistic claims.

8. Full Jacobian rank at sampled known locators does not count all smooth
   points. The algebraic Jacobian theorem states local independence when
   \(g'\neq0\); it is not a global point bound.

9. Endpoint splitting is not random-like. The exact \(q=128\) endpoint is
   enriched by a factor of 96.94. The working heuristic concerns the sixth
   incidence moment, not the complete agreement distribution.

10. The notation \(H_w\) for the common BCH coset spectrum is unrelated to
    the generating function \(H(t)\) and the complete homogeneous functions
    \(h_n(T)\). Context distinguishes them.

11. \(P\), \(Q\), and \(C\) are extended length-256 codes. The locator
    \(A_{37}\) count belongs to the punctured length-255 BCH code.

12. A chosen quotient representative can have weight 128 even when its coset
    contains weight-38 words. Do not confuse representative weight with coset
    minimum weight.

13. The complete BCH folder is untracked. A passing verifier does not mean
    the state is committed or protected by Git.

14. The original attachment bch_spectrum_codex_bundle.zip does not contain
    RATE_HALF_FAMILY_CURVES.md. A search of the current user tree also did
    not find that file.

15. The original SPIN PDF is not retained in this worktree. The extracted
    text and selected rendered pages are under tmp\pdfs\spin_review. Treat
    that directory as scratch and potentially stale.

16. Do not run two benchmarks simultaneously. This is an explicit project
    constraint.

## 9. Open questions

### P0: final application transfer

What exact \(A_{38}(C)\), \(L_{37}\), or \(R_6\) threshold does the final
RM2Sub transfer permit?

Known: the user reports RM2Sub is at least as good as random convolution near
\(m\approx20\).

Unknown: the authoritative RM2Sub coefficients, full shell functional, and
directed certificate are absent from this folder.

### P0: smooth complement bound

Can the smooth incidence locus outside the GF(16)-coefficient family be
bounded by at most six times the random-rank reference?

Known: the direct sample, size ladder, and GF(16) invariant sample are near
one. The derivative-zero family is cheap.

Unknown: no rigorous signed character-sum or point-count bound controls the
remaining component.

### P0: endpoint control

Can the 37-agreement endpoint be bounded separately?

Known: endpoint concentration dominates the exact GF(2) and GF(4) subfield
contributions. At least 5,244 normalized locators exist.

Unknown: no useful upper bound on all endpoint locators is known. Universal
moments allow a catastrophically endpoint-heavy distribution.

### P1: Frobenius decomposition

Can coefficient Frobenius orbits produce an exact decomposition into
GF(16)-fixed and nonfixed contributions that interacts usefully with the
six-dimensional kernel?

Known: nested subfield strata are harmless on their measured incidence scale.

Unknown: the right projector or signed cancellation identity for the
complement.

### P1: LP feedback

Do the new locator and orbit facts imply an exact linear constraint that
tightens the scaled BCH-sandwich LP?

Known: new affine orbits only improve lower bounds unless a structural
divisibility or design relation is extracted.

Unknown: whether the incidence equations yield a usable moment, congruence,
or shell relation.

### P1: occupations beyond one

Does the fixed-BCH spectrum model survive the final RM2Sub high-occupation
diagnostics?

Known: current work evaluates the occupation-one functional.

Unknown: the all-occupation result needed for a full fixed-code claim.

## 10. Recommended next steps

### 1. Re-establish and preserve the current state

Run:

    cd bch_spectrum_work\bch_spectrum_codex_bundle
    python code\check_all.py

Expected final lines include:

    PASS: exact binary/GF(4) and sampled GF(16) Frobenius strata
    PASS: twenty-three new certified affine weight-38 orbits

Then decide whether to add the untracked BCH folder and this handoff to a
branch. Do not assume a commit already exists.

### 2. Import the final RM2Sub transfer

Locate the other agent's RM2Sub artifact and record its exact parameter set.
Recompute the shell transfer at fixed row counts. Then adapt
code\certify_random_inner_threshold_outward.py to produce an RM2Sub-directed
certificate.

This step determines whether \(R_6\leq6\) remains the right theorem target or
whether RM2Sub gives additional slack.

### 3. Write the exact Frobenius partition before more sampling

Define coefficient Frobenius under \(c\mapsto c^{16}\). Partition the
polynomial family and the character sum into the GF(16)-fixed locus and its
complement. Keep the 37-agreement endpoint as a separate term.

The desired output is an exact identity, not an estimate. It should explain
which signed sums remain after the fixed locus is removed.

### 4. Attack the endpoint term

Use the locator equation

\[
z^{37}+g(z^2)
\]

and the derivative identity for its 37 nonzero roots. Search for an
orbit-count, divisibility, or low-state certificate that upper-bounds the
endpoint family. Avoid relying only on the first 18 Reed--Solomon moments.

A useful intermediate result would be any upper bound below
\(2^{34.093}\) for \(L_{37}\), even without a full \(R_6\) theorem.

### 5. Use additional GF(16) runs only as targeted evidence

Independent seeds can measure orbit-discovery saturation and produce more
exact lower witnesses. Record sample counts and seeds. Canonicalize every
endpoint against the existing 138-orbit inventory.

Do not interpret discovery saturation as completeness. Never run this
benchmark concurrently with another benchmark.

### 6. Feed genuinely new identities into the exact LP

If Steps 3 or 4 produce a linear equality, inequality, congruence, or
divisibility condition, add it to code\export_exact_constraints.py. Regenerate
the scaled LP and verify both primal and dual certificates.

### 7. Return to the full SPIN application

After the RM2Sub occupation-one transfer is fixed, evaluate occupations beyond
one with the same fixed-code spectrum model. Keep Expanded Accumulant as the
formal fallback until these checks pass.

## 11. User preferences and project constraints

- Pursue heuristic evidence and a narrow proof in parallel. Do not put the
  entire case on one technique.
- The user is willing to use the random-inner model temporarily, but the
  intended final transfer is RM2Sub.
- The user wants evidence specifically for the BCH code actually used by
  SPIN, not only generic coding-theory bounds.
- A narrow proof of the application functional is sufficient. The whole
  spectrum is not required.
- Keep Expanded Accumulant as the proof-backed fallback. Its performance cost
  is acceptable but undesirable.
- Do not present BA-3 as proof-backed without a variance or concentration
  calculation.
- Performance is a first-class constraint. Preserve fixed-width loops, AVX2
  batching, stack storage, explicit layouts, and predictable memory access in
  hot kernels.
- Never run two benchmarks simultaneously.
- Keep exact claims, experimentally supported claims, and hypotheses visibly
  separate.
- For long mathematical prose, use the Controlled Writing for Cryptography
  skill.

## 12. Astra Bootstrap

GOAL:

- Justify or certify BCH \([256,128,38]\) as the fast rate-half SPIN outer.
- Temporary target: prove \(R_6\leq6\).
- Final threshold must be recomputed for RM2Sub.

CURRENT BEST APPROACH:

- Reduce \(A_{38}(C)\) to normalized 37-agreement locators.
- Split derivative-zero, proper-subfield, endpoint, and smooth-complement
  contributions.
- Use the symmetric equations and six-dimensional Vandermonde kernel.

ESTABLISHED FACTS:

- \(Q\subset C\subset P\), with \(A_w(C)=A_w(Q)+31H_w\).
- Directed RandomStepConv cap:
  \(A_{38}(C)\leq3{,}827{,}351{,}840{,}403\).
- \(R_6\leq6\) suffices.
- Main sample: \(R_6=1.03172\), normal upper 1.13608.
- GF(2) and GF(4) families are exactly exhausted and globally negligible.
- GF(16) invariant sample: 1.00624, normal upper 1.01432.
- Exact current floor:
  \(A_{38}(C)\geq1{,}095{,}168\).
- 138 certified affine orbits; 5,244 normalized locators.

CURRENT PARAMETERS:

- Field modulus: 0x14d.
- Target monomial exponent: 146.
- Base size: 18.
- Endpoint agreement count: 37.
- Incidence extension order: 6.
- Domain size: 255 nonzero field elements.

IMPORTANT FILES:

- bch_spectrum_work\bch_spectrum_codex_bundle\PROGRESS.md
- bch_spectrum_work\bch_spectrum_codex_bundle\LOCATOR_INCIDENCE_DUAL_TRACK.md
- bch_spectrum_work\bch_spectrum_codex_bundle\SPIN_HEURISTIC_EVIDENCE.md
- bch_spectrum_work\bch_spectrum_codex_bundle\generated\random_inner_threshold_outward.json
- bch_spectrum_work\bch_spectrum_codex_bundle\generated\frobenius4_invariant_incidence.json
- bch_spectrum_work\bch_spectrum_codex_bundle\generated\frobenius16_locator_extension_monte_carlo.json
- bch_spectrum_work\bch_spectrum_codex_bundle\code\check_all.py

DO NOT REPEAT:

- generic Delsarte/Johnson-only attempts;
- universal moment LP without monomial constraints;
- naive one-variable Walsh factorization;
- brute-force GF(16) invariant-set enumeration;
- unscaled exact LP solves;
- BA-3 proof claims without variance.

BIGGEST UNCERTAINTIES:

- final RM2Sub shell threshold;
- endpoint upper bound;
- smooth complement outside GF(16);
- high-occupation fixed-BCH behavior;
- missing RATE_HALF_FAMILY_CURVES.md source.

NEXT 3 ACTIONS:

1. Run code\check_all.py and preserve the untracked state.
2. Import the final RM2Sub transfer and regenerate a directed threshold.
3. Derive an exact GF(16)-fixed/complement character-sum decomposition with a
   separate endpoint term.
