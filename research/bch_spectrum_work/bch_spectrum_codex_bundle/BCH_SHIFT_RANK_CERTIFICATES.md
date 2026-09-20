# Independent certificates for the stronger BCH moments

Updated: 2026-09-05.

We independently certify d(Q-perp) >= 24 and d(P-perp) >= 26 for the
length-256 extensions of the dimension-123 and dimension-131 BCH codes.
The certificates do not use the inaccessible SchaubPlus table. They justify
the three additional moment equations previously treated as provisional.
An additional case refinement now proves d(Q-perp) >= 30. The resulting
strength-29 model certifies A38(C) <= 24,866,368,872,377 and gives a full
M22 margin of 39.72179723376 bits with the saved separate other-shell bounds.
Joint optimization of the dominant shells now closes the original target
at 40.72497637925699 bits. See
[RANDOM_INNER_M22_CLOSURE.md](RANDOM_INNER_M22_CLOSURE.md) for the final theorem.

## Code identification

Let E0 be either primitive binary narrow-sense BCH code of length 255,
with designed distance 39 or 37. The generator construction uses the fixed
primitive element alpha=2 of F256 and field modulus 0x14D. The checker
reconstructs the binary generator g and the reciprocal of (X^255+1)/g.
The latter polynomial generates E0-perp: its shifts are orthogonal to all
shifts of g, and their dimensions sum to 255. Direct polynomial evaluations
verify its full defining root set Z.

Every dual generator word has even weight. The root set contains exponent
zero, which also gives this parity constraint directly.

## A finite case cover

Fix a nonzero binary word c in E0-perp. Define its Fourier values by

\[
F(e)=\sum_{j=0}^{254}c_j\alpha^{ej},\qquad e\in\mathbb Z/255\mathbb Z.
\]

The Fourier transform is invertible, so at least one value is nonzero.
Moreover F(2e)=F(e)^2. Thus nonzero Fourier values occur in complete binary
cyclotomic cosets. Order the cosets outside Z by their least representatives.
Assign c to the first coset in this order on which its Fourier values are
nonzero. If b represents that coset, then F(b) is nonzero. All exponents in
Z and all earlier cosets have Fourier value zero.

The checker reconstructs this case partition. It checks 17 cases for Q0
and 16 cases for P0; their union includes every nonzero dual word.

## Independent-set witnesses

For a word c in one case, let S be its support and w=|S|. For each exponent
e define v_e=(alpha^(ej)) indexed by j in S, as a vector over F256.
The coordinate-sum functional sends v_e to F(e).

The witness starts from the independent singleton {v_b}. Each step uses
two operations that preserve linear independence. Multiplying every
exponent by 2^k applies an invertible field automorphism to all entries.
Adding a common exponent t multiplies the columns by nonzero scalars.

The checker verifies that all transformed exponents have Fourier value zero
under the current case assumptions. It then appends b. The new vector v_b
is outside their span because its coordinate sum is nonzero. Hence every
successful step increases the number of independent vectors by one.

A witness ending with r exponents proves w >= r, since these vectors lie
in F256^w. The search stores relative exponents, centered at b. A saved step
(k,s) therefore corresponds to the actual exponent transformation

\[
e\longmapsto2^k e+s+(1-2^k)b\pmod{255}.
\]

The independent checker reconstructs this transformation and every zero-set
containment. It does not rerun or trust the beam search. Search truncation
can miss a witness, but cannot make a checked witness unsound.

## Root-run witnesses for the other cases

Most later cases have long arithmetic progressions of Fourier zeros. The
checker records and verifies a start, a step coprime to 255, and a run length
r. Such a run proves w >= r+1 by the BCH Vandermonde argument. If w <= r,
the first w zero equations form an invertible matrix on the support S,
contradicting that c is nonzero.

These witnesses also avoid an inefficiency in the first search. Dense zero
sets created too many candidate independent sets, and four late Q cases
hit search time limits. Each of those cases instead has a much stronger,
short root-run certificate. The timed-out searches are retained separately.

The completed Q cover has minimum rank 23, and the P cover has minimum rank
25. Since both punctured duals are even, their minimum distances are at
least 24 and 26, respectively. Toy regression exhausts all 31 nonzero binary
words of length five and checks 188 valid witness transitions.

## Transfer to the length-256 duals

Let E be the even extension of E0, with its extra coordinate indexed by
field element zero. The checker verifies that translation x -> x+1 and
scaling x -> 2x preserve E on its full generator basis. These maps generate
all affine coordinate maps, so translations also preserve E-perp.

Suppose a nonzero word in E-perp had weight below the certified bound for
E0-perp. It has a zero coordinate. Translate that coordinate to the extra
position and puncture it. The punctured word is nonzero, has the same weight,
and is orthogonal to E0. This contradicts the punctured dual-distance bound.
The same lower bound therefore holds for E-perp.

## Consequences for the moment model

Projection of Q onto any 23 coordinates is surjective. A proper image would
give a nonzero dual word supported on those coordinates. The linear fibers
have equal size, and the same statement holds for every coset of Q. Thus
Q and its cosets have orthogonal-array strength at least 23. Likewise P has
strength at least 25.

Retain q=A(Q), h for the common nonzero Q-coset spectrum in P, and the
length-256 Krawtchouk transform K. Since A(P)=q+255h, the new equations are

\[
(Kq)_{22}=(Kh)_{22}=0,\qquad(K(q+255h))_{24}=0.
\]

At this intermediate stage, no separate degree-24 zero equation was asserted
for q or h. The stronger
moment model uses these certified identities, not a sampled spectrum or
an unverified literature bound.

## Reproduction

The independent audit is

    python -B code/certify_bch_shift_rank.py --verify

Its receipt is generated/bch256_shift_rank_dual_distances.json. Complete
witnesses are in generated/shift_rank_q24_complete and
generated/shift_rank_p26_complete. The checker reconstructs all generators,
root sets, case assumptions, witness steps, parity, and affine invariance.

The new model is generated/shift_rank_hull_probe. Its exported LP is exactly
identical to the earlier provisional LP; the saved rational solution is
reused only after checking that identity. A new independent primal/dual
audit attaches the now-certified mathematical inputs. No frozen provisional
file is altered, and no new optimizer execution is needed for this model.

## Refinement to dual distance 30

A higher-threshold search gave rank at least 29 in every Q case except the
case represented by 15, which reached rank 28. Independent auditing of that
complete cover therefore certifies distance 28, not the search target 30.

For the remaining case, split according to whether F(23) is zero. If it is
zero, its entire binary cyclotomic coset is zero, and a new rank-29 witness
uses the original nonzero value F(15). If F(23) is nonzero, another rank-29
witness uses exponent 23 as its pivot and retains the parent's zero set.
These two alternatives exhaust the parent case. Their independently checked
witnesses, together with the other 16 cases, prove rank at least 29 for
every nonzero Q0-dual word. Evenness yields distance at least 30.

The extension argument above gives d(Q-perp) >= 30. Since Q is contained in
P, P-perp is contained in Q-perp and inherits the same lower bound. No
stronger P-specific bound is needed here. The new Q moment equations are
(Kq)j=(Kh)j=0 for j=24,26,28. The resulting model has 196 variables and
1163 normalized rows. Its exact primal/dual certificate gives the cap and
39.72179723376-bit separate-cap M22 bound stated above. That calculation
alone did not close the 40-bit goal.

The split checker reconstructs both children and checks that their zero sets
and nonzero pivots agree with the exhaustive split. Its proof does not depend
on how the search selected exponent 23. Replay sequentially:

    python -B code/certify_bch_shift_rank_threshold.py 39 shift_rank_q30_complete bch256_shift_rank_q28 --verify
    python -B code/certify_bch_shift_rank_split.py --verify
    python -B code/audit_bch_extended_hull_cap.py prepare_bch_shift_rank_oa29_probe shift_rank_hull_probe --verify

The refined receipt is generated/bch256_shift_rank_q30_refined.json; the
model and exact solution are in generated/shift_rank_oa29_hull_probe.
The older Q search folder named `q30_complete` certifies only 28 without
this separate split proof. Higher-threshold P searches also stopped short
of their targets; no distance-32 P bound is claimed.

## Joint objective closes M22

The individual worst-case weight counts need not occur in the same feasible
spectrum. The same 1163-row model now bounds the paired contributions of
weights 38, 40, and 42 jointly. Upward-rounded objective coefficients and
an independently checked exact rational dual witness give a full bound with
40.72497637925699 bits. Every other shell and occupation remains included.

The witness is in generated/shift_rank_oa29_joint. It is a weighted-objective
certificate, not an h38-cap certificate; do not run the h38 checker on it.
The final audit is generated/bch256_m22_unconditional_closure_audit.json.
Run `python -B code/verify_bch_m22_closure.py` for the final reconstruction,
independent one-row transfer checks, evidence hashes, and exact full sum.

The ideal-M22 goal is complete. Next: independent review and paper integration
of the finite theorem, with practical-inner transfer kept separate.
