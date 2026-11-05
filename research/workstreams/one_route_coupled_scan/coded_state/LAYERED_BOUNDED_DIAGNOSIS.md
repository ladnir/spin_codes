# Bounded diagnosis: finish the layered proof or strengthen the inner?

## Decision

Recommend an incremental design change, not further neighboring-suffix tuning.
Preserve the present two-shear candidate and its partial certificates. Use
almost-uniform invertible pair maps as the next research model, retaining the
RS outer, regional permutations, 512-symbol state, and reset recurrence.

This recommendation concerns proof tractability and explicit operation counts.
It is not evidence that the present construction has bad distance. No native
encoder was changed or timed. The stronger pair model is not fully certified.

The remaining question after this change is empty-pair geometry. Distinguish
whole-state inactive intervals from holes within a live state. The former
already have a successful counting method in the uniform-feedback control.
The latter may need more local mixing. Test those families before committing
to a second layer or a periodic expensive refresh.

## Scope and preservation

The diagnostic uses K=2^20, N=2^21, rate 1/2, and GF(2^32) histories.
The layered candidate is tested at W=838860, the exploratory 40% threshold.
The certified whole-code control reaches 46%; the long-term goal remains
near-half distance. No new full certificate or scalar-extension claim follows
from a passing selected profile.

The existing 79 Python sources in coded_state were unchanged. The checkpoint
hash is SHA256 of their sorted filename, a NUL byte, and exact file bytes,
excluding new filenames beginning with diagnose_:

    7861e9defd8a29de20a557fae822c482a4a1ffeadda66aa37bb1113fc49dada6

Both frozen theory-helper hashes also match the previous checkpoint.
The suffix-193/194/195 certificates are retained, not superseded.
Only diagnostic source, this note, and index/ledger entries were added.

## What the coarse map actually measures

Let R be the number of regions in a suffix, M=8192, b=512, and
P=(RM-b)/2. The post-initial pair counts are (z,l,r,f) for types
(00,10,01,11). The exceptional-full-pair count is e. B counts internal
births of nonzero state coordinates. Thus l=P-z-r-f and V=P-z+f counts
nonzero incoming-state coordinates. These are aggregate histories, not
sampled encoder runs or observed codewords.

The grid uses:

- R in {129,144,160,176,192,200,208,224,240,256};
- z/P in {0,1/16,1/4,1/2}, rounded down;
- f/(P-z) in {0,1/4,1/2,3/4,1};
- r/(P-z-f) in {0,1/2,1};
- e/f in {0,1/4,1/2,3/4,1};
- B in {0,floor(P/32),floor(P/8),min(floor(P/2),floor((W+b)/2))}.

Duplicate tuples are removed. Necessary wiring/weight tests exclude 2132
of 8480 selected tuples. These tests do not prove the remaining tuples
are realizable. For the retained tuples, the diagnostic evaluates bounds
over every active-row count, not just a sample of those counts. The grid
does not cover the gaps between tuples or omitted suffix lengths.

Bounds combine the retained location-aware transfer, coarse forward and
inverse counts, and the older stronger bound when z=0. Three bounded
weight proposals are replayed with rational local weights and directed
logarithms. New weights live only in the diagnostic process; the old bank
is not edited. Active-row counts below 5901 may use the location bound too,
with d_S <= min(d_bar,(R-128)s), instead of always using the old fixed-weight
small-source fallback.

| R | Retained tuples | Current bound passes | Current bound fails | Of those failures, also fail the ideal-pair comparison |
| ---: | ---: | ---: | ---: | ---: |
|129|636|620|16|7|
|144|636|542|94|8|
|160|636|459|177|12|
|176|636|419|217|15|
|192|636|407|229|30|
|200|636|403|233|60|
|208|634|403|231|60|
|224|634|405|229|75|
|240|632|414|218|95|
|256|632|434|198|95|

Passing means log2 upper bound below -100 on the applicable good-route
event. This is not an estimated fraction of the full domain covered.
There are 1842 unresolved selected tuples. The ideal-pair comparison passes
1385 of them and leaves 457 unresolved. Since the ideal model has no
exceptional label, several old e-values project onto the same ideal class;
these numbers are not counts of distinct ideal classes.

All remaining failures in this selected map have z>0. Most have z/P equal
to 1/4 or 1/2. The small-empty-pair slice has far fewer failures. Thus the
remaining work is not confined to R=196 or a narrow continuation of the
previously studied neighborhood.

## Weight artifacts versus persistent failures

The initial frozen-bank map was misleadingly pessimistic away from its
tuning region. Two bounded probes demonstrate the problem:

| R | (z,r,f,e,B) | Tested s | Retained location slice | Reoptimized slice |
| ---: | --- | ---: | ---: | ---: |
|192|(393088,0,196544,98272,24568)|5901|+8895271|-5568812|
|160|(327552,163776,0,0,81888)|5901|+6576130|-10091250|

The second retained value is before importing the first probe's weight.
The first optimized interpolation exponent is about 1.1184; the second is
zero. This is very different from the approximately 24--26 exponents tuned
near the old crossover. Each result concerns one source-count slice, not
the entire profile. The final map incorporates all three proposals.

One important failure persists under a dedicated probe. For

    R=200, (z,r,f,e,B)=(49136,24568,589632,368520,49136), s=8192,

the location slice improves from +92857 to +53315 at exponent about 25.7285.
This remains a failed bound, not a lower bound on the code's failure probability
or a proof that no better weights exist.

Suffix length is not a monotone progress coordinate. With the old central
counts fixed, the retained bound is negative at R=192, positive at R=200,
and negative again at R=224 and R=256. Following only adjacent R values
therefore obscured the shape of the unresolved domain.

## Comparison A: uniform feedback on the entire live state

Keep the same outer, route, b=512, and

    Q_(t+1) = delta Q_t + Y_t.

Replace feedback by an independent uniform b-by-b matrix over GF(2^32)
at each batch. For each fixed nonzero incoming state, its image is uniform
over all b-symbol vectors. The geometric reset weights do not change the
reset-support ranks used in the running-sum proof.

The existing complete 46% proof replays at both 80 and 120 decimal digits:
459 rectangles cover 4196352 feasible episode-count pairs. The worst directed
class bound is below -65.77481065587761, giving conservative integer margin
42 before the small sandwich-density correction. The partition hash remains

    faf9d54b427678b812fae5cfbdbcb20060b65385c66b4c8af6a5718bd00bef6e.

Consequently the outer, routing, and amount of state are compatible with
more than the current 40% exploratory target. This does not prove that a
cheaper local inner inherits the certificate, or that the same components
already close 49% over this coefficient field.

## Comparison B: nearly uniform randomness within each pair

This is a different setup distribution, not a new assumption about the old
two-shear keys. Fix q=2^32 and delta in GF256 other than 0 or 1. Independently
sample A uniformly from SL(2,q), rejecting unless both A and A+delta I have
all entries and determinants nonzero. Keep the same cycling pair schedule.
The recurrence remains

    Y = X + (A+delta I)Q,       Q_next = delta Q + Y.

SL(2,q) acts transitively on the q^2-1 nonzero pairs. Each specified nonzero
image has exactly q preimages in the group. Rejection removes at most a
fraction

    4/(q+1) + 2q/(q^2-1) + 1/(q-1) = (7q-3)/(q^2-1).

The first two terms exclude zero entries of A and zero diagonal entries of
A+delta I. The last term bounds the trace condition det(A+delta I)=0:
for fixed trace, the determinant equation permits at most q(q+1) matrices.
Hence each image of a nonzero pair has probability at most

    1/(q^2-7q+2).

Relative to q^-2, the product of these density losses over at most N/2
pairs is below 4096/4089, by Bernoulli's inequality. The old one-sided and
exceptional q^-1 point masses are absent. A zero pair still has deterministic
feedback. This is the important remaining difference from full-state uniformity.

Every accepted A=[[a,b_key],[c,d]] has b_key!=0 and factors as

    A = L(v) U(b_key) L(w),
    v=(d+1)/b_key, w=(a+1)/b_key,
    L(t)=[[1,0],[t,1]], U(t)=[[1,t],[0,1]].

Thus the proposal uses three variable products per pair rather than two.
Divisions and rejection are setup work. The variable-product count rises
by 50%; that is not a 50% end-to-end timing prediction. No new timing exists.

Exact GF4/GF8/GF16 enumeration verifies 333 nonzero-input image laws,
retained point masses, and 2763 three-shear factorizations.

## Bounds used for the pair comparison

Drop the old exceptional label and project a profile to (z,r,f,B). Every
nonzero incoming pair now supplies two random inverse-input coordinates,
so D_pair=2(P-z). They are complete edges of the regional matching.
Using the existing actual-moment envelope E_R(u), a simultaneous shortening
cap is

    d_pair <= min(D_pair/2, (R-128)M,
                  max(0,max_u[min(D_pair,2E_R(u))-128u])).

Let A_pat bound the log2 count of support histories with these counts.
Use either the multinomial pair-support count or the initial-support plus
temporal-change count. There are at most 2^V choices of internal zero-output
edges. Their independent path equations give

    C >= max(0,V-W,V+B-b-W).

Three conservative log2 bounds, including one bit for density inflation, are

    RM + 32(W-D_pair+d_pair) + 1,
    A_pat + V + 32(V+b-C+d_pair-D_pair) + 1,
    A_pat + 32((R-128)M-max(0,D_pair-V-b)) + 1.

The last bound charges next-state zeros among the free coordinates. There
are at most V+b nonzero next-state positions over those batches. The first
two bounds count outputs and state amplitudes, respectively. No independence
between outer membership and temporal constraints is assumed.

For the stubborn R=200 central profile, the minimum is -5277460, versus
the reoptimized positive bound +53315 for one slice of the original design.
This is a comparison of bounds with different strengths, not a measured
distance difference. It nevertheless identifies a concrete way to remove
the exceptional-branch accounting burden.

The ideal-pair bound remains positive for, among others,

    R=192, (z,r,f,B)=(393088,0,294816,98272): +478536;
    R=160, (z,r,f,B)=(327552,0,327552,81888): +1326593.

These are empty-pair-heavy classes. Their aggregate counts do not say whether
the empty pairs are simultaneous inactive batches or persistent holes within
live batches. Full-state randomness handles the latter; local pair randomness
does not. The former need episode placement and reset-rank accounting, not
stronger randomness applied to a zero state.

## Next discriminating test

Do not enlarge the old local certificate boxes again as the main task.
Use the rejected-uniform SL2 pair distribution as a simple intermediate
model. Derive an episode-aware bound for its remaining empty-heavy classes.
First separate zero whole states from partly active states; preserve the
successful uniform-control argument on the former. If the partly active
families still fail, compare a second crossing layer or periodic stronger
mixing before doing another large coverage search.

The diagnostic supports this direction, but does not promise that one extra
shear closes 40%, much less the near-half objective. The current candidate
remains a possible faster fallback if a better proof later appears.

## Reproduction and resource limits

Run the new scripts serially from coded_state:

    python -B diagnose_layered_families.py --test
    python -B diagnose_layered_families.py --retuned-map --seconds 120
    python -B diagnose_layered_families.py --ideal-examples
    python -B diagnose_uniform_pairs.py

Each weight proposal has a 20-second cap. The map has a 120-second cap;
the complete selected screen took about five seconds after preparing weights.
One initial run timed out because it recomputed large binomials for each
marker. That partial run was rejected; caching those constants fixed the
diagnostic cost without changing its bounds. Process memory was capped at
1.5 GiB, with one BLAS thread. No native benchmark ran.

C: had only about 10.5 GiB free, below the research safety floor. Numerical
runs therefore used Peach's existing research connection and a temporary
directory with over 60 GiB free. No local numerical jobs, raw dumps, or
archives were created. No unrelated files were deleted. Existing proof
replays and the frozen-helper hashes were checked after the diagnostic.
