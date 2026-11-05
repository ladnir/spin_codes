# Defects enlarge an input envelope and add zero constraints

The unchanged three-shear inner has a closed-episode identity that extends
the clean input-containment argument to births and one-sided pairs. It
also closes a first selected profile with internal defects. This is not
whole-code distance closure, and it does not complete the clean shared-region
families left open by the preceding pass.

## State history and defect charge

Fix a message and route before exposing the inner maps. Consider one closed
episode with incoming supports S_1,...,S_l and S_(l+1)=empty. Its initial
support U=S_1 is the activation input support; write w=|U|. The actual
construction and its probability space are specified in UNIFORM_PAIR_PINS.md.

At update t, let L_t be the union of matching pairs intersecting S_t.
Count the following coordinate events, with multiplicity over time:

- H: births in S_(t+1) from pairs empty in S_t;
- P: one-sided incoming pairs;
- B: all births from inactive to active coordinates;
- D: all deaths from active to inactive coordinates.

Let C=sum_t |L_t minus S_(t+1)| count prescribed zero next-state coordinates
in live pairs. Put e=H+P. A birth in a live pair can only occur at the
inactive partner of a one-sided pair. If B_live counts those births, then

    B=H+B_live,    B_live<=P,    D=w+B.

In a one-sided pair, its inactive partner either becomes active or remains
zero. Therefore

    C=D+P-B_live = w+H+P = w+e.

This identity includes all temporal correlations in the support sequence.
It does not assume uniformly distributed states. In particular B<=e.

## All inputs fit in an envelope of at most C coordinates

Define V as the union of U, all coordinates born from empty pairs, and
the inactive partner coordinates of every one-sided pair. Every active
state coordinate lies in V: its first activation is initial, an empty-pair
birth, or a birth at an inactive partner. Every live pair also lies in V.

The update gives X_t=Q_(t+1)-A_t Q_t. Its support is contained in
S_(t+1) union L_t, hence in V. Including the activation input, every input
in the episode is supported on V, and

    |V| <= w+H+P = C.

Repeated events may name the same coordinate, so this is an upper bound,
not an assertion that all named coordinates are distinct. The clean case
has e=0 and V=U. With defects, initial-support containment itself can fail;
the enlarged envelope is the replacement.

## Probability and history multiplicity

Use the independent uniform pair-image comparison from the retained proof.
For a fixed support history and fixed input, each prescribed zero in a live
pair costs q^-1, with q=2^32. Empty pairs have deterministic next state X;
their compatibility may be discarded for an upper bound. Nonzero conditions
can only reduce probability. The fixed-history contribution is at most
q^(-C)=q^(-w-e). The accepted maps incur one global density factor below two,
not a factor per support history or per episode.

For fixed U, a support history is specified by its birth events and the
death time of each initial or newly born active component. If there are
B births, their coordinate/time labels have at most binom(bl,B) choices,
where b=512. Every component has at most l choices of death time. Since
B<=e, the number of histories of exact defect charge e is at most

    T(l,w,e) = l^w * sum_(B=0)^e binom(bl,B) l^B.

Invalid event descriptions only enlarge this count. They need not satisfy
nonempty intermediate states, cyclic matching compatibility, or properness.
Consequently the bound applies to the subset of valid proper episodes.
The formula counts support histories, not field-valued amplitudes: those
were already summed in the forward pair-image probability.

For a fixed activation support, one extra charge at l=8 contributes at
least 16 probability bits after this history bound, before route and
activation choices. At l=17 it contributes at least 14 bits. These figures
alone are not code certificates because the larger envelope also makes
routing easier.

## Joint route bound for a defective crossing episode

The initial matching has at most e one-sided pairs. The number of possible
initial supports of weight w is therefore at most

    A_b(w,e) = sum_a binom(b/2,a) 2^a binom(b/2-a,(w-a)/2),

where 0<=a<=min(e,w,b/2), w-a is even, and invalid terms are zero.
Fix a split u+v=l+1 of the input batches across two regions. Let n_u,n_v
be their nonzero input counts. Set v_w=min(b,w+e). For each fixed initial
support and support history, the remaining left-region inputs have at
most (u-1)v_w available positions; the right region has at most v*v_w.

The route/history contribution is at most

    sum_(w>=1) A_b(w,e) T(l,w,e)
      * binom((u-1)v_w,n_u-w) binom(v*v_w,n_v)
      / [binom(M,n_u) binom(M,n_v) q^(w+e)].

Here w<=min(n_u,b), M=8192, and invalid binomial coefficients are zero.
The envelope can depend on the support history. The bound fixes that
history first and uses its maximum envelope size; it does not assume
independence of the envelope and routing. Summing over histories and initial
supports permits overlap between the counted events.

For e=0 the verifier uses the sharper previous clean-episode formula,
including w>=2^l, rather than the relaxed general history count. For an
actual defective episode, that clean lower bound on w is not imported.

## First restricted certificate with internal defects

Keep five closed crossing episodes of length three, followed by an
unrestricted live tail, with disjoint touched regions and R=131. Four
episodes remain clean. One designated episode has exact defect charge
e=1 or e=2. The actual RS zero-count sum and endpoint union are retained.
The verifier includes all five choices of the designated episode and all
source counts s=8,...,28 allowed by the route cap.

| Extra charge in one episode | Good-route upper | Route exception | Combined upper |
| --- | ---: | ---: | ---: |
| 0, clean control | 2^-813 | 2^-100 | 2^-99 |
| 1 | 2^-747 | 2^-100 | 2^-99 |
| 2 | 2^-750 | 2^-100 | 2^-99 |

These are separate restricted families. Their union needs its own small
union factor; do not treat the table as a whole-code result. The clean
control improves on its earlier bound because this calculation retains
the actual RS zero count. The e=2 bound being slightly stronger than e=1
is not a monotonicity claim: these are exact-charge classes and conservative
upper bounds. No output-weight restriction was used.

Exhaustive tests check 115893 support histories at widths2,4,8, including
the exact charge identity, envelope containment, and birth/death balance.
They also check 245 history-count bounds and 66 initial-mask enumerators.
The combinatorial arguments above establish the general claims; the toy
checks validate the implementation and can catch indexing mistakes.

The three profile runs take under two seconds each locally. Tests take
about one second. The source has 20-second test and 35-second profile caps,
bounded caches, and no numerical file output. Runs were serial with more
than23GiB free at the preceding checks. No encoder or earlier proof helper
was changed.

## Next discriminating step

The construction already supplies a useful defect tradeoff: the input
envelope cannot grow without increasing prescribed live-pair zeros. This
removes the need to assume that defective inputs stay in the initial
support. It does not settle whether the complete history count is small
enough in the remaining finite regimes.

Next combine this envelope with the shared-region and actual-zero count,
and test the first residual family containing an internal defect. Track
e=H+P explicitly: it counts coordinate births and pair occurrences, not
the number of defective batches. Broad history classes may have many such
events in one batch. That distinction matters when mapping this interface
back into the full certificate partition.

Keep the construction fixed for that test. No structural failure demanding
redesign has been established, but the full defect-rich and extension-field
distance proofs remain open. The certified control is preserved.
