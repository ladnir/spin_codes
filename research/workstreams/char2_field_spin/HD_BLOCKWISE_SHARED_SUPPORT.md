# Joint support constraints for the shared-transform inner

The shared FFT8 map has a useful local MDS description. It also has an exact
cross-tile cancellation symmetry when the output shift is8. Moving the shift
to64 removes that symmetry in the tested class, but does not make state
constraints independent.

This note gives exact rank interfaces and bounded exhaustive comparisons.
It does not give a whole-code distance certificate or a timing claim.

## Variables and local code

Use the coefficient basis and width64/state4 construction from
[HD_BLOCKWISE_SHARED_MAPS.md](HD_BLOCKWISE_SHARED_MAPS.md). Let beta be8 or64.
For a tile at alpha, write its input polynomial as p_alpha, of degree below8.
The four-coordinate state represents a shared polynomial Q of degree below4.
The tile input and output values are

    X_alpha = evaluations of p_alpha on alpha+[0,8),
    Y_alpha = evaluations of p_alpha+Q on (alpha XOR beta)+[0,8).

Here the two evaluation cosets are disjoint. Let L(X) be the monic subspace
polynomial with roots [0,8), and set gamma=L(beta), which is nonzero. Define

    P_alpha(X) = p_alpha(X) + (L(X)+L(alpha))*Q(X)/gamma.

This polynomial has degree below12. On the input coset it evaluates to
p_alpha. On the output coset it evaluates to p_alpha+Q. Conversely, division
by L(X)+L(alpha) recovers p_alpha and Q from every polynomial of degree below12.
Thus the joint 16 observed coordinates form an MDS[16,12,5] code.

This statement includes the four state coordinates as independent local
variables. In a global history they may depend on the input. The rank
interface below retains those dependencies explicitly.

## Output constraints and the shared state

Fix a tile input support S of size w and a set Z of z output positions required
to be zero. Coordinates outside S are fixed to zero. The unknown vector is
the w supported input values followed by the four shared state coefficients.
The output-zero matrix is

    B_alpha(S,Z) = [ (F_8)_(Z,S) | (A_alpha)_Z ].

The MDS[16,12] description proves

    rank B_alpha(S,Z) = min(z,w+4).

Indeed, it imposes 8-w input zeros and z output zeros on the 12-dimensional
polynomial space. Any at most12 evaluation constraints are independent.

The input block is a submatrix of a superregular matrix, so it has rank
min(w,z). Eliminate those input pivots. The remaining equations involve only
the four state coefficients; call their matrix R_alpha(S,Z). Then

    rank R_alpha(S,Z) = min(4,max(0,z-w)).

For a packet's eight tiles, the exact output rank is

    sum_alpha min(w_alpha,z_alpha)
      + rank(stack_alpha R_alpha(S_alpha,Z_alpha)).

The second term is one joint four-column rank. Replacing it by the sum of
the eight local ranks is invalid.

## A paired cancellation family at shift8

Consider adjacent tiles alpha and alpha XOR8. Their input and output cosets
are exchanged. Choose an input support S across their16 positions, and force
output zeros on the complement of S after swapping the two tiles.

Both polynomials P_alpha must then vanish on the same16-|S| evaluation nodes.
Their high four coefficients encode the same Q. Their difference has degree
below8. If |S|<=8, the common zero set has at least8 nodes, so the polynomials
are equal.

Consequently this constraint space is precisely a shortened MDS[16,12] code.
For 5<=w=|S|<=8 it has dimension w-4. Its state image also has dimension w-4:
a vector in the kernel of the state map would have degree below8 and support
at most8, which an MDS[16,8,9] code excludes.

On this space, the output tiles swap the input tiles. The input polynomials
differ by the degree-below4 polynomial Q, so their high four coefficients
cancel. The pair's feedback contribution is zero.

For example, take packet input positions

    S = {0,1,2,8,9,10}

and output-zero positions

    Z = {3,4,5,6,7,11,12,13,14,15}.

Each tile separately leaves two state equations. Together they still have
rank2. The ten-variable output matrix has rank8, leaving a two-dimensional
space with live state, full input support, and zero feedback.

This is not a whole-code counterexample. The other six tiles, later packets,
outer parity constraints, and the initial zero state still matter. In
particular, with an invertible mixer, a live state and zero feedback cannot
return to zero at this packet.

## What shift64 changes

At beta64, the packet's output evaluation nodes64 through127 are disjoint from
its input nodes0 through63. The preceding example now has rank10 and no
nonzero solution. The local MDS rank theorem remains unchanged.

A bounded exhaustive comparison uses the first two tiles. For each tile,
enumerate all C(8,3)C(8,5)=3136 masks with three input positions and five
output-zero positions. These masks produce1408 distinct residual-state
planes per tile.

| Shift | Identical residual planes across the two tiles |
| --- | ---: |
| 8 | 1408 |
| 64 | 0 |

Shift64 therefore removes this exact duplication. It does not imply that
every pair of planes spans all four state constraints. Among the1,982,464
pairs of distinct-plane catalog entries,4704 have stacked rank3. A concrete
remaining example is

    S = {0,1,2,9,12,15},
    Z = {0,1,2,3,4,8,9,13,14,15}.

Its output matrix has rank9. The surviving one-dimensional space has full
input support, a one-dimensional state image, and a one-dimensional feedback
image. Adding the four feedback-zero equations kills that space.

These counts describe a catalog of subspaces, not the distribution of a
protocol or of source words. Different masks can produce the same plane.
No probability estimate is inferred from the fraction4704/1,982,464.

## Feedback, returns, and existing source dependencies

Let S be a packet input support, and order its ambient variables as
(X_S,Q). Write B for its output-zero matrix, and define

    D = [C_S | 0],      J_M = [C_S | M].

The matrices D and J_M give the feedback and the successor state. For the
unrestricted ambient input/state variables, let K have columns spanning ker B.
Then

    feedback image rank = rank(DK),
    successor-state image rank = rank(J_M K),
    rank([B;D]) = rank(B)+rank(DK),
    rank([B;J_M]) = rank(B)+rank(J_M K).

These are exact joint ranks for a fixed mixer M. They do not assign fresh
randomness to input amplitudes or charge the same source constraint twice.

To retain dependencies from earlier packets and the outer code, represent
their joint input/state image as the column space of a full-column-rank matrix
U. Replace B by BU and take K as a basis of ker(BU). The surviving ambient
vectors are UK, so the feedback and successor images are DUK and J_M UK.
The helper accepts this U directly.

It also checks whether the support class contains a feasible vector. Each
required input coordinate must remain a nonzero linear form on UK. A required
live state must have nonzero image. These conditions are sufficient over
GF256 and GF65536 here: the union of at most |S|+1 proper hyperplanes cannot
cover the surviving vector space. For a return-to-zero query, the same check
is applied after adjoining J_M U to BU.

An abstract bound depending only on the rank of U can lose these gains.
For example, setting U to the two-dimensional solution space of the shift8
defect makes BU identically zero. The ambient output rank8 cannot be charged
again on that source class. The tests include this embedding explicitly.

## Further freedom without an independence claim

The tile output cosets can be distinct aligned cosets delta_alpha instead of
alpha XOR one common beta. Define

    F_alpha=E_delta E_alpha^(-1),
    A_alpha=E_delta[:,0:4],
    C_alpha=E_alpha^(-1)[4:8,:].

Interpolation under translation preserves degree below4, so each tile has
C_alpha A_alpha=0. Also F_alpha^2=I,

    F_alpha A_alpha=E_alpha[:,0:4],
    C_alpha F_alpha=E_delta^(-1)[4:8,:].

Thus inversion swaps the input and output coset families. Distinct evaluation
nodes retain the MDS emission and feedback maps. This is an algebraic design
option, not an assertion that different offsets make residual-state planes
independent. No offset search was performed for this note.

## Scope and next proof step

The exact tile theorem reduces output constraints to eight disjoint input
eliminations and one four-column state calculation. The feedback image then
measures which cancellations persist, inject fresh rank, or permit a return.
This is a usable interface for a source-aware proof.

It is not yet a composable probabilistic transfer. Such a transfer must carry
or safely bound the relevant source image U, expose the current mixer after
the output constraints, and retain the outer MDS equations. The uncapped
whole-code first-moment problem is not resolved by this local improvement.

Shift8 remains the recorded construction. Shift64 is a working option with
an exactly verified improvement against one structured cancellation family.
Neither has a whole-code40% distance certificate at this checkpoint.

## Reproduction

[theory/hd_blockwise_shared_support.py](theory/hd_blockwise_shared_support.py)
implements the local elimination, source embedding, feedback and return
ranks, support feasibility, and bounded residual-plane catalog.

Its six tests pass. They verify all1820 twelve-coordinate information sets of
the tile code at both shifts, local rank formulas, direct-versus-eliminated
packet ranks, the paired cancellation family, a dependent source embedding,
the remaining shift64 counterexample, and the full stated plane comparison.

    python -B -m unittest hd_blockwise_shared_support_tests -v
    python -B hd_blockwise_shared_support.py

No native files, benchmarks, or existing proof receipts were changed.
