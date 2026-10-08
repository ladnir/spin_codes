# Two-active-row certificate for the Cauchy diagonal-state construction

The two-active-row class has directed failure bounds below `2^-65`, `2^-72`, and `2^-77` at the padded K16, K18, and K20 sizes. These bounds include every pair of active outer rows and every nonzero word in each row. They do not certify the whole code.

The construction uses fixed GRS[32,22] rows over GF256, sixteen-coordinate packets, four state coordinates, independent nonzero input and state diagonals, and one regional routing stage. This note retains the fixed-word packet envelope from [the first-moment note](CAUCHY_FIRSTMOMENT_GATE.md). It sharpens the route and source counting, not the construction.

## Packet and route interfaces

Fix both outer words before exposing the route and packet diagonals. Let state 0 mean the actual zero state and state 1 mean an arbitrary nonzero state. For input weight `c`, let `T_c(y)` be the positive two-state upper matrix for the output-weight moment, where `0 < y <= 1`. The bound holds uniformly over incoming state values and nonzero input values. Thus its products remain upper bounds after conditioning on previous packets. The input diagonal is exposed only in this packet bound; the outer word count below is deterministic.

Only `c=0,1,2` occurs with two active outer rows. In the notation of the first-moment note,

\[
 T_0=\begin{pmatrix}1&0\\0&M_0\end{pmatrix},\qquad
 T_1=\begin{pmatrix}0&y\\y^{16}/255^4&M_1\end{pmatrix},\qquad
 T_2=\begin{pmatrix}0&y^2\\R_2&M_2\end{pmatrix}.
\]

The one-input return has exactly sixteen nonzero outputs. The two-input return uses the original conservative positive bound

\[
 R_2=\sum_{i=0}^{1}\sum_{z=0}^{3}
 \binom2i\binom{14}z
 \frac{y^{16-i-z}(1-y)^{i+z}}{255^{i+z+2-i}}.
\]

The entries `M_1,M_2` are full marginal moments, so they also count returns. This deliberate overcount is safe. No sharper joint two-input return result is needed for this certificate.

Write `L` for the number of outer rows and `H=L/16` for the packets in a region. If exactly `c` of the two outer coordinates are nonzero, their route slots form a uniform `c`-subset of the `L` positions. Define

\[
 P(z)=T_0+16zT_1+\binom{16}2z^2T_2,
 \qquad
 U_c=\binom Lc^{-1}[z^c]P(z)^H,\quad c\in\{0,1,2\}.
\]

The coefficient sums all packet occupancies in chronological order. Its binomial factors count lane choices. Dividing by `binom(L,c)` gives the exact route average of the conditional upper matrices. With two labeled rows, every unordered slot pair has two equally likely assignments, and the packet bounds apply to either assignment. In particular, the same-packet probability is exactly `15/(L-1)`; it is not replaced by an independent-placement approximation.

The regional permutations are independent. A regional matrix remains a conditional upper bound when previous regions have been exposed because it is uniform over incoming values. All four phase maps obey the same packet bounds.

## Exact source counting

For one specified support of size `w >= 11`, the number of projective GRS[32,22] words is the integer

\[
 n_w=\sum_{j=0}^{w-11}(-1)^j\binom{w-1}j256^{w-11-j};
 \qquad n_w=0\quad(w\le10).
\]

Use separate support markers `u,v` for the two labeled active rows, and define

\[
 F(u,v)=e_0\bigl(U_0+uU_1+vU_1+uvU_2\bigr)^{32}\mathbf1.
\]

The initial state is zero; the final state is discarded, so both terminal states are included. The coefficient `[u^a v^b]F` sums all ordered support pairs with weights `a,b`. It retains the position of every active region while composing the two-state matrices.

Each row has `255 n_w` ordinary words on its specified support. A pair of nonzero rows therefore contributes `255 n_a n_b` after dividing by one global scalar action. Dividing separately by a scalar for each row would be incorrect. For `e=ceil(32L/10)-1`, a complete two-active-row failure upper bound is

\[
 \boxed{\binom L2\,255\,y^{-e}
 \sum_{a=11}^{32}\sum_{b=11}^{32}
 n_a n_b\,[u^a v^b]F(u,v).}
\]

This is a first-moment bound over global projective source words, followed by the union over the two active row locations. Padding enlarges the counted source space and is conservative.

## Frozen directed replay

All polynomial arithmetic is positive. The regional exponentiation truncates only degrees greater than two, which cannot contribute to the retained coefficients. The 32-region recursion uses exact integer support multiplicities. Every noninteger arithmetic operation is enclosed by the existing 80-digit outward-rounded decimal interval type.

| Requested source size | Padded rows `L` | Frozen rational `y` | Failure upper, rounded upward here | Certified margin |
|---:|---:|---:|---:|---:|
| `2^16` | 3008 | 4979/5000 | `1.3600e-20` | 65 bits |
| `2^18` | 11968 | 49947/50000 | `1.3368e-22` | 72 bits |
| `2^20` | 47680 | 24993/25000 | `4.9877e-24` | 77 bits |

The preceding floating search gave margins 66.0047, 72.6637, and 77.4257 bits. Those search values select witnesses; they are not the certificate. The directed witness loses less than one bit at each size. The dominant floating support pair is `(13,13)` at K16 and `(11,11)` at K18/K20.

```sh
python -B research/workstreams/char2_field_spin/theory/cauchy_q2_transfer.py --certify
python -B research/workstreams/char2_field_spin/theory/cauchy_q2_transfer_tests.py
```

Eleven tests pass. Eight checks cover the certificate: exhaustive distinct-slot routing on small regions, exact same-packet probabilities, two-marker support enumeration, the single global projective quotient on a GF4 MDS code, directed packet enclosures, an independent complete small Fraction expression, the floating log implementation against positive arithmetic, and all three frozen witnesses. Three further checks cover the exploratory extension below.

## Scope and scalable continuation

This certificate covers exactly two active outer rows. Together with the separate one-row certificate, it removes the first two row-activity classes. It does not fill the sparse multirow gap, and it is not evidence that an ordinary word first moment will close every class.

For `A` active rows, exact route averaging still has a simple form:

\[
 U_c=\binom Lc^{-1}[z^c]
 \left(\sum_{j=0}^{16}\binom{16}jz^jT_j\right)^H,
 \qquad 0\le c\le A.
\]

An exact source transfer could assign a marker to every row and retain each row's degree through all 32 regions. Its dense storage grows as `33^A`, so enumerating successive small values of `A` is not a scalable whole-class proof. The following executed calculation retains the union of active regions and the exact nonempty MDS row polynomial instead.

## Executed source-union envelope

Fix `A` active row labels and a scalar `theta > 0`. Let `S_A(j)` be the weighted sum over ordered `A`-tuples of row supports whose union is one specified `j`-element region set. A row of support size `w` has weight `n_w theta^w`. In particular, each row separately has support size at least eleven. Set `S_0(0)=1` and all other entries of `S_0` to zero. Then

\[
 S_{a+1}(j)=\sum_{\ell=0}^j K_{j,\ell}S_a(\ell),
 \qquad
 K_{j,\ell}=\binom j\ell
 \sum_{i=0}^{\ell}\binom\ell i
 n_{j-\ell+i}\theta^{j-\ell+i}.
\]

To verify the recurrence, choose the previous support union as an `ell`-subset of the target `j`-set. The new row must contain all `j-ell` new regions, and may contain any `i` regions in the previous union. This is a positive 33-state recurrence. Evaluating `K^A S_0` does not enumerate row-support profiles or assume independent fresh flags.

Keep `E=U_0` unchanged and define the entrywise upper matrix

\[
 B_{r,s}=\max_{1\le c\le A}\theta^{-c}(U_c)_{r,s}.
\]

For a fixed support tuple, every occupied region contributes `theta^c B` and every empty region contributes `E`. Positivity allows entrywise substitution. Summing the exact source weights by their union gives the finite upper expression

\[
 \binom LA255^{A-1}y^{-e}
 \sum_{j=0}^{32}S_A(j)
 [z^j]e_0(E+zB)^{32}\mathbf1.
\]

The factor `255^(A-1)` again comes from one global projective quotient, not one quotient per row. The route matrices `U_c` are extracted exactly in their algebraic definition. Only the entrywise monomial bound loses their dependence on `c`. The implementation evaluates this expression in floating logarithmic arithmetic; none of the following values is a directed certificate.

| Active rows `A`, K20 | Best executed margin | Dominant union size | Count maximizing every entry of `B` |
|---:|---:|---:|---:|
| 2 | 77.2611 bits | 11 | 2 |
| 3 | 41.3873 bits | 27 | 3 |
| 4 | -49.0041 bits | 31 | 4 |
| 8 | -424.0726 bits | 32 | 8 |
| 16 | -1141.1325 bits | 32 | 16 |

This envelope loses only about 0.165 bits against the exact two-row K20 screen, and gives a prospective three-row witness. It does not resolve the larger sparse classes. A negative margin means this upper expression is unusable, not that the construction fails. Each entry of `B` is selected by the largest allowed regional count; exact row nonemptiness and union chronology alone do not control that relaxation. These screens use the original conservative multi-input return bounds. Improved joint returns can be substituted after checking their placement conditioning, but they cannot remove the first-moment obstruction below.

```sh
python -B research/workstreams/char2_field_spin/theory/cauchy_q2_transfer.py --k-log 20 --union-active 2 3 4 8 16
```

The additional tests enumerate the source-union recurrence on small supports, check that its monomial envelope dominates the complete two-row expression, and verify the following obstruction by exact rational comparison.

## Why the complete sparse proof needs more than a word first moment

At K20, fix sixteen outer row labels and take an aligned suffix window of width `W=4752` in every region. The event that all sixteen routed coordinates lie in the window in every region has probability

\[
 p=\left(\frac{\binom{4752}{16}}{\binom{47680}{16}}\right)^{32}.
\]

There are `22*16` source dimensions. Imposing zero terminal state after each of the first 31 regions gives at most `4*31` linear equations. For every choice of packet diagonals and route satisfying the event, the remaining source kernel therefore has dimension at least 228. Its outputs vanish outside the chosen windows: every region begins in zero state, and every input before its suffix window is zero. The last state is discarded. The encoder is injective, and `32W=152064 <= e=152575`, so each nonzero kernel source gives a nonzero bad output word.

Consequently, the expected number of bad global projective words supported on these sixteen row labels is at least

\[
 p\,\frac{256^{228}-1}{255}>2^{111}.
\]

The exact rational comparison passes; the floating logarithm is approximately 111.6447. These kernel words may use fewer than sixteen rows. Thus the witness concerns the union of classes `1 <= A <= 16`, not exactly sixteen active rows. Most importantly, this is a lower bound on an expected word count, **not** a lower bound on bad-code probability. A rare route may contribute many words from one subspace.

The next useful step is to combine the exact regional route/source counting with an event- or rank-based bound that counts such a subspace once. Certifying three or four more individual row counts does not by itself bridge the remaining sparse class.
