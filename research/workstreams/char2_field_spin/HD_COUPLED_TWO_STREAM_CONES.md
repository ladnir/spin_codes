# A five-cone interface for coupled interleaved states

Pure interleaving keeps two independent state trajectories. A source can
still place its inputs in only one trajectory. The alternative below
couples the states without adding another emission, projection, or refresh.
Isolated packet and regional proof helpers now implement this proposal.
Only exact small tests have run. No production regional bank, native
kernel, timing result, or distance receipt is supplied.

## Packet map and implementation cost

Fix the existing scalar field `F=GF(256)`, state dimension 16, packet width
64, and maps `A,C` with `CA=0`. Let `M` be the existing fresh invertible
MDS diagonal refresh. Replace the one-state packet recurrence by

\[
Y=X+AU,\qquad (U',V')=(V,M(U+V+CX)).
\]

The initial state is `(0,0)`. Both states persist across region boundaries;
there is no reset or terminal flush. Each packet still emits 64 symbols.

The native schedule can retain two named physical `State` objects and
alternate which one emits. Before refreshing that object, XOR the other
state into `active+feedback`. Swapping the logical roles then implements
the displayed recurrence. This needs four additional ZMM XORs per packet,
with no additional variable GF product, refresh key, route load, or output
pass. At K20 there are 32768 packets, so the static addition is 131072
vector XORs. Persistent storage rises from four to eight ZMM vectors.

This does not establish unchanged runtime. The extra state can cause spills,
and the coupling restores a dependence between consecutive refreshes.
The current approximately 3.398-ms result leaves little literal 3.400-ms
headroom. A later implementation should first use two explicit inlined
packet calls per loop iteration, retaining contiguous stores. A separate
pass or dynamically indexed state array is unnecessary. Code generation
and serial A/B timing would be required only after a useful proof gate.

For empty input the state transition is invertible: its inverse recovers
`V=U'` and `U=M^-1 V'+U'`. Two consecutive empty emissions are `AU` and
`AV`. Since `A` is injective, they cannot both vanish for a nonzero state
pair. This is an observability fact, not a universal two-return probability.
Input words with `CX=0` at every packet still leave the zero state unchanged;
the previously identified quiet subspace has not disappeared.

## Why the old four product cones are insufficient for the old caps

Write `D=256^16-1`, and let `mu` be the uniform measure on nonzero states.
Fix a packet input with `c=CX!=0`. From input state `(0,V)` with `V` drawn
from `mu`, the event `V=c` has probability `1/D`. The next state is then
`(c,0)`. Its mass is small, but its first coordinate is a point mass.
Relative to `mu` its density coefficient is one, not `1/D`.

Thus the previous flat-return factor cannot be inserted unchanged into a
transition from a one-live product cone to the other. Averaging a one-symbol
nonzero input does not repair this: `CX` then ranges over only 255 vectors,
not over all `D` nonzero states.

The fifth cone below keeps the small total mass without pretending it is
flat. A four-cone comparison could instead inflate the density coefficient,
but that discards the intended advantage.

## Positive measures and conditioning

Fix the complete scalar outer support word and all constituent row sets.
At the current packet let `S` be its support, `w=|S|`, and let `X` have
independent uniform nonzero entries on `S`, with zero entries elsewhere.
The fresh outer output diagonals give this law independently of the past,
even after the past is weighted by its emitted-weight marks. The current
inner refresh is also fresh. These are the conditioning assumptions needed
by the interface.

For `0<z<=1`, propagate positive state measures with weight `z^wt(Y)`.
Use five normalized cone types:

| Type | Reference or family of measures on `(U,V)` |
|---|---|
| `Z` | the point mass at `(0,0)` |
| `L` | `delta0 tensor mu` |
| `R` | `mu tensor delta0` |
| `B` | `mu tensor mu` |
| `P` | any probability measure on `{(u,0):u!=0}` |

A coefficient of `P` bounds total mass, not density relative to `mu`.
The other coefficients are the established pointwise density bounds
relative to their normalized references. The cone family is closed under
positive mixtures, so a newly created `P` measure may depend on the current
input. At the next packet it is part of the past. Fresh current amplitudes
then justify taking a supremum over its first coordinate. There is no extra
union over possible coordinates or row sets.

Let `beta=D/255^16`. The existing MDS refresh bound says that, for every
fixed nonzero refresh input, each nonzero output atom is at most `255^-16`.
All transitions below use this pointwise bound, not an assumption that the
two actual states are independent.

## A sound packet transition interface

For the fixed support `S`, put `c=CX` and `q_S=Pr[c=0]`. Define

\[
J_S=\mathbb E_{U\gets\mu,X}z^{\mathrm{wt}(X+AU)},
\]

\[
K_S=\sup_{v\ne0}\mathbb E_X
 [1_{v+c\ne0}z^{\mathrm{wt}(X+A(v+c))}],
\]

\[
H_S=\sup_{u\ne0}\mathbb E_X
 [1_{c=u}z^{\mathrm{wt}(X+Au)}],\qquad
L_S=\sup_{u\ne0}\mathbb E_X
 [1_{c\ne u}z^{\mathrm{wt}(X+Au)}].
\]

Let `T10(S),T11(S)` be valid existing one-state raw packet caps from an
incoming flat nonzero state. Those two transitions are literally the same
map when the new recurrence starts in cone `R`. With row/column order
`Z,L,R,B,P`, the following nonnegative matrix is valid:

\[
\begin{pmatrix}
q_Sz^w & \beta(1-q_S)z^w & 0 & 0 & 0\\
0 & 0 & 0 & \beta z^w & (1-q_S)z^w/D\\
T10(S) & T11(S) & 0 & 0 & 0\\
0 & 0 & K_S/D & \beta J_S & 0\\
H_S & \beta L_S & 0 & 0 & 0
\end{pmatrix}.
\]

For `L`, a cancellation leaves the atom `c` in the unrefreshed first half;
its total mass is `(1-q_S)z^w/D`. Otherwise the uniform first half and the
pointwise refresh cap yield the stated `B` density. This bound need not
subtract the cancellation mass from its `B` entry.

For `B`, condition on its reference value `V=v`. Cancellation requires the
single value `U=v+c`, producing density `K_S/D` relative to the remaining
uniform first coordinate. On the noncancelling branch, discard the exclusion
inside the marked `U` average and use the fresh refresh cap. This gives
`beta J_S` uniformly in `v`.

For `P`, a fixed input coordinate `u` goes either to zero or to cone `L`.
The displayed suprema dominate every positive mixture defining `P`.
Thus `P` lasts only one packet; it need not become another persistent
geometric family.

The current empty-input inverse-refresh sharpening of `T11` is valid for
`R -> L`, where the refresh input is `U`. It cannot simply be reused for
`B -> B`, whose refresh input is `U+V`. The latter needs `beta J_S` or a
separately proved stronger bound. Likewise `B -> R` needs the displayed
supremum, not the old averaged return entry.

## Implemented occupancy-only bounds

The packet helper uses the following explicit majorants. They depend on
the four block occupancies `l_i`, not on a resampled support. Set
`w=sum_i l_i`, `b0=256^-4`, and `h=z+(1-z)/255`.

First, for `t=16-l`, define the positive sum

\[
E_l=\sum_{j=0}^{\min(3,t)}
\binom tj z^{t-j}(1-z)^j(256^{-j}-b_0).
\]

For a uniform four-dimensional state block, the marked expectation on its
`t` coordinates outside the input support is exactly `b0+E_l`. This follows
by expanding zero-coordinate indicators and using the MDS rank of each
set of selected rows. Conditioning the inputs in `S` to be nonzero costs
at most `(256/255)^w`. Excluding the all-zero sixteen-dimensional state
therefore gives

\[
J_S\le\overline J_S:=\min\left\{1,
(z+1/255)^w
\frac{\prod_i(b_0+E_{l_i})-b_0^4}{1-b_0^4}\right\}.
\]

The implementation evaluates the numerator with a positive recurrence,
not floating-point subtraction. Starting with `T=1,N=0`, each block updates
`N <- N*b0+T*E_l`, then `T <- T*(b0+E_l)`. This is the unsharpened marked
flat-state bound. The production helper does not import the old
inverse-refresh empty-input sharpening into `B -> B`.

For a fixed nonzero state block, an active input coordinate has expected
mark `h` if that state emission is nonzero, and `z` otherwise. Combining
this identity with emission distance 13 bounds that block by

\[
g_l=h^l z^{\max(0,13-l)}.
\]

A zero state block contributes `f_l=z^l`. Consequently

\[
L_S\le\overline L_S:=
\min\left\{1,\max_{\varnothing\ne I\subseteq\{0,1,2,3\}}
\prod_{i\in I}g_{l_i}\prod_{i\notin I}f_{l_i}\right\}.
\]

The exclusion `CX!=u` is discarded in this bound. The maximization also
discards compatibility between the selected codeword supports. Both are
explicit slack, not claims of attainable emission patterns.

For any nonzero incoming state, define the deterministic emission floor

\[
d(S)=\min_{\varnothing\ne I\subseteq\{0,1,2,3\}}
\left(\sum_{i\notin I}l_i+
       \sum_{i\in I}\max(0,13-l_i)\right).
\]

On the `B -> R` branch, the incoming first state is `v+CX!=0`; hence this
floor applies even though that state depends on the input. Also `I+AC` is
invertible because `(AC)^2=0`. Its columns restricted to `S` have rank `w`.
Selecting `w` independent output rows, discarding the other output marks,
and relaxing the nonzero input restriction gives the affine-image bound
`(z+1/255)^w`. Thus the implemented supremum majorant is

\[
K_S\le\overline K_S:=\min\{1,z^{d(S)},(z+1/255)^w\}.
\]

For `P -> Z`, put `rho=rank(C_S)=sum_i min(4,l_i)`. An information set for
the affine system `C_S X=u` leaves `w-rho` free input coordinates, each
nonzero. Every fixed target consequently has probability at most
`255^-rho`. On this return branch, the nonzero state blocks form a set
`I` satisfying

\[
\{i:1\le l_i\le4\}\subseteq I\subseteq\{i:l_i>0\},\qquad I\ne\varnothing.
\]

Every block in `I` also has `CY_i=CX_i=u_i!=0`, so its output weight is at
least `max(1,13-l_i)`. Let

\[
d_{\rm ret}(S)=\min_I\left(
\sum_{i\notin I}l_i+\sum_{i\in I}\max(1,13-l_i)\right)
\]

over these sets. For nonempty input the implemented return majorant is

\[
H_S\le\overline H_S:=\min\{1-q_S,255^{-\rho}\}z^{d_{\rm ret}(S)}.
\]

The kernel probability `q_S` is the product of the exact full-support MDS
kernel probabilities for the four blocks. For one nonzero input,
`H_S<=z^12/255`. This factor follows the `1/D` already carried into `P`;
it is not another unconditional `1/D` factor. Empty input has no feasible
return set `I`, so `H_empty=0` and `L_empty<=z^13`.

The generic rational helper also uses
`T10<=(1-q_S)z^d_ret/D` and `T11<=beta*Jbar`. The production FLOAT table
takes their entrywise minimum with the established old raw `T10,T11` caps,
only in the `R` row. All new entries remain the displayed explicit bounds.

These bounds establish a closed five-cone packet interface. They do not
establish a favorable regional or whole-code bound. Any later provider
must average the original paired placement jointly, retain ordered
five-state products, and preserve the fixed outer row sets through the
whole-word support floor. The two-stream Kronecker factorization does not
apply to this coupled recurrence.

## Implementation and review status

`theory/hd_coupled_state_packet.py` exports:

- `packet_caps(occupancy,z,...)`: exact rational flat 25-entry caps for
  generic block parameters;
- `packet_details(occupancy,z,...)`: the exact intermediate bounds;
- `raw_packet_table(z)`: production FLOAT 5-by-5 caps keyed by four-block
  occupancy, with no state scaling or fractional power.

The state order is `Z,L,R,B,P`. Every cone has terminal mass functional one,
but a useful potential need not assign `P` a value at most one. Its mass
normalization differs from the flat references. No product-potential
restriction is imposed by this packet interface.

Seven small tests pass in 1.071 seconds, including import overhead. The
GF(4) fixture uses scalar state, width-two repetition emission/projection,
`CA=0`, and a uniform nonzero scalar refresh. It enumerates every state,
input support, nonzero input assignment, and refresh. Checks cover each
reference cone, every point defining `P`, weighted mixtures of those
points, and the fixed-input cancellation atom. Additional tests check the
positive emission bound, rank/marked return floors, and that old-row reuse
does not change `B -> B`. No production packet table was built by these
tests; the table-wiring test is mocked.

Independent source review cleared the packet bounds. The root subsequently
passed all 15 packet and regional tests in 0.043 seconds. The independent
regional helper `hd_coupled_state_regions.py` retains ordered five-state
polynomial products and the single paired-route normalizer. It does not
reuse the interleaved tensor factorization.

| File in `theory` | SHA-256 |
|---|---|
| `hd_coupled_state_packet.py` | `faab3a80e582f0c1864b79f63b256788d99a5034f15bfbd9e781940b43bc2a43` |
| `hd_coupled_state_packet_tests.py` | `0d9f37b6f9f52a990691f7fca6d117f15e776fc452dd6eb50cb588616b1ff4bf` |

## Next decision

The fixed pure-interleaving comparison did not close. The coupled recurrence
now has exact small transition tests and merits one bounded regional
proof control before native work. The decisive questions are the
size of `P -> L` and the unsharpened `B -> B` term, not merely the extra
four-XOR operation count. Existing constructions and receipts stay frozen.
