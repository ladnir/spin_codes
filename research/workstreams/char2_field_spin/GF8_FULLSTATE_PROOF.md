# Wider persistent state: proof screen at rate three quarters

The width-eight cubic inner did not close a three-quarter-rate outer.
This checkpoint tests whether a larger persistent state repairs that gap without adding another global route.
It also isolates losses from return counting and from the outer's minimum weight.
No candidate in this note has a whole-code certificate.
Negative margins describe unsuccessful upper bounds, not demonstrated low-distance codewords.

## Construction and exact map checks

All scalar operations are in the AES representation of GF(256), with reduction polynomial `0x11b`.
One record holds sixteen scalar coordinates; the following equations act independently on each coordinate.
The packet width is eight and the persistent state has size \(s\in\{7,8,12,16\}\).
Set

\[
A_{u,v}:=(u+v)^{-1},\qquad C_{v,j}:=(v+j)^{-1},
\]

where \(0\le u<8\), \(128\le v<128+s\), and \(64\le j<72\).
Field addition here is XOR.
Both matrices are superregular by the Cauchy determinant formula.

For packet input \(v\) and incoming state \(S\), sample independent diagonals \(D_x,D_s\) with independent entries in GF(256)\(^*\).
Compute

\[
X:=D_xv,\quad R:=D_sS,\qquad Y:=X+2AR,\quad S':=R+CX.
\]

The setup fixes these diagonals separately at every packet position.
The encoder starts with state zero, emits \(Y\), and retains \(S'\).
The outer and single regional route are unchanged from the preceding screens.

The maps do not satisfy \(CA=0\).
Instead, exact checks establish invertibility of both \(I_s+2CA\) and \(I_8+2AC\).
These are equivalent Schur-complement conditions for invertibility of the full packet map.
The exact minimum weights of \((I_8+2AC)X\), for nonzero inputs of weights one through eight, are:

| State size | Returned-output minimum weights |
|---:|---|
| 7 | 8, 6, 5, 4, 3, 2, 2, 1 |
| 8 | 8, 7, 6, 5, 4, 3, 2, 1 |
| 12 | 5, 5, 5, 4, 3, 3, 2, 1 |
| 16 | 8, 7, 6, 5, 4, 3, 2, 1 |

The helper exhaustively checks support feasibility for these returned-output profiles.
It also enumerates every Cauchy minor for the state-seven and state-eight maps.
The general Cauchy determinant identity covers the rectangular state-twelve and state-sixteen maps.
See `theory/gf8_fullstate_screen.py` and `theory/gf8_fullstate_tests.py`.

For state eight, `GF8_FULLSTATE_COST.md` gives an exactly equivalent normalized implementation.
Its operation counts are proposals, not native timings.

## Why the prior transfer argument remains applicable

Fix an exact input support of size \(w\).
For a fixed vector before the current diagonals, let \(h\) be its state-support size.
A nonzero return requires \(h+w>s\), by superregularity of \(C\).
If \(z\) inactive output coordinates are zero, superregularity of \(A\) requires \(h>z\).
Select also \(i\) active output coordinates that must be zero.
The joint return/output equation matrix has rank at least

\[
r\ge\max\{s,\ i+z,\ w-8+i+2z+1\}.
\]

The first term comes from the return equations.
The second uses active-input pivots and inactive Cauchy rows.
The third uses invertibility of the full packet map: only \(8-i-z\) output coordinates remain unconstrained.
These arguments do not require \(CA=0\).

The nonreturn transfer must include all-output-zero events.
In particular, when \(s>8\), even a quiet packet can emit eight zeros without returning to zero state.
The implementation includes root counts through \(\min(s-1,8)\), including eight when applicable.

`theory/gf8_fullstate_rank.py` also retains a cap on the rank of the incoming state forms.
Its coefficient products and regional mixtures preserve the existing cone of nonnegative continuations that increase with that cap.
Directed evaluations use rational marks and outward Arb arithmetic.
They certify only the stated active-row class, not the unexamined classes of a whole code.

## Gates at three-quarter rate

The outer is an MDS \([128,96]\) code over GF(256).
At requested \(K=2^{20}\), natural padding gives 10,944 rows and \(N=1,400,832\) output records.
The excluded output weights are at most 140,083.
The following floating searches use the two-state transfer; \(A\) counts active outer rows.

| Inner | Sparse A=64 | Sparse A=128 | Middle density .03 | Middle density .10 | Middle density .14 |
|---|---:|---:|---:|---:|---:|
| width8/state8 | -901.70 | -2370.66 | -15272.58 | -2228.78 | 32106.81 |
| width8/state16 | 332.97 | -1102.87 | -14888.27 | 375.65 | 38190.55 |

Entries are margins in bits from individual-class upper bounds.
The state-sixteen searches at densities .03 and .10 reached their iteration limit; the displayed witnesses remain valid floating evaluations.
The larger state substantially helps some classes, but leaves both an intermediate sparse gap and a low-density middle gap.
State-eight middle densities .20, .24, and .28 give margins approximately 120,599, 205,639, and 312,276 bits.
The corresponding state-sixteen margins are approximately 131,103, 215,509, and 319,294 bits.

Rank history does not close the sparse gap.
At the state-sixteen A=128 coarse-optimal marks, the seventeen-state transfer with live-kernel return bounds gives a directed margin greater than -1065 bits.
Its approximate margin is -1064.19 bits, versus -1102.87 for the coarse transfer.

## Return bounds without a state-zero-flat sum

The incoming state coordinates need not all remain nonzero for a returned source vector.
Only the \(w\) active input coordinates must remain nonzero.
This permits a direct live-vector moment bound before exposing either current diagonal.

Let \(U\) be the joint space of state/input forms, of dimension \(d\).
For pure returns, each fixed input-live vector survives with probability at most \(255^{-s}\).
The ambient return kernel has dimension at most \(w\).
If a surviving \(j\)-dimensional kernel contains an input-live vector, it contains at least

\[
L_j:=256^j-1-w(256^{j-1}-1)
\]

such vectors.
The expected number of surviving input-live vectors is at most \((256^d-1)255^{-s}\).
An exact two-point moment optimization over \((0,0)\) and \((L_j,\tau^{j-d})\) bounds the return rank potential.
The helper maximizes over every actual dimension up to the known cap.
It does not assume monotonicity in that dimension or concavity of adjacent interpolation segments.

For selected output zeros, the preceding joint-rank bound supplies probability \(255^{-r}\).
Full packet invertibility gives ambient kernel dimension at most \(\min(w,8-i-z)\).
This extension is applied here only when \(s\ge8\), so a zero-state fixed vector cannot return with nonzero input.
The fixed-vector support conditions are used before the diagonals are sampled; they are not conditioning events on those diagonals.

The resulting transfer takes the minimum of the old valid bound and these new pure/joint bounds.
At the old optimized marks, margins become:

| Inner | A=64 | A=128 |
|---|---:|---:|
| width8/state8 | -892.34 | -2348.01 |
| width8/state16 | 343.56 | -1070.60 |

Thus this improvement removes a real counting loss but not the main sparse obstruction.
An exhaustive GF(4) regression checks the joint probability and kernel-cap hypotheses in 960 fixed-vector/output-selection cases, over every nonzero diagonal.
See `theory/partial_live_rank.py`, `theory/gf8_fullstate_partial_live.py`, and their tests.

## Outer-weight mark and fixed-word controls

For \(A\) active MDS rows, every admissible source support has weight \(W\ge(129-96)A\).
For \(\lambda\ge1\), multiplying its nonnegative contribution by \(\lambda^{W-33A}\) gives an upper bound.
The regional degree-\(h\) coefficient therefore receives \(\lambda^h\), and the final bound receives \(\lambda^{-33A}\).
The rank mark \(\tau\) is unchanged.
This introduces no randomness or independence assumption.

At all four optimized state-eight/state-sixteen A=64/A=128 witnesses, the best added mark is \(\lambda=1\).
The derivatives in bits per unit \(\log\lambda\) are approximately 1743, 3051, 1812, and 3386.
The relaxed tilted support already lies above the total minimum weight, so this particular restriction does not help.
See `theory/gf8_fullstate_sparse_tilt.py`; its tests check the positive-coefficient inequality exactly.

A separate direct fixed-word first moment uses source coefficient \(256^W\), parity penalty \(256^{-32A}\), and one projective normalization by 255.
Its local zero probabilities use powers of \(255^{-1}\), without rank-potential or zero-flat losses.
This control is much weaker: state-eight A=64/A=128 margins are -19177/-41226 bits, and state-sixteen margins are -15084/-37132 bits.
The optimized small rank mark is therefore providing essential first-moment savings.
See `theory/gf8_fullstate_fixedword.py`.

## A larger outer is insufficient by itself

For an MDS \([256,192]\) outer, natural padding at K20 gives \(N=1,409,024\).
The state-sixteen partial-live screen gives sparse margins 1406.97, -743.64, and -4303.74 bits at A=64, 128, and 256.
Middle densities .10 and .14 give 4021.48 and 42511.39 bits; the .10 optimization reached its iteration limit.
This larger outer improves some gates, but does not remove the intermediate sparse gap.

The next construction change should strengthen birth/output mixing, or replace the remaining nonreturn estimate.
Merely increasing state size, removing the state-zero-flat sum, or enforcing total outer minimum weight has not closed this family.
