# Uniform-entry state refresh: a local backup, not a certificate

A fresh uniform-entry four-by-four state matrix has a useful conditional return-rank bound. Its canonical branching factor is `1 + tau^-1 + tau^-2 + tau^-3`, rather than four. The complete product estimate rises from `169/88` to `193/88` GFNI512 instructions per source record, about 14.20%.

This note preserves a backup construction and a local upper action. It does not establish a whole-code bound or select the backup over diagonal refresh. The half-rate, 10% target is already supplied by the scalar-extension control and is not a new success criterion.

The later unequal-window audit in [SPARSE_MULTIROW_REVIEW.md](SPARSE_MULTIROW_REVIEW.md) shows that this four-state backup cannot repair the two saved bad routes. Their dimension witness applies to every linear choice of the state matrices. The local lemmas below remain valid, but further work on this backup is paused pending a change that addresses the route-capacity obstruction.

## Packet interface and a changed empty-packet rule

Keep the fixed Cauchy maps `A,C`, scalar `eta=2`, input diagonals, outer code, and one regional routing stage. Replace the four nonzero state diagonals by a matrix `G` with sixteen independent uniform GF256 entries. Different packets use independent matrices. Apply

\[
 X=D_{\rm in}\ell,\qquad R=GS,\qquad
 Y=X+\eta AR,\qquad S'=R+CX.
\]

The matrix may be singular. The stream encoder remains invertible because its current-input coefficient is the identity: given the earlier inputs and current state, recover `X` from `Y`. In fact `CA=0` still makes the zero-initial-state inner map an involution with the same fixed setup. No inverse of `G` is used.

Unlike nonzero diagonal refresh, this map can erase a nonzero incoming state in an empty-input packet. For each fixed nonzero `S`, the vector `GS` is uniform on GF256^4, so this event has probability `256^-4`. The diagonal construction's rule forbidding empty-input returns must not be carried over.

## Common-history rank potential

Fix the previous homogeneous equation space on the shortened source and all previously exposed randomness. Work on its quotient. A branch requires every prescribed active input form to be nonzero on a witnessing source. A nonzero-state branch requires at least one incoming state coordinate to be nonzero on that witness. The potential for an additional equation rank `r` is `tau^-r`, with `1 <= tau <= 255`.

The independent-diagonal pivot lemma remains valid when the active source forms are linearly dependent. Its proof and the necessary alive predicate are given in the audit section of [SPARSE_MULTIROW_REVIEW.md](SPARSE_MULTIROW_REVIEW.md). For a scalar sampled from all of GF256, the analogous one-pivot bound holds up to `tau <= 256`. We retain the common range `tau <= 255` because input diagonals remain nonzero.

### Canonical nonzero-state cover

Suppose the incoming state forms span an `r`-dimensional space in the current quotient, with `1 <= r <= 4`. Select a deterministic ordered basis `s_1,...,s_r` from these forms. Cover every witness with nonzero state by its first nonzero basis form:

\[
 s_1=\cdots=s_j=0,\qquad s_{j+1}\ne0,
 \qquad 0\le j<r.
\]

The `j` imposed zero forms are independent before this split. They cost `tau^-j`. Branches that kill a required input are discarded. Hence the total deterministic branching factor is

\[
 \Pi_r(\tau):=\sum_{j=0}^{r-1}\tau^{-j}
 \le\Pi_4(\tau).
\]

This is a cover of alive equation spaces, not an assumed distribution on incoming states. At `tau=35`, `Pi_4` is less than 1.03.

On one branch, express the original state coordinates in the remaining basis. Each row of the uniform matrix `G` induces a uniform coefficient vector in that basis. Expose every coefficient except the coefficient of the required nonzero form `s_(j+1)`. The four remaining coefficients are independent and uniform. Each return equation consumes one such scalar pivot. They jointly contribute `tau^-4`.

## Selected-zero return action for at most three inputs

Fix the input placement `I`, of size `w` in `{1,2,3}`. Let `E_I` embed those input coordinates into the packet, and define the returned-output matrix

\[
 B_I:=E_I+\eta AC_I.
\]

On a return, `R=C_I X_I` and therefore `Y=B_I X_I`. For a selected output-zero set `Z`, first impose

\[
 (B_I)_Z D_I\ell_I=0.
\]

These equations consume the current input diagonals and give the alive factor `tau^-rank((B_I)_Z)`. A selected set is impossible if this numerical matrix has no kernel vector with every input coordinate nonzero.

After imposing those equations, choose the canonical basis of the remaining incoming state forms. Apply the preceding first-nonzero cover and four fresh matrix-row pivots to impose the return. The state-form rank may have decreased, so use `Pi_4` as an upper bound. This exposure order is essential: it does not claim that the input diagonals remain fresh after conditioning on return equations.

For `0 < y <= 1`, the resulting conditional upper action is

\[
 \boxed{\Pi_4(\tau)\tau^{-4}
 \sum_{\substack{Z\subseteq[16]\\
       \ker (B_I)_Z\ \text{contains an all-nonzero vector}}}
 y^{16-|Z|}(1-y)^{|Z|}\tau^{-\operatorname{rank}(B_I)_Z}.}
\]

The sum is the positive selected-zero expansion, not a count of field-valued words. The existing exact projective arrangement is used only to determine which selected sets are feasible. Each feasible `Z` appears once, regardless of the number of projective inputs or their injection-support masks.

For one input, every returned output is nonzero. The bound reduces to `Pi_4*y^16/tau^4`. For zero inputs, a return makes every output zero, giving `Pi_4/tau^4`. Averaging placements is valid only when the enclosing route average supplies a fresh uniform placement conditional on the packet weight.

## Marginal action and safe empty continuation

A companion marginal bound follows from the same canonical split. If `i` selected zeros occur at active-input positions and `z` at zero-input positions, then

\[
 \operatorname{rank}[\eta A_Z,E_{Z,I}]=i+\min(z,4).
\]

The active-coordinate columns supply `i` independent pivots. Cauchy superregularity supplies the remaining `min(z,4)` state pivots. Joint exposure of fresh matrix-row coefficients and input diagonals therefore gives

\[
 \Pi_4\sum_{i=0}^{w}\sum_{z=0}^{16-w}
 \binom wi\binom{16-w}z
 y^{16-i-z}(1-y)^{i+z}\tau^{-i-\min(z,4)}.
\]

This marginal includes returns. It must not be multiplied by the separate return action, since both use the same randomness. In an empty packet whose next state is required nonzero, four selected output zeros are impossible. Its nonreturn action may therefore restrict the sum to `z < 4`, with charge `tau^-z`. The separate empty-return entry remains positive.

Repeated use of the conservative `Pi_4` cover can itself cause substantial loss. These formulas establish local interfaces, not evidence that their global composition will close.

## Cost and bounded verification

Pack the four state records in one ZMM register. A uniform four-by-four matrix-vector map uses four GFNI512 products: broadcast one input record per product, apply the four column coefficients, then combine the four results. The current diagonal map uses one product. The change therefore adds three products and three vector XORs per sixteen-record packet, plus the required broadcasts.

The full estimate for GRS[32,22] becomes

\[
 \frac{41}{88}+\frac{19}{16}\frac{32}{22}
 =\frac{193}{88}\approx2.19318.
\]

State coefficients increase from four to sixteen field bytes per packet, before SIMD expansion. Sampling requires no invertibility test or rejection. Register scheduling, coefficient expansion, and complete-call runtime have not been measured.

The scoped reference is [full_state_refresh.py](theory/full_state_refresh.py). Its five tests pass. They check the canonical factor and empty-return change, independently enumerate feasible two-input zero sets, verify the marginal coefficient ranks, and exhaust all 256 two-state GF4 refresh matrices with all current nonzero input diagonals for dependent source-form examples. The last test distinguishes empty nonreturn from the full marginal.

```sh
python -B research/workstreams/char2_field_spin/theory/full_state_refresh_tests.py
```

No global screen, directed distance certificate, native implementation, benchmark, or remote run accompanies this backup. The diagonal-mask composition was subsequently evaluated and failed to close. The window obstruction then showed that a change to outer/state/routing geometry is required before revisiting this local law.
