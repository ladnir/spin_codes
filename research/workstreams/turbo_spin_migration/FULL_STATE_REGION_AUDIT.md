# Full-state regional comparison audit

This note audits Lemma `lem:shared-inner-fixed-region` in
`research/paper/structured_uniform_appendix.tex`. It checks the comparison
before the sparse and fixed-occupation numerical gates. It is a conventional
argument, not a Lean verification.

## Exact process and boundary states

Fix the number `a` of bits in one region, their distinct row identities,
and arbitrary entering and exiting states. The region shuffle puts those
bits at uniformly sampled distinct slots. The state update is
`Q' = M Q + C X`, the output is `Y = X + A Q`, and each step has an
independent uniform `M` in `GL(22,2)`. The assumptions used below are
`C e_j != 0` for every input coordinate and that every output coordinate
of `Aq` is a nonzero linear functional of `q`.

1. **Empty gap.** For a live entering state `q` and `g>=1` empty steps,
   the terminal state is uniform live and its exact weighted entry is
   `z_L^wt(Aq) phi_L^(g-1)/D`, with
   `phi_L = D^-1 sum_{v!=0} z_L^wt(Av)` and `D=2^22-1`.
   Zero stays zero. The uniform-live mean of `wt(Av)` is `64 p0`, where
   `p0=2^21/D`. Hoeffding's lemma bounds the logarithmic error in replacing
   the live gap by `D_gamma(64g/L)`, including a one-step gap and any
   entering live state.
2. **One-bit impulse.** At tilt one, zero becomes live because `C e_j!=0`.
   A uniform live input returns to zero with probability exactly `1/D`:
   `M q` is uniform live and one value cancels `C e_j`. Thus
   `R V_j(1) J = [[0,1],[1/D,1-1/D]]` for every bit coordinate `j`.
   At `z_L=exp(-theta/L)`, replacing `V_j(z_L)` by `V_j(1)` costs at most
   `exp(theta*64/L)` multiplicatively.
3. **Separated placements.** A boundary impulse or two impulses in the
   same or neighboring steps has probability at most
   `2a*64/L + 3*64*binom(a,2)/(L-1)`. Otherwise all `a+1` gaps are
   nonempty. Their exact full-state product is bounded by
   `exp(+-E_a/L)` times the lifted two-state product, with
   `E_a=(a+1)c_theta+a*theta*64` and
   `c_theta=theta*64+theta^2*64/8`. The first and last gaps justify the
   uniform lift even for arbitrary boundary states.
4. **Continuum coupling.** Quantizing `a` labeled independent uniform
   points in `[0,1]` gives the exact uniform distinct-slot law when no
   two land in one slot. Collision probability is at most
   `binom(a,2)/L`. On the coupled good event, the spacing-vector error
   is at most `2(a+1)(64+1)/L` in `ell_1`. The lifted product is
   `gamma`-Lipschitz in that vector, since its other factors are
   stochastic or substochastic. The paper now states the resulting
   explicit additive entrywise bound. Every term is `O_Q(1/L)` for fixed
   `Q` and `theta`.
5. **Multiplicative form.** The only exact zeros are cross-class entries
   for `a=0` and zero-to-zero for `a=1`; the exact regional kernel has
   the same zeros. Every other limiting entry is positive. For fixed
   `Q`, the minimum positive entry over `a<=Q` and both boundary states
   is a positive constant. Dividing the additive error by it gives the
   entrywise factor `1+C_{Q,theta}/L`.
6. **Composition.** The full-state lift satisfies `R J=I_2`. Conditional
   on the row-to-region incidence pattern, regional slot shuffles and
   mixers are independent. Multiplying the entrywise bounds therefore
   costs `(1+C_{Q,theta}/L)^b=1+o(1)` because `b=O(log L)`. This is a
   comparison for each incidence pattern before the nonnegative
   row-weight coefficient sums, so those sums preserve it.

No gap was found in these six steps under the stated setup law. The
comparison does not by itself prove the distance theorem; it is combined
with the selected outer, routing, and fixed-occupation inequalities.

## Lean status

`research/lean/SpinCodes/UniformMixer/Kernel.lean` now verifies the uniform
action law for finite-dimensional invertible linear maps, the exact
full-state empty-gap kernel, and the one-bit impulse quotient. The proof
also gives a 22-bit theorem parameterized by a finite enumeration of
`GL(22,2)`; a separate theorem proves that such an enumeration exists.
No group elements are materialized. The value `D=2^22-1` is checked.

`research/lean/SpinCodes/UniformMixer/Regional.lean` verifies kernel
composition, `RJ=I_2`, the exact nonempty-gap factorization, the literal
matrix product `RV_c(1)J=P_0`, and the weighted impulse sandwich. It
proves both directions of the full-state gap comparison conditional on
the scalar live-weight estimate. Entrywise scalar bounds compose through
nonnegative kernels.

The scalar estimate (including the live mean and Hoeffding bound), the
discrete-to-continuous slot coupling, and the support/positivity argument
remain to be formalized. The finite numerical gates are separate obligations.
The revised Scalable theorem and the finite Turbo certificates are not yet
Lean-verified.
