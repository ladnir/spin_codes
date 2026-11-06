# Turbo inner with the scalable bit-level route

## Exact candidate

Keep the paper's growing, shared Golay--BAA outer and its independent
bit-level group and region shuffles. Thus the regional SPIN architecture
has packet width `p=1`, `L=128m`, and `b` growing with the output length.
Replace only its IMT inner by the certified Turbo K18/K20 inner:

- step width `t=64` and state dimension `s=22`;
- the fixed expansion `A` whose image is `RM(2,6)` in the retained basis;
- feedback `C=A^T`;
- independent uniform `M_i` in `GL(22,2)` at every step;
- zero initial state, state carried across regions, and no final flush.

The recurrence remains `Y_i=X_i+A Q_i`,
`Q_{i+1}=M_i Q_i+A^T X_i`. Since `64` divides `L`, each bit-level region
contains whole steps. The outer, route, and inner setup are independent.
This is an exact inner replacement within the paper's regional architecture,
not a claim that the existing Scalable SPIN theorem or its Lean proof applies.
The native rate remains one half. With fixed `t,s`, the outer, route, and
inner still have `O(N)` ordinary and transposed work.

## Local facts that favor this candidate

The retained s22 expansion spans `RM(2,6)`, so every nonzero state has
`wt(Aq)>=16`. The kernel of `A^T` has minimum weight 8 (the dual
`RM(3,6)`). In particular, a step input with weight 1 through 7 activates
a zero state, and a live state with such an input emits at least
`16-7=9` of the 64 output bits. These are local statements, not a
whole-code distance bound.

For fixed `X` and nonzero `Q`, uniform `M_i` makes `M_i Q` uniform over
the nonzero states. Consequently `Q_{i+1}` is uniform over all states
except `A^T X`. If `A^T X` is nonzero, the probability of returning to
zero is exactly `1/(2^22-1)`; otherwise it is zero. This removes IMT's
lazy branch and may simplify the state transfer, but output cancellation
`X+A Q=0` still requires a proof under the routed input law.

The proved four-bit packet-window obstruction at the paper's `39/4`
growth constant does not directly apply to this one-bit route. In its
20% window witness, the required 52% regional concentration has
large-deviation exponent `D(0.52||0.2)=0.25166965...` per outer bit for
`p=1`, versus one quarter of that exponent for `p=4`. The growth
budget is `ln(2)/(39/4)=0.071092...` per outer bit. This defeats that
particular concentration witness for `p=1`; it does not prove the
candidate's distance.

## Proof gates

1. Recompute the one-active-row first moment with the shared BAA
   spectrum and the exact s22 live-state output law. Check the original
   `b >= (39/4) log_2 N` schedule at distance `0.11N`, allowing a larger
   logarithmic constant if needed.
2. Replace the IMT fixed-occupation transfer. Uniform GL refresh may
   reduce its state description, but the proof must retain the routed
   input support and possible `X+A Q` cancellations.
3. Reprove growing sparse and positive occupation bounds, including the
   shared outer selection event. The old seven-coordinate IMT witness
   is not a certificate for this state law.
4. Only after a complete asymptotic proof, test finite BCH or RS outer
   certificates and measure transposed/ordinary encoding. The Turbo
   K16 profile uses paired s15, so an exact all-size s22 claim would also
   need a new K16 certificate and timing.

At a fixed output length, the t64 inner takes twice as many steps as
the paper's t128 IMT inner, and a uniform GL22 update can cost more than
a transvection. The bit-level route also moves four times as many routing
units as Turbo's four-bit route. The asymptotic work remains linear, but
the finite performance constant must be measured. The retained K18
Turbo stage screen gives 256.569 us for *inner plus packet routing* out
of 461.422 us total; it does not isolate the inner or predict bit-route
performance.

## Source bindings

- `research/paper/structured_proof.tex`: current p1/IMT asymptotic proof.
- `research/paper/turbo_appendix.tex`: exact Turbo s22 map and finite proof.
- `research/workstreams/k16_design/README.md`: RM(2,6) and feedback-kernel
  checks, including the s22 map's minimum weights.
- `research/workstreams/turbo_spin_migration/PACKET_WINDOW_OBSTRUCTION.md`:
  four-bit route obstruction at growth constant 39/4.
- `spin/experiments/rs8_sizes/README.md`: retained stage measurements.
