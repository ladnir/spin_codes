# Cheaper block-diagonal feedback and the alternating certificate

The global-MDS inner applies an RS transform both to emit the state and
to collect input feedback. Only emission needs a global MDS distance
bound. Replacing feedback by smaller independent RS transforms saves
arithmetic, but requires retaining their input occupancies in the proof.
Alternating four/eight-check feedback now has a complete 256-bit directed
certificate: 40% GF65536-symbol distance with more than 51.04 bits of
whole-code setup-failure margin at K20 and rate 1/2. Its validated native
circuit now takes approximately 5.5 ms with streaming output; the earlier
5.7--5.8-ms cached-output implementation is preserved. The pure four-check, approximately
5.6-ms candidate remains uncertified.

## Construction

Keep the outer, uniform regional route, global RS[64,16] emission A, and
joint packet refresh from [the shared-diagonal construction](HD_SHARED_DIAGONAL.md).
All setup distributions and their independence requirements are unchanged.
Put q=65536 and choose h=4 or h=8. Partition each packet's 64 positions
into 16/h consecutive blocks, each of width 4h. Within every block, C
evaluates the transpose of the normalized RS[4h,h] generator. Thus C is
block diagonal, has sixteen outputs, and differs from A^T. Compute

    Y = X + A Q,
    Q_next = Q + D2 (H D1 C X + 1 (v^T Q)).

The incoming state starts at zero and the terminal state is unrestricted.
Identity feedthrough still makes the inner invertible. The state and
feedback refresh share D2 exactly as in the certified reference.
Packet-local sorting is not adopted: this proof needs the uniform input
positions supplied by the original regional permutations.

## Conditional local bounds

Fix a potential input set U. Let l_b be its occupancy in feedback block b,
and let j=sum_b l_b. In the outer comparison experiment, X is uniform
on F^U, including zero coordinate values. For 0<z<1, define

    k_h(l,z) = sum_(c=0)^l binom(l,c) z^(l-c)(1-z)^c
                                  q^(-c-min(h,l-c)),
    K_l = product_b k_h(l_b,z),
    rho_l = sum_b min(h,l_b).

Then K_l=E[1[CX=0] z^wt(X)]. Expand the mark as a positive sum over
output-zero constraints. In each block, c fixed input zeros leave a
rank-min(h,l-c) RS restriction. Different blocks use disjoint variables,
which proves the product formula. T_j and V_j are unchanged from the
[global-MDS moments](HD_SYSTEMATIC_MDS.md), because they do not use C.

For every fixed nonzero incoming state Q, put

    Gamma_l = min(K_l, q^(-rho_l) z^max(0,49-j), (q^16-1)V_j).

This bounds E[1[CX=0] z^wt(X+A Q)]. For the first bound, each
output-zero constraint becomes an affine system in X with the same rank
as the homogeneous quiet system. It has either no solutions or the same
number of solutions. The second bound uses Pr[CX=0]=q^(-rho_l) and
the RS emission's weight floor. The third follows by summing the marked
moment over all nonzero incoming states and retaining one summand.

With d=q^16-1, the ideal two-component positive operator is therefore

    P_l = [[K_l, T_j-K_l],
           [1[j>0] V_j/d, (1[j>0]+1[j=0](q-1)/q)V_j+Gamma_l/q]].

The shared-diagonal pointwise argument applies to every syndrome CX;
it does not require C=A^T or CA=0. It bounds the actual operator by
a*b*P_l, with a=q/(q-1) and b=d/(q-1)^16. The same total density
loss applies after composition. This local comparison alone does not
establish a whole-code margin.

## Occupancy averaging and validation

Uniform regional permutations give a uniform potential-position set
conditional on its cardinality. A labelled occupancy vector has
multiplicity product_b binom(4h,l_b). Grouping permutations of its entries
also multiplies by the corresponding multinomial coefficient. Consecutive
packet operators are multiplied before applying the fractional power;
the regional coefficient then enforces exactly r positions per region.

Merging all dense profiles of one cardinality was much too loose. Keeping
their feedback ranks helps, but still loses information. The tightest
screen retains every occupancy profile, merging only equal rounded
matrix majorants. No profile or tail is omitted.

The exact module and tests check profile multiplicities, concrete mask
ranks, every input subset in a GF4 example, and its full weighted
shared-refresh operator. The native four/eight-check maps match separate
Python digests and dense scalar oracles, including arbitrary input states,
routed offsets, natural lengths, and aliasing. All 27 native variants pass
Release and ASAN/UBSAN checks.

At K20 with precomputed setup and 128-bit records, a sequential 301-call
run gives 6.175/5.810/5.610/5.929 ms for reference/eight-check/four-check/
reference. A reversed 501-call run gives 5.629/6.020/6.005/5.597 ms for
four-check/reference/eight-check/four-check. The four-check reduction is
repeatable; the eight-check advantage is not yet stable across runs.
Use about 5.6 ms for the four-check candidate, not a certified result.

The corresponding modes are `bucket4` and `bucket8` in the isolated
`spin_char2_hd_sumrank_forward` research target. Their checksums are
17504743431244161852 and 7820491484200485028, respectively. The validated
Release executable SHA256 is
`2eec238d0e03cec273ba870c06335088ab7c45cc89ca06414e5a56a7113e8ff3`.
The ~6.1-ms, 49.338-bit reference remains the fully certified checkpoint.

Next, check the medium-density classes with directed arithmetic and
complete coverage before treating either faster candidate as a winner.

## Alternating checks and a common circuit

The four-check screen still fails at r=2560: its optimized paired-regional
log2 upper bound is about +664. This is a limitation of the current proof,
not a counterexample. A directed four-check witness does certify
r=2048,...,2061 before the shared-refresh loss, with union logarithm
below -72.3896. It does not fill the medium-density gap.

The next candidate alternates four-check and eight-check packets, starting
with four checks at packet zero. K20 regions contain an even number of
packets, so each pair has one of each type. Compose P_l^(4) P_l'^(8)
before taking the fractional power. Profile multiplicities and coefficient
extraction remain as above. The directed driver explicitly records the
different left and right feedback maps; it does not reuse a same-map
certificate. Every active-row count and the final union must still pass.

The first native dispatch duplicated two entire packet kernels around a
branch and measured 7.29--7.32 ms. An explicit packet-pair loop measured
6.76--6.78 ms. Separately compiled packet functions reduced that to
6.29--6.30 ms, still slower than the reference. These are negative controls.

A common-circuit variant instead evaluates each local feedback polynomial
at its original packet coordinate, rather than translating every block
back to zero. Each block still evaluates all polynomials of degree below h
at distinct nodes, so every rank used above is unchanged. This changes
the encoded word, but preserves the local bounds and the witness verifier.
It lets both packet types share emission, the first two feedback FFT
levels, and refresh. Only the short final feedback calculation branches.
The native mode is `bucket48-fused`; its performance is checked separately.

The old ~6.1-ms certificate has also been replayed at 384-bit precision
after generalizing the verifier. It still covers all 8192 row counts and
gives 49.33892764811194 bits. The 405-test high-distance regression run
passes; two additional alternating-composition tests check direct ordered
matrix products against the FLOAT and directed coefficient builders.

The common-circuit variant passes all 29 Release and ASAN/UBSAN native
variants. Its two 501-call medians are 5.745 and 5.752 ms around a
5.953-ms frozen control. Checksum: 16994121132285933612. Its validated
executable is frozen as `bucket48-anchor-ca5148b1`, SHA256
`ca5148b15a1c67063e7fcca408fb2cda8bf29a7b939254b1020ad91bd41af72d`.
The complete directed coverage described below now accompanies this
performance/correctness checkpoint. A separate 384-bit replay also passes.

Two same-code follow-ups did not improve it. The unstriped affine circuit
(`bucket48-affine`) takes 5.951/5.974 ms around a 5.691-ms frozen control.
Padding the outer coefficient planes by one cache line (`bucket48-padded`)
takes 5.827 ms between frozen controls of 5.721/5.724 ms. In one executable,
the corresponding A/B/A comparison is 5.727/5.809/5.710 ms. All checksums
match. Neither variant is selected.

The frozen common-circuit phase measurements are 5.721 ms full, 2.456 ms
outer, 3.174 ms routed inner, and 2.921 ms pre-routed inner. Separate phase
medians need not sum to the full median. These controls suggest that
arithmetic, rather than random access alone, is the main remaining cost.

## Complete alternating certificate

The 256-bit directed run checks all 8192 active-row counts, including the
unchanged inverse-output dense tail. Its 104 rational witnesses are pinned
in [the replay module](theory/hd_bucket_feedback_witnesses.py). No numerical
optimization is needed for replay. Each witness rebuilds the left four-check
and right eight-check operators, including all occupancy profiles or their
explicit upward-rounded majorants. Their order is not interchangeable.

The ideal union logarithm is at most -63.30561793743408. The joint-refresh
density loss is at most 12.26300140707000 bits. Thus the final logarithmic
bound is below -51.04261653036408. Over all setup coins, the probability
that any nonzero message has output weight at most 838860 is therefore
below 2^-51.04. Here N=2097152, K=1048576, and weight counts GF65536
symbols. The eight field lanes in one 128-bit record use the same sampled
scalar code; no additional eight-fold union is required.

The completed 384-bit replay independently validates all 8192 row counts
and gives log2(failure) < -51.04261653036409. It is a separate validation
run, not a new witness search.
Run `python hd_bucket_feedback_witnesses.py --precision 384` from the theory
directory. The previously certified global-feedback checkpoint remains
preserved, and no production profile or paper number has been changed.
