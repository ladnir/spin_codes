# A higher-rate K20 certificate

The first complete higher-rate checkpoint has rate `11/16` and relative
symbol distance at least `10%`. At requested source length `K=2^20`, its
ideal setup fails with probability below `2^-119`. The complete replay
joins three bounds for the same construction; it does not extrapolate
from sampled weights. The [first native implementation](FIELD256_NATIVE.md)
is correct but slower than the binary control.

This addresses the corrected objective: half-rate/10% already follows by
scalar extension of binary SPIN. The new result improves rate, not merely
the setup-failure margin. It does not yet establish a performance win or
the separate high-distance target.

## Construction and guarantee

Work over AES `F=GF(256)`. The fixed outer is a systematic GRS[256,176]
code, repeated in `L=6016` rows. Its padded source length is `1,058,816`
symbols and its output length is `N=1,540,096`. Restricting the source to
the requested `1,048,576` symbols gives effective rate `32/47`, still above
`2/3`. The full padded code is the object certified below, so zero padding
does not require another proof.

Collect the same coordinate from each outer row into one region. In each
of the 256 regions, independently sample a uniform permutation of the
`L` rows. This is one regional routing stage, not a uniformly sampled
permutation of all `N` symbols. Multiply each routed coordinate by an
independent uniform element of `F*`.

The inner processes sixteen symbols per packet and carries eight state
symbols across packets and region boundaries. Use the four fixed symmetric
Cauchy phases and `eta=2` from [the map construction](STATE_CAPACITY_REPAIR.md).
For each phase, `A` is `16 x 8`, `C=A^T`, and `CA=0`. For packet input `X`
and incoming state `S`, sample an independent diagonal `D_s` with eight
uniform nonzero entries and apply

\[
R=D_sS,\qquad Y=X+2AR,\qquad S'=R+CX.
\]

Sample all permutations and diagonals once during setup, independently of
one another. Fix them for all messages. Initialize the state to zero and
discard its final value. The diagonal identity block on the current input
makes the complete inner invertible.

Let `E` be the resulting random injective encoder. The proved event is

\[
\Pr_{\mathrm{setup}}\left[
  \exists m\ne0:\ \operatorname{wt}(E(m))\le154009
\right]<2^{-119}.
\]

Outside this event, the minimum distance is at least
`154010=ceil(N/10)`. Weight counts nonzero field symbols. The 119 bits
measure setup failure, not the computational security of an application.

## Complete coverage and numerical replay

Let `a` be the number of nonzero outer rows and let `W` be the weight of
the input to the inner. Each nonzero GRS row has at least 81 nonzero
coordinates. Hence `W>=81a`.

| Bound | Covered class | Certified integer margin |
|---|---|---:|
| Sparse | Every `1<=a<=600`, all `W` | 119 bits |
| Intermediate | Every `46203<=W<=369623` | 344 bits |
| Dense | Every `369624<=W<=N` | 12072 bits |

If `a>=601`, then `W>=48681`, above the intermediate lower endpoint.
Every nonzero word is therefore covered. Overlap is allowed, and the
three upper bounds must be added.

The [sparse replay](theory/state8_sparse_certificate.py) uses 80 disjoint
active-row intervals. It evaluates the positive regional coefficient
formula from [the sparse screen](LARGER_OUTER_STATE8_SCREEN.md) using
rational marks and directed Arb arithmetic. Polynomial truncation removes
only degrees above the largest required regional occupancy. It checks
the source-factor maximum on each interval and sums interval
multiplicities. Its exact dyadic upper endpoint is approximately
`1.156409247462385e-36`.

The [intermediate replay and proof](STATE8_DENSITY_CERTIFICATE.md) cover
all 323421 integer weights in 216 contiguous cells. The improvement is a
joint rank bound for output-zero and state-return equations. It avoids
charging the same random diagonal twice. Each cell fixes rational marks;
all subsequent arithmetic is outward-rounded. Its upper bound is
approximately `1.608974643900775e-104`.

The [dense proof and exact replay](FIELD256_COVERAGE_AUDIT.md) condition on
the independently sampled invertible inner. They use the regional
outer-membership bound and a union over the GF256 low-weight output ball.
All 152 density intervals are checked with rational inequalities, giving
the conservative upper bound `2^-12072`.

[The aggregate replay](theory/field256_whole_certificate.py) executes all
three components, verifies their common geometry and integer coverage,
and adds their upper endpoints as exact fractions. It returns

\[
B_{\rm all}\le
1.1564092474623849969182311231281288098474677334282973911734868933659920640858008
\cdot10^{-36}<2^{-119}.
\]

This is a summed whole-code bound, not the minimum of three reported
margins. The frozen intermediate partition hash is
`01c230076663d1bff17808ee26577d6cc1304f2761a0e599a2e124e072abecff`.

Run from the repository root:

```sh
python -B research/workstreams/char2_field_spin/theory/field256_whole_certificate.py
python -B research/workstreams/char2_field_spin/theory/field256_whole_certificate_tests.py
python -B research/workstreams/char2_field_spin/theory/state8_sparse_certificate_tests.py
python -B research/workstreams/char2_field_spin/theory/state8_density_certificate_tests.py
```

The sparse and intermediate replays have independent audits and matching
replays. Their five and eight tests pass; the aggregate has five passing
regression tests. Exact scalar-map, rank-cover, and finite-field checks
are linked in the component notes. These are auditable proof scripts and
mathematical arguments, not proof-assistant verification.

## Extension fields and implementation status

All realized generator entries lie in GF256. For each fixed setup,
scalar extension preserves the minimum symbol distance exactly. Expand
any extension-field word in a GF256 basis: each component is a base-field
codeword whose support is contained in the original support. Conversely,
base-field codewords remain available after extension. Thus the same
setup event and failure bound apply simultaneously to compatible
`GF(2^m)` representations for `m=8,32,64,128`, without a lane union.

Only GF256 scalar actions are needed, not general large-field products.
The [circuit-cost study](LARGE_OUTER_COST.md) nevertheless identifies a
real implementation risk: the GRS[256,176] outer has substantial Boolean
work and scratch traffic. Its product count alone does not predict time.

At K20 the padded output has 26.5625% fewer symbols than a half-rate
encoder on the requested input. The first matched native comparison does
not reduce total encoding time: arithmetic and scratch costs dominate.
Neither K16/K18 certificates nor width-thirty-two certificates follow
from this checkpoint. The width-thirty-two density screen still has a
gap in its current bound; it remains an unproved alternative.

The exact-map optimization pass preserves this certificate but does not
yet match binary SPIN's speed. Next, use this proof as a reference for a
cheaper compatible construction. Do not spend the next pass increasing the already sufficient
119-bit margin. Preserve the old proofs and production implementations.
