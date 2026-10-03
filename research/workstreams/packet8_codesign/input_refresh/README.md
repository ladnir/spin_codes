# Input-refresh order: a negative bounded gate

Neither tested update change closes the current uniform-input comparison.
These are distinct research constructions, not implementation replacements.
The outer, route, and byte-native maps A and C stay unchanged. The emission is
`y=X+A*s` in both variants; the initial state is zero and the final state is discarded.

In the independent-refresh variant, update to `a*s+b*CX`, where a and b are
independent uniform nonzero GF(2^16) scalars. In the shared-refresh variant,
update to `M*(s+CX)`, where M is independently uniform in GL(2,GF(256)).
Each nonzero updated state is exactly uniform, so the weighted state measure
has only two components: a point mass at zero and the uniform nonzero law.

Let L=65,535, let W_j(0)=E[z^wt(X) 1{CX=0}], and let T_j be the emitted-weight
moment averaged over a uniform nonzero entering state. The zero-source row is
`[W_j(0), E[z^wt(X)]-W_j(0)]` for both variants.
An independent-refresh nonzero-source row is bounded by `[T_j/L,T_j]`, with
zero return mass at j=0. For shared refresh, the exact return mass is

```
J_j = E[z^wt((I+AC)X) 1{CX!=0}]/L.
```

Its exact row is `[J_j,T_j-J_j]`. The evaluator uses
`[J_upper,T_j-J_lower]`; it never subtracts an upper estimate of J_j.
Positive integer censuses provide exact finite formulas for j=1 and j=2,
covering 2,040 and 1,820,700 inputs, respectively. For larger j, it bounds the
full feedforward moment by the minimum of a four-block Hölder bound, an
information-set bound, and one. The exact weighted CX=0 contribution is
subtracted from those moment upper bounds before dividing by L.

All numerical endpoints are floating proposals, not certificates.
The four exhaustive GF(4) tests pass: they check the return census, both
operator dominations, and the blockwise Hölder majorant.

At q=119 and tilt 0.4, the independent-refresh bound gives -4,148.305000 bits,
only 2.63 bits better than the old update at that same point. The shared-refresh
bound gives -4,433.600119 bits. Deleting all shared-refresh returns at j>=2
still gives -3,978.961196 bits; retaining exact j=1/j=2 returns and deleting
the rest gives -3,989.902809 bits. These deletions are deliberately optimistic
diagnostics, not valid bounds. They show that a large higher-j census is not
justified to repair this particular point.

The exact shared-refresh J_1 is 20.24 times T_1/L at this tilt.
For example, a unit byte at position zero has CX=(x,0), and (I+AC)X has
seven copies of that unit byte, of total binary weight seven. Thus moving the
same random map after feedback introduces a cheap one-packet return pattern.
This is a local event, not a whole-code distance counterexample.

Run the bounded gate and tests from the repository root:

```text
python -B -m unittest discover -s research/workstreams/packet8_codesign/input_refresh -p test_refresh.py -v
python -B research/workstreams/packet8_codesign/input_refresh/refresh_gate.py --output <fresh-result.json>
python -B research/workstreams/packet8_codesign/input_refresh/diagnostic.py --q 119 --tilt .4
```

`refresh_v1.json` records the maps, geometry, full tilt grid, and source hashes,
all verified at completion. The next analysis instead separates shared routing
failures from message counting; no further census or implementation work is
recommended for these update variants on the present evidence.
