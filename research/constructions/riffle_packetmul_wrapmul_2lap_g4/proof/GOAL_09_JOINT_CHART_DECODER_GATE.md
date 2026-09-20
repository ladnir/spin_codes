# Goal 09: joint-chart decoder gate

## Result

The joint-chart reduction is correct, but the two tested decoder mechanisms
do not yet scale to the 38-node target.

An exact width-16 implementation used the systematic chart

\[
\Phi(s)=(q_1(s),q_2(s)).
\]

This map has rank 32. The decoder exhausted every pair

\[
(u,v)\in\mathbb F_2^{16}\times\mathbb F_2^{16},
\qquad
\operatorname{wt}(u),\operatorname{wt}(v)\le5.
\]

It recovered minimum weight 121, minimum count one, and state
`0xcd7e70d6`. These values agree with the independent exhaustive census of all
$2^{32}-1$ nonzero width-16 states from Goal 07.

This agreement validates the chart change of basis and the joint objective.
It does not provide a scalable decoder. The exact implementation forms the
entire Cartesian product of the two sparse lists.

A second implementation exposed the chart coordinates directly to a
native-XOR SAT solver. It returned `UNKNOWN` after 30 seconds for both weight
bounds 121 and 120. The first instance contains the known minimum. The second
instance has no solution by the independent census. Thus the solver failed
on both sides of the exact threshold before the width-64 instance was tried.

## Width-16 systematic chart

Let $Q_{24}(s)$ denote the 24-node observation word for the extended BCH
$[32,16,8]$ recurrence. The selected chart consists of two complete output
packets:

\[
u=q_1(s),
\qquad
v=q_2(s).
\]

The linear map $s\mapsto(u,v)$ is bijective. Therefore each pair $(u,v)$
determines exactly one lifted state and exactly one observation word.

The known minimum has packet weights

\[
\operatorname{wt}(q_1(s))=4,
\qquad
\operatorname{wt}(q_2(s))=4.
\]

It lies inside the tested chart domain. Each radius-five list contains

\[
\sum_{a=0}^{5}\binom{16}{a}=6{,}885
\]

vectors. The exact join checked

\[
6{,}885^2=47{,}403{,}225
\]

pairs. It excluded the zero state and minimized the weight of the full
384-bit observation word.

The run took 0.048 seconds on the current machine. This timing is useful as a
small-instance implementation check. It is not evidence that the same join
is practical at width 64.

## Native-XOR solver experiment

The SAT encoding first changed from the original lifted-state basis to the
systematic chart basis. The first 16 information variables were the bits of
$q_1(s)$. The next 16 variables were the bits of $q_2(s)$. The encoding then
imposed three cardinality bounds:

\[
\operatorname{wt}(Q_{24}(s))\le B,
\qquad
\operatorname{wt}(q_1(s))\le5,
\qquad
\operatorname{wt}(q_2(s))\le5.
\]

Only the 352 nonsystematic output coordinates required native XOR
constraints. The encoding used 32 information variables, 352 output
variables, and 64,021 CNF clauses before the XOR constraints.

For $B=121$, the solver did not find the known satisfying assignment in 30
seconds. For $B=120$, it did not establish the known unsatisfiability in 30
seconds. Both results are timeouts. They establish no distance claim.

The small-instance gate therefore rejects this SAT encoding as the immediate
width-64 route. A C38 timeout would add little information.

## Width-64 scaling obstruction

One C38 chart has sparse-coordinate list sizes

\[
L_5=sum_{a=0}^{5}\binom{64}{a}=8{,}303{,}633
\]

and

\[
L_6=sum_{b=0}^{6}\binom{64}{b}=83{,}278{,}001.
\]

The direct join contains

\[
L_5L_6=691{,}509{,}957{,}277{,}633
\]

pairs per chart. This is 14.6 million times the width-16 pair count. The
1,330-chart cover from Goal 08 contains approximately
$9.20\cdot10^{17}$ pairs.

Even an optimistic extrapolation that preserves the measured width-16 pair
rate gives about 8.2 days per chart and 29.7 years for all charts. This
estimate ignores the longer codewords, memory traffic, and list-construction
costs. It is therefore a lower, not upper, estimate of the practical cost.

The required improvement is now precise. A useful decoder must avoid forming
almost all pairs in $L_5\times L_6$. Faster evaluation of individual pairs is
not sufficient.

## Assessment

Goal 09 is positive about the proof model and negative about the first two
algorithms.

The positive conclusion is that the sparse systematic coordinates are the
right interface. On the exact smaller recurrence, their joint enumeration
recovers the independently known global minimum without loss.

The negative conclusion is that neither a Cartesian join nor the tested
native-XOR cardinality encoding closes the width-64 problem. The present
obstruction is algorithmic rather than a new failure of the construction.
No C38 counterexample was found, and no C38 lower bound was proved.

The next goal should implement an exact threshold join with blockwise
Hamming lower bounds on the width-16 chart. A suitable gate is to recover the
same minimum while evaluating only a small fraction of the 47.4 million
pairs. Multi-index Hamming bucketing or a Stern-style collision split are the
most direct candidates. The implementation should report the number of
generated buckets, candidate pairs, duplicate pairs, and exact final
comparisons. Only an algorithm that materially compresses the width-16 join
should be transferred to a sampled width-64 chart.

## Evidence and scope

The exact decoder source is
`scripts/analyze_riffle_packetmul_wrapmul_2lap_g4_goal09_chart_exact_width16.cpp`.
The solver source is
`scripts/probe_riffle_packetmul_wrapmul_2lap_g4_goal09_chart_sat.py`.
The audit source is
`scripts/audit_riffle_packetmul_wrapmul_2lap_g4_goal09_joint_chart.py`.

The authenticated summary is
`constructions/riffle_packetmul_wrapmul_2lap_g4/receipts/goal09_joint_chart_audit.json`.
The raw exact-enumeration receipt and the two SAT receipts are in the same
directory.

The audit reconstructs the width-16 BCH instance, checks that
$s\mapsto(q_1(s),q_2(s))$ has rank 32, replays the minimum witness, and
compares the result with the independent Goal 07 census. The exact decoder
proves only the minimum within the stated width-16 chart. The agreement with
the full census comes from the separate exhaustive result. No artifact in
this goal decides any width-64 chart.
