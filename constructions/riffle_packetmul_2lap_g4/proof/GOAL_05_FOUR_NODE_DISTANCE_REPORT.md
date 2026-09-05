# Goal 05: four-node distance report

## Result

Goal 05 does not determine the exact minimum distance. The certified statement
remains

\[
36\le d(C_4)\le42.
\]

The lower bound is the Goal 04 four-node certificate. The upper bound is an
independently replayed word with block weights \((8,14,6,14)\). No solver found
a word of weight at most 41, but neither solver proved that no such word exists.

## Object and decision predicate

Let \(U:\mathbb F_2^{64}\to\mathbb F_2^{64}\) be the autonomous output map.
The observation code is

\[
C_4=\{(o,U(o),U^2(o),U^3(o)):o\in\mathbb F_2^{64}\}.
\]

It is a systematic binary \([256,64]\) linear code. The unresolved exact
predicate is

\[
\exists o\ne0:\quad
\operatorname{wt}(o,U(o),U^2(o),U^3(o))\le41.
\]

Unsatisfiability would prove \(d(C_4)=42\). A satisfying assignment would give
an exact smaller word.

## Authenticated upper bound

The primary and independent constructions replay the same four states:

| block | state | weight |
|---:|---:|---:|
| 0 | `0x8951a0000` | 8 |
| 1 | `0x1622240e0462` | 14 |
| 2 | `0x142290000` | 6 |
| 3 | `0x25050ab00a5` | 14 |

The audit reconstructs the systematic state map from the generator constant,
applies the accumulator recurrence three times, and obtains the same states.
Their total weight is 42.

## Exact solver attempts

The primary encoding uses 64 information variables, 256 output variables,
native XOR equations, and sequential cardinality counters. CryptoMiniSat
returned `UNKNOWN` after its 60-second budget.

The independent encoding reconstructs the linear map separately and expands
each parity equation into a pure-CNF Tseitin chain. Its instance has 19,591
variables and 50,923 clauses. CaDiCaL returned `UNKNOWN` after exactly 100,000
conflicts, 180,461 decisions, and 103,716,959 propagations. The run took 31.92
seconds. `UNKNOWN` is only a bounded diagnostic; it is not evidence that the
predicate is false.

Both instances use the following valid strengthening. Goal 04 excludes any
total-weight-at-most-41 word containing a block of weight at most 8. Therefore
all four block weights are at least 9. Their total is at most 41, so each block
weight is at most 14.

## Why direct anchor enumeration stalls

At least one block of a putative word of total weight at most 41 has weight 9
or 10. There are

\[
\binom{64}{9}+\binom{64}{10}=179{,}013{,}799{,}328
\]

such anchors in one information set. Checking four possible anchor alignments
gives 716,055,197,312 naive visits. This is about 140 times the 5,130,659,560
low-weight anchors exhausted in Goal 04. Directly extending that enumeration is
therefore not a suitable medium-sized proof step.

## Construction-level consequence

Goal 05 does not improve the four-node lower bound, so it does not move the
sufficient zero-prefix cutoff. The cutoff remains node 20,976. The support-33
closed partial bound remains in the log-base-two interval

\[
[-40.972830672896,-40.972830672895],
\]

and the remaining numerical budget remains

\[
[-41.027690825931,-41.027690825930].
\]

These statements inherit the Goal 04 pointwise argument. They do not assume
that packet multipliers are independent after conditioning on a nonzero
first-lap terminal state.

## Recommended next proof step

Return to the first open construction-level strata rather than spending the
next bounded session on generic minimum-distance SAT. At first occupied nodes
20,972 through 20,975, the completed autonomous blocks already contribute
188,748 output bits. A bad word can therefore have tail weight at most 17.
The next goal should isolate this finite boundary band and prove or refute that
small-tail event using the packet-multiplier structure.

If exact knowledge of \(d(C_4)\) later becomes necessary, use a dedicated
minimum-distance method with information-set lower bounds or meet-in-the-middle
syndrome tables. Repeating generic SAT or enumerating all weight-9 and
weight-10 anchors is not the recommended next mechanism.

## Artifacts

- `receipts/goal05_distance_primary.json`: native-XOR bounded result.
- `receipts/goal05_distance_independent.json`: independent pure-CNF bounded
  result.
- `receipts/goal05_distance_audit.json`: source authentication, recurrence
  replay, search-scale calculation, and retained-ledger audit.

The status is
`GOAL_05_BOUNDED_SOLVER_OBSTRUCTION_EXACT_DISTANCE_OPEN`. This report does not
claim the full construction.
