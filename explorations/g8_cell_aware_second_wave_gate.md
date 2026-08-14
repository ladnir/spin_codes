# Audit of the first cell-aware wave and the next finite gate

The first cell-aware experiment justifies one additional pricing wave. It does
not justify an indefinite sequence of tuning rounds.

## Evidence from the first wave

The input checkpoint has 189 active leaves. The first wave selected eight
distinct `h2` parents and tuned one inner and one outer component per target.
All sixteen candidates passed the `64`-bit reduced-cost gate.

The lane results are asymmetric:

| lane | reduced-cost range | mean reduced cost |
|---|---:|---:|
| inner | `-260639` to `-135039` | `-177204` |
| outer | `-767636` to `-628440` | `-734084` |

The rational replay lowered the expanded and worst contribution diagnostics
from `1085370.7304764104` to `685003.483464075`. Thus the realized improvement
is `400367.24701233534` bits. Every active leaf gives positive rational weight
to a new outer component. A new inner component has positive weight on 142 of
189 leaves. All leaves improve, and no old selector is retained.

These facts show broad component reuse. They do not isolate the marginal
effect of each lane, because the artifact contains no outer-only or inner-only
ablation. The reduced costs also measure pricing at the old dual, not the
realized improvement after both banks change. Nevertheless, the outer signal
is about four times larger and has complete leaf usage. The next bounded wave
should therefore spend its jobs on the outer lane.

## Second-wave targets

Recompute every component LP after loading all first-wave rows. Rank leaves by

```text
E_l=log2(count_l)+U_l.
```

Use canonical leaf identifiers to break ties. Retain at most one leaf per
original `h2` parent. The current first eight are:

```text
h2:146, h2:102, h2:072, h2:063,
h2:071, h2:101, h2:108, h2:022.
```

Their contributions range from about `685003` to `652324` bits. None was a
first-wave target. This change in the active frontier is evidence against
retuning the old eight barycenters.

Run one outer pricing job at each new dual barycenter. Freeze at most one
candidate per target. Reprice every candidate on all 189 leaves before solving
the updated rational marginal LPs. Do not tune inner components in this wave.

An outer candidate is accepted only if:

1. its affine barycenter price agrees with its complete vertex column within
   `1e-6` bit;
2. its measured reduced cost is below `-64` bits;
3. all frozen parameters are finite and support-eligible;
4. its canonical parameter digest is new; and
5. a replayed rational selector gives it positive weight somewhere.

Conditions 1--4 admit the row to the bank. Condition 5 distinguishes a useful
column from one made redundant by other additions.

## Finite continuation and stop policy

Let `A_t` be the expanded rational diagnostic after wave `t`, and define

```text
G_t=max(0,A_t+40),                Delta_t=A_(t-1)-A_t.
```

Wave two is successful only if all of these conditions hold:

1. at least four of eight targets improve by `16384` bits;
2. at least four new outer rows receive positive rational weight;
3. at least one new row is used outside its source `h2` parent;
4. both the expanded diagnostic and worst contribution decrease by `65536`
   bits; and
5. no leaf regresses after retaining the old selector as a fallback.

If any condition fails, stop continuous component pricing. Retain useful
rows, then compare exact geometry refinement with a broader outer family.

Authorize at most one third wave. It is permitted only when wave two succeeds
and

```text
ceil(G_2/Delta_2) <= 2.
```

The third wave uses the same eight-target, one-lane, one-column limits after a
fresh global replay. Stop component pricing after that wave regardless of its
outcome. This cap prevents a local nonlinear oracle from becoming an implicit
completeness assumption.

Diagnostic closure occurs when the expanded rational replay is at most
`-256` bits and every active leaf has a fixed rational selector. The buffer is
only a handoff gate. The theorem closes only after independent outward replay
proves the expanded union upper bound at most `-40` bits.

## Proof artifact path

The second wave should produce four immutable layers.

1. A diagnostic pricing report binds the input checkpoint, old component
   sources, target vertices, dual barycenters, complete columns, and reduced
   costs.
2. A content-addressed outer source freezes only accepted parameter rows. It
   records the conditioned-row variables, BL coefficients, Cauchy split, and
   all source digests required to reconstruct the outer affine intervals.
3. A replay artifact stores two sparse rational marginals per leaf. Every
   reference resolves to a separately hashed inner or outer row. The combined
   evaluator subtracts the packet-profile normalization exactly once.
4. A new manifest version binds the component sources, recombination contract,
   normalization rule, proof code closure, geometry, and transitive data.

The present manifest binds combined atlas rows. It must not be reinterpreted
as a factorized manifest. The independent verifier must reconstruct component
intervals from parameters and ignore producer constants, charges, scores, and
dual data.

## Theorem work in parallel

Pricing should not block the proof lane. The following obligations are already
on the critical path:

1. State and prove that the inner transfer bound is uniform after conditioning
   on the complete outer word and all information retained by the outer lemma.
   This excludes an invalid product of averages over shared randomness.
2. Freeze the factorized selector schema and its component-level support rule.
3. Generalize independent outward replay to the `g=8` conditioned-row outer
   parameters used by the new rows.
4. Verify every exact leaf vertex, count, split, and union-accounting record.
5. Produce a receipt that reports both collapsed and expanded aggregation. The
   expanded upper bound is the theorem gate.

The recommended next action is therefore one eight-job outer-only wave while
the uniform conditional inner lemma and factorized manifest/verifier work
continue independently.
