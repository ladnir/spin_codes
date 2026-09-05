# Workstream brief: paper architecture

## Objective

Design the SPIN paper's dependency-ordered narrative and theorem spine. Do not
rewrite the manuscript during this first wave.

## Inputs

- `SPIN_NAMING.md`;
- `PAPER_RESTRUCTURE_PLAN.md` and `PAPER_INVENTORY.md`;
- the compiled TeX spine;
- the frozen Structured SPIN record;
- the Controlled Writing for Cryptography skill.

## Owned scope

Write only within `workstreams/paper_architecture/`.

## Tasks

1. Define the reader-level question and the common SPIN interface.
2. Produce a section outline for:
   - framework;
   - Accumulator SPIN;
   - Random SPIN;
   - Structured SPIN;
   - finite-length certificates;
   - asymptotic scaling;
   - implementation and performance.
3. Inventory every desired theorem, lemma, definition, and certificate.
4. Map reusable passages from existing TeX into the new outline.
5. Record terminology and notation, including the single composed interleaver.
6. Identify claims blocked by the modeled outer spectrum or non-outward
   arithmetic.
7. Recommend which proof details belong in the main text or appendix.

## Deliverables

- `OUTLINE.md`;
- `THEOREM_INVENTORY.md`;
- `SOURCE_MIGRATION_MAP.md`;
- `TERMINOLOGY_LEDGER.md`;
- `MERGE_SUMMARY.md`.

Do not claim that Structured SPIN has a final unconditional distance theorem.
