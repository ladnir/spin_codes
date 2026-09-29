# Maintaining the Lean formalization

Start from [PROOF_GUIDE.md](PROOF_GUIDE.md) and the public import
[SpinCodes.Native](SpinCodes/Native.lean). Keep the canonical declarations in
`Spin.Structured.ConcreteNativeFamily`; the public import supplies a reading
boundary without duplicating their definitions.

## Proof changes

Use the pinned `lean-toolchain` and `lake-manifest.json`. A toolchain upgrade
is a separate change because it affects elaboration and compiled artifacts.

For a proof change, first identify the theorem's actual conclusion and
remaining hypotheses. A successful axiom audit does not discharge those
hypotheses. Preserve the native setup law, shared outer seed, exact threshold,
and code interpretation when modifying an interface.

Keep hand-written lemmas small enough to check independently. Import the
lowest suitable dependency rather than the public native entry point inside
an implementation module. Document a module's purpose, probability law or
matrix convention where relevant, and the theorem consumed downstream.

Check edited modules and their affected dependents before accepting a change.
Inspect full final statements as well as recursive axiom dependencies. The
original `SpinCodes/Pin.lean` and concrete statement pins guard different
claims; the historical pin alone does not check the native result.

## Generated certificates

The thousands of numerical modules are proof artifacts, not ordinary source
files to reformat or rename in bulk. Generators, exact witnesses, global box
indices, and family aggregates must stay aligned. Regenerate through the
owning script and check both local arithmetic and the semantic association
with the intended map, matrix, majorant, or rectangle.

Splitting data and checks across modules bounds kernel memory and permits
independent replay. Preserve that structure. Keep numerical certificates
kernel-checked, without proof holes, new project axioms, or compiled-trust
shortcuts. Run independent compilers only within the documented memory
budget. Never run two benchmarks simultaneously.

## Verification records

A report binds a particular source/object snapshot and command result.
Comments, whitespace, package metadata, and source filenames can change
recorded hashes even when the proposition stays the same. Do not relabel an
old PASS or replace its hashes to fit edited files. Preserve the old evidence
and obtain a new compiler record for changed artifacts.

A dependency inventory is not an individual source-recompilation record.
The provenance categories in [FINAL_EVIDENCE_REVIEW.md](FINAL_EVIDENCE_REVIEW.md)
keep those claims separate. A clean replay produces its own evidence and may
invalidate the old assembly snapshot by replacing compiled objects.

Use [FINAL_REPRODUCTION.md](FINAL_REPRODUCTION.md) to choose between checking
an existing assembly and replaying every project source. Stop editing the
relevant sources while their checks run. Keep generated reports, compiler
logs, and script hashes with the result they describe.

## Documentation changes

`README.md`, `NATIVE_RESULT.md`, and `STATUS.md` describe the current result.
`PROOF_GUIDE.md` explains its dependency structure. Historical checkpoint
notes and work logs retain their original scope and dates; mark them as
historical rather than rewriting the development history.

When relating Lean to the paper, state the verified scope precisely. The
native theorem establishes distance and rate. Separate statements about
sharper intermediate constants, encoder complexity, and arbitrary requested
lengths need their own correspondence checks. Use the Controlled Writing for
Cryptography skill when drafting that mathematical exposition.
