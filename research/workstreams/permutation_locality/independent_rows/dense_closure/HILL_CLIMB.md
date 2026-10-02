# Hill-climbing after the first four-bit closure

The [0.5% certificate](FIRST_CLOSURE.md) is complete and replayed at
256-bit precision. The first 1% attempt is not complete. Neither its
floating proposals nor its passing subset establish a 1% distance bound.

The first retarget uses cutoff 20971, without changing the code. Rechecking
the 812 old leaves at this cutoff retained 537 and requeued 275. The
bounded exploratory run then processed 100 additional cells. It reported
576 accepted leaves, 277 pending regions, and no depth-limit failures at
that point. The run was interrupted; its intermediate partition was not
saved. The complete 0.5% witness file remains available and unchanged.

Some old witnesses cannot simply retain their parameters. For example,
retuning only the output tilt on the old cell

    x1 in [1/2048, 3/4096],
    x2, x3, x4 in [0, 1/4096]

still gave a floating log2 upper around +2034 at 1%. This is a failure of
that whole-cell bound, not a bad-codeword example, a lower bound on the
failure probability, or an upper bound on the construction's distance.
It includes slack from the cell width and from separate weight majorants.

## Coupled comparison refinement

`atlas.py --coupled-alternates` now couples alternate input tilts to the
outer count through affine moment bounds. The derivation is in
[README.md](README.md#coupling-alternate-comparison-weights). The code
and setup distribution are unchanged. Both the earlier witnesses and
the new coupled witnesses can be replayed by the same verifier.

On the displayed boundary cell at cutoff 20971, the base proposal gives
log2 upper approximately +1544.1123. The coupled proposal gives +12.1644,
confirmed by 192-bit outward arithmetic. The earlier +2034 number used
only output-tilt retuning, not the full base proposal. These comparisons
show substantial removable proof slack. None of these positive bounds
certifies the cell or demonstrates a low-weight codeword.

The CDF support-cover driver also now accepts an explicit output cutoff.
Previously it always used 209715, even when planning a lower-distance
dense proof. The cutoff changes only the final Chernoff factor; it does
not change the encoder, inner operators, or outer counting bounds.
Shell-based driver modes still require their original cutoff and reject
an incompatible option rather than silently ignoring it.

The resulting [sparse extension](SPARSE_HANDOFF.md) is now verified
through q=64 at cutoff 104857 (5%). The dense cover can start at q=65
at any smaller cutoff, including the current 1% target. This does not
establish a whole-code certificate at either new distance.

The first coupled run saved a complete partition at cutoff 20971 and
q_min=59 in `tmp/four-bit-dense-one-coupled-batch1.json`. After rechecking
the original atlas, it processed 57 new cells and saved 569 accepted
leaves and 268 remaining regions. There were no depth-limit failures.
Of the accepted witnesses, 27 use the new coupled alternate bounds.
The accepted-region aggregate has more than 88.244 bits of margin,
but this excludes all remaining regions and is not a whole-code bound.

The run originally requested 100 cells. It was shortened to a 57-cell
batch at a debugger checkpoint by changing only the loop's cell budget.
The cover then finished its current cell and saved normally; no bound,
arithmetic state, witness, or remaining region was changed or discarded.
The saved file has 273663 bytes and SHA-256

    e0c0e9c9beb430bbdde154bbda45d9f8e22a8e5081bcc94f0e2cd27e54d6ecce

Repeated affine-weight fits are now cached during proposal generation.
Replay still reconstructs every supplied dual's exact majorant. This
search optimization does not alter the inequalities or saved witnesses.

## Next work

The [selected-composition diagnosis](BOUNDARY_DIAGNOSIS.md) now identifies
an avoidable shuffle-comparison charge. All nine tested compositions pass
at 1% and 2%; four fail at 5% under the previous bound, but all four pass
under the new exact special-composition comparison, replayed at 256-bit
precision. This is not a full certificate at any new distance.

The [posterior extension](POSTERIOR_COMPARISON.md) now accepts
heterogeneous mixtures and is integrated into complete enumeration of
small cells. Prioritize extending its coverage to the broader pending
cells, rather than further selected-point tests. Keep the existing fixed
construction. The original cover refinements remain useful alongside
this sharper inequality:

- Apply the coupled alternate witnesses to the remaining cells, splitting
  cells where their range remains too wide.
- Use the verified q=65 handoff. If it remains expensive, extend the sparse proof at the
  lower output cutoff so that the dense cover can begin at a larger
  occupancy. There is no requirement to keep the split at q=59.

Do not change the inner or claim that it needs changing merely because a
comparison bound fails. Preserve the 0.5% certificate and the stronger
two-bit certificate while testing either refinement.

Use smaller saved batches for the next search, rather than interrupting
a long batch before it writes its partition:

```sh
python -B research/workstreams/permutation_locality/independent_rows/dense_closure/atlas.py --retarget tmp/four-bit-dense-one-coupled-batch1.json --threshold 20971 --minimum-groups 65 --coupled-alternates --posterior-shuffle --max-cells 10 --max-depth 64 --max-unresolved 2 --output tmp/four-bit-dense-one-handoff65-batch1.json
python -B research/workstreams/permutation_locality/independent_rows/dense_closure/atlas.py --resume tmp/four-bit-dense-one-handoff65-batch1.json --coupled-alternates --posterior-shuffle --max-cells 10 --max-depth 64 --max-unresolved 2 --output tmp/four-bit-dense-one-handoff65-batch2.json
```

Search scores and old stored uppers are never accepted in place of outward
replay. The driver restores the requested precision after authenticating
the inputs; its output prints that precision explicitly. Final acceptance
still requires no unresolved cell and a fresh successful assembly.
