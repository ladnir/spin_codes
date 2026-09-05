# One-stage sparse-EA versus optimized BCH

## Question

The degree-33 sparse-EA outer has a complete finite distance certificate.
Its raw XOR count appeared too high for the current implementation budget.
This experiment compares its optimized transpose with the existing optimized
BCH outer at \(k=2^{20}\).

## Matched outer maps

The sparse constituent is

\[
 E:\mathbb F_2^{256}\longrightarrow\mathbb F_2^{512},
 \qquad
 C:=AE,
\]

where each row of \(E\) has degree 33. The matrix was sampled once with seed
20260903. Its row hash is
42508e8b7d1d74f52ef3b4b9d884655d113fec81408bc095e0908da995ae2fac.
The benchmark computes

\[
 C^{\mathsf T}=E^{\mathsf T}A^{\mathsf T}.
\]

The implementation first performs the suffix XOR for \(A^{\mathsf T}\).
It then evaluates the verified common-subexpression circuit for
\(E^{\mathsf T}\).

The comparator applies two copies of ExtendedBch256x128-Eq3. Thus both
paths map 512 input coordinates to 256 output coordinates. The BCH generator
hash is
ff5f6eae3ab0e9b974e38f20334f25dc50e412018439bcbead0d53e9985dcb3d.

Both kernels pack two independent constituent calls into AVX2 vectors. The
complete workload maps \(2^{21}\) 128-bit input blocks to \(2^{20}\)
128-bit output blocks.

## Exact circuit comparison

| Path | Scalar transposed XORs per \(512\to256\) capacity | Peak live packed values | Compiled function | Stack frame |
|---|---:|---:|---:|---:|
| sparse-EA | 9,733 | 1,305 | 164,189 bytes | 53,312 bytes |
| two BCH constituents | 5,868 | 448 | 51,914 bytes | 19,328 bytes |

The sparse circuit uses 1.659 times as many XORs. Its register-pressure
surrogate is also substantially larger.

The generator verifies both linear circuits exactly. Before each benchmark
process starts, the executable compares the generated packed kernel with an
independent dense reference.

## Peach result

Peach has an AMD Ryzen 9 7950X. GCC 15.2 compiled the source with
-O3, -march=native, -std=c++20, and -pthread. Each process was pinned to
CPU 0. Every run used five warmups and 31 measured trials. The variants ran
in separate processes in the order sparse, BCH, BCH, sparse.

| Path | First median | Second median | Mean of medians |
|---|---:|---:|---:|
| sparse-EA | 8.658 ms | 8.747 ms | 8.703 ms |
| optimized BCH | 4.891 ms | 4.901 ms | 4.896 ms |

The sparse-EA outer is 1.777 times slower. It adds 3.807 ms to the isolated
outer transform, an overhead of 77.7% relative to BCH.

The runtime ratio is 7.2% worse than the static XOR ratio. The sparse
function is 3.16 times larger and uses a 2.76-times larger stack frame.
Instruction-cache pressure and spills are plausible causes. This experiment
does not separate those effects with hardware counters.

## Interpretation

The optimized circuit removes the 16,895-XOR objection to one-stage
sparse-EA, but it does not make the outer competitive with BCH. A direct
replacement in an approximately 11-ms BCH-based encoder would probably
exceed the current performance target. This statement is an inference from
the isolated 3.807-ms difference, not an end-to-end measurement.

The comparison does not weaken the one-stage distance theorem. That theorem
uses the sparse matrix ensemble and gives 41.362 bits of failure margin at
10.9% distance with RandomStepConv-M22. The benchmark fixes one sampled
matrix for code generation. It does not prove that this realized matrix
satisfies the theorem's spectrum caps.

## Artifacts

- pure_expander_accumulate/build_one_stage_vs_bch_benchmark.py reconstructs
  and verifies both circuits.
- pure_expander_accumulate/one_stage_sparse_ea_vs_bch_bench.cpp is the
  generated self-contained benchmark.
- pure_expander_accumulate/one_stage_sparse_ea_vs_bch_build.json binds the
  circuit hashes, schedules, and XOR counts.
- pure_expander_accumulate/one_stage_sparse_ea_vs_bch_peach.json records
  the environment and all timing samples.

## Recommendation

Do not replace the optimized BCH outer with the present degree-33
sparse-EA implementation. The next proof-compatible performance experiment
should reduce the sparse circuit's live set, not merely its XOR count. A
pressure-aware circuit search is worthwhile only if preserving the certified
one-stage ensemble is preferred over returning to the depth-two proof.
