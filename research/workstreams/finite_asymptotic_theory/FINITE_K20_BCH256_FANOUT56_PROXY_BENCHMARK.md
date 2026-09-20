# Fanout-56 cost in the frozen BCH256 encoder structure

The optimized proxy takes 10.650813 ms with zero fanout layers and
20.454284 ms with 56 layers. Thus, Fanout-56 adds 9.803471 ms by the
difference of medians. The total time increases by a factor of 1.920443.

This measurement is a performance proxy. It is not a benchmark of the exact
BCH250 construction in the finite 11% certificate.

## Measured encoder

The proxy starts from the frozen Structured SPIN transpose encoder. It keeps
the following components unchanged:

- the extended BCH256-128 outer transpose;
- the RM2Sub-S19 reverse emitter;
- the packed 24-bit route;
- four tiles of 2,048 outer rows; and
- the fused emit, route, fanout, and outer-transpose pipeline.

The proxy samples 56 independent ParityFanout-31x33 schedules for each of
8,192 rows. Each schedule selects 33 target coordinates and 31 disjoint
source coordinates. The transposed encoder applies layers 55 through 0.

The zero-layer control uses the same binary, input, setup object, workspace,
and measurement loop. It skips only the fanout action.

## Result

The benchmark pins one process to logical processor 15 on an AMD Ryzen 9
7950X. It performs three warmup pairs and 21 measured calls per mode.
Measurement order alternates between the two modes.

| mode | median | mean | minimum | maximum |
|---|---:|---:|---:|---:|
| zero layers | 10.650813 ms | 10.667073 ms | 10.551919 ms | 11.003533 ms |
| 56 layers | 20.454284 ms | 20.461394 ms | 20.295016 ms | 20.758693 ms |

The median paired increment is 9.799123 ms. The difference of mode medians is
9.803471 ms. Setup takes 310.983616 ms and is excluded from all online times.
The explicit 56-layer schedule occupies 28 MiB.

The zero-layer median is 1.60% below the frozen one-layer receipt of
10.823872 ms. That proximity confirms that the proxy preserves the relevant
production structure.

## Correctness checks

The benchmark first checks the RM2Sub-S19 reverse emitter against its dense
transpose oracle. It then constructs two staged paths outside the fused method.
The zero-layer and 56-layer fused outputs agree block-for-block with their
respective staged outputs. The two candidate outputs are distinct.

These checks cover the measured transposed implementation. They do not prove
equivalence to an ordinary forward encoder.

## Relation to the finite 11% construction

The finite certificate uses 8,576 rows of a fixed BCH250-124 subcode. The
proxy uses 8,192 BCH256-128 rows. Therefore, the exact proof-model fanout has
4.6875% more row-layer applications:

\[
 \frac{8576}{8192}=1.046875.
\]

Linear scaling by row count estimates a 10.263009-ms fanout increment for the
exact row count. This estimate is not a measurement. It does not model the
BCH250 circuit or the final shortening projection.

The benchmark generator uses a fixed `std::mt19937_64` seed. It samples the
correct combinatorial schedule shape, but it is not the proof-compatible
setup sampler required by the finite certificate.

## Interpretation

The optimized path rejects the earlier suggestion that Fanout-56 costs only
about 20% end to end. In the preserved production structure, 56 layers add
about 92% to the zero-layer median. The discarded direct-scatter prototype
and its receipts have been removed.

The result still supplies a useful optimization target. The present proof
requires 30,256,128 row-local block XORs for the exact BCH250 row count. A
smaller certified layer count now has direct performance value.

## Reproduction

The modified sources are under `proxy256_fanout56/`. On Peach, the benchmark
was compiled with

```text
g++ -std=c++20 -O3 -DNDEBUG -march=native -pthread -I. \
  RiffleFieldCheckpoint_Bench.cpp -o RiffleFanout56Proxy_Bench
```

The single benchmark command was

```text
taskset -c 15 env RIFFLE_BENCH_CPU=15 \
  ./RiffleFanout56Proxy_Bench 21 fanout56-proxy
```

The complete samples, hashes, and environment are in
`bch256_fanout56_proxy_peach_7950x.json`.
