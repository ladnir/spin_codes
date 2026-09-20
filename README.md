# SPIN codes

Fast binary linear encoders for cryptographic applications. The standalone
C++20 library lives in [`spin/`](spin/README.md); the paper and supporting
research live alongside it and are not build dependencies.

## Build the library

```sh
cmake -S spin -B out/spin -DCMAKE_BUILD_TYPE=Release
cmake --build out/spin --config Release -j2
ctest --test-dir out/spin -C Release --output-on-failure
cmake --install out/spin --config Release --prefix /path/to/install
```

The library supports Linux/GCC and Windows/MSVC on x86-64 with AVX2.
AVX-512 kernels are selected at runtime where available. ARM is not yet supported.
No Python, TeX, libOTe, or Hypercat dependency is required to build the library.

```cmake
find_package(spin 0.1 CONFIG REQUIRED)
target_link_libraries(my_target PRIVATE spin::spin)
```

See the [library guide](spin/README.md) for API examples, supported sizes,
generic XOR-element support, forward and transposed encoding, and tests.
The library exposes explicit code parameters; protocol integrations choose
their own parameter and security policies.

## Repository layout

| Path | Purpose |
|---|---|
| [`spin/`](spin/README.md) | Self-contained encoder library, public API, kernels, and correctness tests. |
| [`paper/`](paper/README.md) | SPIN manuscript, figures, and paper build instructions. |
| [`artifact/`](artifact/README.md) | Paper artifact guide and historical packaging/reproduction tools. |
| [`workstreams/`](workstreams/) | Research implementations, selected results, and proof-development records. |
| `constructions/`, `explorations/`, `scripts/`, `bch_spectrum_work/` | Supporting research, not library dependencies. |
| `BA_paper/`, `enumerator_paper/`, `expander_codes/` | Related manuscripts and research. |

Research paths remain stable because scripts and proof manifests refer to them.
Library consumers need only `spin/`. SPIN is [MIT licensed](spin/LICENSE).
The [generation record](spin/PROVENANCE.md) records the kernel snapshots and hashes.
The library license does not change the terms of the paper or research archives.

## Paper

The manuscript can be built independently of the library:

```sh
cd paper
latexmk -pdf -outdir=../output/pdf -jobname=spin_codes_draft main.tex
```

See the [paper guide](paper/README.md) for dependencies and supporting material.
Compiling the manuscript does not require replaying numerical proofs.
The [selected finite results](workstreams/inner_design/finite_migration/PAPER_RESULTS.md)
and [reproduction guide](artifact/REPRODUCING.md) describe the author-side evidence.

## Contributions and data

Commit code, documentation, compact selected results, and manifests—not build
products or experiment archives. The repository hygiene check rejects tracked
files over 5 MiB and bulk-output formats. See [publication policy](GITHUB_PUBLISH_POLICY.md).
Run performance benchmarks serially.

GitHub is the authoritative workspace. Overleaf is an optional paper-only mirror.
