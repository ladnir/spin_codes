# Kernel generation record

`PreparedEncoder.cpp`, `BankState.h`, and `BankKernel.cpp` are package-owned
implementations of the explicit banked heuristic mode. They generalize the
`row-rotate1` experiment without changing immutable `Code` sampling.
Generic routing accessors are statically dispatched; the import script preserves them.

SPIN is an independently licensed MIT library; see [LICENSE](LICENSE).
The packaged kernels were generated from the SPIN research workspace at commit
`3ced6dc0`, using `ForwardRecommended.cmake`. The wide-kernel development snapshot
was staged at commit `2582a9fb366d28750ef92afdd9b2dbe9538552c6`.
These references record SPIN's development and packaging history, not attribution
of the library to its consumer projects.

The import changes namespaces, internal symbol prefixes, and compiler
portability spellings, and adds a private setup accessor. It also ports
the standalone transpose range-direct dispatch through K=458752, reusing
the forward direct route table where present. It does not
resynthesize circuits or change the full/partial-tile kernel schedules.
The packed-bit state update is expressed without nested conditional lambdas.
Its fixed lookup tables are initialized at compile time for MSVC portability.
The setup-only overlay batches random words, uses exact reciprocal modulo
sampling, and avoids generating discarded index formats. Seeded maps are unchanged.
SetupRandom.h is package-owned; tools/setup_overlay.patch records the overlay.
A private prepared-route hook supports opt-in experiments; public Code sampling is unchanged.
Block storage and capability detection are package-owned files.
The generic transpose circuits come from the same configured build.

The package has no runtime or build dependency on its consumer projects.

| Source file | Import input SHA-256 | Imported SHA-256 |
|---|---|---|
| `Spin.h` | `b657e995bc766598588d103cbb8af572eed48db62dfbd4b90d139f64bef3fe1b` | `a2ef84ad14113443ec5ec075abec4c76132cbe79f7107c20696c9d54ffc06f19` |
| `Spin.cpp` | `3a9be018b19e50209cd456159aae89ffc34ee7bdac040a749c214c6b11112e01` | `202322fb19e55289a4dbd924042555f7715f59f7f8c5f51e1850e0e60029e66b` |
| `Inner.h` | `51c19481981366a689126b8275db98f190446026fc0a450ec8677c5e7e8b3360` | `5ae62182b4b4a2692882384e8de602c5e7d170df862bec235d0b6d2dad3c45dd` |
| `ImtRounds.h` | `d7aad4c445e734d755a5bb6df405ddbe30ea471dc6b26f0c8cd068d92eed2fbc` | `34d4b7cc072bc74f785bd4749bd29803ea01b4e8776f6a14754793686978ddbc` |
| `K16Inner.h` | `fb1556fb7425c2624792cc8fa6df5d3f80f8ac546e1f61e9f8b70b6b31615a52` | `3af4280c9097a8e2f51420da199c01cbdcf8cf98c1384822aa30b4fb9061025b` |
| `K16R2Inner.h` | `2b4ac2e7bc86e4c135c13a4f69a557afee1b18acf50e11f60c9d55b3038d6178` | `1036062f6d29678b601a698df09f5faab89292b1244b87bc11b4214aeb6fb497` |
| `Map64S12.h` | `e2b7254bc21edf0b59df60c24333935a4e7f9ed760bd0ebe778ea2f3db9628b7` | `a77129fdb2b9d9f9c10c57e04196c4fbe51f24c1616c1b727ec22c9ca360e713` |
| `LengthGeometry.h` | `057515f7eaacf9df69fea0e7dd3a4cc2ff3d9218ae8bbe366d4231524030c324` | `58fb31892a536ad6161f77536daf16facd6689b923d55a3288ee38412ce889b6` |
| `WorkspaceRouting.h` | `f14798a8fcf41ed94ceacb879a1e77d94ef13e97c960cad5c3da5d5acacff270` | `2f461862b5e6fdc5d774e77931dc984be13de7d24b005c920fa1fd4298d1206f` |
| `Fast.cpp` | `681f90e605d608e900f32582d132ce0ecea4bef0c5f868e807a90893044e9a7c` | `65ba84f9c0ea604f44e294ceeadaf4ae03c28cb2d6d1f1c13ee29727b9f55398` |
| `ForwardFast.cpp` | `a74358b009f9efa7f42aec3361c7ee5bff5da0eb1782c941cb0ee0e160dcfadf` | `4123787fe9a6bbc66fb0b98a00cd35f53a5366470d472b43603553a99c874779` |
| `generated/BchCircuit.h` | `07bec40ea14e6a71a06c28cfacdc0a087ce7a135cd525536bfa49e3f92c74b41` | `c13131ca9d02861ac5b4688849a4550196120abdce5c5aa0bec7df71eac02424` |
| `generated/BchCircuit.cpp` | `daa9d938f3bcf87a9ecaa0eefec0afa05a304a73e7f55ac7e0cb509a0cea843c` | `a86a0ff94ee74e502a8b9e734a982a71fbdd8fe9e9848963ee6a579872ccef93` |
| `generated/BchForward.h` | `bbc0ac40f1d9a99f3048e18cb26926f45ed52afe54ca5cbfd1425b144f0e3574` | `a4577b1bc6d82c4dc60afc69c2484819082a50b46b72a7988be62d3f40e20182` |
| `generated/BchForward.cpp` | `dc08f4fd9025dc18ff5df8518528a300382731365829823168014fba54505217` | `40bd32889e01f2c87ebcb453b3b79d8ab7b1b2d7101f2375fdd8a0b0d71db99c` |
| `generated/BchAvx512.cpp` | `1858e95f46a23b375e5055faa14f2d1b4f43ad6e51a0ac165c554969c94d770a` | `9439bee22dda7556b1ad0c960f8ceb2f6962b41c4c84177c5fa7b0b929e47eef` |
| `generated/BchForward512.cpp` | `5ab9159387e71dc4fdd574baf75c12fc5473ee78dacb5942fc7304a069aa5018` | `fec76bb7153fe63c68484135aebbb4ed2825bc216cbf60f98e8021105d97565d` |
| `generated/SelectedMaps.h` | `c9eaeeec8bb5b00e038d1564182a06d478e225bafc41fa542b5bfe2e1f541c72` | `f3ac0e3724edf06aa06525ec39064fc857e2fef77e72c477f803bae6a7f0e1c3` |
| `wide/Wide256.cpp` | `b9c22cd843aa6b05603a54b2eac544837e6fac0963d7acb20a1eaf844339bf97` | `b9c22cd843aa6b05603a54b2eac544837e6fac0963d7acb20a1eaf844339bf97` |
| `wide/Wide512.cpp` | `fd288bbae2c59f0a083b79b80391128a38971d6652d764eae21dc139596c7422` | `fd288bbae2c59f0a083b79b80391128a38971d6652d764eae21dc139596c7422` |
| `wide/WideKernel.h` | `2f7525672cb453240ed55332daf5907041d5dcadb4aee1dc69fed1cf23efcb70` | `8b430f53cbbd068f22f09035fd859c63663408c6e3ff2ea0e6278749ef21eb2f` |
| `wide/generated/WideCircuit.h` | `16a3eff87148262d3a16ebfef384ede19ececddf8bd5ea87b34949b33865b6f1` | `ef0dcdf3cba39e2cf13dec78c8b745dc658db66f7a0453011e6a64a46a6c585c` |
