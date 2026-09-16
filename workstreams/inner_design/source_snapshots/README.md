# Historical generator inputs

The digest-named directories preserve the exact `bare_bch_rm2sub/Spin.cpp`
and `CMakeLists.txt` used by the mixer and quarter-outer manifests, before
the optional workspace optimization.
Its SHA-256 remains `053812578e3f5038eb067d3b9d80b1bb1eb67326650b3cbf382a670aa48ba229`.

The CMake snapshot has SHA-256
`10de11c137a16ab717b4c93df36b25a0d36940879c1f71bd07021ea0859c8a91`.

The balanced certificate verifier uses these snapshots only for those historical
build inputs when the supported implementation changes. It still checks
the original recorded digest. Live generated candidates, numerical inputs,
and certificate artifacts must match their existing hashes without substitution.
This snapshot is provenance data, not the supported build input.
