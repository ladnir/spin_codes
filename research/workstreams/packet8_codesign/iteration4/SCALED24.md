# A-only scaling of the 24-bit expansion

This gate tests whether cheap byte scalings improve the remaining low-occupancy bound. The selected scalings give less than 0.4 additional bits at the q16 test point. They do not justify an implementation change on this evidence.

The state is `(a,b,c)` in the AES polynomial-basis field `GF(256)`. For physical packet index `h=0,...,7`, replace only the expansion by

\[
 (A_\sigma(a,b,c))_h=\sigma_h(a+hb+h^2c).
\]

The feedback remains

\[
 C(X)=\left(\sum_hX_h,\sum_hhX_h,\sum_hh^2X_h\right).
\]

The recurrence remains `y=X+A_sigma*s`, followed by `s'=M*s+C*X`. Each step has a fresh independent uniform `M` in `GL(3,GF(256))`. The outer code and the 32-region byte route are unchanged. The initial state is zero, and there is no final flush.

Every tested expansion has rank 24. The scaled maps need not satisfy `C*A_sigma=0`; the original-order local comparison does not require that identity. The character tables used for the weighted feedback distribution are those of `C^T`, not those of `A_sigma^T`. The implementation preserves the former and separately tests the latter.

## Selected gate

All rows below use `K=65536`, `N=131072`, weight cutoff 13107, occupancy `q=16`, tilt `theta=.06`, and fractional exponent `.4`. The comparison first forms the potential-input operators, changes to the potential-birth basis, and then groups four physical steps. The subset union factor stays outside the fractional power.

| Scales | Full expansion minimum | Rank of `C*A_sigma` | Estimated class margin | Gain |
|---|---:|---:|---:|---:|
| All one | 8 | 0 | 17.354182792 bits | — |
| First four one; last four `0x17` | 10 | 16 | 17.580442186 bits | 0.226259394 bits |
| First four one; last four `0x2e` | 10 | 16 | 17.571650491 bits | 0.217467699 bits |
| `sigma_h=3^h` | 11 | 24 | 17.728858222 bits | 0.374675430 bits |

The two-band choices maximize `min_{a!=0}(wt(a)+wt(d*a))` at four, with a favorable low-weight histogram. Consequently their constant-state words have weight at least 16. Their full expansions still contain weight-10 words. The proxy therefore does not establish the full expansion minimum or the distance of the complete code.

These are floating proposal values for one occupancy, not outward endpoints or a whole-code certificate. They neither establish the 40-bit target nor provide a low-distance counterexample. No encoder timing was run.

## Validation and replay

Each variant has a fresh census of all `2^24` states and fresh local operators. The all-one local operators reproduce the existing `.06` cache bit-for-bit. Five tests check literal field evaluation, binary adjoints, unchanged feedback, full expansion rank, and exhaustive small-field profiles. The largest alpha-one macro relative error is `2.94e-14`; the largest regional log discrepancy is `9.95e-13`. All saved source hashes authenticate.

`scaled24_v1/` contains one local cache and one G4 receipt per variant. Each receipt identifies the literal map and pins its sources. The selected fractional gate uses the all-log global backend. Reproduce into a fresh directory:

```powershell
python -m unittest discover -s research/workstreams/packet8_codesign/iteration4 -p test_scaled24.py -v
python research/workstreams/packet8_codesign/iteration4/scaled24.py --output-dir tmp/scaled24-replay
```

The kernel analysis estimates three additional GFNI operations per payload half for two bands, versus seven for `3^h`. This is an operation count, not a measured slowdown. The transpose emission through `C^T` stays unchanged; the extra work enters the feedback through `A_sigma^T`.

The next priority is the wider outer with a positive all-occupancy floating bound. Further scaling searches would need a stronger reason than the constant-state proxy.
