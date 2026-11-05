# Marker-specific regional gates

The exact three-state obstruction at marker `y=2` does not rule out a
smaller marker. This checkpoint compares a few fixed markers before
authorizing another full-population control. The encoder, random setup,
and regional label quotas are unchanged. Network capacities change only
the proof comparison.

For capacities `(a,b,2a)`, let `K_x` be the positive adaptive population
operator and restrict every intermediate population to
`O=(0,0), P=(512,0), H=(512,512)`. Write this restriction as `L_x` and
put `C0=C(8192,4096)`. For a positive three-state vector `v`, the tested
regional comparison is

```
F_3(v)_s = inf_(x>0) x^-4096 (L_x^16 v)_s / C0.
```

A positive subeigenvector of this restriction extends by zero to the
full population space. A certified lower inequality for `F_3` can
therefore obstruct the full *sourcewise Cauchy comparison*. An upper
bound for `F_3` controls only paths confined to these three states; it
does not bound the full population operator. Neither direction is an
encoder failure claim. The calculations below use floating arithmetic
and have not been converted to directed certificates.

The free endpoint is general-network specific:

```
q_x(u,w)=(1+x)^(512-u)
          (y^(a-b)+x)^(u-w) (y^a+x)^w.
```

Thus the restricted endpoint is `[x^4096]L_x^15 q_x/C0`. The complete
tail comparison includes `C0*y^-delta`, with one row-set union and no
additional terminal `y^512` factor. The full exact endpoint is at least
one because its cover cost is nonnegative. This justifies using
`log2 C0-delta*log2 y+255*log2 lambda` as a lower comparison floor when
the restricted regional lower inequality has been established.

## Parent-run bounded gates

[quartet_general_three_state.py](quartet_general_three_state.py) builds
its local coefficients directly from the general gain, loss, and cost
rules. Its independent checks passed: 60 local coefficients, 108
physical restriction entries, 36 powers, 36 quotas, 36 endpoint entries,
three minimization controls, and eight gauge checks. Independent static
source review also passed for the networks used here.

The parent ran the following bounded controls serially. The budget
column is the regional logarithmic budget used for the comparison gate;
passing it only clears this restricted obstruction test.

| Capacities | Marker | Integer deficit threshold | Restricted continuous lower diagnostic | Gate budget |
|---|---:|---:|---:|---:|
| `(2,1,4)` | `3/2` | 184116 | 689.7667 | 390.258 |
| `(1,1,2)` | `3/2` | 60518 | -91.4031 | 106.7278 |
| `(1,1,2)` | `3` | 60518 | about 1342 | about 344 |
| `(1,1,2)` | `4` | 60518 | about 2199 | about 442 |
| `(2,1,4)` | `5/4` | 184116 | -80.7427 | about 200 |

At `(2,1,4),y=3/2`, six half-updates gave pure-state log potential
`(0,379.79695,586.37073)` and a full-sourcewise lower-tail diagnostic
of `+76374.7`. The control took 0.163 seconds. A full physical run at
that marker is parked; this result does not exclude the partial network
at a different marker.

At the original network and `y=3/2`, the same bounded procedure gave
potential `(0,236.036407,467.698388)` and restricted upper growth
`-74.9974`. The three continuous row minima occurred near
`1.38772772,1.63743799,1.80978644`. The coarse grid `{1,3/2,9/4}` instead
gave upper growth `116.37`, so it would falsely reject this promising
gate relative to the useful budget. The selected fixed rational grid
for the full control is `{7/5,33/20,9/5}`.

At `(2,1,4),y=5/4`, the restricted upper growth was `-74.1227`, with
potential `(0,236.977,350.2606)` and row minima near
`1.38908,1.58899,1.67121`. This is a second admissible gate, not a
full-population result. Its source-normalized marker is `y^a=25/16`,
close to `3/2` for the original network. Comparing markers through
`y^a` avoids confusing a capacity rescaling with a new control.

## Next full control and scope

[quartet_general_horner_region.py](quartet_general_horner_region.py)
prepares one full-population witness for the original network at
`y=3/2`, potential log `(236u+232w)/512`, and the three rational quota
evaluations above. It reuses the independently checked positive Horner
action and its input-independent power cache. It keeps all 33,153
physical populations, composes sixteen updates before taking a quota
minimum, and uses the general endpoint formula. Its source includes
small exact physical quota and endpoint checks for both networks.

The caller passed independent static review and exact small checks:
4410 physical coefficients, 42 quota rows, 42 endpoint rows, 84 cached-chain
rows, 84 unit-endpoint rows, 42 scale comparisons, and four composed
comparisons. The intended first production comparison is one bounded
witness, followed by a scope review of its result. It was launched before
the user paused the goal; no completed result was retained. The full
proof remains open. No failed-marker
diagnostic here rules out exact regional coefficient extraction, a
different marker, or the fixed construction.
