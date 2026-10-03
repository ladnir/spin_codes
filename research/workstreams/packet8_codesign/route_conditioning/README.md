# Route-conditioning proof checkpoint

Separating rare routing failures from message counting substantially improves
the proof screen without changing the baseline encoder. The selected q=119
baseline contribution improves from -4,150.93 to -617.49 bits at the same
weight tilt 0.4. It remains negative. No whole-code certificate is claimed.

The best separately tested actual 24-bit-state contribution is -258.39 bits.
That result belongs to a different inner, not to the measured 16-bit-state kernel.

## The event belongs to setup, not to individual messages

Fix q of the 512 outer groups. Each region assigns them a uniform q-subset
of its 512 byte slots. A physical step contains eight such slots.
Call these positions potential slots, whether or not their later byte labels are zero.
For cap c, define

```
H_c(S) = sum over all 2,048 physical steps of min(potential occupancy, c).
```

The good-route event requires H_c(S)>=h_q for every q-group subset S.
Its failure union counts C(512,q) subsets, with no message-count factor beta^q.
On the good event, the expected low-weight message count is bounded using

```
1_good <= exp(nu*(H_c(S)-h_q)),   nu >= 0.
```

The original outer comparison then contributes beta^q. It makes every byte
at a potential slot uniform; activity has probability 255/256.
If T_k is the existing local comparison operator for k nonzero packets, use

```
T_marked[j] = exp(nu*min(j,c))
              * sum_{k=0}^j Binomial(j,k;255/256) T_k.
```

Apply ordered, without-replacement placement with exactly q potential slots
per region. Do not add a second regional binomial mixture. State remains
continuous across the 32 regions. The full expression includes exp(-nu*h_q).
The parent [proof note](../PROOF.md) gives the probability split and its premises.

For c=1, `route_counts.py` counts occupied steps by exact integer coefficients.
For c=2, `cap_counts.py` evaluates the exact weighted polynomial

```
[t^q] (sum_{j=0}^8 C(8,j) x^min(j,2) t^j)^64 / C(512,q).
```

Independent regions and the subset union give rational Chernoff witnesses.
The strict event H_c<h uses h-1 in its exponent.
These exact witnesses certify only the routing event, not the complete code.

## Exact route thresholds used

Each rational witness below passes an exact integer comparison against 2^-60.

| q | H_1 threshold | H_1 Chernoff x | H_2 threshold | H_2 Chernoff x |
|---:|---:|---:|---:|---:|
| 64 | 1,057 | 1/5 | 1,615 | 3/14 |
| 119 | 1,483 | 3/20 | 2,559 | 9/35 |
| 128 | 1,531 | 1/7 | 2,680 | 9/35 |

The first H_1 drivers evaluate floating Chernoff tilts near these rational choices.
The H_2 drivers embed the exact rational witness and use its certified integer
bit floor when adding the route-failure contribution to the message contribution.
Every message-contribution evaluation remains floating, with no outward guarantee.

## Selected results

The H_1 baseline grids checked q=64,119,128 over weight tilts
0.2,0.3,0.4,0.5,0.6,0.8,1.0 and the recorded nonnegative route tilts.
Their best combined margins were respectively -1,521.82, -2,255.30, and -2,407.42 bits.
The q=119 choice was weight tilt 0.4 and route tilt 3.

Actual A-scaled16 and actual24 operators were also checked under H_1 at weight
tilt 0.5. Their best q=119 margins were -2,466.08 and -2,225.73 bits.
These single-weight-tilt checks are not globally optimized comparisons of the maps.

The stronger H_2 gate was restricted to q=119:

| Actual map | Weight tilt | Route tilt | Combined margin, bits |
|---|---:|---:|---:|
| Baseline16 | 0.4 | 3.5 | -617.492878 |
| A-scaled16 | 0.4 | 3.5 | -532.108194 |
| Actual24 | 0.5 | 3.5 | -258.391144 |

The baseline grid covered weight tilts 0.3,0.4,0.5,0.6 and route tilts
0.5,1,1.5,2,3,3.5,4,5 where recorded in the two receipts.
The A-scaled map received one matched point, not another broad search.
The actual24 check reused its authenticated weight-tilt-0.5 local matrices
and tested route tilts 3,3.5,4. It performed no new 2^24-state census.
Its route-tilt-4 result worsens to -459.58 bits, so the tested maximum is bracketed.

A final joint gate requires both H_1>=1483 and H_2>=2559. Its bad-route
bound adds the two exact rational witnesses, without assuming independence;
the sum has a certified 60-bit floor. The message bound uses a joint local
marker on one route, not a product of separately averaged moments.
At weight tilt 0.4, nine positive multiplier pairs gave no improvement.
The best was -688.275596 bits at (nu_1,nu_2)=(0.5,3.5).
The H_2-only control reproduced -617.492878 bits exactly.
This bounded grid does not rule out other statistics or better multiplier choices.

All selected values remain negative. A failed upper-bound screen does not
establish that the sampled code has low distance. These results instead
identify which shared routing events made the earlier first moment expensive.

## Files and reproduction

| Driver | Authenticated receipts | Scope |
|---|---|---|
| `gate.py` | `conditioned_v1.json`, `conditioned_v2.json` | Baseline16, H_1, three q values |
| `variant_gate.py` | `scaled16_v1.json`, `saved24_v1.json` | Two actual alternative maps, H_1 |
| `cap_gate.py` | `cap2_v1.json`, `cap2_boundary_v1.json`, `cap2_scaled_point_v1.json` | H_2, q=119 |
| `cap24_gate.py` | `cap2_saved24_v1.json` | Saved actual24 operators, H_2, q=119 |
| `joint_gate.py` | `joint_v1.json` | Baseline16, joint H_1/H_2, q=119 |

Each receipt records map identity, complete geometry, all tested tilts, and
source hashes verified at completion. The saved24 runs additionally verify
every source pin in `../larger_state/trajectory_v1.json` and pin that receipt itself.

From the repository root:

```text
python -B -m unittest discover -s research/workstreams/packet8_codesign/route_conditioning -p "test_*.py" -v
python -B research/workstreams/packet8_codesign/route_conditioning/gate.py --tilts .4 .5 .6 --nus 0 1 2 3 4 6 --output <fresh-h1.json>
python -B research/workstreams/packet8_codesign/route_conditioning/cap_gate.py --tilts .4 .5 .6 --nus 3.5 4 5 --output <fresh-h2.json>
python -B research/workstreams/packet8_codesign/route_conditioning/cap24_gate.py --output <fresh-saved24-h2.json>
```

The 18 passing tests independently check exact small route distributions and tails,
rational witnesses, noncommuting ordered products, potential-slot thinning,
and the good-event indicator majorant. Setting route tilt zero reproduces
the earlier active-slot comparison; the numeric drivers check this identity too.

No production implementation, prior authenticated source, or encoder parameter
was modified. No encoder benchmark was run in this workstream.
The joint gate did not remove the remaining gap on its tested grid.
Next, inspect the H_2-conditioned occupancy and state trajectories before
adding more constraints. For the actual24 fallback, an H_3 condition is also
a concrete candidate because any three-byte feedback restriction has full rank.
Neither a broad implementation campaign nor a complete outward replay is
justified by these still-negative message bounds.
