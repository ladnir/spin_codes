# Independent audit of the state-eight density comparison

The finite source/route comparison in `state8_density_screen.py` is valid with the symmetric Cauchy packet cover. The original sixteen-wide, pure-return screen nevertheless leaves a large middle-weight gap. The numbers below are floating screening results for isolated exact input weights, not an interval certificate or a construction-failure lower bound.

## Formula and conditioning

Fix outer length `b`, dimension `k`, redundancy `p=b-k`, and `L` rows. Put `N=bL`, output threshold `e=ceil(N/10)-1`, and exact nonzero inner-input weight `W`. Set `alpha=W/N`. The screen evaluates

\[
\frac{(L+1)^b}{\tau(x\lambda)^W y^e}
\left[(1-\alpha)^b+
\tau^{-p}(1-\alpha+\tau\lambda\alpha)^b\right]^L
e_0T(x,y,\tau)^{N/t}\mathbf1.
\]

For each fixed inner-input support, inverse regional routing assigns independent fixed-size subsets of row labels. Compare these subsets with independent Bernoulli columns at their respective mean densities. Conditioning a column on its modal count costs at most `L+1`. There is no additional union over regional profiles.

The exact feasible row factor is one for an empty row and `tau^(r-p)` for an occupied row of size `r>p`. Shorter occupied supports are infeasible. The displayed row polynomial relaxes the second case to every row size, including an extra empty-row term. This only increases the bound. Insert `lambda^W` before conditioning, then use arithmetic--geometric mean separately on the two nonnegative products across regions.

The packet bound is uniform over routes and common source-form dependencies. Its remaining route dependence is therefore only through the shortened source dimension, which is handled by this row comparison. The factor `tau^-1` is the nonzero-kernel rank threshold, not a second projective-word normalization. The transfer starts in actual zero state and permits either terminal state.

The implementation's `binom(t,w)` coefficient counts packet input placements once. The external `x^-W` and `lambda^-W` factors extract the same fixed total input weight from the packet and source bounds. The packet return correction in the original sixteen-wide file replaces the coarse weight lower bound; it does not add another independent return event.

An exact small regression enumerated every regional profile for `b=L=3,p=1`, every routed subset assignment, three rank tilts, and three source tilts. All 558 finite source-profile comparisons passed. This supplements the existing proof; it is not the proof for the large parameters.

## Executed baseline points

For GRS[256,176] at requested `K=2^20`, the padded geometry is `L=6016`, `N=1540096`, and `e=154009`. The following baseline used width sixteen, state eight, weighted state flats, improved pure-return output minima, and all finite source-conditioning factors.

| Requested density | Exact input weight | Floating margin in bits |
|---|---:|---:|
| .025 | 38502 | 377.347960 |
| .03 | 46203 | 1939.954484 |
| .12 | 184812 | 1567.982461 |
| .15 | 231014 | -17533.205813 |
| .20 | 308019 | -67826.796603 |
| .24 | 369623 | -121923.729778 |

The three negative points were independently restarted from both default optimizer seeds. Their best objective values agreed within `1.1e-7` bits. This confirms the executed numerical gap, not global optimality of the witness search.

The six-point replay command was:

```sh
python -B research/workstreams/char2_field_spin/theory/state8_density_screen.py --outer 256 176 --alpha .025 .03 .12 .15 .20 .24 --iterations 80
```

These are historical **pure-return** results. Later joint-return options change the local upper bound and must be labeled separately. Neither interpolation between positive points nor the sign of a failed upper bound supplies a certificate.

A bounded symmetric four-state comparison with the same pure-return method was worse: margins at `.025,.12,.15,.20,.24` were approximately `-9360,-117135,-178099,-298341,-408077` bits. That changed four-state inner is not the older asymmetric cubic construction. Its weak pure-return cover does not decide whether a different four-state proof or design is preferable. Bounds for distinct inners cannot be spliced into one certificate.

## Dense-class scope

The existing GRS[32,22] membership exponent `1.28` transfers exactly to GRS[256,176] for inverse-input density at least `.24`. For the same rational tilt, write the old row bracket as `X+Y`. The new bracket is `X^8+Y^8`, which is at most `(X+Y)^8`.

The finite conditioning factor changes with the outer length. The conservative bound `floor(N/100)-256*bit_length(L)` gives 12072 bits at K20 and 1116 bits at K18. It gives **-1321** at K16, so that coarse finite estimate does not establish forty bits there. The old GRS[32,22] size margins must not be copied to the larger outer.

## Next gate

The baseline motivates a wider packet and a joint return/output equation cover. Their local statements and scalar tests are recorded in [STATE_CAPACITY_REPAIR.md](STATE_CAPACITY_REPAIR.md). New floating points still require directed interval coverage and explicit matching to the sparse and dense classes. No whole-code or runtime result follows from this audit.
