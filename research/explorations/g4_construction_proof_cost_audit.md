# Why the g=4 construction grew

## Scope

The historical encoder runs in 11.085 ms. The current proof-oriented encoder runs in 26.001 ms under the same 31-trial protocol.

This note identifies why each construction feature entered the proof. A feature can be necessary for the current proof without being necessary for distance.

## Proof dependencies

| Construction feature | Role in the current proof | Classification |
| --- | --- | --- |
| Four-element packet permutation | Supplies the multinomial factor in the packet-profile orbit size `Q_4(a)` | Required by the chosen `g=4` profile proof |
| Independent 64-lane bijection before packetization | Makes each packet support uniform conditional on packet weight; supplies the factors `4^(a1+a3) 6^a2` | Required by the current orbit lemma |
| Independent order inside each packet | Completes the same conditional support uniformity | Required by the current orbit lemma |
| Systematic BCH basis and accumulator | Removes an artificial low-weight trajectory in the earlier inner relaxation and raises the certified live rate | Proof repair for the current inner argument |
| Independent state permutation at each recursive step | Makes the incoming state support uniform conditional on its weight | Required by the current 65-state transfer bound |
| Three sloped outer bands | Creates the incidence structure used by the one-conditioned-row outer bound | Required by the current outer lemma |
| Uniform assignment of data blocks to `(t,l)` | Supplies exchangeability across the sloped incidence structure | Required by the current outer averaging argument |
| Independent coordinate bijection for each data block | Makes band occupancy follow the committed projection spectra | Required by the current outer witnesses |
| Independent lane order in each physical group | Makes occupied lanes uniform conditional on group weight | Required by the current packet-profile and outer arguments |
| One global puncture lane and 128 distinct tiles | Places every punctured block in one conditioned matching | Soundness repair for the conditioned-row lemma |
| Dense 24-bit graph map | Supplies the declared graph-syndrome distribution | Bound to the current graph argument |

The table does not show that a cheaper construction lacks distance. It shows that removing a listed feature invalidates at least one current proof step.

## Where the measured loss occurs

The historical inner takes 5.103 ms. The current packet-permuted accumulator takes 10.861 ms before state permutations and 12.267 ms after them.

The present measurements do not separate accumulator cost from packet-gather cost. The generated circuits differ by only 17 XORs per node. The 5.758 ms increase before state permutations therefore cannot be attributed to the circuit count alone.

The historical sequential outer takes 6.594 ms. The current fixed three-band outer takes 7.766 ms. The randomized outer takes 11.259 ms, and the graph stage takes another 1.291 ms.

The current measurements expose two large regions:

1. The packetized inner adds about 5.8 ms before state permutations.
2. The randomized outer plus graph path adds about 6.0 ms over the historical outer.

State permutations add another 1.4 ms in the isolated current-inner comparison.

## Process conclusion

The additions were proof-motivated, but the construction change was not governed correctly. The proof branch treated its evolving random ensemble as frozen before the encoder interface and performance budget were frozen.

The correct rule is stricter:

1. Freeze a pseudocode construction and a transposed cost model.
2. Benchmark every proposed construction change before adopting it.
3. State which proof obligation the change resolves.
4. Require explicit approval before the proof changes the construction.
5. Keep the previous construction and certificate branch available for comparison.

## Next decision

The current certificate proves the current randomized ensemble. It does not prove that every source of randomness is necessary for distance.

The next implementation experiment should separate four effects in the inner path:

1. cyclic split with the historical 64-group permutation;
2. systematic accumulator with the same 64-group permutation;
3. cyclic split with a four-element packet permutation;
4. systematic accumulator with a four-element packet permutation;
5. the fourth variant plus state permutations.

The outer experiment should retain the existing three measurements:

1. historical sequential outer;
2. fixed sloped three-band outer;
3. sloped outer with lane and coordinate maps;
4. graph stage.

These ablations will identify which proof requirements are expensive. We can then choose between preserving the current certificate and proving a cheaper ensemble.
