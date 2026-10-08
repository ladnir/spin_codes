# Local output mixing: exact gains, costs, and proof obstruction

A local mixer can spread sparse packet outputs without another global
route. It is not a free improvement: invertible mixing can also cancel
coordinates, and it does not repair the old return-rank obstruction.
The current screen therefore keeps local mixing as a possible ingredient
of a cheaper inner, not as a certificate-compatible replacement.

## Concrete candidate

Keep the current recurrence through its raw packet output `Y in GF256^16`.
Partition the sixteen coordinates into eight pairs. Setup independently
samples one uniform invertible two-by-two GF256 matrix for each pair and
packet. Apply those matrices only to the emitted output; do not feed the
mixed output back into the state. All setup choices remain fixed across
messages. This candidate uses no multiplication in a larger field.

For adjacent pairs, a ZMM register holds four 128-bit records. One packed
pair swap, two GFNI products, and one XOR implement the two independent
pair matrices. The sixteen-record packet therefore adds eight GFNI
products, four lane shuffles, and four XORs. These counts exclude loads,
coefficient materialization, and any random matching of packet positions.
The mixer adds 32 fresh field coefficients per packet. Their expanded
loads or compact expansion instructions are additional costs.

The dense state-eight inner rises from 70 to 78 GFNI instructions per
packet. Combining this mixer with the 62-product state-six layout merely
returns to 70 products, with more shuffles and setup data. Neither count
establishes an end-to-end speedup. No native kernel was added or measured.

## Exact fixed-word law

Fix a pair before sampling its independent uniform GL(2,q) matrix. A
nonzero pair becomes uniform over the `q^2-1` nonzero pairs. For `0<y<=1`,
its output-weight moment is

    g(y) = (2y+(q-1)y^2)/(q+1).

Zero pairs stay zero. This statement applies to a fixed word, not to the
rank of a whole source subspace conditioned on previous zero equations.

For an optional uniformly sampled perfect matching within a packet of
even width `b`, fix an input word of weight `t`. Let `h_t(y)` average over
the matching and all pair matrices. A positive selected-zero expansion is

    h_t(y) = sum_{j=0}^{b/2} c_{2j}(y) binom(b-t,2j),
    c_{2j}(y) = binom(b/2,j)/binom(b,2j) g(y)^(b/2-j) (1-g(y))^j.

All odd coefficients are zero. To derive this formula, write the moment
for a fixed matching as `g^(b/2)` times `(1/g)` for every pair whose two
coordinates are zero. Expand over selected empty pairs. In a uniform
matching, the endpoints of `j` specified pairs form a uniform subset of
size `2j`. Thus their probability of all being zero is
`binom(b-t,2j)/binom(b,2j)`.

For `b=16`, `q=256`, and `y=.9`, the fixed-word moments are:

| Raw output weight | Without mixer | Uniform matching and GL2 |
|---:|---:|---:|
| 1 | .900000 | .810700 |
| 2 | .810000 | .667466 |
| 4 | .656100 | .472665 |
| 8 | .430467 | .279340 |
| 12 | .282430 | .204172 |
| 16 | .185302 | .186588 |

Smaller moments help a low-weight tail bound. The mixer improves sparse
fixed words but slightly worsens fully occupied packets through new
cancellations. The uniformly sampled matching is part of this formula;
an adjacent fixed pairing does not have a weight-only law.

A cheaper fixed all-nonzero MDS pair matrix after fresh input diagonals
does not even give the same positive expansion. Its moments at pair
weights one and two are

    h_1=y^2,
    h_2=y^2+2y(1-y)/(q-1).

For `0<y<1`, the coefficient `h_1-h_2` is negative. Inserting these moments
into an existing positive selected-zero proof is invalid without a new
dominating envelope.

## What a compatible rank proof must retain

The existing SPIN certificate bounds surviving source rank, not just the
expected weight of a fixed word. The fixed-word expansion above is not an
automatic substitution for that transfer.

There is a small exact rank interface for further work. Fix two output
forms whose joint rank on the current quotient is `r in {0,1,2}`. Sample
a fresh uniform GL(2,q) matrix. For `tau>=1`, selecting one mixed output
zero gives

    E[tau^(-rank)] = 1,                              r=0;
                    1/(q+1)+q/((q+1)tau),            r=1;
                    1/tau,                          r=2.

Selecting both outputs gives `tau^(-r)`, since the pair matrix is
invertible. Selecting neither gives one. These formulas are exact scalar
moments. A following state-cap continuation also depends on how the new
equation intersects the incoming state forms. That intersection needs a
joint bound; the scalar moments do not establish one.

A valid whole-code replay must therefore integrate the pair ranks and
state continuation. It must also retain the relation between mixed output
zeros and return equations. Charging the fresh pair matrix twice would
not give a valid bound.

## Cheap output mixing preserves a known bad coordinate flat

The archived coefficient-state emission map has this exact witness in
phase zero: state coordinates `(0,4)` with values `(1,80)` emit zero in
positions `(0,1,2,3)`. The four corresponding output forms have joint rank
one on those two state coordinates. With fresh nonzero state diagonals,
their simultaneous-zero event has probability `1/255`, not `1/255^4`.

Any block-diagonal invertible output mixer whose blocks are contained in
`{0,1,2,3}` preserves all four zeros. This includes adjacent GL2 pairs and
an adjacent GL4 block. Therefore those mixers do not repair the archived
coordinate-rank counterexample. Mixing across the flat's boundary can
change this particular witness, but then requires a new rank analysis.

Postprocessing the output also leaves the state-return equations
unchanged. The coefficient-state design's rank-five two-input return does
not become an eight-pivot return after an output mixer. Better return
output weight could still help, but the original return charge cannot be
imported.

Finally, an invertible pair matrix can take a weight-two vector to weight
one. For example, `[[1,1],[1,2]]` maps `(2,1)` to `(3,0)` in AES GF256.
The old distance certificate therefore transfers deterministically only
with a factor-two loss in its distance bound, not unchanged. A stochastic
distance improvement requires a new proof over the new setup distribution.

## Decision and replay

Do not append this mixer to the dense inner as a performance optimization.
Its arithmetic and setup data increase, while its whole-code advantage is
unproved. A more plausible use is inside a cheaper packet map, with a
small joint rank transfer derived before native implementation. Even there,
adjacent postmixing alone does not fix the known cheap-basis failures.

The exact screen and seven tests cover every GL2 matrix and fixed word
over GF4, every local matching and support at width eight, both pair laws,
the positive zero expansion, the scalar rank moments, and the preserved
GF256 coordinate flat.

```text
python -B research/workstreams/char2_field_spin/theory/gf8_local_mixing_tests.py
python -B research/workstreams/char2_field_spin/theory/gf8_local_mixing_screen.py
```

These are algebraic checks, not benchmarks or a whole-code certificate.
The next useful experiment is a cheaper width/state geometry with a valid
packet transfer; it offers an actual reduction in work before adding a
local mixer.
