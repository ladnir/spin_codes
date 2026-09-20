# Riffle MultiBlockStripe FieldCheckpoint

This folder records the structured route that balances the two large
permutation dimensions.  At (q=32), both permutations have size 256.

The initial 9% proof attempt finds no loss in the dominant one-active-block
layer.  A proof-safe bound for reused group shifts becomes too loose at 47
active blocks.  Sampling the small group shift independently per outer
coordinate restores 55.9507 bits of modeled margin through 64 distinct active
groups.

That distinct-group calculation does not dominate shared groups.  An exact
rational counterexample uses two adjacent active lanes in one group.  The
counterexample invalidates the proposed compression lemma but still has about
138 modeled bits of margin.  The construction therefore remains unproved at
9%.

See `CONSTRUCTION.md` and `PROOF_ATTEMPT_09.md`.

`COLLISION_SURCHARGE_ATTEMPT.md` records the subsequent bounded attempt to
charge same-group dependencies.  The observed scalar surcharge is tiny and
approximately additive, but no rigorous arbitrary-background composition
bound is known.
