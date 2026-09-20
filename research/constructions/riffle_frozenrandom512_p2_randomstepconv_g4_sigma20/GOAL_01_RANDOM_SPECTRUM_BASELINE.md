# Goal 01: random-spectrum baseline

## Question

Does the one-lap random step convolution still exhibit a low-distance family
when every active 256-bit data group enters a random rate-one-half constituent?

## Model

For each data-group position \(i\), setup samples an independent uniform
linear injection

\[
E_i:\mathbb F_2^{256}\to\mathbb F_2^{512}
\]

without conditioning. Setup then freezes the maps. All probability over the
outer maps refers to this setup experiment.
The packet permutation and the random convolution maps use their existing
independent probability spaces.

The baseline uses two field parity symbols. Each nonzero parity symbol is
encoded by extended BCH \([128,64,22]\).

## Deliverables

1. Replace the RM spectrum by the exact random-linear ensemble spectrum.
2. Optimize the inner bound pointwise at every candidate support.
3. Report the dominant outer occupation, local binary weight, packet support,
   inner episode count, and expected-count exponent.
4. Give a concrete bad family whenever the exponent is nonnegative.
5. Repeat the sparse analysis for two, three, four, and five field parities.

## Decision rule

If the random-spectrum model closes at the target distance, treat the inner
and parity structure as viable. The next phase may search for a structured
constituent with the required spectrum envelope.

If the model fails, attribute the failure to a canonical placement and state
trajectory. Do not weaken the constituent merely to simplify enumeration.

## Scope

The first calculation proves or refutes an ensemble statement. It does not
yet replace the independently sampled local maps by one reused efficient
constituent.
