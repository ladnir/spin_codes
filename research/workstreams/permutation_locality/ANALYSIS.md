# Grouped routing: first proof checks

These checks concern the existing BCH [256,128] outer and IMT (t,s)=(128,19).
They do not certify the new distribution's distance.

## Distribution

Partition the outer rows into fixed groups of g rows, where g divides 128 and
the number of outer rows. Each row sends one coordinate to each of 256 regions.
The independent variant samples a separate coordinate permutation for each row.
The shared variant samples one coordinate permutation for all rows of a group.
In each region, independently permute the groups and independently permute the
g row labels within each group. Random draws in this mathematical model are
independent and uniform over their specified permutation sets.

The implementation realizes this model with the existing seeded sampler. The
following counting identities concern the mathematical model, not pseudorandom
generator security. The production distribution is recovered at g=1.

## The one-row marginal is unchanged

Fix a row and its outer word before sampling the route. Its coordinate
permutation is uniform in either variant. In each region, the row's group has
a uniform group position, and the row has a uniform position within its group.
These two choices give a uniform position over the whole region. Across regions,
the position choices are independent. Thus the full routing marginal of this
fixed row matches the original construction.

This permits reuse of a bound whose only active row is that row. It does not
permit reusing the original multi-row occupancy distribution. Shared coordinate
permutations also correlate which regions are active for different rows.

## Local feedback enumeration

Let B be the actual 19-by-128 feedback matrix, with columns B_p. An aligned
group window starts at j*g for 0 <= j < 128/g. For a binary support S within
that window, its feedback is XOR_{p in S} B_p.

The experiment enumerates every subset in every aligned window. Conditional
on a fixed active-lane count a, a uniform window and lane permutation give
each (window, a-subset) equal probability. The total count is
(128/g)*binomial(g,a). The census records zero-feedback counts and the largest
count of any particular feedback syndrome, separately for each a.

Exact enumeration finds:

- For g=2,4,8, every nonempty subset has nonzero feedback.
- For g=4, for each fixed a>0, every (window, a-subset) has a distinct syndrome.
- For g=8, at fixed a, the largest syndrome multiplicity is at most two.
- For g=16, there are two zero-feedback subsets, both of weight eight, among
  102,960 (window, weight-eight subset) choices.

Consequently, when the IMT state is zero and just one group supplies a nonzero
input in an epoch, g<=8 forces the next state to be nonzero, unless that epoch
is the final one and the next state is not used. This statement follows from
q_next = F*q + B*x; it does not depend on the sampled transvection F.

It says nothing by itself about cancellation from a nonzero state or between
several active groups. Window positions are also dependent when several groups
occupy the same region. A full proof must retain those dependencies.

For g<=8 this also excludes an identically zero state trajectory for every
nonzero message confined to a single group, for every route in these families.
Indeed, a nonzero constituent word occupies at least 38 distinct regions, so
its group's first active epoch cannot be the final epoch. At that first epoch
the state is zero and the group's nonzero input has nonzero feedback.

## Full-group zero-state test

For a realized route, restrict messages to the first group of g outer rows.
There are 128*g message bits. For every nonfinal epoch, form the 19 linear
constraints B*x_epoch=0, using the actual BCH generator and route. A nonzero
solution keeps the state zero throughout, so the output equals the routed
outer word and has weight at most 256*g.

At K=2^20, route seed 1, both variants have full constraint rank for
g=1,2,4,8,16. This excludes this particular failure for those sampled routes
and that group only. It is not a bound on the probability of failure.

There is also a deterministic upper limit on grouping. A group with g<=128
and g dividing 128 occupies at most one epoch per region. Hence these
constraints have rank at most 256*19=4864. If 128*g>4864, a nonzero solution
must exist. In particular g=64 and g=128 cannot meet the 10% target at K=2^20:
their output-weight upper bounds are 16384 and 32768, respectively, versus
the required 209715.2. This obstruction holds for every sampled route in
these families, not merely an unlucky seed.

## Next proof interface

For independent row-coordinate shuffles, keep the existing row-weight
description but replace independent row positions by group occupancy counts.
For shared shuffles, track column-vector types within a group: the binary
g-vector in each outer column. Its Hamming weight determines active-lane count;
a shared permutation samples without replacement from these column types.

The shared variant therefore needs bounds on joint supports of tuples of BCH
words, not only the ordinary single-word weight spectrum. One option is to
bound tuples through the dimension and support of the subcode they span.
This is an analysis direction, not a completed reduction or certificate.

## A first bound on joint support

The span dimension gives an inexpensive rigorous starting point. Suppose the
g BCH words in a group span a binary subcode of dimension h>0. Let U be the
union of their coordinate supports. Each coordinate in U is one in exactly
2^(h-1) words of this subcode, whereas every nonzero word has weight at least
38. Counting total Hamming weight over the subcode in these two ways gives

    2^(h-1)*|U| >= 38*(2^h-1).

Thus |U| is at least 38, 57, 67, and 72 for h=1,2,3,4, respectively. A shared
coordinate permutation turns U into a uniformly sampled subset of that size
among the 256 regions. This is a bound on the number of active regions, not
on the active-lane counts within them or on the final SPIN distance.

For h=1, all nonzero rows of the group equal the same BCH word. There are
(2^g-1)*A_w such tuples of support size w, where A_w counts BCH words of
weight w. This makes the rank-one class accessible using the existing BCH
spectrum bounds. Higher-dimensional classes need multiplicity bounds as well
as the support lower bounds above. Start with g=4 and these classes rather
than attempting the complete joint enumerator immediately.
