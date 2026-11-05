# Random pairings with input constraints

The marginal matching moment in [RANDOM_PAIRING.md](RANDOM_PAIRING.md)
is dominated by rare matchings that do not cross the current support.
Such matchings also fix more inverse inputs. This continuation counts
those fixed inputs together with the next state and output weight.

The construction is unchanged from that note. These are local identities
and small exact checks, not a whole-code distance certificate.

## Why births matter

Fix the current support H and next support H'. A live pair intersects H;
let F contain the coordinates of all live pairs. On its complement,
the inverse recurrence gives X=Q_next. Consequently:

- A coordinate in H' outside F fixes a nonzero input value.
- A coordinate outside both H' and F fixes a zero input value.
- Coordinates in F incur the inverse field-point cost.

We call the first constraint a pin. A matching with few crossings reduces
the field-point cost, but can increase the number of pins. Counting these
effects separately can select incompatible worst cases.

## The pin identity does not require fixed matchings

Fix a complete state-support history starting at Q_0=0. Select the
coordinates with zero output between successive active states. Each such
coordinate imposes Q_next=delta Q. These equalities form disjoint paths
along individual state coordinates, regardless of the matchings.

Let A count active state coordinates across time, and let Z count the
selected equalities. The resulting amplitude space has dimension a=A-Z.
Ignoring the required nonzero inequalities gives q^a assignments.

Every pin lies at the first vertex of a distinct path component. Indeed,
its incoming pair is empty, so its preceding state coordinate is zero.
Initial active coordinates are pins for the same reason. If H_pin is the
number of pins, projection to the pinned values has rank H_pin and
kernel dimension a-H_pin. The terminal state is free throughout.

This extends the identity in [UNIFORM_PAIR_PINS.md](UNIFORM_PAIR_PINS.md)
to independently changing matchings. It does not assert independence of
the outer evaluations at these coordinates.

## A weight-indexed transfer count

Let n be the state width, q the field size, and w the incoming weight.
Let P_w(k) be the exact probability of k crossing pairs from the previous
note. Use z to mark next-state weight and y to mark output weight.
Assign a positive price alpha to each pin and beta to each fixed zero.

The following polynomial averages a relaxed amplitude count over the
uniform matching:

    K_w(z;y) = sum_k P_w(k) q^(-w-k)
               [y+(1+qy)z]^w
               [1+qyz]^k
               [beta+qy alpha z]^(n-w-k).

The three coordinate classes explain each factor:

1. An incoming active coordinate either becomes zero, emitting a nonzero
   output, or remains active. In the latter case a zero output identifies
   its amplitude with the preceding amplitude; a nonzero output retains
   a free amplitude. The factor is y+(1+qy)z.
2. An inactive coordinate in a live pair either remains zero or becomes
   active with a nonzero output. The factor is 1+qyz.
3. A coordinate in an empty pair either fixes a zero input or becomes a
   pin with a nonzero output. The factor is beta+qy alpha z.

There are w+k live coordinates, giving q^(-w-k). The accepted SL2 density
loss is excluded here; its global factor less than 2 is paid once.

The coefficient of z^v sums next supports of weight v. This is an exact
identity for the relaxed count, not a transition probability. In particular,
q-valued amplitudes overcount nonzero assignments, and rows need not sum
to one. Independent uniform matchings make the local count depend only
on w; they do not remove correlations among RS evaluations across regions.

## Connecting the prices to the outer

Fix a suffix of R regions and an active set S of s outer rows. For
R>128, write k_outer=R-128 for the shortened RS dimension per row.
Let J_zero count fixed-zero inputs routed into S. The earlier counting
argument bounds the number of compatible amplitude/message pairs by
both q^(a+d_s) and q^(a-H_pin+k_outer*s-J_zero), where d_s bounds the
remaining outer dimension on random inverse coordinates. Impossible
nonzero-row patterns contribute zero.

For 0<=theta<=1, interpolation gives

    q^[a+theta*d_s+(1-theta)(k_outer*s-J_zero-H_pin)].

Under an independent Bernoulli-p surrogate for membership in S, the
local prices are therefore

    alpha = p q^(-(1-theta)),
    beta  = 1-p+p q^(-(1-theta)).

Pins must belong to S. A fixed zero costs q^(-(1-theta)) only when its
row belongs to S. Actual routing has exactly s selected positions per
region, not independent membership. Conditioning each region costs at
most the reciprocal of

    [binom(8192,s) p^s (1-p)^(8192-s)]^R.

The outer dimension term, row-set choices, conditioning cost, output
threshold, and accepted-map density are not included in K. Moreover,
d_s can depend on the route and matching history. It cannot be pulled
outside the average without a uniform bound or a corresponding partition.
These are remaining obligations for the full-size calculation.

## A diagnostic that distinguishes the joint count

Take n=512, w=256, q=2^32, p=1/2 and theta=0. Fix b inactive coordinates
that must become active. For a matching with k crossings, let j of those
births be covered by live pairs. Apart from factors common to all these
matchings, the weighted term is

    P_w(k) Pr[j|k] q^-k alpha^(b-j)
        beta^(n-w-k-b+j).

For b=0 the weighted mean crossing count is approximately 1.42e-14:
the no-crossing term still dominates. For b=256 the weighted mean is
approximately 171.001, and the no-crossing fraction is about 8.30e-123.
Here every crossing replaces a pin, so its field-point cost is offset
by releasing that pin. The surviving tilt favors crossings.

These are weighted sums for two fixed support patterns, not actual
conditional trajectory probabilities or distance margins. They show why
birth constraints can change the dominant contribution. They do not
resolve histories that remain confined to the same support.

## Exact checks and next test

`diagnose_random_pairing_pins.py --test` compares 30 polynomial rows with
literal enumeration of matchings and next supports through width 6.
It also checks the pin identity on 128943 three-state histories with
changing width-4 matchings and selected zero-output constraints. All pass.

`diagnose_random_pairing_outer_toy.py` tests composition over GF(4), with
four independent repetition rows, two independently routed regions, and
two width-4 updates. It uses uniform pair images as the comparison model,
not accepted SL2 matrices. For a fixed active row set, it enumerates all
nonzero messages, the relative route, matchings, and live pair images.
The joint count upper-bounds every output-weight coefficient.

For one active row, the expected number of outputs of weight at most 3
is 69/32, versus the bound 115/16. For two active rows these are 83/128
and 617/48. Expectations may exceed one. The slack includes relaxed
amplitudes and outer constraints; these checks validate accounting, not
tightness at production parameters. The runs enumerate 3456 and 114048
images respectively. All runs are serial, under 30 seconds, and write
no data files; a pre-run guard requires at least 20 GiB free.

The [projection continuation](RANDOM_PAIRING_PROJECTION.md) combines actual
RS shortening with a forward count and closes two constant-support families.
For the full transfer, next combine the count with actual RS shortening and per-region
conditioning, starting with confined-support and activation/reset
histories. Exclude the all-zero message explicitly. Test these classes
before attempting a complete 4096-step count or native pairing benchmark.
