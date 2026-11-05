# Destination-conditioned regional quota bounds

2026-11-01. The continuum obstruction in
[QUARTET_REGIONAL_CAUCHY_OBSTRUCTION.md](QUARTET_REGIONAL_CAUCHY_OBSTRUCTION.md)
applies when the quota evaluation is chosen after summing destinations.
This comparison changes that order: partition the destination states,
bound each part with its own evaluation, then sum. It retains every
physical source state and all sixteen updates within each region.
No construction parameter changes.

Use the original `(1,1,2)` adaptive operator `K_x`, marker `y=2`,
`C0=binom(8192,4096)`, and exact regional operator
`B=[x^4096]K_x^16/C0`. The fixed groups and potential are

```
J0={(a,w): w<256},        X0={1.6,1.86,2.04},
J1={(a,w): w>=256},       X1={2.04,2.45,2.76},
phi(a,w)=2^((315a+372w)/512).
```

The grids are selected for this one bounded control. In the parent's
three-state floating diagnostic, the partition `O | {P,H}` gave growth
log `258.69596`, above the roughly `205.2268`-bit bulk budget. The
partition `{O,P} | H` gave `177.10420`. The displayed `w` split agrees
with the latter partition on those three states. These diagnostics do
not bound the full operator or predict that this control will close.

## Positive regularization and the bulk bound

Put `epsilon=2^-8192`, `kappa=1+epsilon`, and

```
r_j=1_(Jj)+epsilon*1_(Jj complement),
f_j=phi*r_j.
```

The vectors are strictly positive and `f_0+f_1=kappa*phi` exactly.
Each coefficient-Cauchy bound applies separately, giving

```
kappa*Bphi
 <=F(phi):=sum_(j=0)^1 min_(x in Xj) x^-4096*K_x^16(f_j)/C0. (1)
```

The caller drops the denominator `kappa`; hence `Bphi<=F(phi)` remains
valid. The group minima are physical pointwise minima after sixteen
actions at a fixed evaluation. They are summed only afterward. No
quota evaluation changes within one chain.

## The grouped endpoint

The exact first-region endpoint is

```
b=[x^4096]K_x^15 q_x/C0,
q_x(a,w)=(2+x)^w(1+x)^(512-w).
```

For each evaluation define the common multiplier
`mu_x=max_s q_x(s)/phi(s)`. The group endpoint polynomial is `q_x*r_j`.
Since `q_x*r_j<=mu_x*f_j`, positivity gives

```
kappa*b
 <=E:=sum_(j=0)^1 min_(x in Xj)
            mu_x*x^-4096*K_x^15(f_j)/C0.            (2)
```

Again the caller drops division by `kappa`. The multiplier is computed
against `phi`, not the much smaller regularized vector `f_j`. The
fifteenth vector is saved from the same chain used in (1), so (2)
requires no additional action chain.

For this one positive potential set

```
rho=max_s F(phi)(s)/phi(s),
nu=max_s E(s)/phi(s).
```

The whole geometry comparison at threshold 60518 is

```
C0*2^-60518*phi(512,512)*rho^255*nu.                 (3)
```

The row-set union appears once. There is no additional terminal
`2^512` factor. No iteration or convergence assumption enters (3).

## The regularization loss is uniformly negligible

For this raw potential, `1<=phi<=2^687`. Moreover

```
K_x 1(a,w)=(2+x)^w(1+x)^(512-w)<=(2+x)^512.
```

Every evaluation lies strictly between 1 and 4, so
`(2+x)^2/x<9`. Using `9^4<2^13` and `C0>2^8178` gives

```
x^-4096*K_x^16(phi)/C0
 <2^687*9^4096/2^8178 <2^5821.
```

Thus the added complement contribution is less than `2^-2371` per
group, both absolutely and relative to `phi>=1`. The same bound holds
for the endpoint shortcut: `mu_x<=(2+x)^512` and
`K_x^15(phi)<=2^687*(2+x)^7680`. Taking the minimum cannot enlarge
the error beyond a uniform per-evaluation bound. Summing the two
groups adds less than `2^-2370` relative to the corresponding
unregularized group comparison. This is an algebraic estimate; it
does not certify floating rounding errors.

## Implementation and verification scope

[quartet_destination_quota.py](quartet_destination_quota.py) supplies
one fixed witness. It uses the exact-algebra Horner action and the
audited optional T-power cache. One cache remains live during each
sixteen-step chain and is then released; no cache from another
evaluation or group is retained. Each physical vector has its own
explicit common scale. Positive log sums combine the group results.

The production command is `--control`, with a 300-second runtime
budget, 295-second watchdog, 2-GiB address-space cap, 20-GiB free-space
guard, and one numerical thread. It prints compact per-evaluation,
group, and completed-witness summaries, with no raw vectors or dumps.

Independent tiny checks are in `--checks`. They use the separately
enumerated physical width-four and width-eight kernels to verify
two-step quota coefficients, regularized group identities, direct and
shortcut endpoints, dropping `kappa`, actual matrix composition,
scale invariance, and fixed-evaluation ordering. Independent static
review passed the grouped composition, common endpoint multiplier,
positive group sum, normalizers, scale handling, and single-cache
lifetime. The parent ran the tiny checks serially: 261 physical
coefficients; 168 each for grouped quota, direct endpoint, and shortcut
endpoint; 84 each for regularization identities and dropping `kappa`;
four composition and four fixed-evaluation-order checks; 84 floating
row comparisons; and 42 scale comparisons all passed.

The single full control completed in 197.291 seconds over all 33,153
physical rows. Its growth log was `642.6366722863` and endpoint-ratio
log was `642.2330577667`, both maximized at `O=(0,0)`. The full tail
upper log was `+112868.75869868135`. The minimum physical growth log
was `518.93702738364`, at `(460,1)`.

For group zero, the evaluations `1.6,1.86,2.04` were selected on
`0,10044,23109` source rows, respectively. For group one, the
evaluations `2.04,2.45,2.76` were selected on `29586,3567,0` rows.
These counts concern the comparison minima, not a physical history
distribution.

The control did not close. If its minimum growth ratio were certified,
monotonicity, homogeneity, and positive-vector comparability would
exclude every scalar supersolution below that floor for this same
regularized two-group comparison. The reported floating minimum is
not such a certificate. It does not prove an obstruction for the exact
quota operator or a construction defect.

Further full controls are parked. The next bounded mathematical test
retains explicit return paths in the exact quota operator, rather than
another potential or grid. No numerical job was run by the author;
the full-code certificate remains open.
