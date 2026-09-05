# Correcting the weight-38 application claim

Updated: 2026-09-04

Current status: the separate weight-40/42 tests have now passed full verification.
Using all three tested caps, the full outward Q1+Q2+Q3 bound has 40.0476888685
bits (display approximation). The sufficient remaining target,
sum_(q=4)^8192 F_q<=2^-45, has now been met by a bound below 2^-83.
See RANDOM_MODEL_FULL_DISTANCE_BOUND.md for the all-occupation result.
This progress does not change the historical correction below: weight 38 alone
did not close the application bound. The completed random-inner distance bound
is conditional on all three tested shell caps; it is not an unconditional
deterministic BCH spectrum theorem or an RM2Sub transfer.

The weight-38 sampling test addresses a valid BCH shell-count target. That
target does not, by itself, certify the full RandomStepConv occupation-one
functional at 2^-40. Earlier notes and the Sol handoff incorrectly said it did.

The source of the error is visible in
`code/certify_random_inner_threshold_outward.py`. The script first reserves
a rigorous upper bound for weights 52 and above. It then distributes the
remaining budget across weights 38 through 50 using one common multiplier
of a binomial reference spectrum. Its reported weight-38 cap is that shell's
allocation within this joint envelope. The computation never substituted
the other shells' actual BCH-sandwich LP caps when declaring this allocation.

Let c_w be the directed coefficient for the paired weights w and 256-w.
Let U_w be the current multiplicity cap for the length-256 code C.
The audited interior-shell upper bound is

\[
 F^{\rm ub}=F_{\ge52}^{\rm ub}+\sum_{w\in\{38,40,42,44,46,48,50\}}U_wc_w.
\]

Substituting U_38=3,827,351,840,403 and the existing exact LP caps for the
other six weights yields -log2(F^ub)=36.8687886, as a display approximation.
The exact rational comparison gives F^ub>2^-40. The same comparison fails
even if the weight-38 term is removed: the remaining margin is 36.9955341.

The following contributions are expressed as fractions of the entire 2^-40
budget. A sum below one would close this upper-bound check.

| Weight | Contribution divided by 2^-40, approximate |
|---:|---:|
| 38, at the tested cap | 0.736901 |
| 40, at the LP cap | 7.540808 |
| 42, at the LP cap | 0.394471 |
| 44, at the LP cap | 0.069597 |
| 46, at the LP cap | 0.001919 |
| 48, at the LP cap | 0.000052 |
| 50, at the LP cap | 0.000002 |

These are contributions of upper bounds, not estimates of the actual shell
multiplicities or evidence that the code fails its distance target.
Weight 40 is the largest remaining obstruction in this bound. At the tested
weight-38 cap, eliminating weight 40 alone would still not suffice; another
improvement, for example at weight 42 or weight 38, would also be needed.

Reproduce the audit with:

    python code/audit_endpoint_application_budget.py

The script requires python-flint. It converts the already-directed transfer
logarithms into upper dyadic coefficients using Arb, and then aggregates exact
rational numbers. It stores the coefficients, terms, comparisons, and source
hashes in `generated/endpoint_application_budget_audit.json`.

The endpoint experiment remains valid as a test of U_38. It uses a fixed
sample count and a fixed shell threshold, so this application correction does
not alter its sampling law or acceptance event. Its randomness qualification
and its limitation to a shell bound remain explicit in the retained certificate.

The next application step is to improve and jointly allocate the other low
shells, starting with weight 40. One exact conditional budget calculation uses

\[
 A_{40}(C)\le50{,}000{,}000{,}000{,}000,\qquad
 A_{42}(C)\le5{,}000{,}000{,}000{,}000{,}000.
\]

With these two additional assumptions, the tested weight-38 cap, and the
existing caps at weights 44 and above, the audited interior-shell sum has
40.0476889 bits. The script checks its exact rational upper bound against
2^-40. The two displayed multiplicity bounds are unproved targets, not
newly certified counts. This calculation concerns the existing interior-shell
functional; it does not replace an audit of every application interface class.

The final RM2Sub functional and occupations
beyond one remain separate work.
