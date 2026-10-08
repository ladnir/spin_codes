# Fresh feedback without state refresh

This proposal has an exact small transfer for each fixed message, but its
current whole-code first-moment screen does not close at K18 or K20.
It is parked as a negative screen. No implementation or timing result is
claimed, and the retained GL construction is unchanged.

Let F=GF(16), K=8L, and N=16L. Use the same GRS[16,8] outer,
sixteen independent uniform regional permutations, and independent uniform
nonzero input diagonals. Process the region-major stream in blocks of sixteen.
For state dimension s, define A_j=(1,j,...,j^(s-1)) for all j in F.
Initialize Q=0 and discard the final state.

Independently sample every entry of each s-by-16 matrix C_i uniformly from F,
including zero. Setup fixes these matrices for every message. Let X in F^16
be the current routed block after its diagonal. At step i use

    Y=X+A Q,       Q'=Q+C_i X.

The intended cheap variant uses s=4. It removes the separate state refresh
and replaces the fixed feedback map with fresh unrestricted feedback.
The input/output map remains invertible because its direct block is the identity.

## Exact state-weight transfer

Fix a base-field source word before sampling setup. Conditional on its routed
support, the diagonal makes all nonzero entries of X independent uniform
elements of F minus zero. For every fixed nonzero X, C_i X is uniform in F^s.
Consequently Q' is uniform even after conditioning on Q, X, and Y.
For X=0, the exact state persists.

Retain v=wt(AQ), with categories

    V={0,17-s,...,16}.

The s=4 variant has five categories: 0,13,14,15,16. Let p_v be the probability
that a uniform state has category v. Then p_0=16^(-s). For v>0,

    p_v = binom(16,v) n_v / 16^s,
    n_v = sum_(i=0)^v (-1)^i binom(v,i) 16^max(0,v-(16-s)-i).

Here n_v counts MDS words with one prescribed support of size v.

Fix z in (0,1] and input weight w. Conditional on the step weights in a
uniform regional support, each step mask is uniform among its w-subsets.
For a fixed incoming state of category v, let a count active input lanes
where AQ vanishes. Define f(z)=z+(1-z)/15. Averaging over the mask and diagonal
gives the exact output moment

    G_(v,w)(z) = 1/binom(16,w)
      * sum_a binom(16-v,a) binom(v,w-a)
              z^(v-w+2a) f(z)^(w-a),

where max(0,w-v)<=a<=min(w,16-v).
The marked transfer, with rows and columns indexed by V, is

    T_0[v,v'] = 1[v=v'] z^v,
    T_w[v,v'] = G_(v,w)(z) p_(v')       for w>0.

Thus inactive blocks preserve the weight category. Every active block resets
its distribution independently of the current emission. In particular, an
active block resets to zero with probability 16^(-s).

For a region containing h active positions, its averaged transfer is

    R_h(z) = [u^h] (sum_(w=0)^16 binom(16,w) T_w(z) u^w)^(L/16)
             / binom(L,h).

The source-counting argument from the [GL certificate](GF16_PROOF.md) applies
to these larger matrices. The script uses exact outer multiplicities through four
active rows and a positive multiplicity envelope above four.
An optional source tilt rho replaces n_w<=16^(w-8) by

    n_w <= eta(rho) 16^(w-8) rho^w,
    eta(rho) = max_(9<=w<=16) n_w / (16^(w-8) rho^w).

## Floating screen

The following values are margins for individual active-row classes at 10%
distance, using s=4. They are not directed certificates or whole-code margins.

| Size | Active outer rows | Untilted source | Source rho=0.871 |
|---|---:|---:|---:|
| K16 | 24 | 44.79 bits | Not evaluated |
| K18 | 16 | 14.33 bits | 17.98 bits |
| K18 | 32 | -19.88 bits | -11.18 bits |
| K18 | 64 | -48.41 bits | -31.23 bits |
| K20 | 16 | -17.87 bits | -14.20 bits |
| K20 | 32 | -84.78 bits | -75.98 bits |
| K20 | 64 | -180.06 bits | -162.52 bits |

A negative margin means the computed first-moment bound exceeds one.
It does not prove that the construction has a bad codeword.
These screens give no basis for replacing the retained construction to avoid
larger state at K18/K20. An untilted s=5 check also fails at K18 with 64
active rows, giving -5.07 bits.

Exhaustive GF4 checks enumerate feedback rows, every incoming state, and every
input support and value. They verify the reset law and averaged emission.
The script also checks the unmarked GF16 transition row sums.

From the theory directory:

```sh
python -B reset_field_screen.py --check
python -B reset_field_screen.py --k-log 16 --rows 24
python -B reset_field_screen.py --k-log 18 20 --rows 16 32 64
python -B reset_field_screen.py --k-log 18 20 --rows 16 32 64 --source-rho 0.871
```

The four-state arithmetic count is comparable to the original direct kernel:
fresh feedback replaces its fixed feedback and state-matrix products.
Stored embedded feedback coefficients grow from 16 to 64 bytes per step,
and coefficient preparation also changes. Neither observation establishes
a performance gain. The next proof audit remains the GL construction's K18
six-state certificate.
