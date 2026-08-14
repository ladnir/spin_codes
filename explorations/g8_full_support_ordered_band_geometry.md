# Full-support ordered-band diagnostic

The full-support profile simplex has nine positive integer coordinates with
sum `M`. General affine BSP cuts can accumulate many unrelated constraints.
That growth makes exact vertex enumeration expensive.

Lane A instead chooses one class order from competing affine witnesses. It
fixes one prefix of that order as a band `B`. Dominance data selects only the
order, prefix width, and integer thresholds. Every proof cell has the form

```text
L <= sum(a[j] for j in B) <= U,  a[j] >= 1,  sum(a[j]) = M.
```

The integer thresholds are copied from supplementary target profiles. No
binary64 value is rounded into a proof threshold. The BSP stores each cut in
the verifier's support-relative primitive form, including the required child
swap when class 8 belongs to the band.

Let `s=9`, `q=|B|`, and `r=s-q`. A slab has at most

```text
s + 2*q*r <= 49
```

vertices. An internal boundary contributes exactly `q*r` candidates. Endpoint
slabs also retain the applicable simplex corners. All candidates are integral:
at boundary mass `m`, one band coordinate equals `m-q+1`, one complement
coordinate equals `M-m-r+1`, and every other coordinate equals one.

The number of positive profiles in a slab is exactly

```text
sum(m=L..U) C(m-1,q-1) * C(M-m-1,r-1).
```

One scan over feasible `m` records prefix totals at every node endpoint. Thus
all node counts cost `O(M)` binomial-integer operations and `O(nodes)` extra
big integers. Child counts are exact, disjoint, and sum to the parent count.

This decomposition is deliberately one-dimensional. It keeps proof geometry
small, but it cannot localize residuals transverse to the chosen band mass.
The producer emits only `UNRESOLVED` leaves and makes no coverage claim.

Static and synthetic checks:

```powershell
python scripts/produce_packet_group_g8_full_support_ordered_band_slabs.py --self-test
python -m py_compile scripts/produce_packet_group_g8_full_support_ordered_band_slabs.py
```

First bounded root-run command:

```powershell
python scripts/produce_packet_group_g8_full_support_ordered_band_slabs.py --manifest G8_SUPPORT_MANIFEST.json --atlas out/g8_support_seed_atlas.json --supplementary-atlas out/g8_highdim_supplementary_final32.json --max-boundaries 15 --hardening-reserve-bits 8 --output out/g8_full_support_ordered_band_diagnostic.json
```
