# Full-Split Interval Recompute Audit

This artifact records an opt-in regeneration pass for the heavier interval families in the
finite RM/EBCH full-split certificate. The default verifier re-sums committed artifacts and
manifest constants; this report records a deeper run of the command-reproducible intervals.

Tolerance: `5e-06` bits.
Rows recomputed: `45`.
Overall status: `PASS`.

A row passes when its recomputed log2 contribution is at most the stored
upper bound plus the verifier tolerance. Negative deltas are safer than
the stored value and are not treated as failures.

## Family Summary

| family | rows | stored total log2 | recomputed total log2 | max upward delta bits | min delta bits |
| --- | ---: | ---: | ---: | ---: | ---: |
| `complement_high` | 7 | -35150.097366 | -35150.097366 | 0 | 0 |
| `early_accelerated` | 10 | -1056.187188 | -1056.187188 | 0 | 0 |
| `early_postprefix` | 6 | -380.829185 | -380.829185 | 0 | 0 |
| `high` | 15 | -586.597103 | -586.597103 | 0 | -1.48 |
| `late_postprefix` | 7 | -545.690526 | -545.690526 | 0 | 0 |

## Rows

| family | interval | stored log2 | recomputed log2 | delta bits | status |
| --- | --- | ---: | ---: | ---: | --- |
| `high` | `2001--7858` | -586.597103 | -586.597103 | 0 | `PASS` |
| `high` | `7859--20550` | -2655.944828 | -2655.944829 | -1e-06 | `PASS` |
| `high` | `20551--75000` | -6025.574047 | -6025.574048 | -1e-06 | `PASS` |
| `high` | `75001--250000` | -12844.868837 | -12844.868837 | 0 | `PASS` |
| `high` | `250001--350000` | -28476.640541 | -28477.517697 | -0.877 | `PASS` |
| `high` | `350001--400000` | -81555.594538 | -81555.594539 | -1e-06 | `PASS` |
| `high` | `400001--450000` | -91421.772986 | -91422.398395 | -0.625 | `PASS` |
| `high` | `450001--550000` | -82081.961328 | -82081.961329 | -1e-06 | `PASS` |
| `high` | `550001--650000` | -1518.493461 | -1518.493462 | -1e-06 | `PASS` |
| `high` | `650001--725000` | -116328.133946 | -116328.133947 | -1e-06 | `PASS` |
| `high` | `725001--750000` | -169887.902921 | -169887.902922 | -1e-06 | `PASS` |
| `high` | `750001--850000` | -134523.583646 | -134523.583646 | -2.91e-11 | `PASS` |
| `high` | `850001--950000` | -117534.843530 | -117535.750802 | -0.907 | `PASS` |
| `high` | `950001--1050000` | -285052.550668 | -285054.035442 | -1.48 | `PASS` |
| `high` | `1050001--1148736` | -562860.716638 | -562862.201412 | -1.48 | `PASS` |
| `complement_high` | `1048577--1148736` | -128915.315746 | -128915.315746 | 0 | `PASS` |
| `complement_high` | `1148737--1300000` | -124889.651755 | -124889.651755 | 0 | `PASS` |
| `complement_high` | `1300001--1500000` | -113710.357653 | -113710.357653 | 0 | `PASS` |
| `complement_high` | `1500001--1700000` | -35150.097366 | -35150.097366 | 0 | `PASS` |
| `complement_high` | `1700001--1900000` | -524628.366199 | -524628.366199 | 0 | `PASS` |
| `complement_high` | `1900001--2097089` | -894140.350713 | -894140.350713 | 0 | `PASS` |
| `complement_high` | `2097090--2097152` | -893792.285010 | -893792.285010 | 0 | `PASS` |
| `early_postprefix` | `501--2000` | -380.829185 | -380.829185 | 0 | `PASS` |
| `early_postprefix` | `2001--7858` | -852.998772 | -852.998772 | 0 | `PASS` |
| `early_postprefix` | `7859--20550` | -3712.365087 | -3712.365087 | 0 | `PASS` |
| `early_postprefix` | `20551--30000` | -8614.382763 | -8614.382763 | 0 | `PASS` |
| `early_postprefix` | `30001--50000` | -21332.131033 | -21332.131033 | 0 | `PASS` |
| `early_postprefix` | `50001--75000` | -34847.929043 | -34847.929043 | 0 | `PASS` |
| `early_accelerated` | `75001--90000` | -28574.869802 | -28574.869802 | 0 | `PASS` |
| `early_accelerated` | `90001--100000` | -22149.918612 | -22149.918612 | 0 | `PASS` |
| `early_accelerated` | `100001--110000` | -1056.187188 | -1056.187188 | 0 | `PASS` |
| `early_accelerated` | `110001--125000` | -10362.025949 | -10362.025949 | 0 | `PASS` |
| `early_accelerated` | `125001--160000` | -45589.980049 | -45589.980049 | 0 | `PASS` |
| `early_accelerated` | `160001--250000` | -60293.751685 | -60293.751685 | 0 | `PASS` |
| `early_accelerated` | `250001--350000` | -94682.265141 | -94682.265141 | 0 | `PASS` |
| `early_accelerated` | `350001--524288` | -181215.853042 | -181215.853042 | 0 | `PASS` |
| `early_accelerated` | `524289--750000` | -194565.849984 | -194565.849984 | 0 | `PASS` |
| `early_accelerated` | `750001--1048576` | -98386.449256 | -98386.449256 | 0 | `PASS` |
| `late_postprefix` | `501--2000` | -545.690526 | -545.690526 | 0 | `PASS` |
| `late_postprefix` | `2001--7858` | -2237.593414 | -2237.593414 | 0 | `PASS` |
| `late_postprefix` | `7859--20550` | -8919.160723 | -8919.160723 | 0 | `PASS` |
| `late_postprefix` | `20551--75000` | -23772.761808 | -23772.761808 | 0 | `PASS` |
| `late_postprefix` | `75001--150000` | -93856.032305 | -93856.032305 | 0 | `PASS` |
| `late_postprefix` | `150001--250000` | -210536.528296 | -210536.528296 | 0 | `PASS` |
| `late_postprefix` | `250001--380736` | -417970.896784 | -417970.896784 | 0 | `PASS` |

## Command

```powershell
python scripts\verify_fullsplit_finite_ledger.py --check-high-intervals all --check-complement-high-intervals all --check-early-postprefix-intervals all --check-early-accelerated-intervals all --check-late-postprefix-intervals all --write-interval-recompute-audit-md scripts\fullsplit_interval_recompute_audit.md
```
