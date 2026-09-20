# Current numerical tables

Collected from existing receipts; no new numerical search is performed.
Read [the engineering overview](CURRENT_ENGINEERING_RESULTS.md) for interpretation.

## Complete BCH-64/128 grid

Positive entries are full margin in bits, followed by the full-bound loss relative to Q1 in parentheses.
W gives a weak full upper bound and its negative margin. O gives a positive lower exponent a,
meaning the bad-word first moment is at least 2^a; it does not assert actual code failure.
Every original geometry appears in these tables.

Audit: 77 useful full bounds, 7 weak full bounds, 46 first-moment obstructions; no missing geometries.

### Message size at T=64, S=20

| log2 K | BCH-64: margin (loss) | BCH-128: margin (loss) |
| ---: | ---: | ---: |
| 12 | 16.639117 (0.000320032) | 37.661913 (6.12408e-08) |
| 13 | 16.473634 (0.000292576) | 38.142752 (2.66561e-10) |
| 14 | 15.924509 (0.000376143) | 37.943933 (2.26398e-10) |
| 15 | 15.153949 (0.000597125) | 37.353947 (2.87181e-10) |
| 16 | 14.271217 (0.00105968) | 36.564526 (4.53795e-10) |
| 17 | 13.331531 (0.0019915) | 35.673924 (8.03486e-10) |
| 18 | 12.360479 (0.00385805) | 34.729010 (1.51108e-09) |
| 19 | 11.371835 (0.00758885) | 33.756062 (2.93159e-09) |
| 20 | 10.371924 (0.0150234) | 32.769598 (5.77436e-09) |
| 21 | 9.361065 (0.0297824) | 31.776503 (1.14592e-08) |
| 22 | 8.333885 (0.0588738) | 30.780002 (2.28302e-08) |
| 23 | 7.278242 (0.115485) | 29.781637 (4.55713e-08) |
| 24 | 6.171264 (0.222935) | 28.782539 (9.10589e-08) |
| 25 | 4.972020 (0.422422) | 27.782972 (1.82029e-07) |
| 26 | 3.574425 (0.820138) | 26.783189 (3.63971e-07) |

### Epoch/state sweep at K=2^20

| T | S | BCH-64: margin (loss) | BCH-128: margin (loss) |
| ---: | ---: | ---: | ---: |
| 64 | 7 | O (35733.8) | O (40709.9) |
| 64 | 8 | O (6201.17) | O (11168.3) |
| 64 | 9 | W (-17824.3) | W (-5127.07) |
| 64 | 10 | 6.674265 (0.584376) | 27.041016 (2.1406e-06) |
| 64 | 11 | 8.004586 (0.264722) | 28.850828 (4.92853e-07) |
| 64 | 12 | 8.880604 (0.133217) | 30.197059 (1.53382e-07) |
| 64 | 13 | 9.464553 (0.0746147) | 31.166105 (6.03804e-08) |
| 64 | 14 | 9.842649 (0.0464641) | 31.824891 (2.90776e-08) |
| 64 | 15 | 10.077023 (0.0319564) | 32.244724 (1.66027e-08) |
| 64 | 16 | 10.215546 (0.0240858) | 32.494912 (1.09872e-08) |
| 64 | 17 | 10.293941 (0.0196824) | 32.635244 (8.25071e-09) |
| 64 | 18 | 10.336863 (0.017193) | 32.710487 (6.86003e-09) |
| 64 | 19 | 10.359841 (0.0157952) | 32.749575 (6.14293e-09) |
| 64 | 20 | 10.371924 (0.0150234) | 32.769598 (5.77436e-09) |
| 128 | 8 | O (130660) | O (135655) |
| 128 | 9 | O (114283) | O (119279) |
| 128 | 10 | O (97913.8) | O (102909) |
| 128 | 11 | O (81559.1) | O (86553.6) |
| 128 | 12 | O (65235.9) | O (70228.8) |
| 128 | 13 | O (48963.5) | O (53953.8) |
| 128 | 14 | O (32786.7) | O (37772.7) |
| 128 | 15 | O (16834.5) | O (21814.2) |
| 128 | 16 | O (1050.59) | O (6022.97) |
| 128 | 17 | W (-25469.7) | W (-12756.5) |
| 128 | 18 | W (-10331.2) | 32.699951 (6.91928e-09) |
| 128 | 19 | 10.353300 (0.0158728) | 32.739117 (6.19266e-09) |
| 128 | 20 | 10.365380 (0.0150961) | 32.759209 (5.81929e-09) |
| 256 | 9 | O (188000) | O (192996) |
| 256 | 10 | O (179808) | O (184804) |
| 256 | 11 | O (171616) | O (176612) |
| 256 | 12 | O (163424) | O (168420) |
| 256 | 13 | O (155232) | O (160228) |
| 256 | 14 | O (147040) | O (152036) |
| 256 | 15 | O (138848) | O (143844) |
| 256 | 16 | O (130656) | O (135652) |
| 256 | 17 | O (122464) | O (127460) |
| 256 | 18 | O (114272) | O (119268) |
| 256 | 19 | O (106080) | O (111076) |
| 256 | 20 | O (97887.9) | O (102884) |

### Additional message/state slices at T=64

The exponent-20 slice is included above.

| BCH block | S | log2 K=16 | log2 K=18 | log2 K=22 | log2 K=24 |
| ---: | ---: | ---: | ---: | ---: | ---: |
| 64 | 10 | 11.134982 (0.0458772) | 9.073813 (0.167723) | 3.575407 (1.69492) | W (-139.638) |
| 64 | 12 | 12.891468 (0.00972202) | 10.956394 (0.0351428) | 6.543416 (0.476092) | W (-1.64401) |
| 64 | 16 | 14.123995 (0.00169926) | 12.210955 (0.00619524) | 8.151754 (0.0936488) | 5.901222 (0.345587) |
| 128 | 10 | 30.878456 (1.59988e-07) | 29.008666 (5.54399e-07) | 25.065063 (8.47911e-06) | 23.073433 (3.34241e-05) |
| 128 | 12 | 33.963512 (1.28646e-08) | 32.150517 (4.0669e-08) | 28.208952 (6.04478e-07) | 26.213317 (2.41143e-06) |
| 128 | 16 | 36.285079 (8.80112e-10) | 34.453344 (2.8863e-09) | 30.505582 (4.34001e-08) | 28.508122 (1.73054e-07) |

## Matched family Q1 references at T=64, S=20

Every entry in this section is Q1-only, including the random setup charge where shown.
A random mean is an ensemble reference. A caps60 row is conditional on one shared spectrum event;
its parentheses give the margin after adding that event failure once. None is a full-tail claim.

| Family | B | log2 K=16 | log2 K=18 | log2 K=20 |
| --- | ---: | ---: | ---: | ---: |
| bch | 8 | 6.476549 | 4.531074 | 2.544576 |
| bch | 32 | 8.837122 | 6.986231 | 5.023266 |
| bch | 64 | 14.272276 | 12.364337 | 10.386948 |
| bch | 128 | 36.564526 | 34.729010 | 32.769598 |
| rm | 8 | 6.476549 | 4.531074 | 2.544576 |
| rm | 32 | 8.837122 | 6.986231 | 5.023266 |
| rm | 128 | 14.527139 | 12.600669 | 10.618636 |
| rm | 512 | 43.042551 | 41.169498 | 39.200621 |
| random_mean | 8 | -12.559667 | -14.555148 | -16.554050 |
| random_mean | 32 | -4.619372 | -6.606068 | -8.602842 |
| random_mean | 64 | 5.464231 | 3.490271 | 1.496587 |
| random_mean | 128 | 25.128964 | 23.182962 | 21.196112 |
| random_mean | 256 | 63.937645 | 62.057399 | 60.086657 |
| random_mean | 512 | 140.987777 | 139.279054 | 137.350609 |
| random_caps60 | 8 | -16.422839 (-16.422839) | -18.419105 (-18.419105) | -20.418197 (-20.418197) |
| random_caps60 | 32 | -20.379022 (-20.379022) | -22.368049 (-22.368049) | -24.365387 (-24.365387) |
| random_caps60 | 64 | -26.313709 (-26.313709) | -28.290947 (-28.290947) | -30.285425 (-30.285425) |
| random_caps60 | 128 | -18.522840 (-18.522840) | -20.508082 (-20.508082) | -22.504505 (-22.504505) |
| random_caps60 | 256 | 11.699936 (11.699936) | 9.741180 (9.741180) | 7.751216 (7.751216) |
| random_caps60 | 512 | 81.368857 (59.999999) | 79.535048 (59.999998) | 77.575858 (59.999993) |

## Positive full diagnostics in the earlier broad ledger

This ledger is a separate snapshot; it does not contain the newer BCH grid or the separate RM certificate.
Rounded random margins of 60 do not prove a strict greater-than-60-bit bound.

| Series | log2 K | T | S | Full margin |
| --- | ---: | ---: | ---: | ---: |
| random full-rank [512,256] reused | 16 | 64 | 12 | 60.0 |
| RM(4,9) [512,256,32] exact | 16 | 64 | 18 | 42.96946389269532 |
| random full-rank [512,256] reused | 16 | 64 | 18 | 60.0 |
| random full-rank [512,256] reused | 16 | 64 | 20 | 60.0 |

## Input fingerprints

These identify the local snapshots used for this document, not a new proof replay.

| Local input | SHA-256 |
| --- | --- |
| bch_evidence_grid_audit_v1.json | 935c87c8e5f18f197cc69ed90f1800f3a8365dc44c2efc1e90f6a01f5d0f3114 |
| bch_engineering_evidence_v8.csv | e42f3362f9550b85ab2d73367b3c0a971f91701c2b5c8460992351671402f177 |
| engineering_surfaces.json | 9542e2e052e52d2c6cc88a47baa2787e80a8160f4bc324bcb752995805afd754 |
| engineering_surfaces.csv | 120f6e868c2984d9eaf6fbaf79d3645d6b91b13f0d1859c325e03ab8f35e2cba |
| complete_grid_unions.json | f79b81fcce5130f845b03b3c7e68cc14c85a57c98649d485020b1ce5e2c8954b |
| complete_grid_unions.csv | ee6951bb7a2ea209b88a374c928485073e178fb78b6ea0e15bec8be044ca3d3a |
