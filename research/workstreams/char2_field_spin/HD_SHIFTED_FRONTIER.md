# General support profiles: FLOAT frontier

This note records isolated stress tests for the alternating adjacent-pair
route. Each test uses the exact joint support count and the improved
two-step bound from `HD_SHIFTED_GOOD_PAIRS.md`. None establishes coverage
between the displayed profiles.

The driver `theory/hd_shifted_frontier_screen.py` stores reusable witnesses.
It first searches the inexpensive global coefficient bound, then evaluates
the good-pair alias bound. The `--replay` option evaluates a saved witness
without repeating the search. Output remains explicitly nondirected.

## Initial stress tests

The first eight profiles sample small and medium weights, several levels
of pair fullness, and unequal parity weights. Each global search used
180 iterations. All displayed alias evaluations use grid 64.

| (E,O,a,b) | Active groups | log2 upper screen |
| --- | ---: | ---: |
| (128,128,100,100) | 256 | -591.117968 |
| (256,256,200,200) | 512 | -801.963707 |
| (512,512,400,400) | 1024 | -291.716788 |
| (512,512,480,480) | 1024 | -738.573170 |
| (512,512,500,400) | 1024 | -546.743895 |
| (600,424,400,380) | 1024 | +257.431557 |
| (1024,1024,800,800) | 2048 | -2028.525510 |
| (1100,948,900,850) | 2048 | -1004.710303 |

The unequal-parity profile (600,424,400,380) is the only positive result
in this initial sample. Its two histograms are

    H0 = (3472,200,24,400),
    H1 = (3452,44,220,380).

It is a candidate for further analysis, not evidence that the construction
fails. The initial witness was optimized for the ordinary global bound,
and a finite alias grid adds positive terms. Refining that witness and
increasing the grid can distinguish these losses from a limitation of
the regional bound.

The initial witness had

    z     = 0.4993504868994289,
    alpha = 0.09004029815447658,
    even marker log2 values
          = (-3.7895746104396224,-7.066683441897507,-2.474115241015154),
    odd marker log2 values
          = (-6.190557247192267,-3.68421922993356,-2.5402735481047323).

An additional 450 global iterations found a substantially better witness.
Its good-pair bound is -392.510113931 at grid 64 and -557.366455156
at grid 128. The saved witness is

    z     = 0.48346682697675003,
    alpha = 0.10354771628792783,
    even marker log2 values
          = (-3.6670491924654645,-6.711178743827144,-2.304473696143282),
    odd marker log2 values
          = (-5.828374644650113,-3.521695518880656,-2.37029381976133).

Thus the initial positive result did not require a stronger operator.
The bounded follow-up batch samples substantially larger parity imbalances
and lower pair fullness at weight 1024. Each new global search uses
320 iterations, followed by good-pair evaluation on grid 64.

| (E,O,a,b) | log2 upper screen |
| --- | ---: |
| (768,256,200,200) | -3212.071387 |
| (896,128,100,100) | -3963.500476 |
| (512,512,256,256) | -3078.295368 |
| (512,512,128,128) | -7664.714932 |
| (512,512,450,128) | -2357.718847 |

All five follow-up profiles have negative bounds. After refining the
initial unequal-parity witness, all thirteen sampled profiles have
grid-64 bounds below -40. The weakest displayed bound is -291.716788
at (512,512,400,400). The driver stores all thirteen witnesses,
including the refined witness for (600,424,400,380).

No unresolved class remains in this bounded sample. This conclusion
does not imply a bound for unsampled classes. The next turn should
build a disjoint adaptive cover, using failed boxes to select further
profiles. If point probes are needed first, the sampled trend suggests
checking balanced weight-1024 profiles with a,b between 400 and 480;
that suggestion is an interpolation heuristic, not a proven location
of the worst class.

For coefficient validation, run
`python -m unittest hd_shifted_good_pair_screen_tests -v` from the
theory directory. The six tests are described in
`HD_SHIFTED_GOOD_PAIRS.md`. The support-count module has seven
additional tests in `hd_shifted_general_count_tests.py`.
