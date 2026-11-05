"""Bounded FLOAT probe at seven retained output-marker anchors.

Each source size uses only its own retained z, three nearby alpha values,
four origin penalties, and the five existing regional coefficient tilts.
The retained markers are proposals, not evidence that an old certificate
applies to the new whole-state kernel or paired route. All arithmetic here
is a diagnostic; representative profiles do not cover a whole source size.

This file imports the current bridge helpers without changing them. It
writes no files. Run with Python -B; the root owns every numerical run.
"""

import argparse
from fractions import Fraction
import json
import math
from pathlib import Path
import shutil
import time

import screen_small_state_bridge as bridge


ANCHORS = {
    7: ("9634151/9689600", "5470045/5470051"),
    16: ("2238127/2263291", "5162271/5451841"),
    64: ("1607513/1676799", "2062633/9650309"),
    256: ("6945437/8121969", "1012405/6913782"),
    1024: ("2825730/5340629", "14853/101062"),
    4096: ("370391/5732639", "734878/2992721"),
    8192: ("53362/9139873", "2384018/4197155"),
}
ALPHA_SCALES = (Fraction(4, 5), Fraction(1), Fraction(6, 5))
PENALTIES = (0., 8., 16., 32.)


def alpha_variants(anchor):
    value = Fraction(anchor)
    if not 0 < value <= 1:
        raise ValueError("Invalid retained alpha anchor")
    return tuple(sorted({min(Fraction(1), scale * value) for scale in ALPHA_SCALES}))


def probe(groups=tuple(ANCHORS), fractional="packet", seconds=60., memory_mib=1024,
          source_mode="cutoff"):
    if (not groups or len(set(groups)) != len(groups)
            or any(g not in ANCHORS for g in groups)
            or fractional not in ("packet", "region")
            or source_mode not in ("cutoff", "joint")
            or not 0 < seconds <= 60 or not 32 <= memory_mib <= 1024):
        raise ValueError("Use retained group anchors, at most 60 seconds and 1024 MiB")
    if shutil.disk_usage(Path(__file__).resolve().parent).free < 20 * 1024 ** 3:
        raise RuntimeError("At least 20 GiB free space is required")
    started = time.monotonic()
    deadline = started + seconds

    def guard():
        if time.monotonic() >= deadline:
            raise TimeoutError("Anchor probe reached its time cap")
        if bridge._resident_bytes() > memory_mib * 1024 ** 2:
            raise MemoryError("Anchor probe reached its memory cap")

    checks = bridge.selfcheck()
    best = {}
    evaluated = 0
    finished_groups = []
    status = "complete-selected-anchor-grid"
    try:
        for g in groups:
            guard()
            z = Fraction(ANCHORS[g][0])
            source = bridge.source_factor(g, source_mode)
            # Joint mode uses eta=1/256, not a claim that a cutoff retained
            # every auxiliary source. No retained fraction is reported there.
            penalty_fraction = (1. if source_mode == "joint" else
                                float(Fraction(source["retained_fraction"])))
            profiles = bridge.selected_profiles(g)
            full_counts = sorted({f for profile in profiles for f in profile})
            table = bridge.local_table(z)
            guard()
            for exact_alpha in alpha_variants(ANCHORS[g][1]):
                alpha = float(exact_alpha)
                prefix = (bridge.log2_binomial(bridge.GROUPS, g)
                          + alpha * (source["log2_prefactor"]
                                     - bridge.LIMIT * math.log2(float(z))))
                for penalty in PENALTIES:
                    guard()
                    weighted = bridge.penalized_table(table, penalty, penalty_fraction)
                    regions = {f: bridge.regional_log_matrix(
                        weighted, g, f, alpha, fractional=fractional, tilt=True)
                        for f in full_counts}
                    for f0, f1 in profiles:
                        pair = bridge.log_matrix_product(regions[f0], regions[f1])
                        run = bridge.log_matrix_power(pair, 128)
                        operator = bridge.log2_sum(run[0])
                        value = prefix + operator
                        evaluated += 1
                        key = (g, f0, f1)
                        if key not in best or value < best[key]["log2_uncapped_bound"]:
                            best[key] = dict(source_mode=source_mode, actual_groups=g,
                                pair_profile=[f0, f1],
                                z=str(z), alpha=str(exact_alpha), alpha_float=alpha,
                                origin_log2_penalty=penalty, log2_uncapped_bound=value,
                                log2_operator=operator,
                                log2_prefactor_and_source_choices=prefix)
                    guard()
            finished_groups.append(g)
    except (TimeoutError, MemoryError) as error:
        status = str(error)

    rows = []
    for g in groups:
        selected = [row for key, row in best.items() if key[0] == g]
        if not selected:
            continue
        low = min(selected, key=lambda row: row["log2_uncapped_bound"])
        high = max(selected, key=lambda row: row["log2_uncapped_bound"])
        rows.append(dict(source_mode=source_mode, actual_groups=g,
            requested_profiles=len(bridge.selected_profiles(g)),
            evaluated_profiles=len(selected), completed_anchor_grid=g in finished_groups,
            best_profile_log2=low["log2_uncapped_bound"],
            worst_profile_log2=high["log2_uncapped_bound"],
            best_profile_witness=low, worst_profile_witness=high,
            retained_anchor=dict(z=ANCHORS[g][0], alpha=ANCHORS[g][1]),
            source_factor=bridge.source_factor(g, source_mode), fractional_grouping=fractional,
            directed=False, whole_code_certificate=False,
            scope="best marker per selected actual pair profile; no source-size cover"))
    summary = dict(schema="small-state-bridge-anchor-probe-v2-striped", status=status,
        source_mode=source_mode,
        source_groups=list(groups), completed_groups=finished_groups,
        profile_marker_evaluations=evaluated, returned_profiles=len(best),
        requested_profiles=sum(len(bridge.selected_profiles(g)) for g in groups),
        local_marker_tables=len({key[0] for key in best}),
        alpha_scales=list(map(str, ALPHA_SCALES)), origin_log2_penalties=list(PENALTIES),
        iid_tilt_exponents=[0., .5, 1., 1.5, 2.], fractional_grouping=fractional,
        seconds=time.monotonic() - started, peak_resident_bytes=bridge._resident_bytes(),
        output_limit=bridge.LIMIT, selfchecks=checks,
        directed=False, whole_code_certificate=False, generated_file_bytes=0)
    return rows, summary


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--groups", type=int, nargs="+", choices=tuple(ANCHORS),
                        default=list(ANCHORS))
    parser.add_argument("--fractional", choices=("packet", "region"), default="packet")
    parser.add_argument("--seconds", type=float, default=60.)
    parser.add_argument("--memory-mib", type=int, default=1024)
    parser.add_argument("--source-mode", choices=("cutoff", "joint"), default="cutoff",
                        help="Use the missing-group cutoff or joint-Jensen source bound")
    args = parser.parse_args()
    rows, summary = probe(tuple(args.groups), args.fractional, args.seconds, args.memory_mib,
                          args.source_mode)
    for row in rows:
        print(json.dumps(row, allow_nan=False), flush=True)
    print(json.dumps(summary, allow_nan=False), flush=True)


if __name__ == "__main__":
    main()
