"""Summarize serial process medians, including paired ratios to fixed encoding."""
import csv
import pathlib
import statistics
import sys

root = pathlib.Path(sys.argv[1])
groups = {}
for path in sorted(root.glob("*.csv")):
    trial = int(path.stem.rsplit("-", 1)[1])
    with path.open() as stream:
        for row in csv.DictReader(stream):
            groups.setdefault(row["mode"], {})[trial] = row
fixed = groups["fixed"]
base = statistics.median(float(row["total_ms"]) for row in fixed.values())
print("mode,processes,setup_ms,encode_ms,total_ms,ratio_of_medians,median_paired_ratio,min_paired_ratio,max_paired_ratio")
for mode, trials in groups.items():
    if len({row["checksum"] for row in trials.values()}) != 1:
        raise ValueError(f"checksum changed for {mode}")
    values = [statistics.median(float(row[key]) for row in trials.values())
              for key in ("setup_ms", "encode_ms", "total_ms")]
    paired = [float(row["total_ms"]) / float(fixed[t]["total_ms"])
              for t, row in trials.items() if t in fixed]
    print(f"{mode},{len(trials)}," + ",".join(f"{v:.6f}" for v in values +
          [values[2] / base, statistics.median(paired), min(paired), max(paired)]))
