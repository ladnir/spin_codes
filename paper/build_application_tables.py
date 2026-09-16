"""Render compact, pinned application measurements; never runs benchmarks."""
import argparse
import hashlib
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent
DATA = ROOT / "data" / "application_results.json"
BOLT_DATA = ROOT / "data" / "bolt_opening_projection.json"
LIGERITO_DATA = ROOT / "data" / "ligerito_standalone.json"
PATHS = {
    "encoding": "results/spin-brakedown/peach-paired/summary.json",
    "pcs": "results/spin-brakedown/peach-security/summary.json",
    "bolt": "results/bolt-one-thread/summary.json",
    "flock": "results/flock-spin/native-current-comparison/summary.json",
}


def import_results(source):
    records, receipts = {}, {}
    for key, name in PATHS.items():
        raw = (source / name).read_bytes()
        records[key] = json.loads(raw)
        receipts[key] = {"path": name, "sha256": hashlib.sha256(raw).hexdigest()}
    enc = records["encoding"]
    return {
        "source_repository": "hypercat",
        "source_revision": "015a4f8",
        "flock_revision": "2d667ca",
        "receipts": receipts,
        "encoding": {
            "method": enc["method"],
            "source_revision": enc["spin_revision"],
            "runs": {k: v for k, v in enc["runs"].items()
                     if k.startswith("fused-") or k.startswith("confirm-fused-")},
        },
        "pcs": {
            "scope": records["pcs"]["scope"],
            "budget": records["pcs"]["budget"],
            "rows": {k: records["pcs"]["stages"]["pooled_optimized"][k]
                     for k in ("k16-w2", "k18-w2")},
        },
        "bolt": records["bolt"],
        "flock": {
            log: {backend: records["flock"][log]["backends"][backend]["mean_process_medians"]
                  for backend in ("spin", "native-fast100")}
            for log in ("14", "16")
        },
    }


def import_imt(data, source):
    raw = source.read_bytes()
    result = json.loads(raw)
    assert result["verified_proofs"] == 168
    assert set(result["pcs"]["rows"]) == {"k16-w2", "k18-w2"}
    data["source_revision"] = "f98b02a"
    data["flock_revision"] = "e15d047"
    data["historical_receipts"] = data.get("historical_receipts", data["receipts"])
    data["receipts"] = {
        "imt": {"path": "results/imt-20260916/final-summary.json",
                "sha256": hashlib.sha256(raw).hexdigest(),
                "run_receipt_sha256": result["receipt_sha256"]},
        "bolt": data["historical_receipts"]["bolt"],
    }
    data["encoding"] = {"method": result["method"]["encoding_method"],
                        "source_revision": "f98b02a", "prefix": "imt",
                        "runs": result["encoding"]}
    data["pcs"].update(result["pcs"])
    data["flock"] = result["flock"]
    data["measurement_method"] = result["method"]
    return data


def import_bolt(source):
    name = "tools/standalone/opening-results/projection.json"
    raw = (source / name).read_bytes()
    projection = json.loads(raw)
    case = next(c for c in projection["cases"] if c["input_mib"] == 512)
    fields = ("input_mib", "phase_medians_ms",
              "one_opening_amortized_limit_projection_ms",
              "one_opening_non_amortized_projection_ms",
              "mulperm_leading_work_projection_ms")
    return {
        "source_repository": "bolt-rs",
        "source_revision": "9b45089",
        "receipt": {"path": name, "sha256": hashlib.sha256(raw).hexdigest()},
        **{key: projection[key] for key in
           ("label", "scope", "model", "amortization", "omitted", "proxy_scope")},
        "case": {key: case[key] for key in fields},
        "flock_cases": [
            {key: c[key] for key in
             ("input_mib", "phase_medians_ms", "commit_measured_ms",
              "two_openings", "flock")}
            for c in projection["cases"] if c["input_mib"] in (32, 128)
        ],
    }


def flock_projection(data, case):
    log = {32: "14", 128: "16"}[case["input_mib"]]
    spin = data["flock"][log]["spin"]
    outer = spin["total_ms"] - spin["commit_ms"] - spin["open_ms"]
    p = case["phase_medians_ms"]
    shared = (2 * p["row_evaluations"] + p["random_column_fold"]
              + 2 * p["one_sumcheck"] + p["inner_message_proxy"]
              + p["inner_syndrome_proxy"])
    assert abs(shared - case["two_openings"]["shared_amortized_limit_ms"]) < 1e-8
    total = outer + case["commit_measured_ms"] + shared
    # The calibrated Bolt components are retained; the surrounding Flock
    # residual comes from the current SPIN measurement campaign.
    return total


def render(data, bolt, ligerito):
    rows = []
    for field, label in (("forward_ms", "Ordinary"), ("transpose_ms", "Transposed")):
        prefix = data["encoding"].get("prefix", "fused")
        values = [data["encoding"]["runs"][f"{prefix}-k{k}"][field]["median"]
                  for k in (16, 18, 20)]
        rows.append(label + " & " + " & ".join(f"{v:.3f}" for v in values) + r" \\")
    tables = {"ordinary_encoding": rows}
    rows = []
    for key, label in (("k16-w2", "square"), ("k18-w2", "longer rows")):
        r = data["pcs"]["rows"][key]
        rows.append(f"SPIN--Brakedown ({label}) & {r['total_ms']:.0f} & "
                    f"{r['verify_ms']:.0f} & {r['proof_bytes']/2**20:.1f}" + r" \\")
    for key, label in (("fast100", "Fast100"), ("slim100", "Slim100")):
        r = ligerito["rows"][key]
        rows.append(f"Ligerito ({label}) & {r['total_ms']:.0f} & "
                    f"{r['verify_ms']:.0f} & {r['opening_bytes']/2**20:.1f}" + r" \\")
    tables["pcs_standalone"] = rows
    case = bolt["case"]
    rows = []
    for key, label in (("k16-w2", "square"), ("k18-w2", "longer rows")):
        value = data["pcs"]["rows"][key]["open_ms"]
        rows.append(f"SPIN--Brakedown ({label}) & Measured & {value:.0f}" + r" \\")
    for key, label in (("fast100", "Fast100"), ("slim100", "Slim100")):
        value = ligerito["rows"][key]["open_ms"]
        rows.append(f"Ligerito ({label}) & Measured & {value:.0f}" + r" \\")
    for label, key in (("amortized limit", "one_opening_amortized_limit_projection_ms"),
                       ("non-amortized", "one_opening_non_amortized_projection_ms")):
        rows.append(f"Bolt-max ({label}) & Calibrated projection & {case[key]:.0f}" + r" \\")
    tables["pcs_opening"] = rows
    rows = []
    cases = {c["input_mib"]: c for c in bolt["flock_cases"]}
    for log in ("14", "16"):
        for backend, label in (("spin", "SPIN--Brakedown"), ("native-fast100", "Ligerito")):
            r = data["flock"][log][backend]
            rows.append(f"${2**int(log):,}$ & {label} & {r['total_ms']:.0f} & "
                        f"{r['verify_ms']:.0f} & {r['proof_bytes']/2**20:.1f}" + r" \\")
        case = cases[{"14": 32, "16": 128}[log]]
        rows.append(f"${2**int(log):,}$ & Bolt (projection) & "
                    f"{flock_projection(data, case):.0f} & --- & ---" + r" \\")
    tables["flock_comparison"] = rows
    digest = hashlib.sha256(DATA.read_bytes()).hexdigest()
    bolt_digest = hashlib.sha256(BOLT_DATA.read_bytes()).hexdigest()
    ligerito_digest = hashlib.sha256(LIGERITO_DATA.read_bytes()).hexdigest()
    outputs = {ROOT / "tables" / (name + ".tex"):
            f"% Generated by build_application_tables.py; data SHA-256 {digest}\n"
            + (f"% Bolt projection SHA-256 {bolt_digest}\n"
               if name in ("pcs_opening", "flock_comparison") else "")
            + (f"% Ligerito data SHA-256 {ligerito_digest}\n"
               if name in ("pcs_standalone", "pcs_opening") else "")
            + "\n".join(rows) + "\n\\bottomrule\n" for name, rows in tables.items()}
    return outputs


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", type=Path, help="Import the pinned summaries from Hypercat")
    parser.add_argument("--imt-summary", type=Path, help="Refresh SPIN and Flock from the final IMT campaign")
    parser.add_argument("--bolt-source", type=Path, help="Import the pinned Bolt opening projection")
    parser.add_argument("--ligerito-source", type=Path, help="Import the standalone Ligerito summary JSON")
    parser.add_argument("--check", action="store_true", help="Check generated tables without writing")
    args = parser.parse_args()
    if (args.source or args.imt_summary or args.bolt_source or args.ligerito_source) and args.check:
        parser.error("Import and --check are separate operations")
    if args.source:
        DATA.parent.mkdir(parents=True, exist_ok=True)
        DATA.write_text(json.dumps(import_results(args.source), indent=2) + "\n", encoding="utf-8")
    data = json.loads(DATA.read_text(encoding="utf-8"))
    if args.imt_summary:
        data = import_imt(data, args.imt_summary)
        DATA.write_text(json.dumps(data, indent=2) + "\n", encoding="utf-8")
    if args.bolt_source:
        BOLT_DATA.write_text(json.dumps(import_bolt(args.bolt_source), indent=2) + "\n", encoding="utf-8")
    bolt = json.loads(BOLT_DATA.read_text(encoding="utf-8"))
    if args.ligerito_source:
        raw = args.ligerito_source.read_bytes()
        record = json.loads(raw)
        record["source_summary_sha256"] = hashlib.sha256(raw).hexdigest()
        LIGERITO_DATA.write_text(json.dumps(record, indent=2) + "\n", encoding="utf-8")
    ligerito = json.loads(LIGERITO_DATA.read_text(encoding="utf-8"))
    assert ligerito["input_bytes"] == 2**29
    assert all(r["trials"] == 10 for r in ligerito["rows"].values())
    case = bolt["case"]
    phases = case["phase_medians_ms"]
    expected = (phases["row_evaluations"] + phases["random_column_fold"]
                + 2 * phases["one_sumcheck"] + phases["inner_message_proxy"]
                + phases["inner_syndrome_proxy"])
    assert abs(expected - case["one_opening_amortized_limit_projection_ms"]) < 1e-8
    assert abs(expected + case["mulperm_leading_work_projection_ms"]
               - case["one_opening_non_amortized_projection_ms"]) < 1e-8
    for row in data["pcs"]["budget"]:
        a, n, t = row["e"] + 1, row["n"], row["queries"]
        assert 3 * row["e"] < row["assumed_distance"]
        assert (a * n**t + (n-a)**t * 2**128) * 2**100 <= 2**128 * n**t
    row = data["pcs"]["budget"][0]
    a, n, t = row["e"] + 1, row["n"], 2110
    assert 2 * (a * n**t + (n-a)**t * 2**128) * 2**102 < 2**128 * n**t
    for path, content in render(data, bolt, ligerito).items():
        if args.check:
            assert path.read_text(encoding="utf-8") == content, f"Stale table: {path}"
        else:
            path.parent.mkdir(parents=True, exist_ok=True)
            path.write_text(content, encoding="utf-8")
    print("Application tables and exact conditional parameter bounds checked.")


if __name__ == "__main__":
    main()
