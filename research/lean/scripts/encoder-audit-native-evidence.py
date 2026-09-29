"""Read-only evidence coverage audit for the final project's import closure."""
import collections
from datetime import datetime, timezone
import hashlib
import importlib.util
import json
from pathlib import Path
import re
import time

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "scripts/map_data"
TARGET = "SpinCodes.Structured.ConcreteNativeTheoremPin"
spec = importlib.util.spec_from_file_location("replay", ROOT / "scripts/replay-native-closure.py")
replay = importlib.util.module_from_spec(spec)
spec.loader.exec_module(replay)


def sha(path):
    with path.open("rb") as f:
        return hashlib.file_digest(f, "sha256").hexdigest()


def relpath(module):
    return module.replace(".", "/") + ".lean"


def discover():
    graph, pending = {}, [TARGET]
    while pending:
        module = pending.pop()
        if module in graph:
            continue
        path = ROOT / relpath(module)
        if not path.exists():
            graph[module] = None
            continue
        text = replay.uncomment(path.read_text(encoding="utf-8-sig"))
        imports = []
        for line in text.splitlines():
            line = line.strip()
            if not line or line == "prelude":
                continue
            match = re.fullmatch(r"(?:(?:public|meta)\s+)*import\s+(.+)", line)
            if not match:
                break
            imports += [m for m in match[1].split() if m.startswith("SpinCodes.") or m == "SpinCodes"]
        graph[module] = imports
        pending += imports
    return graph


def main():
    started = datetime.now(timezone.utc).isoformat()
    graph = discover()
    current, source_index = {}, collections.defaultdict(list)
    for module in sorted(graph):
        src = ROOT / relpath(module)
        obj = (ROOT / ".lake/build/lib/lean" / relpath(module)).with_suffix(".olean")
        source_hash = sha(src) if src.exists() else None
        object_hash = sha(obj) if obj.exists() else None
        current[module] = {"source_sha256": source_hash, "object_sha256": object_hash}
        if source_hash:
            source_index[source_hash].append(module)
    evidence, source_only = collections.defaultdict(list), collections.defaultdict(list)
    record_files, unresolved = [], []

    def resolve(label, source_hash):
        if isinstance(label, str):
            label = label.replace("\\", "/")
            if label.startswith("SpinCodes/") and label.endswith(".lean"):
                return label[:-5].replace("/", ".")
            if label.startswith("SpinCodes."):
                return label
            if re.fullmatch(r"[A-Za-z_][A-Za-z_0-9]*", label):
                candidate = "SpinCodes.Structured." + label
                if candidate in graph or (ROOT / relpath(candidate)).exists():
                    return candidate
        matches = source_index.get(str(source_hash).lower(), [])
        return matches[0] if len(matches) == 1 else None

    def add_pair(label, sh, oh, origin, location, kind):
        module = resolve(label, sh)
        if module in current:
            evidence[module].append({"record": origin, "location": location, "kind": kind,
                                     "source_sha256": sh.lower(), "object_sha256": oh.lower()})
        elif not module:
            unresolved.append({"record": origin, "location": location, "label": label, "source_sha256": sh})

    def walk(value, origin, location="$", accepted=False, batch=False, label=None):
        if isinstance(value, list):
            for i, child in enumerate(value):
                walk(child, origin, f"{location}[{i}]", accepted, batch, label)
            return
        if not isinstance(value, dict):
            return
        status = value.get("status")
        if status is not None:
            accepted = status == "PASS"
        # These Fourier controllers append boxes only after both compiles return
        # zero and all four theorem axiom audits pass. Retain this distinct class.
        if (origin.startswith("scripts/map_data/dense_fourier_box_") and
                isinstance(value.get("index"), int) and value.get("audits") == 4 and "modules" in value):
            batch = True
        explicit = next((value[k] for k in ["source", "path", "module", "name"]
                         if isinstance(value.get(k), str)), label)
        sh = value.get("source_sha256")
        oh = next((value[k] for k in ["object_sha256", "olean_sha256", "output_sha256"]
                   if isinstance(value.get(k), str)), None)
        if origin.endswith("dense_occupation_fixed_rates_verification.json") and isinstance(value.get("box"), str):
            explicit = "DenseOccupationFixed" + value["box"] + "Rate"
        kind = "pass_record_pair" if accepted else "completed_fourier_batch_pair"
        if origin.endswith("/native_theorem_verification.json") and location.startswith("$.modules."):
            module = resolve(explicit, sh)
            if module not in {"SpinCodes.Structured.ConcreteNativeTheorem", "SpinCodes.Structured.ConcreteNativeTheoremPin"}:
                kind = "final_assembly_dependency_snapshot"
        if isinstance(sh, str) and (accepted or batch):
            if oh:
                add_pair(explicit, sh, oh, origin, location, kind)
            else:
                module = resolve(explicit, sh)
                if module in graph:
                    item = {"record": origin, "location": location, "source_sha256": sh.lower()}
                    if value.get("fresh_object_bytes_match_original") is False:
                        item.update(fresh_object_sha256=value.get("fresh_object_sha256"),
                                    original_object_sha256=value.get("original_object_sha256"),
                                    fresh_object_bytes_match_original=False)
                    source_only[module].append(item)
        if isinstance(sh, dict) and accepted:
            for path, digest in sh.items():
                module = resolve(path, digest)
                if module in graph:
                    source_only[module].append({"record": origin, "location": location + ".source_sha256",
                                                "source_sha256": digest.lower()})
        for key, child in value.items():
            if isinstance(child, (dict, list)) and key != "source_sha256":
                walk(child, origin, location + "." + key, accepted, batch, key)

    paths = list(DATA.glob("*.json"))
    paths += list((ROOT / "scripts/majorant_data").glob("*.json"))
    paths += list((ROOT / "scripts/sparse_data").glob("*.json"))
    for path in sorted(paths):
        # Derived inventories must not become new compiler evidence merely
        # because their overall consistency check says PASS.
        if path.name.startswith("encoder_native_evidence") or path.name == "final_closure_evidence.json":
            continue
        for attempt in range(3):
            try:
                raw = path.read_bytes()
                data = json.loads(raw.decode("utf-8-sig"))
                break
            except (json.JSONDecodeError, PermissionError):
                if attempt == 2:
                    data = None
                time.sleep(0.05)
        if data is None:
            continue
        origin = path.relative_to(ROOT).as_posix()
        record_files.append({"path": origin, "sha256": hashlib.sha256(raw).hexdigest(),
                             "status": data.get("status") if isinstance(data, dict) else None})
        walk(data, origin)
        # The initial dependency archive and mutable transfer cache are snapshots,
        # not independent per-module successful compilation reports.
        if path.name in {"peach_dependencies.json", "peach_extra_files.json"}:
            entries = ({r["path"]: r["sha256"] for r in data["files"]}
                       if path.name == "peach_dependencies.json" else data)
            for rel, digest in entries.items():
                if rel.startswith("SpinCodes/") and rel.endswith(".lean"):
                    obj = ".lake/build/lib/lean/" + rel[:-5] + ".olean"
                    if obj in entries:
                        kind = "initial_peach_dependency_snapshot" if path.name == "peach_dependencies.json" else "mutable_transfer_cache_snapshot"
                        add_pair(rel, digest, entries[obj], origin, "$.files" if "files" in data else "$", kind)

    rows = []
    for module, hashes in current.items():
        # A batch may finish after the initial file scan but before its report
        # was read. Recheck missing objects only when a completed pair appeared.
        if hashes["object_sha256"] is None and evidence[module]:
            obj = (ROOT / ".lake/build/lib/lean" / relpath(module)).with_suffix(".olean")
            if obj.exists():
                hashes["object_sha256"] = sha(obj)
        matches, stale = [], []
        for record in evidence[module]:
            same_source = record["source_sha256"] == hashes["source_sha256"]
            same_object = record["object_sha256"] == hashes["object_sha256"]
            if same_source and same_object:
                matches.append(record)
            else:
                stale.append(dict(record, source_matches=same_source, object_matches=same_object))
        sources = [r for r in source_only[module] if r["source_sha256"] == hashes["source_sha256"]]
        stale_sources = [r for r in source_only[module] if r["source_sha256"] != hashes["source_sha256"]]
        kinds = {r["kind"] for r in matches}
        expected_numeric = bool(re.fullmatch(r"SpinCodes\.Structured\.Dense(?:FourierExact|OccupationFixed|ScalarExact)[BW]\d+(?:Data|Rate)?", module))
        expected_assembly = module.rsplit(".", 1)[-1] in {
            "ConcreteNativeTheorem", "ConcreteNativeTheoremPin", "DenseOccupationAllCertified",
            "DenseOccupationScalarCertified", "DenseOccupationFourierCertified", "DenseOccupationMixedCertified"}
        if not hashes["source_sha256"]:
            category = "expected_missing_fourier_source" if "DenseFourierExactB" in module else "unexpected_missing_source"
        elif "pass_record_pair" in kinds:
            category = "pass_pair_match"
        elif "completed_fourier_batch_pair" in kinds:
            category = "completed_batch_pair_match"
        elif "initial_peach_dependency_snapshot" in kinds:
            category = "initial_dependency_snapshot_match"
        elif any(r.get("fresh_object_bytes_match_original") is False for r in sources):
            category = "fresh_source_check_different_object"
        elif sources:
            category = "source_only_pass_record"
        elif expected_numeric or expected_assembly:
            category = "expected_pending_numeric_or_assembly"
        elif "final_assembly_dependency_snapshot" in kinds:
            category = "final_assembly_snapshot_only"
        elif "mutable_transfer_cache_snapshot" in kinds:
            category = "transfer_cache_only"
        else:
            category = "evidence_gap"
        rows.append({"module": module, "category": category, **hashes,
                     "matching_pairs": matches, "matching_source_only_records": sources,
                     "nonmatching_historical_pairs": stale,
                     "nonmatching_source_only_records": stale_sources})
    counts = collections.Counter(r["category"] for r in rows)
    result = {"status": "REVIEW_SNAPSHOT", "started_utc": started,
              "finished_utc": datetime.now(timezone.utc).isoformat(), "target": TARGET,
              "scope": "Read-only current file hashes versus existing recorded evidence; no compilation, no timestamp inference",
              "source_module_nodes": len(graph), "counts": dict(counts), "records": record_files,
              "modules": rows, "unresolved_record_labels": unresolved,
              "caveat": "Generation/checking may proceed concurrently; rerun after final assembly. Snapshot/cache matches are weaker than compile records."}
    out = DATA / "encoder_native_evidence_coverage.json"
    out.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"nodes": len(graph), "counts": dict(counts), "unresolved_labels": len(unresolved)}, indent=2))
    print(out)


if __name__ == "__main__":
    main()
