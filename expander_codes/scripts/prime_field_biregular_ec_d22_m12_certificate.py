#!/usr/bin/env python3
"""Generate and verify the degree-22, memory-12 prime-field EC result."""

from __future__ import annotations

import argparse
import hashlib
import inspect
import json
import math
from pathlib import Path
import tempfile
from typing import Any

import flint
from flint import arb, ctx

import prime_field_biregular_ec_certificate as ec_certificate

from ea_certificate import decimal_marker
from prime_field_biregular_ec_certificate import (
    VerificationResult,
    balanced_occupancy_prefix_arb,
    full_support_markers_float,
    full_support_term_arb,
    result_json,
    singleton_block_term_arb,
    singleton_exact_band_branch_term_arb,
    singleton_exact_band_term_arb,
)
from prime_field_biregular_ec_diagnostic import (
    candidate_parameters,
    singleton_constraint_balanced_ec_block_logbound,
)
from prime_field_ea_diagnostic import gv_distance


SCHEMA = "prime-field-biregular-ec-low-equation-regions-v1"
CACHE_SCHEMA = "prime-field-biregular-ec-arb-cache-v1"
PRIME = 2**127 - 1
DEGREE = 22
MEMORY = 12
SPARSE_LIMIT = 48
K, N, REGION_LENGTH = candidate_parameters(2**21, DEGREE)
CUTOFF = math.floor(N * gv_distance(PRIME, 0.5))
DEFAULT_CACHE_DIR = (
    Path(__file__).resolve().parent.parent / ".verification_cache" / "d22_m12"
)

SPARSE_FIELD_MARKERS = {
    (1, 4): ("0.9999117386619436", "3.881725035646193e-06"),
    (5, 8): ("0.9998223421378374", "7.761798365941198e-06"),
    (9, 16): ("0.9996439837979334", "1.5548455485401082e-05"),
    (17, 24): ("0.9994660272970728", "2.329672930248755e-05"),
    (25, 32): ("0.9992874014485273", "3.094617295393988e-05"),
    (33, 40): ("0.9991092828116812", "3.8663752108617314e-05"),
    (41, 48): ("0.9989307965037183", "4.647847858811414e-05"),
}

EXACT_MARKERS = {
    (49, 64): (("0.9985761812702978", "6.105589040113493e-05"), ("0.9985739367907277", "6.193856697792455e-05")),
    (65, 96): (("0.9979279810228457", "8.913018126321766e-05"), ("0.997927983742769", "9.017949320586518e-05")),
    (97, 128): (("0.9971851395762397", "0.00012122647675903623"), ("0.9971852098862992", "0.00012210346510064464")),
    (129, 160): (("0.9964626686281732", "0.00015281564588009807"), ("0.9964626716983358", "0.0001534967723690676")),
    (161, 192): (("0.9957552036787093", "0.00018327364729992147"), ("0.9957552039268005", "0.0001842029376419928")),
    (193, 224): (("0.995059652549349", "0.00021383822868757733"), ("0.9950596500480957", "0.0002125994066517015")),
    (225, 256): (("0.9943739715495167", "0.00024293406601854996"), ("0.9943739914908493", "0.00024377712804258457")),
    (257, 288): (("0.9936967768736145", "0.00027229418877859954"), ("0.9936967809584294", "0.0002731827810501416")),
    (289, 320): (("0.9930269602488256", "0.0003012547424649579"), ("0.9930269626843563", "0.00030216683915570443")),
    (321, 360): (("0.9922826352724878", "0.0003334247055087766"), ("0.9922826469121759", "0.000334297609900208")),
    (361, 400): (("0.9914625325563562", "0.000368676031701799"), ("0.9914625334040166", "0.0003696323787173457")),
    (401, 440): (("0.9906505730497129", "0.0004036543704911642"), ("0.990650578819153", "0.00040463567729111553")),
    (441, 480): (("0.9898459688755932", "0.0004383456321527603"), ("0.989845968301158", "0.0004391993477064778")),
    (481, 520): (("0.9890480502362656", "0.00047275431366246897"), ("0.9890480448824667", "0.0004735525738242994")),
    (521, 560): (("0.9882562626928676", "0.0005065894414190316"), ("0.9882562553046488", "0.0005074777788324528")),
    (561, 600): (("0.9874701685119442", "0.0005408804796429879"), ("0.9874701619807598", "0.0005412812274527123")),
    (601, 640): (("0.986682695486827", "0.0005740102973814832"), ("0.9866834822383735", "0.0005751729202397466")),
    (641, 800): (("0.9847846202039392", "0.0006555372409042515"), ("0.984784619612031", "0.0006565587891781301")),
    (801, 1000): (("0.9813803317895631", "0.0008011121090061921"), ("0.9813803351693822", "0.0008019818838110974")),
    (1001, 1200): (("0.9776649699926743", "0.0009594005700555795"), ("0.9776649589463212", "0.0009606280870954773")),
    (1201, 1400): (("0.9735824674275545", "0.0011333753949033225"), ("0.9735673128001461", "0.001134382616081565")),
}


def support_blocks(k: int) -> tuple[tuple[int, int], ...]:
    blocks = [
        (1401, 2000), (2001, 3200), (3201, 6400), (6401, 12800),
        (12801, 25600), (25601, 51200), (51201, 102400),
        (102401, 204800), (204801, 409600), (409601, 819200),
    ]
    tail_widths = [2**power for power in range(16, 2, -1)] + [3]
    start = blocks[-1][1] + 1
    first_width = (k - 1) - start + 1 - sum(tail_widths)
    if first_width <= 0:
        raise ValueError("support range is too short for the fixed partition")
    blocks.append((start, start + first_width - 1))
    start = blocks[-1][1] + 1
    for width in tail_widths:
        blocks.append((start, start + width - 1))
        start += width
    if start != k:
        raise AssertionError("support partition does not end at k - 1")
    return tuple(blocks)


BLOCKS = support_blocks(K)


def sparse_structural_term_arb(
    *, k: int, n: int, cutoff: int, region_count: int, memory: int,
    support_limit: int,
) -> arb:
    """Bound low-equation traces by regions without singleton groups."""
    region_length = n // region_count
    group_size = k // region_length
    laws = balanced_occupancy_prefix_arb(
        region_length=region_length,
        group_size=group_size,
        limit=support_limit,
    )
    total = arb(0)
    for r in range(1, support_limit + 1):
        # A region without a singleton group occupies at most floor(r/2)
        # groups.  The occupancy tail is therefore an upper bound on the
        # probability that this region is bad.
        bad_probability = sum(laws[r][:r // 2 + 1], arb(0))
        zero_intervals = 1 + (r - 1) // memory
        needed_regions = (
            (n - cutoff - (r - 1) + region_length - 1) // region_length
            - 2 * zero_intervals
        )
        if not 1 <= needed_regions <= region_count:
            raise ValueError(f"sparse region bound is vacuous at r={r}")
        total += (
            arb(math.comb(k, r))
            * math.comb(region_count, needed_regions)
            * bad_probability ** needed_regions
        )
    return total


def _canonical_digest(value: Any) -> str:
    encoded = json.dumps(
        value, sort_keys=True, separators=(",", ":"), ensure_ascii=True
    ).encode("utf-8")
    return hashlib.sha256(encoded).hexdigest()


def _mathematical_verifier_digest() -> str:
    """Fingerprint code that determines cached interval values."""
    digest = hashlib.sha256()
    digest.update(Path(ec_certificate.__file__).read_bytes())
    digest.update(inspect.getsource(sparse_structural_term_arb).encode("utf-8"))
    digest.update(flint.__version__.encode("ascii"))
    return digest.hexdigest()


def _verification_result_from_json(value: dict[str, Any]) -> VerificationResult:
    largest_exact = value["largest_exact_band"]
    largest_block = value["largest_block"]
    return VerificationResult(
        success=bool(value["success"]),
        total_bound=arb(value["total_bound"]),
        security_bits=arb(value["security_bits"]),
        exact_bound=arb(value["exact_bound"]),
        block_bound=arb(value["block_bound"]),
        full_support_bound=arb(value["full_support_bound"]),
        largest_exact_band=(
            int(largest_exact["lo"]), int(largest_exact["hi"])
        ),
        largest_exact_term=arb(largest_exact["bound"]),
        largest_block=(int(largest_block["lo"]), int(largest_block["hi"])),
        largest_block_term=arb(largest_block["bound"]),
    )


class VerificationCache:
    """Advisory content-addressed cache for completed Arb calculations."""

    def __init__(
        self, *, certificate: dict[str, Any], directory: Path,
        refresh: bool = False,
    ) -> None:
        self.directory = directory
        self.refresh = refresh
        self.base_key = {
            "schema": CACHE_SCHEMA,
            "certificate_sha256": _canonical_digest(certificate),
            "precision_bits": int(certificate["verification"]["precision_bits"]),
            "verifier_sha256": _mathematical_verifier_digest(),
        }

    def _payload(self, stage: str) -> dict[str, Any]:
        return {**self.base_key, "stage": stage}

    def _path(self, stage: str) -> Path:
        key = _canonical_digest(self._payload(stage))
        return self.directory / f"{key}.json"

    def _read(self, stage: str) -> dict[str, Any] | None:
        if self.refresh:
            return None
        path = self._path(stage)
        if not path.is_file():
            return None
        try:
            entry = json.loads(path.read_text(encoding="utf-8"))
        except (OSError, json.JSONDecodeError):
            return None
        if entry.get("key") != self._payload(stage):
            return None
        return entry

    def _write(self, stage: str, body: dict[str, Any]) -> None:
        self.directory.mkdir(parents=True, exist_ok=True)
        path = self._path(stage)
        entry = {"key": self._payload(stage), **body}
        with tempfile.NamedTemporaryFile(
            mode="w", encoding="utf-8", newline="\n", delete=False,
            dir=self.directory, prefix=f"{path.stem}.", suffix=".tmp",
        ) as stream:
            json.dump(entry, stream, indent=2)
            stream.write("\n")
            temporary = Path(stream.name)
        temporary.replace(path)

    def get_bound(self, stage: str) -> arb | None:
        entry = self._read(stage)
        if entry is None or "bound" not in entry:
            return None
        return arb(entry["bound"])

    def put_bound(self, stage: str, bound: arb) -> None:
        self._write(stage, {"bound": str(bound)})

    def get_result(self) -> VerificationResult | None:
        entry = self._read("complete")
        if entry is None or "result" not in entry:
            return None
        return _verification_result_from_json(entry["result"])

    def put_result(self, result: VerificationResult) -> None:
        self._write("complete", {"result": result_json(result)})


def _cached_bound(
    cache: VerificationCache | None, stage: str, compute,
) -> arb:
    if cache is not None:
        value = cache.get_bound(stage)
        if value is not None:
            print(f"cache hit: {stage}", flush=True)
            return value
    value = compute()
    if cache is not None:
        cache.put_bound(stage, value)
    return value


def validate_coverage(certificate: dict[str, Any]) -> None:
    k = int(certificate["parameters"]["k"])
    if int(certificate["sparse_support"]["hi"]) != SPARSE_LIMIT:
        raise ValueError("unexpected sparse-support limit")
    expected = SPARSE_LIMIT + 1
    for section in ("exact_bands", "blocks"):
        for item in certificate[section]:
            lo, hi = int(item["lo"]), int(item["hi"])
            if lo != expected or hi < lo:
                raise ValueError(f"coverage breaks at r={expected}")
            expected = hi + 1
    if expected != k or int(certificate["full_support"]["r"]) != k:
        raise ValueError("support coverage does not end at full support")


def generate_certificate(*, target_bits: int, precision_bits: int) -> dict[str, Any]:
    exact_bands = [{
        "lo": lo,
        "hi": hi,
        "markers": {"structural": list(markers[0]), "field": list(markers[1])},
    } for (lo, hi), markers in EXACT_MARKERS.items()]
    blocks = []
    for lo, hi in BLOCKS:
        bound = singleton_constraint_balanced_ec_block_logbound(
            prime=PRIME, k=K, n=N, cutoff=CUTOFF, region_count=DEGREE,
            memory=MEMORY, support_start=lo, support_limit=hi,
        )
        blocks.append({
            "lo": lo,
            "hi": hi,
            "markers": {
                "structural": [decimal_marker(x) for x in bound.structural_markers],
                "field": [decimal_marker(x) for x in bound.field_markers],
            },
            "diagnostic_log2": bound.total_log2,
        })
        print(f"optimized block [{lo},{hi}]: {bound.total_log2:.6f}", flush=True)
    return {
        "schema": SCHEMA,
        "parameters": {
            "prime": str(PRIME), "k": K, "n": N, "cutoff": CUTOFF,
            "region_count": DEGREE, "memory": MEMORY,
            "target_bits": target_bits,
        },
        "verification": {"precision_bits": precision_bits},
        "sparse_support": {
            "lo": 1,
            "hi": SPARSE_LIMIT,
            "field_bands": [{"lo": lo, "hi": hi, "markers": list(values)}
                            for (lo, hi), values in SPARSE_FIELD_MARKERS.items()],
        },
        "exact_bands": exact_bands,
        "blocks": blocks,
        "full_support": {
            "r": K,
            "markers": full_support_markers_float(
                prime=PRIME, n=N, cutoff=CUTOFF, memory=MEMORY,
                message_weight=K,
            ),
        },
    }


def verify_certificate(
    certificate: dict[str, Any], *, cache: VerificationCache | None = None,
) -> VerificationResult:
    if certificate.get("schema") != SCHEMA:
        raise ValueError("unsupported certificate schema")
    validate_coverage(certificate)
    parameters = certificate["parameters"]
    prime = int(parameters["prime"])
    k, n = int(parameters["k"]), int(parameters["n"])
    cutoff = int(parameters["cutoff"])
    region_count = int(parameters["region_count"])
    memory = int(parameters["memory"])
    target_bits = int(parameters["target_bits"])
    ctx.prec = int(certificate["verification"]["precision_bits"])
    if cache is not None:
        cached_result = cache.get_result()
        if cached_result is not None:
            print(
                "cache hit: complete certificate "
                "(use --no-cache for a full recomputation)",
                flush=True,
            )
            return cached_result

    sparse = certificate["sparse_support"]
    sparse_structural = _cached_bound(
        cache,
        f"sparse-structural-1-{int(sparse['hi'])}",
        lambda: sparse_structural_term_arb(
            k=k, n=n, cutoff=cutoff, region_count=region_count, memory=memory,
            support_limit=int(sparse["hi"]),
        ),
    )
    sparse_field = arb(0)
    for item in sparse["field_bands"]:
        lo, hi = int(item["lo"]), int(item["hi"])
        sparse_field += _cached_bound(
            cache,
            f"sparse-field-{lo}-{hi}",
            lambda item=item, lo=lo, hi=hi: (
                singleton_exact_band_branch_term_arb(
                    prime=prime, k=k, n=n, cutoff=cutoff,
                    region_count=region_count, memory=memory, lo=lo, hi=hi,
                    values=item["markers"], field=True,
                )
            ),
        )
        print(f"verified sparse field band [{lo},{hi}]", flush=True)
    exact_total = sparse_structural + sparse_field
    largest_exact_band, largest_exact_term = (1, SPARSE_LIMIT), exact_total
    for item in certificate["exact_bands"]:
        lo, hi = int(item["lo"]), int(item["hi"])
        term = _cached_bound(
            cache,
            f"exact-{lo}-{hi}",
            lambda item=item, lo=lo, hi=hi: singleton_exact_band_term_arb(
                prime=prime, k=k, n=n, cutoff=cutoff,
                region_count=region_count, memory=memory, lo=lo, hi=hi,
                markers=item["markers"],
            ),
        )
        exact_total += term
        if term > largest_exact_term:
            largest_exact_band, largest_exact_term = (lo, hi), term
        print(f"verified exact band [{lo},{hi}]", flush=True)

    block_total = arb(0)
    largest_block, largest_block_term = (0, 0), arb(0)
    for item in certificate["blocks"]:
        lo, hi = int(item["lo"]), int(item["hi"])
        term = _cached_bound(
            cache,
            f"block-{lo}-{hi}",
            lambda item=item, lo=lo, hi=hi: singleton_block_term_arb(
                prime=prime, k=k, n=n, cutoff=cutoff,
                region_count=region_count, memory=memory, lo=lo, hi=hi,
                markers=item["markers"],
            ),
        )
        block_total += term
        if largest_block == (0, 0) or term > largest_block_term:
            largest_block, largest_block_term = (lo, hi), term
        print(f"verified block [{lo},{hi}]", flush=True)

    full = certificate["full_support"]
    full_term = _cached_bound(
        cache,
        f"full-{k}",
        lambda: full_support_term_arb(
            prime=prime, n=n, cutoff=cutoff, memory=memory,
            message_weight=k, markers=full["markers"],
        ),
    )
    total = exact_total + block_total + full_term
    result = VerificationResult(
        success=bool(total < arb(2) ** (-target_bits)), total_bound=total,
        security_bits=-total.log() / arb(2).log(), exact_bound=exact_total,
        block_bound=block_total, full_support_bound=full_term,
        largest_exact_band=largest_exact_band,
        largest_exact_term=largest_exact_term,
        largest_block=largest_block, largest_block_term=largest_block_term,
    )
    if cache is not None:
        cache.put_result(result)
    return result


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest="command", required=True)
    generate = commands.add_parser("generate")
    generate.add_argument("--output", type=Path, required=True)
    generate.add_argument("--target-bits", type=int, default=20)
    generate.add_argument("--precision-bits", type=int, default=256)
    verify = commands.add_parser("verify")
    verify.add_argument("certificate", type=Path)
    verify.add_argument("--cache-dir", type=Path, default=DEFAULT_CACHE_DIR)
    cache_mode = verify.add_mutually_exclusive_group()
    cache_mode.add_argument("--no-cache", action="store_true")
    cache_mode.add_argument("--refresh-cache", action="store_true")
    seed = commands.add_parser(
        "seed-cache",
        help="cache the output of a previously completed full verification",
    )
    seed.add_argument("certificate", type=Path)
    seed.add_argument("result", type=Path)
    seed.add_argument("--cache-dir", type=Path, default=DEFAULT_CACHE_DIR)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    if args.command == "generate":
        certificate = generate_certificate(
            target_bits=args.target_bits, precision_bits=args.precision_bits
        )
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(json.dumps(certificate, indent=2) + "\n", encoding="utf-8")
        print(f"wrote {args.output}")
        return
    certificate = json.loads(args.certificate.read_text(encoding="utf-8"))
    if args.command == "seed-cache":
        validate_coverage(certificate)
        ctx.prec = int(certificate["verification"]["precision_bits"])
        result = _verification_result_from_json(
            json.loads(args.result.read_text(encoding="utf-8"))
        )
        target_bits = int(certificate["parameters"]["target_bits"])
        if not result.success or not result.total_bound < arb(2) ** (-target_bits):
            raise ValueError("the supplied result does not meet the target")
        cache = VerificationCache(
            certificate=certificate, directory=args.cache_dir
        )
        cache.put_result(result)
        print(
            "seeded advisory complete-result cache; "
            "verify --no-cache recomputes every interval",
            flush=True,
        )
        return
    cache = None if args.no_cache else VerificationCache(
        certificate=certificate,
        directory=args.cache_dir,
        refresh=args.refresh_cache,
    )
    result = verify_certificate(certificate, cache=cache)
    print(json.dumps(result_json(result), indent=2))
    if not result.success:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
