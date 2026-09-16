"""Audit the BCH-derived [256,64] construction; import a sourced parent spectrum.

No benchmarks, network access, or enumeration of the full code are performed.
Without --parent-spectrum, only the construction and published [256,63] table
are audited. A [256,64] spectrum is never fabricated from distance bounds.
"""
import argparse
import hashlib
import json
from pathlib import Path
import sys

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
BCH_PATH = ROOT / "bch_spectrum_work/bch_spectrum_codex_bundle/code"
sys.path.insert(0, str(BCH_PATH))
import bch_quotient as bch


def require(condition, message):
    if not condition:
        raise ValueError(message)


def read_spectrum(path, n=256):
    values = [0] * (n + 1)
    seen = set()
    for line in Path(path).read_text(encoding="utf-8-sig").splitlines():
        line = line.split("#", 1)[0].strip()
        if not line:
            continue
        w, count = map(int, line.split())
        require(0 <= w <= n and w not in seen and count >= 0,
                "invalid weight, duplicate weight, or negative coefficient")
        seen.add(w)
        values[w] = count
    return values


def macwilliams(values, dimension):
    """Exact integer Krawtchouk recurrence, O(n * support size)."""
    n = len(values) - 1
    totals = [0] * (n + 1)
    for w, count in enumerate(values):
        if not count:
            continue
        previous, current = 1, n - 2*w
        totals[0] += count
        if n:
            totals[1] += count * current
        for j in range(1, n):
            numerator = (n-2*w)*current - (n-j+1)*previous
            require(numerator % (j+1) == 0, "nonintegral Krawtchouk value")
            following = numerator // (j+1)
            totals[j+1] += count * following
            previous, current = current, following
    mass = 1 << dimension
    require(all(v >= 0 and v % mass == 0 for v in totals),
            "MacWilliams transform has negative or nonintegral coefficients")
    return [v // mass for v in totals]


def audit_spectrum(values, dimension, distance_lower_bound, even=True):
    n = len(values)-1
    require(all(type(v) is int and v >= 0 for v in values), "noninteger coefficient")
    require(values[0] == 1 and sum(values) == 1 << dimension, "wrong spectrum mass")
    require(values == values[::-1], "missing complement symmetry")
    require(not even or not any(values[1::2]), "odd weight in parity extension")
    distance = next(w for w in range(1, n+1) if values[w])
    require(distance >= distance_lower_bound, "distance below BCH guarantee")
    dual = macwilliams(values, dimension)
    require(dual[0] == 1 and sum(dual) == 1 << (n-dimension), "wrong dual mass")
    return dict(dimension=dimension, length=n, minimum_distance=distance,
                nonzero_coefficients=sum(bool(v) for v in values),
                mass=str(sum(values)), macwilliams_checked=True)


def rank(rows):
    pivots = {}
    for row in rows:
        while row:
            top = row.bit_length()-1
            if top not in pivots:
                pivots[top] = row
                break
            row ^= pivots[top]
    return len(pivots)


def extend(word):
    return word | ((word.bit_count() & 1) << 255)


def construction():
    p = bch.generator_polynomial(59)
    q = bch.generator_polynomial(61)
    h, remainder = bch.binary_poly_divmod(q, p)
    require(remainder == 0 and p.bit_length()-1 == 184 and q.bit_length()-1 == 192,
            "unexpected BCH dimensions or nesting")
    require(h == 0x187, "unexpected quotient polynomial")
    # Multiplication by X visits every nonzero residue modulo h. This proves
    # the quotient is a field and the cyclic action is transitive, directly.
    orbit, x = set(), 1
    while x not in orbit:
        orbit.add(x)
        x = bch.multiply_by_x_mod(x, h)
    require(orbit == set(range(1, 256)) and x == 1, "quotient not transitive")
    require(bch.binary_poly_divmod((1 << 255) | 1, q)[1] == 0,
            "subcode generator does not divide X^255+1")
    small = [extend(q << i) for i in range(63)]
    rows = small + [extend(p)]
    require(rank(small) == 63 and rank(rows) == 64, "generator rank mismatch")
    require(rank(rows + [(1 << 256)-1]) == 64, "all-one word missing")
    for row in rows:
        require(row.bit_count() % 2 == 0 and row.bit_length() <= 256,
                "invalid parity extension")
        require(bch.binary_poly_divmod(row & ((1 << 255)-1), p)[1] == 0,
                "selected generator not in parent")
    # Independent action check: cycle the actual length-255 word p and reduce
    # its message polynomial modulo h, rather than just iterating residues.
    word, actual = p, set()
    for _ in range(255):
        message, rem = bch.binary_poly_divmod(word, p)
        require(rem == 0, "cyclic shift left parent")
        actual.add(bch.binary_poly_divmod(message, h)[1])
        word = ((word << 1) & ((1 << 255)-1)) | (word >> 254)
    require(actual == orbit and word == p, "actual coordinate action mismatch")
    return dict(length=256, dimension=64, distance_lower_bound=60,
                parent_designed_distance=59, subcode_designed_distance=61,
                parent_generator_hex=hex(p), subcode_generator_hex=hex(q),
                quotient_polynomial_hex=hex(h), nonzero_coset_orbit=255,
                generator_rows_hex=[f"{v:064x}" for v in rows])


def reconstruct(small, parent, cosets=255):
    require(len(small) == len(parent), "spectrum length mismatch")
    differences = [p-s for s, p in zip(small, parent)]
    require(all(v >= 0 and v % cosets == 0 for v in differences),
            "coset difference is negative or not divisible by orbit size")
    return [s + v//cosets for s, v in zip(small, differences)]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--parent-spectrum", type=Path)
    parser.add_argument("--parent-source", help="citation/URL and table identifier")
    parser.add_argument("--output", type=Path, help="optional JSON receipt")
    parser.add_argument("--spectrum-output", type=Path, help="optional reconstructed .wd file")
    args = parser.parse_args()
    if bool(args.parent_spectrum) != bool(args.parent_source):
        parser.error("--parent-spectrum and --parent-source must be supplied together")
    if args.spectrum_output and not args.parent_spectrum:
        parser.error("--spectrum-output requires --parent-spectrum")
    small_path = HERE / "sources/EBCH256_63.wd"
    small = read_spectrum(small_path)
    receipt = dict(schema=1, status="construction_only",
                   construction=construction(),
                   source63_sha256=hashlib.sha256(small_path.read_bytes()).hexdigest(),
                   source63_audit=audit_spectrum(small, 63, 64),
                   source_dependencies_sha256={name: hashlib.sha256((BCH_PATH/name).read_bytes()).hexdigest()
                       for name in ("bch_quotient.py", "affine_wambach.py")})
    if args.parent_spectrum:
        parent = read_spectrum(args.parent_spectrum)
        receipt["parent_audit"] = audit_spectrum(parent, 71, 60)
        selected = reconstruct(small, parent)
        receipt["selected_audit"] = audit_spectrum(selected, 64, 60)
        receipt["parent_source"] = args.parent_source
        receipt["parent_sha256"] = hashlib.sha256(args.parent_spectrum.read_bytes()).hexdigest()
        receipt["spectrum"] = {str(w): str(v) for w, v in enumerate(selected) if v}
        receipt["status"] = "reconstructed_from_supplied_parent_table"
        if args.spectrum_output:
            args.spectrum_output.write_text(
                "# BCH-derived [256,64] intermediate: C63 + span(extended parent generator).\n"
                "# W64 = W63 + (W71-W63)/255; see SPECTRUM_AUDIT.json and README.md.\n"
                + "".join(f"{w} {v}\n" for w, v in enumerate(selected) if v),
                encoding="utf-8", newline="\n")
    content = json.dumps(receipt, indent=2) + "\n"
    if args.output:
        args.output.write_text(content, encoding="utf-8")
    print(content, end="")


if __name__ == "__main__":
    main()
