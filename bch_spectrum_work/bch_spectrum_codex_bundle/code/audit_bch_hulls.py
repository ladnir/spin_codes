"""Construct the two BCH hulls by exact binary nullspaces; no enumeration."""
import json
import sys
from pathlib import Path
sys.dont_write_bytecode = True
from audit_bch_quadratic_sums import rows
from bch_quotient import binary_poly_divmod, generator_polynomial
from affine_wambach import gf_mul, gf_pow
from run_higher_endpoint_preflight import sha, write_new

ROOT = Path(__file__).resolve().parents[1]
OUTPUT = ROOT / 'generated/bch256_hull_structure.json'
MASK = (1 << 255) - 1


def echelon(words):
    pivots = {}
    for x in words:
        while x:
            j = x.bit_length() - 1
            if j not in pivots:
                pivots[j] = x
                break
            x ^= pivots[j]
    return pivots


def contains(basis, word):
    pivots = echelon(basis)
    while word:
        j = word.bit_length() - 1
        if j not in pivots:
            return False
        word ^= pivots[j]
    return True


def nullspace(equations, n):
    pivots = echelon(equations)
    out = []
    for free in range(n):
        if free in pivots:
            continue
        x = 1 << free
        for j in sorted(pivots):
            if (x & pivots[j]).bit_count() % 2:
                x ^= 1 << j
        assert all((x & a).bit_count() % 2 == 0 for a in equations)
        out.append(x)
    assert len(echelon(out)) == len(out) == n - len(pivots)
    return out


def hull(basis):
    gram = [sum(((x & y).bit_count() % 2) << j for j, y in enumerate(basis))
            for x in basis]
    result = []
    for mask in nullspace(gram, len(basis)):
        x = 0
        for j, y in enumerate(basis):
            if (mask >> j) & 1:
                x ^= y
        result.append(x)
    assert len(echelon(result)) == len(result)
    assert all(contains(basis, x) for x in result)
    assert all((x & y).bit_count() % 2 == 0 for x in result for y in basis)
    return result


def cyclic_description(basis):
    # Puncturing is injective on every even code in this audit.
    assert all(x.bit_count() % 2 == 0 for x in basis)
    punctured = [x & MASK for x in basis]
    assert len(echelon(punctured)) == len(basis)
    assert all(contains(punctured, ((x << 1) & MASK) | (x >> 254)) for x in punctured)
    g = (1 << 255) | 1
    for x in punctured:
        while x:
            g, x = x, binary_poly_divmod(g, x)[1]
    assert 255 - (g.bit_length() - 1) == len(basis)
    assert all(binary_poly_divmod(x, g)[1] == 0 for x in punctured)
    roots = []
    for e in range(255):
        a, value = gf_pow(2, e), 0
        for j in range(g.bit_length() - 1, -1, -1):
            value = gf_mul(value, a) ^ ((g >> j) & 1)
        if value == 0:
            roots.append(e)
    assert len(roots) == g.bit_length() - 1
    rootset = set(roots)
    runs = []
    for start in roots:
        length = 0
        while (start + length) % 255 in rootset:
            length += 1
        runs.append((length, start))
    length, start = max(runs)
    bch_bound = length + 1
    even_bound = bch_bound + bch_bound % 2
    return dict(dimension=len(basis), generator_hex=hex(g), roots=roots,
                consecutive_root_start=start, consecutive_root_length=length,
                punctured_BCH_distance_lower=bch_bound,
                extended_even_distance_lower=even_bound,
                narrow_sense_matches=[d for d in range(2, 81) if generator_polynomial(d) == g])


def build():
    p, q = rows(37, 131), rows(39, 123)
    hp, hq = hull(p), hull(q)
    dp, dq = nullspace(hp, 256), nullspace(hq, 256)
    codes = dict(P=p, Q=q, HP=hp, HQ=hq, HPdual=dp, HQdual=dq,
                 L71=rows(59, 71), U187=rows(19, 187))
    assert len(hp) == 85 and len(hq) == 93
    for h in (hp, hq):
        assert all(x.bit_count() % 4 == 0 for x in h)
        assert all((x & y).bit_count() % 2 == 0 for x in h for y in h)
        assert contains(h, (1 << 256) - 1)
    details = {}
    for label in ('HP', 'HQ', 'HPdual', 'HQdual'):
        details[label] = cyclic_description(codes[label])
        details[label]['basis_hex'] = [hex(x) for x in codes[label]]
        if label in ('HP', 'HQ'):
            bound = details[label]['extended_even_distance_lower']
            details[label]['doubly_even_distance_lower'] = 4 * ((bound + 3) // 4)
    containments = {a: {b: all(contains(codes[b], x) for x in codes[a])
                        for b in codes} for a in codes}
    intersections = {a: {b: len(codes[a]) + len(codes[b]) - len(echelon(codes[a] + codes[b]))
                         for b in codes} for a in codes}
    return dict(classification='Exact binary hull construction and punctured cyclic root audit',
                codes=details, containments=containments, intersection_dimensions=intersections,
                new_A38_cap_proved=False,
                source_sha256={str(p.relative_to(ROOT)): sha(p) for p in
                    (Path(__file__), ROOT/'code/audit_bch_quadratic_sums.py',
                     ROOT/'code/bch_quotient.py', ROOT/'code/affine_wambach.py')})


if __name__ == '__main__':
    value = build()
    if '--verify' in sys.argv:
        assert json.loads(OUTPUT.read_text()) == value
    else:
        write_new(OUTPUT, value)
    print(json.dumps({k: {x: y for x, y in v.items() if x not in ('basis_hex', 'roots')}
                      for k, v in value['codes'].items()}, indent=2))
    print(json.dumps(value['containments'], indent=2))
