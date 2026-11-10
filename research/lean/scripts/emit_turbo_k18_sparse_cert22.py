"""Untrusted producer of K18 sparse q1 certificates (regional chain, q1 rows,
z^d power traces) as small generated Lean modules.

Layout follows the old Scalable dense cover: ignored generated sources, a
manifest, one small `decide +kernel` theorem per item. Lean checks every
emitted integer; nothing here is trusted.
"""
from __future__ import annotations
import json
from fractions import Fraction as Q
from math import comb, log2
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / 'SpinCodes/Turbo/Generated/K18'
S_BITS, DEN_BITS = 384, 80
S, DEN = 1 << S_BITS, 1 << DEN_BITS
MULT = S // DEN
R = 128            # K18: L = 2048 = 16 R groups
L = 16 * R
D = (L * 256) // 10  # 52428 = floor(N/10)
CHUNK = 8
TILTS = {1: (1, 1250), 2: (3, 2500), 3: (1, 625)}
CHOICE = [None] * 8 + [1] * 6 + [2] + [3] * 50   # diagnostic's 3-tilt choice


def cdiv(a, b):
    return -((-a) // b)


def rs8_counts():
    M = (1 << 32) - 1
    def A(h):
        return comb(8, h) * sum((-1) ** j * comb(h, j) * ((1 << 32) ** (h - 4 - j) - 1)
                                for j in range(h - 4))
    base = [comb(8, i) * 15 ** i for i in range(9)]
    base[0] -= 1                                   # (1+15X)^8 - 1
    def pmul(a, b):
        c = [0] * (len(a) + len(b) - 1)
        for i, x in enumerate(a):
            if x:
                for j, y in enumerate(b):
                    c[i + j] += x * y
        return c
    total = [0] * 65
    for h in range(5, 9):
        p = [1]
        for _ in range(h):
            p = pmul(p, base)
        for v, c in enumerate(p):
            total[v] += A(h) * M ** (8 - h) * c
    return total, M ** 8


def row_fn(fam):
    births = {0: [DEN] + [0] * 10, 1: fam['birth_one_caps'], 2: fam['birth_two_caps']}
    def get(l, i):
        return l[i] if i < len(l) else 0
    def row(d, i, k):
        if i == 0:
            return get(births[d], k)
        if k == 0:
            return get(fam['return_caps'][d], i)
        if k == 2:
            return get(fam['first_caps'][d], i)
        return 0
    return row


def matmul(A, B):
    return [[sum(A[i][k] * B[k][j] for k in range(11)) for j in range(11)] for i in range(11)]


def madd(A, B):
    return [[A[i][j] + B[i][j] for j in range(11)] for i in range(11)]


def polymul_up(a, b):
    prods = [matmul(a[0], b[0]),
             madd(matmul(a[0], b[1]), matmul(a[1], b[0])),
             madd(madd(matmul(a[0], b[2]), matmul(a[1], b[1])), matmul(a[2], b[0]))]
    return [[[cdiv(x, S) for x in r] for r in P] for P in prods]


def regional(fam):
    row = row_fn(fam)
    binom = [1, 16, 120]
    base = [[[binom[d] * row(d, i, k) * MULT for k in range(11)] for i in range(11)]
            for d in range(3)]
    bits = bin(R)[3:]          # MSB-first after the implicit leading one
    trace, prev, e = [], base, 1
    for b in bits:
        assert b == '0'        # R = 128 is a pure power of two
        sq = polymul_up(prev, prev)
        trace.append(sq)
        prev, e = sq, 2 * e
    assert e == R
    dens = [1, 16 * R, (16 * R) * (16 * R - 1) // 2]
    result = [[[cdiv(prev[d][i][k], dens[d]) for k in range(11)] for i in range(11)]
              for d in range(3)]
    return base, trace, result


def q1_rows(r0, r1):
    rows = [[[S] + [0] * 10] + [[0] * 11 for _ in range(64)]]
    for _ in range(64):
        prev = rows[-1]
        nxt = []
        for v in range(65):
            out = []
            for j in range(11):
                acc = sum(prev[v][i] * r0[i][j] for i in range(11))
                if v:
                    acc += sum(prev[v - 1][i] * r1[i][j] for i in range(11))
                out.append(cdiv(acc, S))
            nxt.append(out)
        rows.append(nxt)
    return rows


def power_trace(znum):
    base = znum * MULT
    prev, trace = S, []
    for b in bin(D)[2:]:
        nxt = (prev * prev * base) // (S * S) if b == '1' else (prev * prev) // S
        assert nxt >= 0
        trace.append((b == '1', nxt))
        prev = nxt
    assert prev > 0
    return base, trace, prev


def lean_mat(M):
    return '[' + ',\n  '.join('[' + ','.join(map(str, r)) + ']' for r in M) + ']'


def lean_poly(P):
    return '[' + ',\n '.join(lean_mat(M) for M in P) + ']'


def ns(k):
    return f'Spin.Turbo.ProfileFamilyK18Tilt{k}22'


def emit_tilt(k, fam, manifest):
    base, trace, result = regional(fam)
    rows = q1_rows(result[0], result[1])
    zbase, ptrace, zlast = power_trace(fam['z_num'])
    mod = f'SpinCodes.Turbo.Generated.K18.Tilt{k}'
    header = (f'import SpinCodes.Turbo.SparseListDefs22\n'
              f'import SpinCodes.Turbo.PhysicalSmallBirth22K18Tilt{k}BudgetCheck\n\n'
              f'namespace {ns(k)}.Cert\n')
    reg = [header,
           f'def s : Int := {S}', f'def mult : Int := {MULT}',
           'def births (d : Nat) : List Int :=\n  match d with\n'
           '  | 0 => birthZeroCaps22\n  | 1 => birthOneCaps\n  | _ => birthTwoCaps',
           'def row (d i k : Nat) : Int :=\n'
           '  profileRowInt22 (births d) (firstCaps.getD d []) (returnCaps.getD d []) i k',
           f'def base : List (List (List Int)) := {lean_poly(base)}']
    for t, sq in enumerate(trace):
        reg.append(f'def sq{t} : List (List (List Int)) := {lean_poly(sq)}')
    reg.append('def trace : List (Bool × List (List (List Int)) × List (List (List Int))) := ['
               + ', '.join(f'(false, sq{t}, [])' for t in range(len(trace))) + ']')
    for d in range(3):
        reg.append(f'def res{d} : List (List Int) := {lean_mat(result[d])}')
    reg.append('def result : List (List (List Int)) := [res0, res1, res2]')
    reg.append(f'def zbaseInt : Int := {zbase}')
    reg.append('def ztrace : List (Bool × Int) := ['
               + ', '.join(f'({str(b).lower()}, {n})' for b, n in ptrace) + ']')
    reg.append('set_option maxHeartbeats 0 in\n'
               f'theorem regionalChecked :\n    sparseRegionalCertificateCheckL22 {R} s row mult base trace result = true := by\n'
               '  decide +kernel')
    reg.append('set_option maxHeartbeats 0 in\n'
               f'theorem powerChecked :\n    sparseIntPowerCertificateCheck22 s zbaseInt {D} ztrace = true := by\n'
               '  decide +kernel')
    reg.append('#print axioms regionalChecked\n#print axioms powerChecked')
    reg.append(f'end {ns(k)}.Cert\n')
    paths = []
    def write(name, text):
        p = OUT / f'Tilt{k}{name}.lean'
        p.parent.mkdir(parents=True, exist_ok=True)
        p.write_text(text, encoding='utf-8', newline='\n')
        paths.append(str(p.relative_to(ROOT)).replace('\\', '/'))
    write('Regional', '\n\n'.join(reg))
    nchunks = (64 + CHUNK) // CHUNK          # rows 0..64
    for c in range(nchunks):
        lo, hi = c * CHUNK, min(64, c * CHUNK + CHUNK - 1)
        body = [f'import SpinCodes.Turbo.SparseListDefs22\n\nnamespace {ns(k)}.Cert']
        for n in range(lo, hi + 1):
            if n == 0:
                body.append('def rows0 : List (List Int) := sparseQ1InitialRowsL22 '
                            f'{S}')
            else:
                body.append(f'def rows{n} : List (List Int) := {lean_mat(rows[n])}')
        body.append(f'end {ns(k)}.Cert\n')
        write(f'Rows{c}', '\n\n'.join(body))
    for c in range(nchunks):
        lo, hi = c * CHUNK, min(63, c * CHUNK + CHUNK - 1)
        if lo > 63:
            continue
        imports = [f'import {mod}Regional', f'import {mod}Rows{c}']
        if (hi + 1) // CHUNK != c:
            imports.append(f'import {mod}Rows{(hi + 1) // CHUNK}')
        body = ['\n'.join(imports) + f'\n\nnamespace {ns(k)}.Cert']
        for n in range(lo, hi + 1):
            body.append('set_option maxHeartbeats 0 in\n'
                        f'theorem step{n} :\n    sparseQ1StepCheckL22 s res0 res1 rows{n} rows{n + 1} = true := by\n'
                        '  decide +kernel')
        body.append(f'end {ns(k)}.Cert\n')
        write(f'Steps{c}', '\n\n'.join(body))
    manifest['tilts'][k] = dict(paths=paths, zlast=zlast,
                               terminal=[sum(rows[64][v]) for v in range(65)])
    return rows, zlast


def main():
    counts_num, counts_den = rs8_counts()
    ref = json.loads(Path('C:/Users/peter/repo/permute_conv-github-bch/tmp/k18-rs8-gf32-s22-v1.json').read_text())
    exact = list(map(Q, ref['exact_expected_shells']))
    assert all(Q(counts_num[v], counts_den) == exact[v] for v in range(65)), 'RS8 count mismatch'
    (OUT).mkdir(parents=True, exist_ok=True)
    (OUT.parent / '.gitignore').write_text('*\n!.gitignore\n', encoding='utf-8')
    manifest = dict(schema='turbo-k18-sparse-cert-manifest-1', R=R, d=D, s_bits=S_BITS,
                    choice=CHOICE, tilts={})
    terms = []
    for k in (1, 2, 3):
        fam = json.loads((ROOT / f'.lake/turbo/profile_family_k18_tilt{k}.json').read_text())
        rows, zlast = emit_tilt(k, fam, manifest)
        manifest['tilts'][k]['rows'] = None
        for v in range(65):
            if CHOICE[v] == k:
                # term = (count_v / choose(64,v)) * (s / zlast) * (sum_j rows/s)
                terms.append((v, Q(counts_num[v], counts_den) / comb(64, v) /
                              Q(zlast) * sum(rows[64][v])))
    for v in range(8):
        terms.append((v, Q(counts_num[v], counts_den)))
    q1 = L * sum(t for _, t in terms)
    manifest['q1_upper'] = dict(num=q1.numerator, den=q1.denominator)
    manifest['q1_bits'] = -log2(q1.numerator) + log2(q1.denominator)
    (OUT / 'manifest.json').write_text(json.dumps(manifest, indent=1, default=str) + '\n')
    print('q1 bits', manifest['q1_bits'])


if __name__ == '__main__':
    main()
