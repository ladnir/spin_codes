#!/usr/bin/env python3
"""Geometry-generic s22 Turbo/RS8 producers for the K-ladder (untrusted; Lean checks all).

    python scripts/ladder22.py E stage [stage ...]
stages: config  families  q1  q2  (in order). E = log2 K, K in {2^14,..,2^26}.

Names: tilt families `Spin.Turbo.ProfileFamilyG{E}Tilt{k}22`
(SpinCodes/Turbo/PhysicalProfileFamilyG{E}Tilt{k}22*.lean), certificates under
SpinCodes/Turbo/Generated/G{E}/, endpoints `Spin.Turbo.G{E}Q1.q1_le`,
`Spin.Turbo.G{E}Q2.q2_le`. The K20 producers are reused by name substitution
(their arithmetic is geometry-parametric through R, L, D).
"""
from __future__ import annotations

import json
import re
import sys
import types
from decimal import Decimal, ROUND_HALF_EVEN, localcontext
from fractions import Fraction as Q
from math import comb, log2
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'scripts'))
sys.set_int_max_str_digits(0)
from emit_turbo_profile_family22 import produce_family, emit_lean  # noqa: E402

T = ROOT / 'SpinCodes/Turbo'
NL = chr(10)
CANDIDATE = None


def geom(e):
    K = 1 << e
    L = K // 128
    return dict(e=e, K=K, L=L, R=L // 16, D=(L * 256) // 10)


def cfg_path(e):
    p = ROOT / f'.lake/turbo/ladder/g{e}'
    p.mkdir(parents=True, exist_ok=True)
    return p / 'config.json'


def load_module(path, subs, name):
    src = (ROOT / path).read_text(encoding='utf-8')
    for a, b in subs:
        assert a in src, (path, a)
        src = src.replace(a, b)
    mod = types.ModuleType(name)
    mod.__file__ = str(ROOT / path)
    exec(compile(src, str(ROOT / path), 'exec'), mod.__dict__)
    sys.modules[name] = mod
    return mod


def cert_module(e, g, tilts, choice):
    m = load_module('scripts/emit_turbo_k20_sparse_cert22.py',
                    [('K20', f'G{e}'), ('k20', f'g{e}'), ('for k in (1, 2, 3, 4):', 'for k in sorted(TILTS):')],
                    f'ladder_cert_{e}')
    m.R, m.L, m.D = g['R'], g['L'], g['D']
    m.TILTS = {k: tuple(v) for k, v in tilts.items()}
    m.CHOICE = choice
    return m


# ---------------------------------------------------------------- config
def family_json(lam):
    global CANDIDATE
    cache = ROOT / '.lake/turbo/q1screen'
    cache.mkdir(parents=True, exist_ok=True)
    path = cache / f'fam_{lam.numerator}_{lam.denominator}.json'
    if path.exists():
        return json.loads(path.read_text())
    if CANDIDATE is None:
        CANDIDATE = json.loads((ROOT / '.lake/turbo/coarse_packet_profiles.json').read_text())
    zden = 1 << 80
    with localcontext() as ctx:
        ctx.prec = 100
        znum = int(((-Decimal(lam.numerator) / Decimal(lam.denominator)).exp() * Decimal(zden))
                   .to_integral_value(rounding=ROUND_HALF_EVEN))
    data = produce_family(CANDIDATE, znum, zden, 'Spin.Turbo.Screen')
    path.write_text(json.dumps(data))
    return data


def stage_config(e):
    g = geom(e)
    m = cert_module(e, g, {}, [None] * 65)
    cnt, den = m.rs8_counts()
    grid = sorted({Q(round(x * (1 << 20) / g['D']), 1 << 20) for x in
                   [2 ** (k / 8) for k in range(24, 69)] if round(x * (1 << 20) / g['D']) > 0})
    terms = {}
    for lam in grid:
        fam = family_json(lam)
        _, _, result = m.regional(fam)
        rows = m.q1_rows(result[0], result[1])
        _, _, zlast = m.power_trace(fam['z_num'])
        terms[lam] = [Q(cnt[v], den) / comb(64, v) / Q(zlast) * sum(rows[64][v]) for v in range(65)]
    # monotone tilt assignment over shells 8..64 (contiguous ranges), DP on exact rationals
    lams = sorted(terms)
    shells = list(range(8, 65))
    INF = None
    best = {li: (terms[lams[li]][8], [li]) for li in range(len(lams))}
    for v in shells[1:]:
        nxt = {}
        run_best = None
        for li in range(len(lams)):
            if run_best is None or best[li][0] < run_best[0]:
                run_best = best[li]
            nxt[li] = (run_best[0] + terms[lams[li]][v], run_best[1] + [li])
        best = nxt
    total, path = min(best.values(), key=lambda x: x[0])
    used = sorted(set(path))
    idx = {li: k + 1 for k, li in enumerate(used)}
    tilts = {idx[li]: [lams[li].numerator, lams[li].denominator] for li in used}
    choice = [None] * 8 + [idx[li] for li in path]
    q1 = g['L'] * (sum(Q(cnt[v], den) for v in range(8)) + total)
    bits = -log2(q1.numerator) + log2(q1.denominator)
    # q2 tilt: the used tilt whose lambda*d is closest to K20's 1/2500 * 209715
    target = Q(209715, 2500)
    q2k = min(tilts, key=lambda k: abs(Q(*tilts[k]) * g['D'] - target))
    cfg = dict(**g, tilts=tilts, choice=choice, q1_bits=bits, q2_tilt=q2k)
    cfg_path(e).write_text(json.dumps(cfg, indent=1, default=str))
    print(f'G{e}: q1 bits {bits:.4f} with {len(tilts)} tilts; q2 tilt {q2k}', flush=True)


def load_cfg(e):
    c = json.loads(cfg_path(e).read_text())
    c['tilts'] = {int(k): v for k, v in c['tilts'].items()}
    return c


# ---------------------------------------------------------------- families
TEMPLATES = ['PhysicalProfileFamilyK18Tilt122Check.lean', 'PhysicalProfileFamilyK18Tilt122.lean',
             'PhysicalProfileBirthFamilyK18Tilt122.lean', 'PhysicalSmallBirth22K18Tilt1BudgetCheck.lean',
             'PhysicalSmallBirth22K18Tilt1BudgetSmoke.lean']


def stage_families(e):
    c = load_cfg(e)
    for k, (num, den) in c['tilts'].items():
        lam = Q(num, den)
        data = family_json(lam)
        data = dict(data)
        ns = f'Spin.Turbo.ProfileFamilyG{e}Tilt{k}22'
        fam = produce_family(json.loads((ROOT / '.lake/turbo/coarse_packet_profiles.json').read_text()),
                             data['z_num'], data['z_den'], ns)
        fam['selection_only_lambda'] = dict(num=num, den=den)
        (ROOT / f'.lake/turbo/profile_family_g{e}_tilt{k}.json').write_text(json.dumps(fam, indent=2) + NL)
        emit_lean(fam, T / f'PhysicalProfileFamilyG{e}Tilt{k}22Data.lean')
        for name in TEMPLATES:
            text = (T / name).read_text(encoding='utf-8').replace('K18Tilt1', f'G{e}Tilt{k}')
            assert 'K18' not in text
            (T / name.replace('K18Tilt1', f'G{e}Tilt{k}')).write_text(text, encoding='utf-8', newline=NL)
    print(f'G{e}: {len(c["tilts"])} families', flush=True)


# ---------------------------------------------------------------- q1
def ranges(choice):
    out = []
    for v in range(8, 65):
        k = choice[v]
        if out and out[-1][0] == k:
            out[-1][2] = v
        else:
            out.append([k, v, v])
    return out


def stage_q1(e):
    c = load_cfg(e)
    g = geom(e)
    m = cert_module(e, g, c['tilts'], c['choice'])
    m.main()
    asm = load_module('scripts/emit_turbo_k20_q1_assembly22.py',
                      [('K20', f'G{e}'), ('D = 209715', f'D = {g["D"]}'),
                       ('R = 512', f'R = {g["R"]}'),
                       ('for k in (1, 2, 3, 4):', f'for k in {tuple(sorted(c["tilts"]))}:')],
                      f'ladder_asm_{e}')
    for k in sorted(c['tilts']):
        asm.emit(k)
    emit_q1_endpoint(e, c, g, m)


def emit_q1_endpoint(e, c, g, P):
    OUT = ROOT / f'SpinCodes/Turbo/Generated/G{e}'
    K = 256
    M8 = ((1 << 32) - 1) ** 8
    ks = sorted(c['tilts'])
    Tn = {k: f'Spin.Turbo.ProfileFamilyG{e}Tilt{k}22' for k in ks}
    rg = ranges(c['choice'])
    cond = {}
    for k, lo, hi in rg:
        cond[k] = f'{lo} ≤ v ∧ v ≤ {hi}' if hi < 64 else f'{lo} ≤ v'
    man = json.loads((OUT / 'manifest.json').read_text())
    cnt, den = P.rs8_counts()
    a = []
    for v in range(65):
        if v < 8:
            a.append(-((-cnt[v] * (1 << K)) // M8))
        else:
            t = man['tilts'][str(c['choice'][v])]
            num = cnt[v] * t['terminal'][v] * (1 << K)
            dd = M8 * comb(64, v) * t['zlast']
            a.append(-((-num) // dd))
    tmpl = (ROOT / 'scripts/templates/k20_q1_endpoint.lean.tmpl').read_text(encoding='utf-8')
    head_end = tmpl.index('import SpinCodes.Turbo.SparseQ1Choice22')
    tmpl = ''.join(f'import SpinCodes.Turbo.Generated.G{e}.Tilt{k}Q1' + NL for k in ks) + tmpl[head_end:]
    z0 = tmpl.index('def zK (v : ℕ) : ℝ :=')
    z1 = tmpl.index('theorem chosen_le')
    def chain(field):
        lines = []
        for i, (k, lo, hi) in enumerate(rg):
            if i < len(rg) - 1:
                lines.append(f'  if v ≤ {hi} then {Tn[k]}.{field} else')
            else:
                lines.append(f'  {Tn[k]}.{field}')
        return NL.join(lines)
    bul = lambda f: NL.join(f'  · exact {f(k)}' for k, _, _ in rg)
    def cases_proof(head, unf, ex, both=False):
        L = [f'theorem {head} := by', f'  unfold {unf}']
        for i, (k, lo, hi) in enumerate(rg):
            if i == len(rg) - 1:
                L.append(f'  exact {ex(k)}')
            else:
                r = f'if_pos h{i}, if_pos h{i}' if both else f'if_pos h{i}'
                n = f'if_neg h{i}, if_neg h{i}' if both else f'if_neg h{i}'
                L.append(f'  by_cases h{i} : v ≤ {hi}')
                L.append(f'  · rw [{r}]; exact {ex(k)}')
                L.append(f'  rw [{n}]')
        return NL.join(L)
    mid = f'''def zK (v : ℕ) : ℝ :=
{chain('birthSmokeZ22')}

def BK (v : ℕ) : Fin 17 → Matrix (Fin 11) (Fin 11) ℝ :=
{chain('physicalRows')}

def useC (v : ℕ) : Bool := decide (8 ≤ v)

@ZPOS@

{cases_proof('hz0 (v : ℕ) : 0 < zK v', 'zK', lambda k: f'zpos{k}')}

{cases_proof('hz1 (v : ℕ) : zK v ≤ 1', 'zK', lambda k: f'zle{k}')}

{cases_proof('hBK (v : ℕ) (j : Fin 17) :' + NL + '    physicalRowSound22 (Matrix.of (turboPhysicalKernel22 j (zK v))) (BK v j)', 'zK BK', lambda k: f'{Tn[k]}.physicalRows_sound j', both=True)}

'''
    tmpl = tmpl[:z0] + mid + tmpl[z1:]
    d0 = tmpl.index('    rcases (by omega')
    d1 = tmpl.index('  · rw [if_neg (by simp [useC, h8])]')
    alts = ' ∨ '.join(f'({cond[k]})' for k, _, _ in rg)
    pats = ' | '.join(['ht'] * len(rg))
    disp = [f'    rcases (by omega : {alts})', f'      with {pats}']
    for i, (k, lo, hi) in enumerate(rg):
        negs = ', '.join(['if_neg (by omega)'] * i)
        last = (i == len(rg) - 1)
        step = negs if last else (negs + (', ' if negs else '') + 'if_pos (by omega)')
        disp.append(f'    · have hz : zK v = {Tn[k]}.birthSmokeZ22 := by')
        disp.append(f'        unfold zK; rw [{step}]')
        disp.append(f'      have hB : BK v = {Tn[k]}.physicalRows := by')
        disp.append(f'        unfold BK; rw [{step}]')
        disp.append('      rw [hz, hB]')
        disp.append(f'      exact tiltTerm{k} points hp v hv\' ht')
    tmpl = tmpl[:d0] + NL.join(disp) + NL + tmpl[d1:]
    tmpl = tmpl.replace('8192', str(g['L'])).replace('K20', f'G{e}')
    # checks / tilt terms / zpos from the K20 emitter helpers, with generic conditions
    ep = load_module('scripts/emit_turbo_k20_q1_endpoint22.py',
                     [('import emit_turbo_k20_sparse_cert22 as P', f'P = sys.modules["ladder_cert_{e}"]'),
                      ('K20', f'G{e}'), ('k20', f'g{e}')], f'ladder_ep_{e}')
    ep.T = Tn
    ep.CONDS = cond
    subs = {'@ALIT@': '[' + ','.join(map(str, a)) + ']', '@ASUM@': str(sum(a)),
            '@CHECKS@': NL.join(ep.check_block(k) for k in ks),
            '@TILTS@': NL.join(ep.tilt_block(k) for k in ks),
            '@ZPOS@': NL.join(ep.zpos_block().split('theorem zpos1')[0:1]) if False else None,
            '@K@': str(K), '@D@': str(g['D']), '@R@': str(g['R'])}
    zp = []
    for k in ks:
        zp.append(f'''theorem zpos{k} : 0 < {Tn[k]}.birthSmokeZ22 := by
  unfold {Tn[k]}.birthSmokeZ22
  exact div_pos (by exact_mod_cast (show (0 : ℤ) < {Tn[k]}.zNum by decide +kernel))
    (by exact_mod_cast (show (0 : ℤ) < {Tn[k]}.zDen by decide +kernel))

theorem zle{k} : {Tn[k]}.birthSmokeZ22 ≤ 1 := by
  unfold {Tn[k]}.birthSmokeZ22
  rw [div_le_one (by exact_mod_cast (show (0 : ℤ) < {Tn[k]}.zDen by decide +kernel))]
  exact_mod_cast (show {Tn[k]}.zNum ≤ {Tn[k]}.zDen by decide +kernel)
''')
    subs['@ZPOS@'] = NL.join(zp)
    for key, val in subs.items():
        tmpl = tmpl.replace(key, val)
    assert '@' not in tmpl.replace('@[', ''), re.findall(r'@\w+@', tmpl)
    (OUT / 'Q1Endpoint.lean').write_text(tmpl, encoding='utf-8', newline=NL)
    q = Q(g['L'] * sum(a), 1 << K)
    man['q1_endpoint'] = dict(K=K, aSum=sum(a), bits=-log2(q.numerator) + log2(q.denominator))
    (OUT / 'manifest.json').write_text(json.dumps(man, indent=1, default=str) + NL)
    print(f'G{e}: q1 endpoint bits {man["q1_endpoint"]["bits"]:.5f}', flush=True)


# ---------------------------------------------------------------- q2
def stage_q2(e):
    c = load_cfg(e)
    g = geom(e)
    k2 = c['q2_tilt']
    if f'ladder_cert_{e}' not in sys.modules:
        cert_module(e, g, c['tilts'], c['choice'])
    m = load_module('scripts/emit_turbo_k20_q2_endpoint22.py',
                    [('import emit_turbo_k20_sparse_cert22 as P', f'P = sys.modules["ladder_cert_{e}"]'),
                     ("'.lake/turbo/profile_family_k20_tilt2.json'", f"'.lake/turbo/profile_family_g{e}_tilt{k2}.json'"),
                     ("zlast = man['tilts']['2']['zlast']", f"zlast = man['tilts']['{k2}']['zlast']"),
                     ('C2 = comb(8192, 2)', f'C2 = comb({g["L"]}, 2)'),
                     ("TEMPLATE = ROOT / 'scripts/templates/k20_q2_endpoint.lean.tmpl'",
                      f"TEMPLATE = ROOT / '.lake/turbo/ladder/g{e}/q2.tmpl'"),
                     ("'@T3@': 'Spin.Turbo.ProfileFamilyK20Tilt222'}", f"'@T3@': 'Spin.Turbo.ProfileFamilyG{e}Tilt{k2}22'}}"),
                     ('K20', f'G{e}')], f'ladder_q2_{e}')
    t = (ROOT / 'scripts/templates/k20_q2_endpoint.lean.tmpl').read_text(encoding='utf-8')
    for a, b in [('def C2 : Int := 33550336', f'def C2 : Int := {comb(g["L"], 2)}'),
                 ('@R@ * 16 = 8192', f'@R@ * 16 = {g["L"]}'),
                 ('zpos2 zle2', f'zpos{k2} zle{k2}'), ('lastT2', f'lastT{k2}'), ('K20', f'G{e}')]:
        assert a in t, a
        t = t.replace(a, b)
    (ROOT / f'.lake/turbo/ladder/g{e}/q2.tmpl').write_text(t, encoding='utf-8', newline=NL)
    m.main()


if __name__ == '__main__':
    e = int(sys.argv[1])
    for st in sys.argv[2:]:
        dict(config=stage_config, families=stage_families, q1=stage_q1, q2=stage_q2)[st](e)
