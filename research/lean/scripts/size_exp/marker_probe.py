#!/usr/bin/env python3
"""Size experiment (A), Lean probe: certify a G22 marker *range* piece over q values
that the current plan covers with exact pieces, reusing the already-built G22 family
and pair modules (imported read-only).  Untrusted producer; Lean checks everything.

Writes only: the scratch cache dir (argv[1]) and SpinCodes/Turbo/SizeExp/<Name>.lean.
Usage: python scripts/size_exp/marker_probe.py SCRATCH qlo qhi [Name]
"""
import json
import shutil
import sys
import types
from fractions import Fraction as Q
from math import log2
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(ROOT / 'scripts'))
sys.path.insert(0, str(HERE))
sys.set_int_max_str_digits(0)
import qstar_corrected as QC  # noqa: E402

E = 22
SCR = Path(sys.argv[1])
QLO, QHI = int(sys.argv[2]), int(sys.argv[3])
NAME = sys.argv[4] if len(sys.argv) > 4 else 'G22MarkerProbe'
REAL = ROOT / '.lake/turbo/ladder/g22'
lch = QC.lch


def load_K():
    L = (1 << E) // 128
    NP, D = L // 32, (L * 256) // 10
    src = (ROOT / 'scripts/emit_turbo_k20_tail22.py').read_text(encoding='utf-8')
    subs = [("OUT = ROOT / 'SpinCodes/Turbo/Generated/K20Tail'", "OUT = ROOT / 'SpinCodes/Turbo/Generated/G22Tail'"),
            ("CACHE = ROOT / '.lake/turbo/k20tail'", f"CACHE = Path(r'{SCR}')"),
            ("GROUPS, NPAIRS, DCUT = 8192, 256, 209715", f"GROUPS, NPAIRS, DCUT = {L}, {NP}, {D}"),
            ("NSQ_EXACT = 8", f"NSQ_EXACT = {NP.bit_length() - 1}"),
            ("NSQ_MARK = 14", f"NSQ_MARK = {(64 * NP).bit_length() - 1}"),
            ("        for q in range(p['qlo'], p['qhi'] + 1):\n            Cq = comb(GROUPS, q)\n            vtr, G, yn = prefactor",
             "        for q in ((p['qlo'], p['qhi']) if p.get('range') else range(p['qlo'], p['qhi'] + 1)):\n            Cq = comb(GROUPS, q)\n            vtr, G, yn = prefactor"),
            ('K20Tail', 'G22Tail')]
    for a, b in subs:
        assert a in src, a[:80]
        src = src.replace(a, b)
    m = types.ModuleType('size_exp_k22')
    m.__file__ = str(ROOT / 'scripts/emit_turbo_k20_tail22.py')
    exec(compile(src, m.__file__, 'exec'), m.__dict__)
    return m, L, NP, D


def main():
    K, L, NP, D = load_K()
    cfg = json.loads((REAL / 'config.json').read_text())
    target = -(cfg['q1_bits'] + 20)
    ps = json.loads((REAL / 'tail_plan.json').read_text())
    fidx = {t: i for i, t in enumerate(K.tilts(ps), start=1)}
    # screen the existing (tilt, alpha) pair keys for this q range (Lean normalizer C^-63)
    tq = [Q(2.0 ** (x / 16)).limit_denominator(4096) for x in range(-192, 16 * 16)]
    ts = np.array([float(x) for x in tq])
    lgt = np.log2(ts)
    best, seen = None, set()
    for p in ps:
        key = (p['tilt'], p['alpha'])
        if key in seen:
            continue
        seen.add(key)
        lam = Q(p['tilt'])
        fpath = ROOT / f'.lake/turbo/q1screen/hi/{lam.numerator}_{lam.denominator}.json'
        if not fpath.exists():
            continue
        B, lz = QC.fam_B(fpath)
        al = float(Q(p['alpha']))
        lr = QC.log2rho(QC.paired(B, al), ts)
        v = np.maximum(*[-63 * lch(L, q) + al * (q * QC.C1 - D * lz) - 64 * q * lgt + 64 * NP * lr
                         for q in (QLO, QHI)])
        i = int(np.argmin(v))
        if best is None or v[i] < best[0]:
            best = (float(v[i]), p, tq[i])
    v, p0, t = best
    print(f'screen: q {QLO}..{QHI} best float term {v:.1f} (target {target:.2f}) tilt {p0["tilt"]} '
          f'alpha {p0["alpha"]} t {t}', flush=True)
    probe = dict(kind='marker', tilt=p0['tilt'], alpha=p0['alpha'], p=None, t=f'{t.numerator}/{t.denominator}',
                 qlo=QLO, qhi=QHI, range=True)
    ps2 = ps + [probe]
    pidx = len(ps2)
    f = fidx[p0['tilt']]
    key = K.pair_key(probe, f)
    SCR.mkdir(parents=True, exist_ok=True)
    for nm in (f'family{f}.json', f'pairs_{key}.json'):
        if not (SCR / nm).exists():
            shutil.copyfile(REAL / 'tail' / nm, SCR / nm)
    if (SCR / f'marker{pidx}.json').exists():
        (SCR / f'marker{pidx}.json').unlink()
    K.stage_marker(ps2, {pidx})
    mk = json.loads((SCR / f'marker{pidx}.json').read_text())
    emit(K, mk, probe, f, key)


def emit(K, mk, p, f, key):
    lm = lambda M: '[' + ','.join('[' + ','.join(map(str, r)) + ']' for r in M) + ']'
    ea, eb = mk['ea'], mk['eb']
    t = key.split('_', 1)[1]
    famns = K.fam_ns(f)
    opens = f'Spin.Turbo Spin.Turbo.Hi {famns} {famns}.A{t}'
    tn0, td0, kp, kn, up, down = mk['tn0'], mk['td0'], mk['kp'], mk['kn'], mk['up'], mk['down']
    st = '[' + ', '.join(f'({dn}, {u}, {lm(U)})' for dn, u, U in mk['st']) + ']'
    ztr = '[' + ', '.join(f'({str(b).lower()}, {sh}, {mm})' for b, sh, mm in mk['ztr']) + ']'
    recs = {int(q): r for q, r in mk['recs'].items()}
    lo, hi = p['qlo'], p['qhi']
    A = max(recs[lo]['a'], recs[hi]['a'])
    S, KT = K.S, K.KT
    tex = f'((({tn0} * 2 ^ {kp} : ℤ) : ℝ) / (({td0} * 2 ^ {kn} : ℤ) : ℝ))'
    al = f'(({ea} : ℕ) : ℝ) / (({eb} : ℕ) : ℝ)'

    def cert(q, nm):
        r = recs[q]
        vt = '[' + ', '.join(f'({str(b).lower()}, {sh}, {mm})' for b, sh, mm in r['vtr']) + ']'
        return f"""def vt_{q} : List (Bool × Nat × Int) := {vt}

set_option maxHeartbeats 0 in
theorem {nm} : TailMarkerRange22.markerF {ea} {eb} {K.NPAIRS} {K.DCUT} selectedRows birthSmokeZ22 {tex} {q} ≤
    (({A} : ℤ) : ℝ) / 2 ^ {KT} :=
  TailMarkerRange22.markerF_le {ea} {eb} (by norm_num) {K.NPAIRS} {K.DCUT} 80 selectedRows
    selectedRows_nonneg {q} (binomQ {q}) (binomQ_eq' {q} (by norm_num)) (binomQ_pos {q} (by norm_num))
    (16 ^ 16 * endpointDen) (2 ^ {S}) (by norm_num) R hR Ys rpowChecked {tn0} {td0} (by norm_num) (by norm_num)
    {kp} {kn} (2 ^ {S}) {S} rfl {down} {up} U0 baseChecked stages (by decide) chainChecked birthSmokeZ22 zNum hz
    ztrace zChecked vt_{q} (by decide +kernel) ({r["yn"]}) ({r["G"]}) (by decide +kernel) ({A}) {KT}
    (by decide +kernel)
"""
    src = f"""import SpinCodes.Turbo.Generated.G22Tail.Tail{f}A{t}Pair
import SpinCodes.Turbo.TailMarkerRange22

/-! Size experiment (generated by scripts/size_exp/marker_probe.py; not part of the
ladder).  One G22 marker range piece q = {lo}..{hi}, which the current G22 plan
covers with exact pieces.  Certified bound per q: A / 2^{KT} with
log2(A) - {KT} = {log2(A) - KT:.1f}. -/

noncomputable section
namespace Spin.Turbo.SizeExp.{NAME}
open {opens}

def U0 : List (List Int) := {lm(mk['U0'])}
def stages : List (Nat × Nat × List (List Int)) := {st}
def ztrace : List (Bool × Nat × Int) := {ztr}

set_option maxHeartbeats 0 in
theorem baseChecked : tailMarkerPairBaseCheck22 Ys ({tn0} * 2 ^ {kp}) ({td0} * 2 ^ {kn}) (2 ^ {S}) (2 ^ {S})
    {down} {up} U0 = true := by
  decide +kernel

set_option maxHeartbeats 0 in
theorem chainChecked : sparseSquareChainCheckL22 (2 ^ {S}) U0 stages = true := by decide +kernel

set_option maxHeartbeats 0 in
theorem zChecked : sparseScaledPowerCheck22 zNum 80 {K.DCUT} ztrace = true := by decide +kernel

{cert(lo, 'Flo')}
{cert(hi, 'Fhi')}

def aP (_q : ℕ) : Int := {A}

theorem zpos : 0 < birthSmokeZ22 := by
  unfold birthSmokeZ22
  exact div_pos (by exact_mod_cast (show (0 : ℤ) < zNum by decide +kernel))
    (by exact_mod_cast (show (0 : ℤ) < zDen by decide +kernel))

theorem pieceTerm (q : ℕ) (hq : q ∈ Finset.Icc {lo} {hi}) :
    (Nat.choose (32 * {K.NPAIRS}) q : ℝ) * ((turboOldUniformFactor22 ^ q / birthSmokeZ22 ^ {K.DCUT}) ^ ({al}) *
      ∑ k : Fin 11, ((pairedRegionDP22 selectedRows ({al}) {K.NPAIRS} q) ^ 64) 0 k) ≤ (aP q : ℝ) / 2 ^ {KT} := by
  obtain ⟨h1, h2⟩ := Finset.mem_Icc.mp hq
  have ht : (0 : ℝ) < {tex} := by norm_num
  exact (TailMarkerRange22.term_le_endpoints {ea} {eb} {K.NPAIRS} {K.DCUT} selectedRows selectedRows_nonneg
    birthSmokeZ22 {tex} zpos ht {lo} {hi} (by norm_num) q h1 h2).trans (max_le Flo Fhi)

#print axioms pieceTerm
end Spin.Turbo.SizeExp.{NAME}
"""
    out = ROOT / f'SpinCodes/Turbo/SizeExp/{NAME}.lean'
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(src, encoding='utf-8', newline='\n')
    print(f'wrote {out} ({out.stat().st_size} bytes); A = {A}, log2(A)-KT = {log2(A) - KT:.2f}', flush=True)


if __name__ == '__main__':
    main()
