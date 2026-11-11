"""Lean emission for the K20 tail (text only; Lean checks every statement).

Module layout under SpinCodes/Turbo/Generated/K20Tail (prefix `Tail`):
  {f}Data/Check/BudgetCheck/BudgetSmoke/BirthFamily/Family   family f (stage families)
  {f}A{t}R, {f}A{t}Y{a}, {f}A{t}Ys                         pair data (core-only)
  {f}A{t}Rp{a}                                             per-a rpow checks (core-only)
  {f}A{t}Pair                                              rows, hR, hz, zChecked, rpowChecked
  P{p}Base, P{p}Sq{s}                                      exact piece polynomials (per-matrix defs)
  P{p}BaseC{c}                                             base-check chunks
  P{p}St{s}R{i}[E{j}], P{p}St{s}                           packed squaring steps
  P{p}Chain                                                baseChecked, chainChecked
  P{p}QD{c}, P{p}T{c}                                      per-q data and term theorems
  P{p}MBase, P{p}MQ{c}                                     marker base/chain, per-q marker terms
  P{p}Agg                                                  aP (balanced lookup) + pieceTerm
"""
from __future__ import annotations

import json
from math import comb
from pathlib import Path

import emit_turbo_k20_tail22 as K

OUT, CACHE = K.OUT, K.CACHE
NL = chr(10)
PRE = 'SpinCodes.Turbo.Generated.' + K.OUT.name + '.Tail'
S, KT = K.S, K.KT
DEN = '(16 ^ 16 * 2 ^ 384)'
# Terms per term module. Lean keeps kernel memory across declarations, so large-L
# rungs (G22+: one more regional squaring per term) need smaller chunks; override
# with TURBO_QCHUNK_EXACT / TURBO_QCHUNK_MARK instead of raising memory limits.
import os
QCHUNK_EXACT = int(os.environ.get('TURBO_QCHUNK_EXACT', 30))
QCHUNK_MARK = int(os.environ.get('TURBO_QCHUNK_MARK', 40))
BASE_CH = 12


def lm(M):
    return '[' + ','.join('[' + ','.join(map(str, r)) + ']' for r in M) + ']'


def w(name, text):
    (OUT / f'Tail{name}.lean').write_text(text, encoding='utf-8', newline=NL)


def poly_defs(name, P):
    parts = [f'def {name}_part{i} : List (List Int) := {lm(M)}' for i, M in enumerate(P)]
    parts.append(f'def {name} : List (List (List Int)) := [' + ', '.join(f'{name}_part{i}' for i in range(len(P))) + ']')
    return NL.join(parts)


def alpha_lean(ea, eb):
    return f'(({ea} : ℕ) : ℝ) / (({eb} : ℕ) : ℝ)'


def ptag(key):
    return key.split('_', 1)[1]


# ---------------------------------------------------------------- pairs
def emit_pairs(key):
    pr = json.loads((CACHE / f'pairs_{key}.json').read_text())
    f, ea, eb = pr['f'], pr['ea'], pr['eb']
    t = ptag(key)
    ns = f'{K.fam_ns(f)}.A{t}'
    base = f'{f}A{t}'
    core = 'import SpinCodes.Turbo.SparseListDefs22'
    w(f'{base}R', f'{core}{NL}{NL}namespace {ns}{NL}{poly_defs("R", pr["R"])}{NL}end {ns}{NL}')
    for a in range(17):
        w(f'{base}Y{a}', f'{core}{NL}{NL}namespace {ns}{NL}{poly_defs(f"Y{a}", pr["Y"][a])}{NL}end {ns}{NL}')
    w(f'{base}Ys', NL.join([f'import {PRE}{base}Y{a}' for a in range(17)]) + f'{NL}{NL}namespace {ns}{NL}'
      f'def Ys : List (List (List (List Int))) := [' + ', '.join(f'Y{a}' for a in range(17)) + f']{NL}end {ns}{NL}')
    for a in range(17):
        w(f'{base}Rp{a}', f'import {PRE}{base}R{NL}import {PRE}{base}Ys{NL}import SpinCodes.Turbo.PairedFracABDefs22'
          f'{NL}{NL}namespace {ns}{NL}open Spin.Turbo{NL}{NL}set_option maxHeartbeats 0 in{NL}'
          f'theorem rpow_{a} : pairedRpowABRowL22 {ea} {eb} {DEN} (2 ^ {S}) R Ys {a} = true := by{NL}'
          f'  decide +kernel{NL}{NL}end {ns}{NL}')
    arms = NL.join(f'    | {a}, _ => rpow_{a}' for a in range(17))
    fam = K.fam_mod(f)
    pair = f'''import {fam}Family
import {PRE}{base}R
import {PRE}{base}Ys
{NL.join(f'import {PRE}{base}Rp{a}' for a in range(17))}
import SpinCodes.Turbo.Hi.RowSelected22
import SpinCodes.Turbo.TailMarkerFracAB22
import SpinCodes.Turbo.BinomFast22

noncomputable section
namespace {ns}
open Spin.Turbo Spin.Turbo.Hi {K.fam_ns(f)}

def births (n : Nat) : List Int :=
  match n with
  | 0 => birthZeroCaps22
  | 1 => birthOneCaps
  | _ => birthTwoCaps

def physInt (j i k : Nat) : Int := profileRowIntAll22 births firstCaps returnCaps zeroHoldCaps j i k

def row (d i k : Nat) : Int := rowSel physInt d i k

theorem hsmall (j : Fin 3) (i : Fin 11) :
    smallBirthSmokeCaps22 j i = ((births j.val).getD i.val 0 : ℝ) / endpointDen := by
  fin_cases j <;> simp [smallBirthSmokeCaps22, births]

theorem hphys (x : Fin 17) (i k : Fin 11) :
    physicalRows x i k = (physInt x.val i.val k.val : ℝ) / endpointDen :=
  rows_eq smallBirthSmokeCaps22 firstCaps returnCaps zeroHoldCaps endpointDen births hsmall x i k

theorem hrow (d : ℕ) (hd : d ≤ 16) (i k : Fin 11) :
    selectedRows ⟨d, by omega⟩ i k ≤ (row d i.val k.val : ℝ) / ((16 ^ 16 * endpointDen : ℤ) : ℝ) := by
  show turboThinnedRow22 ⟨d, by omega⟩ birthSmokeZ22 physicalRows i k ≤ _
  rw [selected_eq birthSmokeZ22 physicalRows physInt endpointDen (by decide +kernel) hphys d
    (by omega) i k]
  unfold row
  norm_num

set_option maxHeartbeats 0 in
theorem hRlit : R = pairedRowMatsL22 row := by decide +kernel

theorem hR (a : Fin 17) (i k : Fin 11) :
    selectedRows a i k ≤ (sparseLMat22 (R.getD a.val []) i.val k.val : ℝ) / ((16 ^ 16 * endpointDen : ℤ) : ℝ) := by
  rw [hRlit]
  exact PairedFracChain22.rowMats_hyp selectedRows row (16 ^ 16 * endpointDen) hrow a i k

theorem hden : (16 ^ 16 * endpointDen : ℤ) = {DEN} := by decide +kernel

theorem rpowChecked : pairedRpowABCheckL22 {ea} {eb} (16 ^ 16 * endpointDen) (2 ^ {S}) R Ys = true := by
  rw [hden]
  exact pairedRpowABCheck_of _ _ _ _ _ _ (by decide) fun a ha =>
    match a, ha with
{arms}
    | _ + 17, h => absurd h (by omega)

theorem hz : birthSmokeZ22 = (zNum : ℝ) / 2 ^ 80 := by
  unfold birthSmokeZ22
  norm_num [zDen]

-- Kernel-cheap C(N, q) (BinomFast22): N! / (q! (N-q)!) costs ~18 GB per check at N = 2^17.
def binomQ (q : ℕ) : ℤ := (binomFast {K.GROUPS} {K.GROUPS.bit_length()} q : ℕ)

theorem binomQ_eq (q : ℕ) (hq : q ≤ {K.GROUPS}) : (binomQ q : ℝ) = (Nat.choose {K.GROUPS} q : ℝ) := by
  unfold binomQ
  rw [binomFast_eq {K.GROUPS} {K.GROUPS.bit_length()} q hq (by norm_num)]
  simp

theorem binomQ_eq' (q : ℕ) (hq : q ≤ {K.GROUPS}) : (binomQ q : ℝ) = (Nat.choose (32 * {K.NPAIRS}) q : ℝ) :=
  binomQ_eq q hq

theorem binomQ_pos (q : ℕ) (hq : q ≤ {K.GROUPS}) : 0 < binomQ q := by
  have h : (0 : ℝ) < binomQ q := by rw [binomQ_eq q hq]; exact_mod_cast Nat.choose_pos hq
  exact_mod_cast h

end {ns}
'''
    w(f'{base}Pair', pair)


# ---------------------------------------------------------------- lookup tree
def tree(vals, ind='  '):
    """Balanced if-tree over sorted (q, a) pairs; leaves are literals."""
    if len(vals) == 1:
        return f'{vals[0][1]}'
    mid = len(vals) // 2
    pivot = vals[mid][0]
    return (f'(if q < {pivot} then {tree(vals[:mid], ind)} else {tree(vals[mid:], ind)})')


def stmt(q, ea, eb, aname):
    al = alpha_lean(ea, eb)
    return (f'(Nat.choose (32 * {K.NPAIRS}) {q} : ℝ) * ((turboOldUniformFactor22 ^ {q} / birthSmokeZ22 ^ {K.DCUT}) ^ ({al}) *' + NL +
            f'      ∑ k : Fin 11, ((pairedRegionDP22 selectedRows ({al}) {K.NPAIRS} {q}) ^ 64) 0 k) ≤ ({aname} : ℝ) / 2 ^ {KT}')


def emit_agg(p, pidx, ns, opens, imports, ea, eb, avals, term_prefix, chunk):
    qs = sorted(avals)
    lo, hi = qs[0], qs[-1]
    al = alpha_lean(ea, eb)
    blocks, names = [], []
    for c0 in range(0, len(qs), 100):
        cq = qs[c0:c0 + 100]
        a, b = cq[0], cq[-1]
        cases = NL.join(f'  · exact {term_prefix}{q}' for q in cq)
        blocks.append(f'''theorem pieceC{c0 // 100} (q : ℕ) (h1 : {a} ≤ q) (h2 : q ≤ {b}) :
    (Nat.choose (32 * {K.NPAIRS}) q : ℝ) * ((turboOldUniformFactor22 ^ q / birthSmokeZ22 ^ {K.DCUT}) ^ ({al}) *
      ∑ k : Fin 11, ((pairedRegionDP22 selectedRows ({al}) {K.NPAIRS} q) ^ 64) 0 k) ≤ (aP q : ℝ) / 2 ^ {KT} := by
  interval_cases q
{cases}
''')
        names.append((a, b, c0 // 100))
    disp = []
    for a, b, c in names[:-1]:
        disp.append(f'  by_cases hc{c} : q ≤ {b}{NL}  · exact pieceC{c} q (by omega) hc{c}')
    a, b, c = names[-1]
    disp.append(f'  exact pieceC{c} q (by omega) (by omega)')
    text = f'''{NL.join(imports)}

noncomputable section
namespace {ns}
open {opens}

def aP (q : ℕ) : Int := {tree([(q, avals[q]) for q in qs])}

{NL.join(blocks)}
theorem pieceTerm (q : ℕ) (hq : q ∈ Finset.Icc {lo} {hi}) :
    (Nat.choose (32 * {K.NPAIRS}) q : ℝ) * ((turboOldUniformFactor22 ^ q / birthSmokeZ22 ^ {K.DCUT}) ^ ({al}) *
      ∑ k : Fin 11, ((pairedRegionDP22 selectedRows ({al}) {K.NPAIRS} q) ^ 64) 0 k) ≤ (aP q : ℝ) / 2 ^ {KT} := by
  simp only [Finset.mem_Icc] at hq
  obtain ⟨hlo, hhi⟩ := hq
{NL.join(disp)}

#print axioms pieceTerm
end {ns}
'''
    w(f'P{pidx}Agg', text)


# ---------------------------------------------------------------- exact pieces
def step_modules(pidx, ns, opens, s, D, v, wbits, prev_mod, prev, sq_mod, per_entry):
    A = f'(pkMatL22 {wbits} {prev})'
    imps = [f'import {PRE}{prev_mod}', f'import {PRE}{sq_mod}', 'import SpinCodes.Turbo.SparseListPackedDefs22']
    rowmods = []
    for i in range(11):
        if per_entry:
            emods = []
            for j in range(11):
                name = f'P{pidx}St{s}R{i}E{j}'
                w(name, NL.join(imps) + f'{NL}{NL}namespace {ns}{NL}open {opens}{NL}{NL}set_option maxHeartbeats 0 in{NL}'
                  f'theorem step{s}_e{i}_{j} : pkEntryCheckL22 {D} {wbits} (2 ^ {S}) {A} {A} sq{s} {i} {j} = true := by{NL}'
                  f'  decide +kernel{NL}{NL}end {ns}{NL}')
                emods.append(f'import {PRE}{name}')
            body = NL.join(imps + emods)
        else:
            body = NL.join(imps)
        L = [body, '', f'namespace {ns}', f'open {opens}', '']
        if not per_entry:
            for j in range(11):
                L += ['set_option maxHeartbeats 0 in',
                      f'theorem step{s}_e{i}_{j} : pkEntryCheckL22 {D} {wbits} (2 ^ {S}) {A} {A} sq{s} {i} {j} = true := by',
                      '  decide +kernel', '']
        L += [f'theorem step{s}_r{i} : pkRowCheckL22 {D} {wbits} (2 ^ {S}) {A} {A} sq{s} {i} = true :=',
              '  SparseListPacked22.pkRowCheck_of_entries _ _ _ _ _ _ _ fun j hj =>', '    match j, hj with']
        L += [f'    | {j}, _ => step{s}_e{i}_{j}' for j in range(11)]
        L += ['    | _ + 11, h => absurd h (by omega)', '', f'end {ns}', '']
        w(f'P{pidx}St{s}R{i}', NL.join(L))
        rowmods.append(f'import {PRE}P{pidx}St{s}R{i}')
    L = imps + rowmods + ['', f'namespace {ns}', f'open {opens}', '',
                          'set_option maxHeartbeats 0 in',
                          f'theorem step{s}_hdr : pkHdrCheckL22 {D} {v} {wbits} {prev} {prev} sq{s} = true := by',
                          '  decide +kernel', '',
                          f'theorem step{s} : pkPolyMulCheckL22 {D} {v} {wbits} (2 ^ {S}) {prev} {prev} sq{s} = true :=',
                          f'  SparseListPacked22.pkPolyMulCheck_of_rows _ _ _ _ _ _ _ step{s}_hdr fun i hi =>',
                          '    match i, hi with']
    L += [f'    | {i}, _ => step{s}_r{i}' for i in range(11)]
    L += ['    | _ + 11, h => absurd h (by omega)', '', f'end {ns}', '']
    w(f'P{pidx}St{s}', NL.join(L))


def packed_vw(polys, D):
    v = max(max(x.bit_length() for M in P for r in M for x in r) for P in polys)
    wb = 1
    while not 11 * (2 * D + 1) * 4 ** v < 2 ** wb:
        wb += 1
    return v, wb


def emit_exact(p, pidx, key):
    ex = json.loads((CACHE / f'exact{pidx}.json').read_text())
    f, ea, eb, D = ex['f'], ex['ea'], ex['eb'], ex['D']
    t = ptag(key)
    pairns = f'{K.fam_ns(f)}.A{t}'
    ns = f'{K.fam_ns(f)}.P{pidx}'
    opens = f'Spin.Turbo Spin.Turbo.Hi {K.fam_ns(f)} {pairns}'
    core = 'import SpinCodes.Turbo.SparseListDefs22'
    w(f'P{pidx}Base', f'{core}{NL}{NL}namespace {ns}{NL}{poly_defs("base", ex["base"])}{NL}end {ns}{NL}')
    for s, P in enumerate(ex['stages']):
        w(f'P{pidx}Sq{s}', f'{core}{NL}{NL}namespace {ns}{NL}{poly_defs(f"sq{s}", P)}{NL}end {ns}{NL}')
    v, wb = packed_vw([ex['base']] + ex['stages'], D)
    per_entry = D > 300
    for s in range(len(ex['stages'])):
        prev_mod = f'P{pidx}Base' if s == 0 else f'P{pidx}Sq{s - 1}'
        prev = 'base' if s == 0 else f'sq{s - 1}'
        step_modules(pidx, ns, 'Spin.Turbo', s, D, v, wb, prev_mod, prev, f'P{pidx}Sq{s}', per_entry)
    nch = (D + 1 + BASE_CH - 1) // BASE_CH
    for c in range(nch):
        parts = ', '.join(f'base_part{i}' for i in range(c * BASE_CH, min(D + 1, (c + 1) * BASE_CH)))
        w(f'P{pidx}BaseC{c}', f'import {PRE}P{pidx}Base{NL}import {PRE}{f}A{t}Ys{NL}import SpinCodes.Turbo.PairedFracChain22'
          f'{NL}{NL}namespace {ns}{NL}open Spin.Turbo {pairns}{NL}{NL}def baseC{c} : List (List (List Int)) := [{parts}]{NL}{NL}'
          f'set_option maxHeartbeats 0 in{NL}theorem baseC{c}_ok : pairedBaseGoL22 Ys {c * BASE_CH} baseC{c} = true := by{NL}'
          f'  decide +kernel{NL}{NL}end {ns}{NL}')
    nst = len(ex['stages'])
    appended = ' ++ '.join(f'baseC{c}' for c in range(nch))
    steps = ', '.join(f'step{s}' for s in range(nst))
    trace = '[' + ', '.join(f'(false, sq{s}, [])' for s in range(nst)) + ']'
    if nch == 1:
        base_tail = '  exact baseC0_ok'
    else:
        base_tail = ('  simp only [List.append_assoc, PairedFracChain22.baseGo_append, Bool.and_eq_true]' + NL +
                     '  exact ⟨' + ', '.join(f'baseC{c}_ok' for c in range(nch)) + '⟩')
    chain = NL.join([f'import {PRE}{f}A{t}Pair', f'import {PRE}P{pidx}Base'] +
                    [f'import {PRE}P{pidx}St{s}' for s in range(nst)] +
                    [f'import {PRE}P{pidx}BaseC{c}' for c in range(nch)] +
                    ['import SpinCodes.Turbo.PairedFracABChain22']) + f'''

noncomputable section
namespace {ns}
open {opens}

def trace : List (Bool × List (List (List Int)) × List (List (List Int))) := {trace}

theorem baseChecked : pairedBaseCheckDL22 {D} Ys base = true := by
  have hb : base = {appended} := by decide +kernel
  simp only [pairedBaseCheckDL22, Bool.and_eq_true, decide_eq_true_eq]
  refine ⟨by decide +kernel, ?_⟩
  rw [hb]
{base_tail}

theorem chainChecked : pkRegionalPowerCheckL22 {D} {v} {wb} (2 ^ {S}) base base trace = true := by
  simp only [trace, pkRegionalPowerCheckL22, {steps}, Bool.false_eq_true,
    ↓reduceIte, Bool.and_self]

end {ns}
'''
    w(f'P{pidx}Chain', chain)
    ztr = '[' + ', '.join(f'({str(b).lower()}, {sh}, {mm})' for b, sh, mm in ex['ztr']) + ']'
    recs = {int(q): r for q, r in ex['recs'].items()}
    qs = sorted(recs)
    termimps = []
    for c0 in range(0, len(qs), QCHUNK_EXACT):
        cq = qs[c0:c0 + QCHUNK_EXACT]
        c = c0 // QCHUNK_EXACT
        qd = []
        for q in cq:
            r = recs[q]
            st = '[' + ', '.join(f'({dn}, {up}, {lm(U)})' for dn, up, U in r['st']) + ']'
            vt = '[' + ', '.join(f'({str(b).lower()}, {sh}, {mm})' for b, sh, mm in r['vtr']) + ']'
            qd.append(f'def U0_{q} : List (List Int) := {lm(r["U0"])}{NL}def st_{q} : List (Nat × Nat × List (List Int)) := {st}{NL}'
                      f'def vt_{q} : List (Bool × Nat × Int) := {vt}{NL}def G_{q} : Int := {r["G"]}{NL}'
                      f'def yn_{q} : Int := {r["yn"]}{NL}def a_{q} : Int := {r["a"]}{NL}')
        w(f'P{pidx}QD{c}', f'{core}{NL}{NL}namespace {ns}{NL}' + NL.join(qd) + f'end {ns}{NL}')
        terms = []
        for q in cq:
            terms.append(f'''set_option maxHeartbeats 0 in
theorem eterm_{q} :
    {stmt(q, ea, eb, f'a_{q}')} :=
  TailTermFrac22.tail_term_fracAB_shift_le {ea} {eb} (by norm_num) {K.NPAIRS} {K.DCUT} 80 {v} {wb} selectedRows
    selectedRows_nonneg {D} {q} (by norm_num) (2 ^ {S}) {S} rfl (16 ^ 16 * endpointDen) R hR Ys rpowChecked base
    baseChecked trace (by decide +kernel) chainChecked (binomQ {q}) (binomQ_eq' {q} (by norm_num))
    (binomQ_pos {q} (by norm_num)) U0_{q} (by decide +kernel) {recs[q]['dn0']} {recs[q]['up0']} (by decide +kernel)
    st_{q} rfl (by decide +kernel) birthSmokeZ22 zNum hz ztrace
    zChecked vt_{q} (by decide +kernel) yn_{q} G_{q} (by decide +kernel) a_{q} {KT} (by decide +kernel)
''')
        w(f'P{pidx}T{c}', NL.join([f'import {PRE}P{pidx}Chain', f'import {PRE}P{pidx}Z', f'import {PRE}P{pidx}QD{c}']) +
          f'{NL}{NL}noncomputable section{NL}namespace {ns}{NL}open {opens}{NL}{NL}' + NL.join(terms) + f'{NL}end {ns}{NL}')
        termimps.append(f'import {PRE}P{pidx}T{c}')
    w(f'P{pidx}Z', f'import {PRE}{f}A{t}Pair{NL}{NL}namespace {ns}{NL}open {opens}{NL}{NL}def ztrace : List (Bool × Nat × Int) := {ztr}'
      f'{NL}{NL}set_option maxHeartbeats 0 in{NL}theorem zChecked : sparseScaledPowerCheck22 zNum 80 {K.DCUT} ztrace = true := by'
      f'{NL}  decide +kernel{NL}{NL}end {ns}{NL}')
    emit_agg(p, pidx, ns, opens, termimps, ea, eb, {q: recs[q]['a'] for q in qs}, 'eterm_', QCHUNK_EXACT)


# ---------------------------------------------------------------- marker pieces
def emit_marker(p, pidx, key):
    mk = json.loads((CACHE / f'marker{pidx}.json').read_text())
    if p.get('range'):
        return emit_marker_range(p, pidx, key, mk)
    f, ea, eb = mk['f'], mk['ea'], mk['eb']
    t = ptag(key)
    pairns = f'{K.fam_ns(f)}.A{t}'
    ns = f'{K.fam_ns(f)}.P{pidx}'
    opens = f'Spin.Turbo Spin.Turbo.Hi {K.fam_ns(f)} {pairns}'
    tn0, td0, kp, kn, up, down = mk['tn0'], mk['td0'], mk['kp'], mk['kn'], mk['up'], mk['down']
    st = '[' + ', '.join(f'({dn}, {u}, {lm(U)})' for dn, u, U in mk['st']) + ']'
    ztr = '[' + ', '.join(f'({str(b).lower()}, {sh}, {mm})' for b, sh, mm in mk['ztr']) + ']'
    w(f'P{pidx}MBase', f'''import {PRE}{f}A{t}Pair

-- Data defs are only consumed by `decide +kernel`; skip codegen, which hits
-- maxRecDepth on the 18-stage chains at G24+ (no limit raised).
noncomputable section
namespace {ns}
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

end {ns}
''')
    recs = {int(q): r for q, r in mk['recs'].items()}
    qs = sorted(recs)
    termimps = []
    for c0 in range(0, len(qs), QCHUNK_MARK):
        cq = qs[c0:c0 + QCHUNK_MARK]
        c = c0 // QCHUNK_MARK
        terms = []
        for q in cq:
            r = recs[q]
            vt = '[' + ', '.join(f'({str(b).lower()}, {sh}, {mm})' for b, sh, mm in r['vtr']) + ']'
            terms.append(f'''def vt_{q} : List (Bool × Nat × Int) := {vt}

set_option maxHeartbeats 0 in
theorem mterm_{q} :
    {stmt(q, ea, eb, f'({r["a"]} : ℤ)')} :=
  TailMarkerFracAB22.tail_marker_fracAB_le {ea} {eb} (by norm_num) {K.NPAIRS} {K.DCUT} 80 selectedRows
    selectedRows_nonneg {q} (binomQ {q}) (binomQ_eq' {q} (by norm_num)) (binomQ_pos {q} (by norm_num))
    (16 ^ 16 * endpointDen) (2 ^ {S}) (by norm_num) R hR Ys rpowChecked {tn0} {td0} (by norm_num) (by norm_num)
    {kp} {kn} (2 ^ {S}) {S} rfl {down} {up} U0 baseChecked stages (by decide) chainChecked birthSmokeZ22 zNum hz
    ztrace zChecked vt_{q} (by decide +kernel) ({r["yn"]}) ({r["G"]}) (by decide +kernel) ({r["a"]}) {KT}
    (by decide +kernel)
''')
        w(f'P{pidx}MQ{c}', f'import {PRE}P{pidx}MBase{NL}{NL}noncomputable section{NL}namespace {ns}{NL}open {opens}{NL}{NL}'
          + NL.join(terms) + f'{NL}end {ns}{NL}')
        termimps.append(f'import {PRE}P{pidx}MQ{c}')
    emit_agg(p, pidx, ns, opens, termimps, ea, eb, {q: recs[q]['a'] for q in qs}, 'mterm_', QCHUNK_MARK)


def emit_marker_range(p, pidx, key, mk):
    """Marker range piece: two endpoint certificates of F (TailMarkerRange22.markerF_le)
    at a common bound A, then term_le_endpoints for every q in [qlo, qhi]."""
    f, ea, eb = mk['f'], mk['ea'], mk['eb']
    t = ptag(key)
    pairns = f'{K.fam_ns(f)}.A{t}'
    ns = f'{K.fam_ns(f)}.P{pidx}'
    opens = f'Spin.Turbo Spin.Turbo.Hi {K.fam_ns(f)} {pairns}'
    tn0, td0, kp, kn, up, down = mk['tn0'], mk['td0'], mk['kp'], mk['kn'], mk['up'], mk['down']
    st = '[' + ', '.join(f'({dn}, {u}, {lm(U)})' for dn, u, U in mk['st']) + ']'
    ztr = '[' + ', '.join(f'({str(b).lower()}, {sh}, {mm})' for b, sh, mm in mk['ztr']) + ']'
    recs = {int(q): r for q, r in mk['recs'].items()}
    lo, hi = p['qlo'], p['qhi']
    A = max(recs[lo]['a'], recs[hi]['a'])
    tex = f'((({tn0} * 2 ^ {kp} : ℤ) : ℝ) / (({td0} * 2 ^ {kn} : ℤ) : ℝ))'
    al = alpha_lean(ea, eb)
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
    w(f'P{pidx}MBase', f"""import {PRE}{f}A{t}Pair
import SpinCodes.Turbo.TailMarkerRange22

-- Data defs are only consumed by `decide +kernel`; skip codegen, which hits
-- maxRecDepth on the 18-stage chains at G24+ (no limit raised).
noncomputable section
namespace {ns}
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

end {ns}
""")
    w(f'P{pidx}MQ0', f"""import {PRE}P{pidx}MBase

noncomputable section
namespace {ns}
open {opens}

{cert(lo, 'Flo')}
{cert(hi, 'Fhi') if hi != lo else ''}
end {ns}
""")
    fhi = 'Fhi' if hi != lo else 'Flo'
    w(f'P{pidx}Agg', f"""import {PRE}P{pidx}MQ0

noncomputable section
namespace {ns}
open {opens}

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
    birthSmokeZ22 {tex} zpos ht {lo} {hi} (by norm_num) q h1 h2).trans (max_le Flo {fhi})

#print axioms pieceTerm
end {ns}
""")


def emit_all(ps):
    fidx = {t: i for i, t in enumerate(K.tilts(ps), start=1)}
    keys = []
    for p in ps:
        key = K.pair_key(p, fidx[p['tilt']])
        if key not in keys:
            keys.append(key)
    for key in keys:
        emit_pairs(key)
    man = []
    for pidx, p in enumerate(ps, start=1):
        key = K.pair_key(p, fidx[p['tilt']])
        if p['kind'] == 'exact':
            if (CACHE / f'exact{pidx}.json').exists():
                emit_exact(p, pidx, key)
        else:
            if (CACHE / f'marker{pidx}.json').exists():
                emit_marker(p, pidx, key)
        ea, eb = K.alpha_ab(p['alpha'])
        man.append(dict(piece=pidx, kind=p['kind'], range=bool(p.get('range')), f=fidx[p['tilt']], ea=ea, eb=eb, qlo=p['qlo'], qhi=p['qhi'],
                        ns=f'{K.fam_ns(fidx[p["tilt"]])}.P{pidx}'))
    (OUT / 'pieces.json').write_text(json.dumps(man, indent=1) + NL)
    print('emitted', len(man), 'pieces')
