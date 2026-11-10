#!/usr/bin/env python3
"""Switch K18 exact-family squaring chains to packed (Kronecker) checks.

The interpreted degree-D chain check (`sparseRegionalPowerCheckDL22`) costs
~2 s of kernel time per 11x11 matrix product and needs (D+1)(D+2)/2 products
per squaring; it is infeasible for D >= 10 (family 2 hit the kernel memory
limit after 70 minutes). `SparseListPacked22.pkPolyMulCheckL22` checks a whole
truncated product with 1331 GMP multiplications of base-2^w packed entries.

For each exact family k (2..6, 8..11) this script
  * computes v = max entry bit length over base and all squaring stages and
    the least w with 11 (2D+1) 4^v < 2^w;
  * writes Tail{k}ExStep{s}.lean: one kernel check per squaring step
    (separate declarations/modules keep kernel caches small and parallel);
  * rewrites `chainChecked` in Tail{k}ExChain.lean as the packed chain check,
    assembled from the step lemmas, and the ExTerm call to `chainPk_sound`.
Untrusted: Lean re-checks every step.
"""
import json
import re
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
GEN = ROOT / 'SpinCodes/Turbo/Generated/K18Tail'
FAMS = [2, 3, 4, 5, 6, 8, 9, 10, 11]
STEPS = 7
INT_RE = re.compile(r'-?\d+')


def max_bits(paths):
    m = 0
    for p in paths:
        for line in p.read_text(encoding='utf-8').splitlines():
            if not line.startswith('def '):
                continue
            body = line.split(':=', 1)[1]
            for t in INT_RE.findall(body):
                x = int(t)
                assert x >= 0, (p, t[:30])
                m = max(m, x.bit_length())
    return m


def step_text(imports, ns, s, D, v, w, prev, rowmods):
    """One squaring step: header check + assembly of 11 row modules. Lean does
    not return kernel memory between declarations of one process, so each row
    (11 per-entry declarations) is its own module/process."""
    L = list(imports) + [f'import {m}' for m in rowmods] + ['', f'namespace {ns}', '']
    L += ['set_option maxHeartbeats 0 in',
          f'theorem step{s}_hdr : pkHdrCheckL22 {D} {v} {w} {prev} {prev} sq{s} = true := by',
          '  decide +kernel', '']
    L += [f'theorem step{s} : pkPolyMulCheckL22 {D} {v} {w} (2 ^ 704) {prev} {prev} sq{s} = true :=',
          f'  SparseListPacked22.pkPolyMulCheck_of_rows _ _ _ _ _ _ _ step{s}_hdr fun i hi =>',
          '    match i, hi with']
    L += [f'    | {i}, _ => step{s}_r{i}' for i in range(11)]
    L += ['    | _ + 11, h => absurd h (by omega)', '',
          f'#print axioms step{s}', f'end {ns}', '']
    return chr(10).join(L)


def row_text(imports, ns, s, D, w, prev, i):
    A = f'(pkMatL22 {w} {prev})'
    L = list(imports) + ['', f'namespace {ns}', '']
    for j in range(11):
        L += ['set_option maxHeartbeats 0 in',
              f'theorem step{s}_e{i}_{j} : pkEntryCheckL22 {D} {w} (2 ^ 704) {A} {A} sq{s} {i} {j} = true := by',
              '  decide +kernel', '']
    L += [f'theorem step{s}_r{i} : pkRowCheckL22 {D} {w} (2 ^ 704) {A} {A} sq{s} {i} = true :=',
          f'  SparseListPacked22.pkRowCheck_of_entries _ _ _ _ _ _ _ fun j hj =>',
          '    match j, hj with']
    L += [f'    | {j}, _ => step{s}_e{i}_{j}' for j in range(11)]
    L += ['    | _ + 11, h => absurd h (by omega)', '', f'end {ns}', '']
    return chr(10).join(L)


def write_step(prefix, k, s, ns, D, v, w, prev_mod, prev):
    base_imports = [f'import SpinCodes.Turbo.Generated.K18Tail.{prev_mod}',
                    f'import SpinCodes.Turbo.Generated.K18Tail.Tail{k}{prefix}Sq{s}',
                    'import SpinCodes.Turbo.SparseListPackedDefs22']
    rowmods = []
    for i in range(11):
        name = f'Tail{k}{prefix}Step{s}R{i}'
        (GEN / f'{name}.lean').write_text(row_text(base_imports, ns, s, D, w, prev, i),
                                          encoding='utf-8', newline=chr(10))
        rowmods.append(f'SpinCodes.Turbo.Generated.K18Tail.{name}')
    (GEN / f'Tail{k}{prefix}Step{s}.lean').write_text(
        step_text(base_imports, ns, s, D, v, w, prev, rowmods), encoding='utf-8', newline=chr(10))


def main():
    man = json.loads((GEN / 'exact_manifest.json').read_text())
    out = {}
    for k in FAMS:
        D = int(man[str(k)]['D'])
        ns = f'Spin.Turbo.Hi.ProfileFamilyK18Tail{k}22.Ex'
        data = [GEN / f'Tail{k}ExBase.lean'] + [GEN / f'Tail{k}ExSq{s}.lean' for s in range(STEPS)]
        v = max_bits(data)
        w = 1
        while not 11 * (2 * D + 1) * 4 ** v < 2 ** w:
            w += 1
        out[k] = {'D': D, 'v': v, 'w': w}
        for s in range(STEPS):
            prev_mod = f'Tail{k}ExBase' if s == 0 else f'Tail{k}ExSq{s - 1}'
            prev = 'base' if s == 0 else f'sq{s - 1}'
            write_step('Ex', k, s, ns, D, v, w, prev_mod, prev)
        # chain module
        cp = GEN / f'Tail{k}ExChain.lean'
        c = cp.read_text(encoding='utf-8')
        old = re.search(r'set_option maxHeartbeats 0 in\ntheorem chainChecked : '
                        r'sparseRegionalPowerCheckDL22 .*?\n  decide \+kernel\n', c)
        if old:
            steps = ', '.join(f'step{s}' for s in range(STEPS))
            new = (f'theorem chainChecked : pkRegionalPowerCheckL22 {D} {v} {w} (2 ^ 704) base base trace = true := by\n'
                   f'  simp only [trace, pkRegionalPowerCheckL22, {steps}, Bool.false_eq_true,\n'
                   f'    ↓reduceIte, Bool.and_self]\n')
            c = c[:old.start()] + new + c[old.end():]
            imports = ''.join(f'import SpinCodes.Turbo.Generated.K18Tail.Tail{k}ExStep{s}\n'
                              for s in range(STEPS))
            c = imports + c
            cp.write_text(c, encoding='utf-8', newline='\n')
        else:
            assert 'pkRegionalPowerCheckL22' in c, k
        tp = GEN / f'Tail{k}ExTerm.lean'
        t = tp.read_text(encoding='utf-8')
        if 'import SpinCodes.Turbo.SparseListPacked22' not in t:
            t = 'import SpinCodes.Turbo.SparseListPacked22' + chr(10) + t
        t2 = t.replace(f'SparseListPolyD22.chainD_sound {D} selectedRows',
                       f'SparseListPacked22.chainPk_sound {D} {v} {w} selectedRows')
        assert t2 != t or 'chainPk_sound' in t, k
        tp.write_text(t2, encoding='utf-8', newline='\n')
        print(k, out[k])
    out['7'] = frac7()
    (GEN / 'packed_manifest.json').write_text(json.dumps(out, indent=1))


def frac7():
    """Fractional family 7 (alpha = 4/5): 6 squarings at D = 143."""
    k, D, steps = 7, 143, 6
    ns = f'Spin.Turbo.Hi.ProfileFamilyK18Tail{k}22.Fr'
    data = [GEN / f'Tail{k}FrBase.lean'] + [GEN / f'Tail{k}FrSq{s}.lean' for s in range(steps)]
    v = max_bits(data)
    w = 1
    while not 11 * (2 * D + 1) * 4 ** v < 2 ** w:
        w += 1
    for s in range(steps):
        prev_mod = f'Tail{k}FrBase' if s == 0 else f'Tail{k}FrSq{s - 1}'
        prev = 'base' if s == 0 else f'sq{s - 1}'
        write_step('Fr', k, s, ns, D, v, w, prev_mod, prev)
    cp = GEN / f'Tail{k}FrChain.lean'
    c = cp.read_text(encoding='utf-8')
    old = (f'theorem chainChecked : sparseRegionalPowerCheckDL22 {D} (2 ^ 704) base base trace = true '
           f':= by decide +kernel\n')
    if old in c:
        stp = ', '.join(f'step{s}' for s in range(steps))
        c = c.replace('set_option maxHeartbeats 0 in\n' + old,
                      f'theorem chainChecked : pkRegionalPowerCheckL22 {D} {v} {w} (2 ^ 704) base base trace = true := by\n'
                      f'  simp only [trace, pkRegionalPowerCheckL22, {stp}, Bool.false_eq_true,\n'
                      f'    ↓reduceIte, Bool.and_self]\n')
        c = ''.join(f'import SpinCodes.Turbo.Generated.K18Tail.Tail{k}FrStep{s}\n' for s in range(steps)) + c
        cp.write_text(c, encoding='utf-8', newline='\n')
    assert 'pkRegionalPowerCheckL22' in cp.read_text(encoding='utf-8')
    tp = GEN / f'Tail{k}FrTerm.lean'
    t = tp.read_text(encoding='utf-8')
    if 'tail_term_frac_pk_le' not in t:
        t = t.replace('TailTermFrac22.tail_term_frac_le selectedRows',
                      f'TailTermFrac22.tail_term_frac_pk_le {v} {w} selectedRows')
        t = 'import SpinCodes.Turbo.TailTermFracPk22\n' + t
        tp.write_text(t, encoding='utf-8', newline='\n')
    print(k, {'D': D, 'v': v, 'w': w})
    return {'D': D, 'v': v, 'w': w}


if __name__ == '__main__':
    main()
