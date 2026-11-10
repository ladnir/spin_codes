"""Emit small, untrusted birth candidates for kernel checking (240+27000 inputs)."""
import json
from pathlib import Path

BASE = Path(__file__).resolve().parents[1]
candidate = json.loads((BASE / '.lake/turbo/coarse_packet_profiles.json').read_text())
rows = [int(v, 16) for v in candidate['expansion_rows_hex']]
def linear(rs, v):
    out = 0
    for k, row in enumerate(rs):
        if v >> k & 1:
            out ^= row
    return out
cols = [sum(((row >> j) & 1) << i for i, row in enumerate(rows)) for j in range(64)]
images = [linear(rows, linear(cols, v << (4*s))) for s in range(16) for v in range(1,16)]
levels = {0:0,16:3,24:4,28:5,32:6,36:7,40:8,48:9,64:10}
def code(image, weight):
    return 9*levels[image.bit_count()] + weight
one = [0]*99
shards = []
pairs = {}
for s in range(16):
    shard = [0]*99
    for u in range(1,16):
        one[code(images[15*s+u-1],u.bit_count())] += 1
    for t in range(s+1,16):
        pair = [0]*99
        for u in range(1,16):
            for v in range(1,16):
                k = code(images[15*s+u-1]^images[15*t+v-1],u.bit_count()+v.bit_count())
                shard[k] += 1
                pair[k] += 1
        pairs[s,t] = pair
    shards.append(shard)
two = [sum(shard[i] for shard in shards) for i in range(99)]
for record, hist in zip(candidate['small_births'], (one,two)):
    for level, counts in zip(record['expansion_weight_levels'],record['coefficients']):
        assert counts == hist[9*levels[level]:9*levels[level]+len(counts)]
def lean_list(xs):
    return '[' + ', '.join(map(str,xs)) + ']'
out = ['import SpinCodes.Structured.MapSpectrumDefs',
       'import Init.Data.List.OfFn', '',
       'namespace Spin.Turbo', 'open Spin.Structured', '',
       'def smallBirthClassWeight22 (w : Nat) : Fin 11 :=',
       '  if w = 0 then 0 else if w = 16 then 3 else if w = 24 then 4 else',
       '  if w = 28 then 5 else if w = 32 then 6 else if w = 36 then 7 else',
       '  if w = 40 then 8 else if w = 48 then 9 else if w = 64 then 10 else 1', '',
       'def smallBirthSingleImages22 : List Nat := ' + lean_list(images), '',
       'def smallBirthOneHistogram22 : List Nat := ' + lean_list(one), '',
       'def smallBirthTwoHistogram22 : List Nat := ' + lean_list(two), '']
for i,hist in enumerate(shards):
    out += [f'def smallBirthTwoShard{i}Histogram22 : List Nat := ' + lean_list(hist), '']
for (s,t),hist in pairs.items():
    out += [f'def smallBirthPair{s}_{t}Histogram22 : List Nat := ' + lean_list(hist), '']
out += ['def smallBirthPairHistograms22 (a b : Fin 16) : List Nat :=', '  (([']
for s in range(16):
    out += ['    [' + ', '.join(f'smallBirthPair{s}_{t}Histogram22' if s<t else '[]' for t in range(16)) + ']' + (',' if s<15 else '')]
out += ['  ] : List (List (List Nat))).getD a.val []).getD b.val []', '']
out += ['end Spin.Turbo','']
(BASE / 'SpinCodes/Turbo/PhysicalSmallBirth22NumericData.lean').write_text('\n'.join(out), encoding='utf-8')
semantic_data = ['import SpinCodes.Turbo.PhysicalSmallBirth22Encoding', '',
                 'namespace Spin.Turbo', 'open Spin.UniformMixer Spin.Structured', '',
                 'set_option maxRecDepth 100000 in', 'set_option maxHeartbeats 1000000 in',
                 'theorem smallBirthSingleImages22_checked : ∀ a : SmallBirthOne22,',
                 '    PackedMap.eval productionPackedRows (PackedMap.eval smallBirthFeedbackRows22',
                 '      (smallBirthOnePacked22 a)) = smallBirthSingleImages22.getD (15*a.1.val+a.2.val) 0 := by',
                 '  decide +kernel', '', 'end Spin.Turbo', '']
(BASE / 'SpinCodes/Turbo/PhysicalSmallBirth22Data.lean').write_text('\n'.join(semantic_data), encoding='utf-8')
for s in range(15):
    chunk = ['import SpinCodes.Turbo.PhysicalSmallBirth22NumericCensusDefs',
             'import Lean.Elab.Tactic.Decide', '',
             'namespace Spin.Turbo', '', 'set_option maxRecDepth 100000',
             'set_option maxHeartbeats 2000000', '']
    if s == 0:
        chunk += ['theorem smallBirthOneHistogram22_checked :',
                  '    smallBirthCodeHistogram22 smallBirthOneCodes22 = smallBirthOneHistogram22 := by',
                  '  decide +kernel', '']
    for t in range(s+1,16):
        chunk += [f'theorem smallBirthPair{s}_{t}Histogram22_checked :',
                  f'    smallBirthCodeHistogram22 (smallBirthPairCodes22 {s} {t}) =',
                  f'      smallBirthPair{s}_{t}Histogram22 := by decide +kernel', '']
    chunk += ['end Spin.Turbo', '']
    (BASE / f'SpinCodes/Turbo/PhysicalSmallBirth22PairChecks{s}.lean').write_text('\n'.join(chunk), encoding='utf-8')
pair_checks = [f'import SpinCodes.Turbo.PhysicalSmallBirth22PairChecks{s}' for s in range(15)] + ['']
(BASE / 'SpinCodes/Turbo/PhysicalSmallBirth22PairChecks.lean').write_text('\n'.join(pair_checks), encoding='utf-8')
checks = ['import SpinCodes.Turbo.PhysicalSmallBirth22PairChecks',
          'import Mathlib.Data.Fin.VecNotation',
          'import Mathlib.Algebra.BigOperators.Fin',
          'import Mathlib.Tactic.FinCases', 'import Mathlib.Tactic.NormNum', '',
          'namespace Spin.Turbo', 'set_option maxRecDepth 100000',
          'set_option maxHeartbeats 2000000', '']
checks += ['theorem smallBirthPairHistograms22_checked (a b : Fin 16) (h : a < b) :',
           '    smallBirthCodeHistogram22 (smallBirthPairCodes22 a b) =',
           '      smallBirthPairHistograms22 a b := by', '  fin_cases a']
for s in range(16):
    checks += ['  · fin_cases b']
    for t in range(16):
        checks += [f'    · exact smallBirthPair{s}_{t}Histogram22_checked' if s<t else '    · norm_num at h']
checks += ['def smallBirthTwoShardHistograms22 : Fin 16 → List Nat := ![',
           ', '.join(f'smallBirthTwoShard{i}Histogram22' for i in range(16)) + ']', '']
for s in range(16):
    checks += [f'theorem smallBirthTwoShard{s}Sum22_checked : ∀ k : Fin 99,',
               f'    smallBirthTwoShard{s}Histogram22.getD k.val 0 =',
               f'      ∑ b : Fin 16, if ({s} : Fin 16) < b then',
               f'        (smallBirthPairHistograms22 {s} b).getD k.val 0 else 0 := by decide +kernel', '']
checks += [
           'theorem smallBirthTwoShardHistograms22_sum_checked :',
           '    ∀ (a : Fin 16) (k : Fin 99),',
           '      (smallBirthTwoShardHistograms22 a).getD k.val 0 =',
           '        ∑ b : Fin 16, if a < b then (smallBirthPairHistograms22 a b).getD k.val 0 else 0 := by',
           '  intro a k', '  fin_cases a']
checks += [f'  · exact smallBirthTwoShard{s}Sum22_checked k' for s in range(16)]
checks += ['', 'theorem smallBirthTwoHistogram22_sum_checked :',
           '    (List.ofFn (fun b : Fin 99 => ∑ a : Fin 16,',
           '      (smallBirthTwoShardHistograms22 a).getD b.val 0)) =',
           '        smallBirthTwoHistogram22 := by decide +kernel', '', 'end Spin.Turbo', '']
(BASE / 'SpinCodes/Turbo/PhysicalSmallBirth22Checks.lean').write_text('\n'.join(checks), encoding='utf-8')
print('Emitted 240 untrusted image candidates and 120 pair histogram candidates.')

