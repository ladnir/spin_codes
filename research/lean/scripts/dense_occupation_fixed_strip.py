from pathlib import Path
import json
r=Path('.')
ns=['First']+[f'W{i:03d}' for i in range(1,9)]
cs=[372]+[json.load(open(f'scripts/map_data/dense_occupation_fixed_{tag}_candidate.json'))['prefactor'] for tag in ns[1:]]
lines=[f'import SpinCodes.Structured.DenseOccupationFixed{tag}Box' for tag in ns]
lines+=['','noncomputable section','namespace Spin.Structured.DenseOccupationFixed.Strip','open Spin.Numeric Spin.Imt Set','', '''/-- A concrete contraction witness for this affine outer segment, with a common finite prefactor. -/
structure PointCertificate (α x : ℝ) where
  p : ℝ
  y : ℝ
  z : ℝ
  radius : ℝ
  w : Coords 5
  p_pos : 0<p
  p_lt_one : p<1
  y_pos : 0<y
  y_lt_one : y<1
  z_pos : 0<z
  z_lt_one : z<1
  radius_pos : 0<radius
  radius_lt_one : radius<1
  zero_eq_one : w.Z=1
  diffuse_floor : (1:ℝ)/568≤w.D
  shell_floor : ∀ i, (1:ℝ)/568≤w.S i
  collatz : ((Occupation.Sparse.numericalMatrix (p*y) z).applyCol w).le (Coords.smul radius w)
  exponent_le : boxExponent (833/500) (-43311/250000) p y radius z α x≤-(4/10000000)
''']
for tag,C in zip(ns,cs):
 box=tag+'Box'
 lines.append(f'''def witness{tag} {{α x : ℝ}}
    (hα : α∈Icc {box}.a0.real {box}.a1.real)
    (hx : x∈Icc {box}.x0.real {box}.x1.real) : PointCertificate α x where
  p := {box}.p.real
  y := {box}.y.real
  z := {box}.z.real
  radius := {box}.radius.real
  w := {tag}.w
  p_pos := by norm_num [{box}.p,QInput.real]
  p_lt_one := by norm_num [{box}.p,QInput.real]
  y_pos := by norm_num [{box}.y,QInput.real]
  y_lt_one := by norm_num [{box}.y,QInput.real]
  z_pos := by rw [{box}.parameters_match.2.2]; exact {tag}.parameters.2.2.1
  z_lt_one := by rw [{box}.parameters_match.2.2]; exact {tag}.parameters.2.2.2.1
  radius_pos := by rw [{box}.parameters_match.2.1]; exact {tag}.parameters.2.2.2.2.1
  radius_lt_one := by rw [{box}.parameters_match.2.1]; exact {tag}.parameters.2.2.2.2.2
  zero_eq_one := by norm_num [{tag}.w,{tag}.v,Fix.sc,scale]
  diffuse_floor := (show (1:ℝ)/568≤1/{C} by norm_num).trans {tag}.witness_floor.2.1
  shell_floor := fun i => (show (1:ℝ)/568≤1/{C} by norm_num).trans ({tag}.witness_floor.2.2 i)
  collatz := {box}.collatz
  exponent_le := by
    simpa only [{box}.m,{box}.c,QInput.real] using {box}.exponent_bound hα hx
''')
lines+=['''/-- Exact coverage of the first nine adjacent boxes, with one of nine checked witnesses at every point. -/
theorem covered {α x : ℝ} (hα : α∈Icc (1/10000) (10031/320000))
    (hx : x∈Icc (13/125) (18296026121/120000000000)) : Nonempty (PointCertificate α x) := by''']
for tag in ns[:-1]:
 box=tag+'Box'
 lines+=[f'  by_cases h{tag} : α≤{box}.a1.real',f'  · apply Nonempty.intro (witness{tag} (α := α) (x := x) ?_ ?_)',f'    · constructor',f'      · norm_num [{box}.a0,QInput.real] at *', '        linarith',f'      · exact h{tag}',f'    · simpa [{box}.x0,{box}.x1,QInput.real] using hx']
# the negative h tags need normalizing a1 definitions before next branch linarith
# insert all a1 normalization in each lower goal because norm_num at * without definitions cannot unfold them
for idx,line in enumerate(lines):
 if line.startswith('      · norm_num [') and 'a0' in line:
  defs=','.join(f'{t}Box.a1' for t in ns[:-1])
  lines[idx]=line.replace(',QInput.real]',','+defs+',QInput.real]')
tag=ns[-1];box=tag+'Box';defs=','.join(f'{t}Box.a1' for t in ns[:-1])
lines += [f'  apply Nonempty.intro (witness{tag} (α := α) (x := x) ?_ ?_)',f'  · constructor',f'    · norm_num [{box}.a0,{defs},QInput.real] at *', '      linarith',f'    · simpa [{box}.a1,QInput.real] using hα.2',f'  · simpa [{box}.x0,{box}.x1,QInput.real] using hx','#print axioms covered','end Spin.Structured.DenseOccupationFixed.Strip']
(r/'SpinCodes/Structured/DenseOccupationFixedStrip.lean').write_text('\n'.join(lines)+'\n',encoding='utf-8')
print('Strip module prepared, awaits batch dependencies')
