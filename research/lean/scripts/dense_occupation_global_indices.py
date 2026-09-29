from pathlib import Path
import json
r=Path(__file__).resolve().parents[1]
d=json.loads((r.parent/'workstreams/inner_design/imt_asymptotic/d11/DENSE_REPLAY.json').read_text())
families={'occupation':'occupationIndices','statefree':'scalarIndices','fourier':'fourierIndices'}
s='import SpinCodes.Structured.DenseGeometryRate\n\nnamespace Spin.Structured.DenseGeometry\nset_option maxRecDepth 1000000\nset_option maxHeartbeats 0\n'
for family,name in families.items():
 ids=[i for i,b in enumerate(d['leaves']) if b['rational_witness']['family']==family]
 s+=f'def {name} : List ℕ := ['+', '.join(map(str,ids))+']\n'
s+='\ntheorem cover_indices : ∀ i : Fin 1023, i.val ∈ occupationIndices ∨ i.val ∈ scalarIndices ∨ i.val ∈ fourierIndices := by decide\n\n'
s+='theorem denseRates_of_families {η C : ℝ}\n    (ho : ∀ i ∈ occupationIndices, CertifiedBox i η C)\n    (hs : ∀ i ∈ scalarIndices, CertifiedBox i η C)\n    (hf : ∀ i ∈ fourierIndices, CertifiedBox i η C) :\n    DenseOccupationFixed.DenseRates η C := by\n  apply denseRates_of_boxes\n  intro i hi\n  rcases cover_indices ⟨i,hi⟩ with h | h | h\n  · exact ho i h\n  · exact hs i h\n  · exact hf i h\n\n#print axioms cover_indices\n#print axioms denseRates_of_families\nend Spin.Structured.DenseGeometry\n'
(r/'SpinCodes/Structured/DenseOccupationGlobalIndices.lean').write_text(s,encoding='utf-8')
print('Prepared exact three-family index cover')
