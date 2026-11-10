"""Emit untrusted, small quadratic-profile zero-set witnesses (no state census)."""
from __future__ import annotations
import argparse
import json
from pathlib import Path

ROWS = [0xffff,0xaaaa,0xcccc,0xf0f0,0xff00,0x8888,0xa0a0,0xaa00,0xc0c0,0xcc00,0xf000]
MONOMIALS = [0,1,2,4,8,3,5,9,6,10,12]

def evaluate(rows, value):
    result = 0
    for i, row in enumerate(rows):
        if value >> i & 1:
            result ^= row
    return result

def quadratic_coefficients(word):
    coefficients = 0
    for i, monomial in enumerate(MONOMIALS):
        value = 0
        submask = monomial
        while True:
            value ^= word >> submask & 1
            if submask == 0:
                break
            submask = (submask - 1) & monomial
        coefficients |= value << i
    if evaluate(ROWS, coefficients) != word:
        raise ValueError('permuted basis is not quadratic')
    return coefficients

def make_witness(zero_mask):
    points = [v for v in range(16) if zero_mask >> v & 1]
    if len(points) not in (4,8,16):
        raise ValueError('nonempty affine zero set has invalid size')
    dimension = len(points).bit_length() - 1
    origin = points[0]
    basis, span = [], {0}
    for direction in [v ^ origin for v in points]:
        if direction not in span:
            basis.append(direction)
            span |= {v ^ direction for v in tuple(span)}
    if len(basis) != dimension or span != {v ^ origin for v in points}:
        raise ValueError('zero set is not an affine subspace')
    for direction in (1,2,4,8):
        if direction not in span:
            basis.append(direction)
            span |= {v ^ direction for v in tuple(span)}
    if len(basis) != 4 or span != set(range(16)):
        raise ValueError('basis completion failed')
    permutation = [origin ^ evaluate(basis, v) for v in range(16)]
    inverse = [0] * 16
    for i, value in enumerate(permutation):
        inverse[value] = i
    if sorted(permutation) != list(range(16)):
        raise ValueError('invalid permutation')
    if any(bool(zero_mask >> permutation[v] & 1) != (v < 1 << dimension) for v in range(16)):
        raise ValueError('canonical zero-subcube alignment failed')
    images = []
    for row in ROWS:
        word = sum(((row >> permutation[v]) & 1) << v for v in range(16))
        images.append(quadratic_coefficients(word))
    return dict(kind=dimension-2, zero_mask=zero_mask,
                permutation=sum(p << (4*i) for i,p in enumerate(permutation)),
                inverse=sum(p << (4*i) for i,p in enumerate(inverse)), images=images)

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, default=Path('.lake/turbo/packet_profile_witnesses.json'))
    parser.add_argument('--lean-output', type=Path)
    parser.add_argument('--table-output', type=Path)
    parser.add_argument('--numeric-table-output', type=Path)
    parser.add_argument('--profile-candidate', type=Path, default=Path('.lake/turbo/coarse_packet_profiles.json'))
    args = parser.parse_args()
    affine_masks = [evaluate(ROWS[:5], a) for a in range(32)]
    masks = [0xffff ^ (affine_masks[p % 32] | affine_masks[p // 32]) for p in range(1024)]
    distinct = sorted(set(masks) - {0})
    bank = [make_witness(mask) for mask in distinct]
    index = {mask: i+1 for i,mask in enumerate(distinct)}
    classifier = [index.get(mask, 0) for mask in masks]
    counts = [0] * 4
    for i in classifier:
        counts[0 if i == 0 else bank[i-1]['kind']+1] += 1
    if len(bank) != 171 or counts != [93,840,90,1]:
        raise ValueError('unexpected finite affine-pair classifier')
    result = dict(schema='map22-packet-profile-zero-witnesses-1',
                  scope='Untrusted finite basis-preserver witnesses; not a checked census.',
                  canonical_rows=ROWS, bank=bank, classifier=classifier, class_counts=counts)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, separators=(',',':'))+'\n', encoding='utf-8')
    if args.lean_output:
        lines = ['import SpinCodes.UniformMixer.Map22PacketProfileWitness',
                 'namespace Spin.UniformMixer', 'open Spin.Structured',
                 'set_option maxRecDepth 10000', 'set_option maxHeartbeats 0']
        for i, w in enumerate(bank):
            images = ','.join(map(str, w['images']))
            lines += [f'def profileZeroWitness{i} : ProfileZeroWitness :=',
                      f'  ⟨{w["kind"]}, {w["zero_mask"]}, {w["permutation"]}, {w["inverse"]}, [{images}]⟩',
                      f'theorem profileZeroWitness{i}_checked : profileZeroWitness{i}.check = true := by decide']
        lines += ['def profileZeroBank : List ProfileZeroWitness := [',
                  ',\n'.join(f'  profileZeroWitness{i}' for i in range(len(bank))), ']']
        lines += ['theorem profileZeroBank_checked : ∀ w ∈ profileZeroBank, w.check = true := by',
                  '  simp only [profileZeroBank, List.mem_cons, List.not_mem_nil, or_false]',
                  '  intro w h', '  rcases h with ' + ' | '.join('h' for _ in bank)]
        for i in range(len(bank)):
            lines += [f'  · subst w; exact profileZeroWitness{i}_checked']
        lines += ['def profileZeroClassifier : List Nat := [',
                  ','.join(map(str,classifier)), ']',
                  'def profilePackedAffine (p : Nat) : Nat := PackedMap.eval (profileQuad4Rows.take 5) p',
                  'def profileClassWitness (p : Fin 1024) : ProfileZeroWitness :=',
                  '  profileZeroBank.getD ((profileZeroClassifier.getD p.val 0)-1) profileZeroWitness0',
                  'def profileZeroClass (p : Fin 1024) : Nat :=',
                  '  if profileZeroClassifier.getD p.val 0 = 0 then 0 else (profileClassWitness p).kind.val+1',
                  'theorem profileZeroClassifier_checked : ∀ p : Fin 1024,',
                  '  if profileZeroClassifier.getD p.val 0 = 0 then',
                  '    ∀ v : Fin 16, ((profilePackedAffine (p.val%32)).testBit v.val ||',
                  '      (profilePackedAffine (p.val/32)).testBit v.val) = true',
                  '  else (profileZeroClassifier.getD p.val 0)-1 < profileZeroBank.length ∧',
                  '    ∀ v : Fin 16, (profileClassWitness p).zeroMask.testBit v.val =',
                  '      !((profilePackedAffine (p.val%32)).testBit v.val ||',
                  '        (profilePackedAffine (p.val/32)).testBit v.val) := by decide',
                  'def profileZeroKinds : List Nat := [' + ','.join(str(0 if i==0 else bank[i-1]['kind']+1) for i in classifier) + ']',
                  'theorem profileZeroKinds_eq : profileZeroKinds = profileZeroClassifier.map (fun i =>',
                  '  if i = 0 then 0 else (profileZeroBank.getD (i-1) profileZeroWitness0).kind.val+1) := by decide',
                  'theorem profileZeroKinds_counts : ∀ k : Fin 4,',
                  '  profileZeroKinds.count k.val = ([93,840,90,1] : List Nat).getD k.val 0 := by decide',
                  'end Spin.UniformMixer']
        split = lines.index('def profileZeroClassifier : List Nat := [')
        args.lean_output.write_text('\n'.join(lines[:split]+['end Spin.UniformMixer'])+'\n', encoding='utf-8')
        classification = args.lean_output.with_name('Map22PacketProfileClassifier.lean')
        classification.write_text('\n'.join([
            'import SpinCodes.UniformMixer.Map22PacketProfileBank',
            'namespace Spin.UniformMixer', 'open Spin.Structured',
            'set_option maxRecDepth 10000', 'set_option maxHeartbeats 0'] + lines[split:])+'\n', encoding='utf-8')
        source = classification.read_text(encoding='utf-8')
        chunked = 'def profileZeroChunks : List (List Nat) := [\n' + ',\n'.join(
            '  [' + ','.join(map(str,classifier[i:i+32])) + ']' for i in range(0,1024,32)) + ']\n'
        chunked += 'def profileZeroIndex (p : Fin 1024) : Nat := (profileZeroChunks.getD (p.val/32) []).getD (p.val%32) 0\n'
        source = source.replace('def profilePackedAffine', chunked+'def profilePackedAffine')
        source = source.replace('profileZeroClassifier.getD p.val 0', 'profileZeroIndex p')
        source = source.replace('theorem profileZeroClassifier_checked : ∀ p : Fin 1024,',
                                'def profileZeroClassifierValid (p : Fin 1024) : Prop :=')
        source = source.replace(')).testBit v.val) := by decide', ')).testBit v.val)')
        marker = 'def profileZeroKinds : List Nat'
        checks = []
        for block in range(32):
            checks += [f'theorem profileZeroBlock{block}_checked : ∀ p : Fin 32,',
                       f'  profileZeroClassifierValid ⟨{block}*32+p.val, by omega⟩ := by decide']
        checks += ['theorem profileZeroClassifier_checked (p : Fin 1024) : profileZeroClassifierValid p := by',
                   '  have hb : p.val / 32 < 32 := by omega',
                   '  let hc := (⟨p.val / 32, hb⟩ : Fin 32)',
                   '  have hr : ∀ b : Fin 32, ∀ r : Fin 32, profileZeroClassifierValid ⟨b.val*32+r.val, by omega⟩ := by',
                   '    intro b; fin_cases b']
        for block in range(32):
            checks += [f'    · exact profileZeroBlock{block}_checked']
        checks += ['  have hh := hr hc ⟨p.val%32, Nat.mod_lt _ (by decide)⟩',
                   '  have he : (⟨hc.val*32+(p.val%32), by omega⟩ : Fin 1024) = p := by',
                   '    apply Fin.ext; dsimp [hc]; omega',
                   '  rw [he] at hh; exact hh']
        checks = ['instance : DecidablePred profileZeroClassifierValid := fun p => by',
                  '  unfold profileZeroClassifierValid; infer_instance'] + checks
        source = source.replace(marker, '\n'.join(checks)+'\n'+marker)
        kinds = [0 if i == 0 else bank[i-1]['kind']+1 for i in classifier]
        kindchunks = [kinds[i:i+32] for i in range(0,1024,32)]
        kindproof = ['def profileZeroKindChunks : List (List Nat) := [' + ','.join(
            '['+','.join(map(str,chunk))+']' for chunk in kindchunks) + ']']
        for block in range(32):
            kindproof += [f'theorem profileZeroKinds_block{block} :',
                f'  List.ofFn (fun p : Fin 32 => profileZeroClass ⟨{block}*32+p.val, by omega⟩) =',
                '  ['+','.join(map(str,kindchunks[block]))+'] := by decide']
        source = source.replace('theorem profileZeroKinds_eq', '\n'.join(kindproof)+'\ntheorem profileZeroKinds_eq')
        # The block ofFn proofs are the semantic route to the class count;
        # avoid re-running the equivalent long flat-list computation.
        start = source.index('theorem profileZeroKinds_eq')
        end = source.index('theorem profileZeroKinds_counts', start)
        source = source[:start] + source[end:]
        classification.write_text(source, encoding='utf-8')
    if args.table_output:
        candidate = json.loads(args.profile_candidate.read_text(encoding='utf-8'))
        entries = ['  ('+str(count)+', !['+','.join(map(str,profile))+'])'
                   for profile,count in zip(candidate['profiles'],candidate['multiplicities'])]
        table = ['import SpinCodes.UniformMixer.Map22PacketProfileCensus',
                 'noncomputable section', 'namespace Spin.UniformMixer',
                 'set_option maxHeartbeats 1000000', 'set_option maxRecDepth 10000',
                 'def productionPacketProfileTable : List (Nat × (Fin 5 → Nat)) := [',
                 ',\n'.join(entries), ']',
                 'theorem productionPacketProfile_table_sum {M : Type*} [AddCommMonoid M]',
                 '    (f : (Fin 5 → Nat) → M) :',
                 '    (∑ q : State22, f (packetProfile64 (fun j => coordMap j q))) =',
                 '      (productionPacketProfileTable.map (fun e => e.1 • f e.2)).sum := by',
                 '  rw [productionPacketProfile_census]',
                 '  simp only [profileEvenClassSum, Fin.sum_univ_succ]',
                 '  simp [profileQuad4WeightCounts, profileQuad3WeightCounts, profileQuad2WeightCounts,',
                 '    productionPacketProfileTable, profileOddShape, profileEvenShape,',
                 '    smul_smul]',
                 '  <;> abel',
                 'theorem productionPacketProfileTable_total :',
                 '    (productionPacketProfileTable.map Prod.fst).sum = 2^22 := by decide',
                 'theorem productionPacketProfileTable_length : productionPacketProfileTable.length = 25 := by decide',
                 'theorem productionPacketProfileTable_counts_checked : ∀ i : Fin 25,',
                 '    (productionPacketProfileTable.map (fun e =>',
                 '      if e.2 = (productionPacketProfileTable.getD i.val (0,0)).2 then e.1 else 0)).sum =',
                 '      (productionPacketProfileTable.getD i.val (0,0)).1 := by decide',
                 'theorem productionPacketProfile_table_count (i : Fin 25) :',
                 '    (Finset.univ.filter (fun q : State22 => packetProfile64 (fun j => coordMap j q) =',
                 '      (productionPacketProfileTable.getD i.val (0,0)).2)).card =',
                 '      (productionPacketProfileTable.getD i.val (0,0)).1 := by',
                 '  rw [Finset.card_filter, productionPacketProfile_table_sum (fun p =>',
                 '    if p = (productionPacketProfileTable.getD i.val (0,0)).2 then (1 : Nat) else 0)]',
                 '  simpa only [nsmul_eq_mul, mul_ite, Nat.mul_one, Nat.mul_zero, Nat.cast_id] using',
                 '    productionPacketProfileTable_counts_checked i',
                 '#print axioms productionPacketProfile_table_sum',
                 '#print axioms productionPacketProfile_table_count',
                 'end Spin.UniformMixer']
        args.table_output.write_text('\n'.join(table)+'\n', encoding='utf-8')
    if args.numeric_table_output:
        candidate = json.loads(args.profile_candidate.read_text(encoding='utf-8'))
        entries = ['  ('+str(count)+', ['+','.join(map(str,profile))+'])'
                   for profile,count in zip(candidate['profiles'],candidate['multiplicities'])]
        lines = ['import SpinCodes.Turbo.PhysicalProfileNumericAggregateDefs22',
                 'namespace Spin.Turbo',
                 '/- Untrusted candidate literals; authenticated in the semantic bridge. -/',
                 'def profileNumericEntries22 : List (Nat × List Nat) := [',
                 ',\n'.join(entries), ']', 'end Spin.Turbo']
        args.numeric_table_output.write_text('\n'.join(lines)+'\n', encoding='utf-8')
    print(f'Emitted {len(bank)} deduplicated zero-set witnesses and 1024 affine-pair indices to {args.output}.')

if __name__ == '__main__':
    main()
