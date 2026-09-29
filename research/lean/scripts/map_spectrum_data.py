"""Generate untrusted block histograms for kernel replay of both map spectra."""
import argparse
import hashlib
import json
from pathlib import Path
import re

from concrete_maps import encode

ROOT = Path(__file__).resolve().parents[1]


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--size', type=int, default=1024)
    args = parser.parse_args()
    assert args.size > 0 and 524288 % args.size == 0
    bits = args.size.bit_length()-1
    assert args.size == 2**bits
    paper = ROOT.parent/'paper/structured_imt_appendix.tex'
    entries = re.findall(r'(?m)^(\d+) & \\texttt\{([0-9a-f]+)\} & \\texttt\{([0-9a-f]+)\}',
                         paper.read_text(encoding='utf-8'))
    assert [int(j) for j, _, _ in entries] == list(range(19))
    maps = {'a': [int(a,16) for _,a,_ in entries],
            'ct': [int(c,16) for _,_,c in entries]}
    directory = ROOT/'SpinCodes/Structured/MapSpectrumData'
    directory.mkdir(exist_ok=True)
    data = ROOT/'scripts/map_data'
    data.mkdir(exist_ok=True)
    basis_lines = ['import SpinCodes.Structured.ConcreteMapData',
                   'namespace Spin.Structured.MapSpectrum',
                   'set_option maxRecDepth 1000000', 'set_option maxHeartbeats 0']
    for name, rowname in [('a','aRows'),('ct','cTransposeRows')]:
        values = ', '.join(hex(encode(maps[name][:bits],q)) for q in range(args.size))
        basis_lines += [f'def {name}Low : List Nat := [{values}]',
            f'theorem {name}Low_checked : (List.range {args.size}).map (PackedMap.eval (ConcreteMaps.{rowname}.take {bits})) = {name}Low := by decide']
    basis_lines += ['end Spin.Structured.MapSpectrum', '']
    basis_path = directory/'Basis.lean'
    basis_path.write_text('\n'.join(basis_lines),encoding='utf-8')
    (data/'basis_manifest.json').write_text(json.dumps(dict(path=basis_path.relative_to(ROOT).as_posix(),
        sha256=hashlib.sha256(basis_path.read_bytes()).hexdigest()),indent=2)+'\n')
    totals = {name: [0]*129 for name in maps}
    manifest = []
    for block in range(524288//args.size):
        start = block * args.size
        histograms = {name: [0]*129 for name in maps}
        weights = {name: [] for name in maps}
        for q in range(start, start+args.size):
            for name, rows in maps.items():
                weight = encode(rows, q).bit_count()
                weights[name].append(weight)
                histograms[name][weight] += 1
        lines = ['import SpinCodes.Structured.MapSpectrumDefs',
                 'import SpinCodes.Structured.MapSpectrumData.Basis',
                 'namespace Spin.Structured.MapSpectrum',
                 'set_option maxRecDepth 1000000', 'set_option maxHeartbeats 0']
        for name, rows in [('a','aRows'),('ct','cTransposeRows')]:
            nameblock = f'{name}{block}'
            high = encode(maps[name][bits:],block)
            values = ', '.join(map(str,histograms[name]))
            wordweights = ', '.join(map(str,weights[name]))
            lines += [f'def {nameblock}High : Nat := {hex(high)}',
                      f'theorem {nameblock}_high_checked : PackedMap.eval (ConcreteMaps.{rows}.drop {bits}) {block} = {nameblock}High := by decide',
                      f'def {nameblock}Weights : List Nat := [{wordweights}]',
                      f'theorem {nameblock}_fast_checked : {name}Low.map (fun q => fastWeight 128 (q ^^^ {nameblock}High)) = {nameblock}Weights := by decide +kernel',
                      f'theorem {nameblock}_weights_checked : {name}Low.map (fun q => PackedMap.weight 128 (q ^^^ {nameblock}High)) = {nameblock}Weights := by',
                      f'  simpa only [fastWeight_eq] using {nameblock}_fast_checked',
                      f'def {nameblock} : List Nat := [{values}]',
                      f'theorem {nameblock}_hist_checked : {nameblock}Weights.foldl (fun h w => bump w h) (List.replicate 129 0) = {nameblock} := by decide +kernel']
            totals[name] = [x+y for x,y in zip(totals[name],histograms[name])]
        lines += ['end Spin.Structured.MapSpectrum', '']
        path = directory/f'Block{block}.lean'
        path.write_text('\n'.join(lines), encoding='utf-8')
        manifest.append(dict(block=block, start=start, size=args.size,
            path=path.relative_to(ROOT).as_posix(), sha256=hashlib.sha256(path.read_bytes()).hexdigest()))
    assert totals['a'][0] == 1
    assert {i:c for i,c in enumerate(totals['a']) if c} == {
        0:1, 48:5166, 56:110288, 64:293455, 72:110128, 80:5250}
    (data/'manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
    (data/'spectra.json').write_text(json.dumps(totals,indent=2)+'\n')
    print(f'Generated {len(manifest)} blocks, {args.size} inputs each, both maps.')


if __name__ == '__main__':
    main()
