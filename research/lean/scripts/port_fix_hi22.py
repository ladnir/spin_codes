"""Port the Fix-based local-family chain to a high-precision Fix (scale 10^120).

Mechanical text transformation only; every ported module is re-checked by Lean.
Originals are untouched. Ported modules: SpinCodes/Turbo/Hi/*, namespace
Spin.Turbo.Hi (enclosing-namespace resolution keeps base references intact),
Fix in namespace Spin.NumericHi.
"""
import re
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
DEP = Path('C:/Users/peter/repo/permute_conv-github-bch/research/lean')
HI = ROOT / 'SpinCodes/Turbo/Hi'
PREC = 120
CHAIN = ['PhysicalProfileNumericDefs22', 'PhysicalProfileNumeric22', 'PhysicalProfileNumericAggregateDefs22',
         'PhysicalProfileNumericAggregate22', 'PhysicalProfileNumericFactors22', 'PhysicalProfileNumericCacheDefs22',
         'PhysicalProfileNumericCache22', 'PhysicalProfileNumericData22', 'PhysicalProfileNumericTable22',
         'PhysicalProfileFamilyDefs22', 'PhysicalProfileFamily22', 'PhysicalProfileFamilyExtract22',
         'PhysicalProfileFamilyFourier22', 'PhysicalProfileMoment22', 'PhysicalProfileWeighted22',
         'PhysicalProfileRows22', 'PhysicalProfileAssembly22', 'PhysicalProfileCompleteRows22',
         'PhysicalProfileBirthFamily22', 'PhysicalSmallBirth22BudgetDefs', 'PhysicalSmallBirth22Budget']


def numeric(text):
    return re.sub(r'\bSpin\.Numeric\b(?!Hi)', 'Spin.NumericHi', text)


def port_fix():
    HI.mkdir(parents=True, exist_ok=True)
    d = (DEP / 'SpinCodes/Numeric/FixedDefs.lean').read_text(encoding='utf-8')
    d = numeric(d)
    d = re.sub(r'def prec : Nat := 30', f'def prec : Nat := {PREC}', d)
    d = re.sub(r'def scale : Int := 1(0{30})\b', 'def scale : Int := 1' + '0' * PREC, d)
    assert f'def prec : Nat := {PREC}' in d and ('1' + '0' * PREC) in d
    (HI / 'FixedDefs.lean').write_text(d, encoding='utf-8', newline='\n')
    f = (DEP / 'SpinCodes/Numeric/Fixed.lean').read_text(encoding='utf-8')
    f = f.replace('import SpinCodes.Numeric.FixedDefs', 'import SpinCodes.Turbo.Hi.FixedDefs')
    f = numeric(f)
    (HI / 'Fixed.lean').write_text(f, encoding='utf-8', newline='\n')


def port_module(name):
    t = (ROOT / f'SpinCodes/Turbo/{name}.lean').read_text(encoding='utf-8')
    t = t.replace('import SpinCodes.Numeric.FixedDefs', 'import SpinCodes.Turbo.Hi.FixedDefs')
    t = t.replace('import SpinCodes.Numeric.Fixed\n', 'import SpinCodes.Turbo.Hi.Fixed\n')
    for other in CHAIN:
        t = re.sub(rf'^import SpinCodes\.Turbo\.{other}$', f'import SpinCodes.Turbo.Hi.{other}', t, flags=re.M)
    t = numeric(t)
    t = re.sub(r'^namespace Spin\.Turbo$', 'namespace Spin.Turbo.Hi', t, flags=re.M)
    t = re.sub(r'^end Spin\.Turbo$', 'end Spin.Turbo.Hi', t, flags=re.M)
    bad = re.findall(r'Spin\.Turbo\.(' + '|'.join(CHAIN) + r')\b', t)
    assert not bad, (name, bad)
    (HI / f'{name}.lean').write_text(t, encoding='utf-8', newline='\n')


if __name__ == '__main__':
    port_fix()
    for m in CHAIN:
        port_module(m)
    print('ported', len(CHAIN), 'modules')
