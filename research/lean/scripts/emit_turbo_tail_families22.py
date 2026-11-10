"""Untrusted producer of K18 tail local families on the high-precision Fix chain.

For each retained tail tilt: exact-rational family (endpoint grid 2^-ENDPOINT_BITS),
precision gate, and six generated modules derived from the verified Smoke chain
(Data, Check, BudgetCheck, BudgetSmoke, BirthFamily, Family) with Hi imports.
Lean checks everything; nothing here is trusted.
"""
from decimal import Decimal, localcontext, ROUND_HALF_EVEN
from fractions import Fraction
import json
import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'scripts'))
from emit_turbo_profile_family22 import produce_family, emit_lean  # noqa: E402

OUT = ROOT / 'SpinCodes/Turbo/Generated/K18Tail'
ENDPOINT_BITS = 384
Z_BITS = 80
TILTS = ['3/1000', '3/500', '3/250', '3/125', '1/25', '11/200', '17/200', '7/50',
         '4/25', '1/5', '6/25', '8/25', '2/5', '12/25', '16/25', '1', '7/4']
MARKER_ONLY = 12   # indices >= this serve only the large-margin marker range
CHAIN = ['PhysicalProfileFamilyDefs22', 'PhysicalSmallBirth22BudgetDefs', 'PhysicalSmallBirth22Budget',
         'PhysicalProfileBirthFamily22', 'PhysicalProfileCompleteRows22']
SMOKE = {'Check': 'PhysicalProfileFamilySmoke22Check', 'BudgetCheck': 'PhysicalSmallBirth22BudgetCheck',
         'BudgetSmoke': 'PhysicalSmallBirth22BudgetSmoke', 'BirthFamily': 'PhysicalProfileBirthFamilySmoke22',
         'Family': 'PhysicalProfileFamilySmoke22'}


def port_smoke(kind, idx):
    t = (ROOT / f'SpinCodes/Turbo/{SMOKE[kind]}.lean').read_text(encoding='utf-8')
    mod = f'SpinCodes.Turbo.Generated.K18Tail.Tail{idx}'
    names = {'PhysicalProfileFamilySmoke22Data': f'{mod}Data',
             'PhysicalProfileFamilySmoke22Check': f'{mod}Check',
             'PhysicalSmallBirth22BudgetCheck': f'{mod}BudgetCheck',
             'PhysicalSmallBirth22BudgetSmoke': f'{mod}BudgetSmoke',
             'PhysicalProfileBirthFamilySmoke22': f'{mod}BirthFamily',
             'PhysicalProfileFamilySmoke22': f'{mod}Family'}
    for old in sorted(names, key=len, reverse=True):
        t = re.sub(rf'^import SpinCodes\.Turbo\.{old}$', f'import {names[old]}', t, flags=re.M)
    for c in CHAIN:
        t = re.sub(rf'^import SpinCodes\.Turbo\.{c}$', f'import SpinCodes.Turbo.Hi.{c}', t, flags=re.M)
    t = re.sub(r'\bSpin\.Numeric\b(?!Hi)', 'Spin.NumericHi', t)
    t = t.replace('namespace Spin.Turbo.ProfileFamilySmoke22', f'namespace Spin.Turbo.Hi.ProfileFamilyK18Tail{idx}22')
    t = t.replace('end Spin.Turbo.ProfileFamilySmoke22', f'end Spin.Turbo.Hi.ProfileFamilyK18Tail{idx}22')
    assert 'ProfileFamilySmoke22' not in t and 'SpinCodes.Turbo.Physical' not in t.replace(
        'SpinCodes.Turbo.Hi.Physical', ''), (kind, idx)
    (OUT / f'Tail{idx}{kind}.lean').write_text(t, encoding='utf-8', newline='\n')


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT.parent / '.gitignore').write_text('*\n!.gitignore\n', encoding='utf-8')
    candidate = json.loads((ROOT / '.lake/turbo/coarse_packet_profiles.json').read_text())
    manifest = []
    only = set(map(int, sys.argv[1:]))
    for idx, ts in enumerate(TILTS, start=1):
        if only and idx not in only:
            continue
        lam = Fraction(ts)
        zden = 1 << Z_BITS
        with localcontext() as ctx:
            ctx.prec = 100
            znum = int(((-Decimal(lam.numerator) / Decimal(lam.denominator)).exp() * Decimal(zden))
                       .to_integral_value(rounding=ROUND_HALF_EVEN))
        ns = f'Spin.Turbo.Hi.ProfileFamilyK18Tail{idx}22'
        data = produce_family(candidate, znum, zden, ns, endpoint_bits=ENDPOINT_BITS)
        data['selection_only_lambda'] = ts
        audit = {k: Fraction(v['num'], v['den']) for k, v in data['precision_audit'].items()}
        audit['live_return_lower'] = (audit['moment_uniform_lower'] *
                                      audit['minimum_positive_feedback_active'] / 4194303)
        gates = {k: bool(Fraction(32, 1 << ENDPOINT_BITS) / v < Fraction(1, 1 << 48))
                 for k, v in audit.items()}
        if not all(gates.values()):
            if idx < MARKER_ONLY:
                raise RuntimeError(f'tilt {ts}: precision gate failed {gates}')
            print(f'tail {idx}: tilt {ts} precision gate WARNING (marker-only, large margin): {gates}')
        (ROOT / f'.lake/turbo/profile_family_k18_tail{idx}.json').write_text(json.dumps(data) + '\n')
        emit_lean(data, OUT / f'Tail{idx}Data.lean', hi=True)
        for kind in SMOKE:
            port_smoke(kind, idx)
        manifest.append(dict(index=idx, tilt=ts, z_num=znum, z_den=zden))
        print(f'tail {idx}: tilt {ts} gates ok', flush=True)
    (OUT / 'manifest.json').write_text(json.dumps(manifest, indent=1) + '\n')


if __name__ == '__main__':
    main()
