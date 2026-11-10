"""Untrusted producer of the K20 q1 tilt families (physical rows are
geometry-independent; only the tilt differs from K18).

K20 q1 retained tilts (q1_coarse_complete65_diagnostic.json):
1/5000, 1/2500, 3/5000, 1/1250. Each family gets the Data module from
`produce_family` and the five fixed proof modules stamped from the green
K18 tilt-1 templates by name substitution. Lean checks everything.
"""
from decimal import Decimal, localcontext, ROUND_HALF_EVEN
from fractions import Fraction
import json
from pathlib import Path

from emit_turbo_profile_family22 import produce_family, emit_lean

T = Path('SpinCodes/Turbo')
TILTS = ((1, 5000), (1, 2500), (3, 5000), (1, 1250))
TEMPLATES = ['PhysicalProfileFamilyK18Tilt122Check.lean', 'PhysicalProfileFamilyK18Tilt122.lean',
             'PhysicalProfileBirthFamilyK18Tilt122.lean', 'PhysicalSmallBirth22K18Tilt1BudgetCheck.lean',
             'PhysicalSmallBirth22K18Tilt1BudgetSmoke.lean']


def stamp(index):
    for name in TEMPLATES:
        text = (T / name).read_text(encoding='utf-8')
        out = name.replace('K18Tilt1', f'K20Tilt{index}')
        text = text.replace('K18Tilt1', f'K20Tilt{index}')
        assert 'K18' not in text, (name, index)
        (T / out).write_text(text, encoding='utf-8', newline='\n')


def main():
    candidate = json.loads(Path('.lake/turbo/coarse_packet_profiles.json').read_text())
    manifest = []
    for index, (lnum, lden) in enumerate(TILTS, start=1):
        zden = 1 << 80
        with localcontext() as ctx:
            ctx.prec = 100
            znum = int(((-Decimal(lnum) / Decimal(lden)).exp() * Decimal(zden))
                       .to_integral_value(rounding=ROUND_HALF_EVEN))
        namespace = f'Spin.Turbo.ProfileFamilyK20Tilt{index}22'
        data = produce_family(candidate, znum, zden, namespace)
        data['selection_only_lambda'] = dict(num=lnum, den=lden)
        json_path = Path(f'.lake/turbo/profile_family_k20_tilt{index}.json')
        lean_path = T / f'PhysicalProfileFamilyK20Tilt{index}22Data.lean'
        json_path.write_text(json.dumps(data, indent=2) + '\n')
        emit_lean(data, lean_path)
        audit = data['precision_audit']
        minima = {key: Fraction(value['num'], value['den']) for key, value in audit.items()}
        minima['live_return_lower'] = (minima['moment_uniform_lower'] *
            minima['minimum_positive_feedback_active'] / 4194303)
        gates = {key: bool(Fraction(32, 1 << 80) / value < Fraction(1, 1 << 48))
                 for key, value in minima.items()}
        if not all(gates.values()):
            raise RuntimeError(f'Insufficient proposed local precision at tilt {index}: {gates}')
        stamp(index)
        manifest.append(dict(index=index, z_num=znum, z_den=zden,
                             selection_only_lambda=dict(num=lnum, den=lden), namespace=namespace,
                             json_path=str(json_path), lean_path=str(lean_path)))
        print(f'K20 tilt{index}: z={znum}/{zden}; gates ok')
    Path('.lake/turbo/profile_family_k20_manifest.json').write_text(json.dumps(manifest, indent=2) + '\n')


if __name__ == '__main__':
    main()
