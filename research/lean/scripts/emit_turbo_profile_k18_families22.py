"""Bounded three-family K18 candidate batch; no DP or exponential proof claim."""
from decimal import Decimal, localcontext, ROUND_HALF_EVEN
from fractions import Fraction
import json
from pathlib import Path

from emit_turbo_profile_family22 import produce_family, emit_lean


def main():
    candidate = json.loads(Path('.lake/turbo/coarse_packet_profiles.json').read_text())
    manifest = []
    for index,(lnum,lden) in enumerate(((1,1250),(3,2500),(1,625)),start=1):
        # Decimal exp only selects an untrusted dyadic candidate. Every proof
        # concerns this rational z; no theorem identifies z with exp(-lambda).
        zden = 1 << 80
        with localcontext() as ctx:
            ctx.prec = 100
            znum = int(((-Decimal(lnum)/Decimal(lden)).exp()*Decimal(zden))
                       .to_integral_value(rounding=ROUND_HALF_EVEN))
        namespace = f'Spin.Turbo.ProfileFamilyK18Tilt{index}22'
        data = produce_family(candidate,znum,zden,namespace)
        data['selection_only_lambda'] = dict(num=lnum,den=lden)
        json_path = Path(f'.lake/turbo/profile_family_k18_tilt{index}.json')
        lean_path = Path(f'SpinCodes/Turbo/PhysicalProfileFamilyK18Tilt{index}22Data.lean')
        json_path.write_text(json.dumps(data,indent=2)+'\n')
        emit_lean(data,lean_path)
        audit = data['precision_audit']
        minima = {key: Fraction(value['num'],value['den']) for key,value in audit.items()}
        # Safe diagnostic lower bound for every positive live-return target:
        # sqrt(S F)>=S F, S>=z^128, and M>=z^64>=z^128 F.
        minima['live_return_lower'] = (minima['moment_uniform_lower']*
            minima['minimum_positive_feedback_active']/4194303)
        # 32/D is a conservative local perturbation allowance (including the
        # sqrt branch and return rounding); report exact powers-of-two gates.
        gates = {key: bool(Fraction(32,1<<80)/value < Fraction(1,1<<48))
                 for key,value in minima.items()}
        if not all(gates.values()):
            raise RuntimeError(f'Insufficient proposed local precision: {gates}')
        manifest.append(dict(index=index,z_num=znum,z_den=zden,
            selection_only_lambda=dict(num=lnum,den=lden),namespace=namespace,
            json_path=str(json_path),lean_path=str(lean_path),
            endpoint_base_int_512=znum*(1<<432),
            relative_rounding_less_than_2_pow_minus_48=gates))
        print(f'K18 tilt{index}: z={znum}/{zden}; precision gates={gates}')
    manifest_path = Path('.lake/turbo/profile_family_k18_manifest.json')
    manifest_path.write_text(json.dumps(manifest,indent=2)+'\n')
    print(f'Finished exactly3 unique local families; q2 reuses tilt3. {manifest_path}')


if __name__ == '__main__':
    main()
