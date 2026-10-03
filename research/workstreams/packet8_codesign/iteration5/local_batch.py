"""Generate the selected outward witnesses using one exact24-state census.

The libm exponential only selects a dyadic proof parameter. Every receipt
identifies that parameter exactly; no exponential approximation is certified.
Existing authenticated outputs are reused, never overwritten.
"""
from __future__ import annotations

import os
for _name in ('OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'OMP_NUM_THREADS'):
    os.environ[_name] = '1'

import argparse
from fractions import Fraction
import hashlib
import json
from pathlib import Path
import platform
import sys
from time import monotonic
from types import SimpleNamespace

import numpy as np
import local_outward24 as local
import macro_outward as macro

HERE = Path(__file__).resolve().parent
POINTS = (('.01', '.7'), ('.06', '.4'), ('.12', '.3'), ('.2', '.35'),
          ('.3', '.35'), ('.4', '.35'), ('.5', '.4'), ('.8', '.4'),
          ('1', '.4'), ('1.4', '.5'), ('1.8', '.6'), ('2.2', '1'))


def paths(directory, theta, alpha):
    suffix = f'{theta.numerator}_{theta.denominator}'
    return (directory/f'local_theta_{suffix}.json',
            directory/f'macro_theta_{suffix}_alpha_{alpha.numerator}_{alpha.denominator}.json')


def runtime_metadata():
    macro.op.check_runtime()
    return dict(python=sys.version, numpy=np.__version__, platform=platform.platform(),
        machine=platform.machine(), binary64_mantissa_bits=np.finfo(float).nmant+1,
        positive_arithmetic_runtime_check_passed=True,
        thread_environment={name: os.environ.get(name) for name in
            ('OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'OMP_NUM_THREADS')})


def run(args):
    args.directory.mkdir(parents=True, exist_ok=True)
    pins = local.sources()
    pins.update(macro.source_pins())
    pins[str(Path(__file__).resolve())] = hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    runtime = runtime_metadata()
    data, census = None, None
    started = monotonic()
    for theta_text, alpha_text in POINTS:
        theta, alpha = Fraction(theta_text), Fraction(alpha_text)
        z = Fraction.from_float(float(np.exp(-float(theta))))
        local_path, macro_path = paths(args.directory, theta, alpha)
        if local_path.exists():
            _, stored_z, saved = macro.load_local(local_path)
            if stored_z != z or saved.get('proposal_theta') != str(theta):
                raise ValueError('existing local receipt has different parameters')
            print(f'reusing authenticated {local_path.name}', flush=True)
        else:
            if data is None:
                data = local.maps24.make_maps()
                census = local.maps24.census(data)
                print(f'shared literal24 census ready in {monotonic()-started:.2f}s', flush=True)
            beginning = monotonic()
            lower, upper, diagnostic = local.local_operators(data, census, z, progress=True)
            if not macro.checked_pins(pins):
                raise ArithmeticError('source changed during local batch')
            saved = dict(schema='packet8-actual24-local-outward-1', z=str(z),
                z_numerator=z.numerator, z_denominator=z.denominator, z_binary64_hex=float(z).hex(),
                proposal_theta=str(theta), selected_alpha=str(alpha), map_record=data['record'],
                local_operator_lower=lower.tolist(), local_operator_upper=upper.tolist(),
                local_operator_lower_hex=[[[float(x).hex() for x in row] for row in matrix] for matrix in lower],
                local_operator_upper_hex=[[[float(x).hex() for x in row] for row in matrix] for matrix in upper],
                source_sha256=pins, source_pins_verified_at_finish=True, diagnostics=diagnostic,
                runtime=runtime, reused_exact_integer_census=True,
                outward_under_stated_arithmetic_contract=True, whole_code_certificate=False,
                elapsed_seconds=monotonic()-beginning)
            local_path.write_text(json.dumps(saved, indent=2)+'\n', encoding='utf-8')
            print(f'saved {local_path.name}; maxrelativewidth={diagnostic["max_relative_width"]:.3g}; '
                  f'elapsed={monotonic()-started:.2f}s', flush=True)
        if macro_path.exists():
            saved_macro = json.loads(macro_path.read_text(encoding='utf-8'))
            if (saved_macro.get('schema') != 'packet8-actual24-G4-outward-1'
                    or Fraction(saved_macro['z']) != z or Fraction(saved_macro['alpha']) != alpha
                    or not saved_macro['source_pins_verified_at_finish']
                    or not macro.checked_pins(saved_macro['source_sha256'])):
                raise ValueError('existing macro receipt cannot be authenticated')
            print(f'reusing authenticated {macro_path.name}', flush=True)
        else:
            macro.run(SimpleNamespace(local=local_path, alpha=str(alpha), bins=256, output=macro_path))
        if not macro.checked_pins(pins):
            raise ArithmeticError('source changed during completed batch point')
    print(f'all{len(POINTS)} local/macro witnesses authenticate; elapsed={monotonic()-started:.2f}s', flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--directory', type=Path, default=HERE/'cache_v1')
    run(parser.parse_args())
