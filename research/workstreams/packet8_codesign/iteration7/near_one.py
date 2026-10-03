"""Four bounded near-one local/G4 witnesses, reusing one literal24 census.

This invokes frozen, audited outward arithmetic without changing its sources.
The libm exponential selects a dyadic z; the selected rational is certified.
"""
from __future__ import annotations
import os
for _key in ('OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','OMP_NUM_THREADS'):
    os.environ[_key]='1'
from fractions import Fraction
import hashlib
import json
from pathlib import Path
import sys
from time import monotonic
from types import SimpleNamespace
import numpy as np
HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parent/'iteration5'))
import local_outward24 as local
import local_batch as batch
import macro_outward as macro

POINTS=((Fraction(1,400),Fraction(7,10)),(Fraction(1,1600),Fraction(7,10)),
        (Fraction(3,400),Fraction(1,2)),(Fraction(3,1600),Fraction(1,2)))


def run():
    destination=HERE/'near_one_v1'
    destination.mkdir(parents=True,exist_ok=True)
    pins=local.sources()
    pins.update(macro.source_pins())
    pins[str(Path(__file__).resolve())]=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    pins[str(Path(batch.__file__).resolve())]=hashlib.sha256(Path(batch.__file__).read_bytes()).hexdigest()
    runtime=batch.runtime_metadata()
    data=census=None
    started=monotonic()
    for theta,alpha in POINTS:
        local_path,macro_path=batch.paths(destination,theta,alpha)
        z=Fraction.from_float(float(np.exp(-float(theta))))
        if not local_path.exists():
            if data is None:
                data=local.maps24.make_maps()
                census=local.maps24.census(data)
                print(f'exact census ready: {monotonic()-started:.2f}s',flush=True)
            begin=monotonic()
            lower,upper,diagnostic=local.local_operators(data,census,z,progress=True)
            if not macro.checked_pins(pins):raise ArithmeticError('source changed')
            saved=dict(schema='packet8-actual24-local-outward-1',z=str(z),
                z_numerator=z.numerator,z_denominator=z.denominator,z_binary64_hex=float(z).hex(),
                proposal_theta=str(theta),selected_alpha=str(alpha),map_record=data['record'],
                local_operator_lower=lower.tolist(),local_operator_upper=upper.tolist(),
                local_operator_lower_hex=[[[float(x).hex() for x in row] for row in m] for m in lower],
                local_operator_upper_hex=[[[float(x).hex() for x in row] for row in m] for m in upper],
                source_sha256=pins,source_pins_verified_at_finish=True,diagnostics=diagnostic,
                runtime=runtime,reused_exact_integer_census=True,
                outward_under_stated_arithmetic_contract=True,whole_code_certificate=False,
                elapsed_seconds=monotonic()-begin)
            local_path.write_text(json.dumps(saved,indent=2)+'\n',encoding='utf-8')
        else:
            _,stored_z,_=macro.load_local(local_path)
            if stored_z!=z:raise ArithmeticError('existing z mismatch')
        if not macro_path.exists():
            macro.run(SimpleNamespace(local=local_path,alpha=str(alpha),bins=256,output=macro_path))
        print(f'completed theta={theta} alpha={alpha}: total={monotonic()-started:.2f}s',flush=True)
    return monotonic()-started


if __name__=='__main__':run()
