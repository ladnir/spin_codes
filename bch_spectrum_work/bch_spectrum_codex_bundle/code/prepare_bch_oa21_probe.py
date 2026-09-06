"""Exploratory exact LP with extra OA21 constraints; not yet a certificate.

The extension of the published length-255 dual-distance bound must be audited
before these constraints are used in an application theorem.
"""
import json
from pathlib import Path
from spectrum_exact import symmetric_kraw_coeff
from export_scaled_rational_lp import export_scaled
from run_higher_endpoint_preflight import write_new

ROOT=Path(__file__).resolve().parents[1]


def main():
    folder=ROOT/'generated/oa21_closure_probe'
    folder.mkdir(exist_ok=False)
    model=json.loads((ROOT/'generated/coupled_lp_exact.json').read_text())
    for row in model['constraints']:
        if row['name'].startswith(('Wambach_','orbit_search_')):
            row['rhs']='0'
    for prefix in ('q','h'):
        for degree in (16,18,20):
            model['constraints'].append(dict(name=f'provisional_OA21_{prefix}_{degree}',
                coeffs={f'{prefix}_{w}':str(symmetric_kraw_coeff(256,degree,w))
                        for w in range(0,129,2) if symmetric_kraw_coeff(256,degree,w)},
                sense='eq',rhs='0'))
    model['metadata']['provisional_extra_assumption']='Dual minimum at least 22 for the extended Q; source and extension proof pending audit'
    path=folder/'model.json'
    write_new(path,model)
    for objective in ('h_38','c_40','c_42'):
        export_scaled(path,objective,folder/(objective+'.lp'))
    print(folder)


if __name__=='__main__':
    main()
