"""Compare the long-region cancellation model with finite surface slices."""
import hashlib
import json
from pathlib import Path

import activation_refresh_native as native
import continuous_state_model as continuous
from report_engineering_surfaces import read,select

HERE=Path(__file__).resolve().parent


def main():
    rows,payload=read();results=[];kernel=native.RefreshKernel()
    for b in (8,32,128,512):
        spectra={family:{sh['weight']:sh['count'] for sh in value['shells']}
                 for key,value in payload['models'].items()
                 for family in (key.rsplit('_',1)[0],)
                 if value['block_bits']==b and family in ('bch','rm','random_mean')}
        for s in (10,12,16,20):
            answers=continuous.screen(kernel,b,s,spectra)
            for family,answer in answers.items():
                errors={str(e):answer['intercept_bits']-select(rows,family=family,block_bits=b,
                    step_bits=64,state_bits=s,message_exponent=e)[0]['intercept_bits'] for e in (16,20,24)}
                row=dict(family=family,block_bits=b,state_bits=s,intercept_bits=answer['intercept_bits'],
                         model_minus_finite_bits_by_exponent=errors)
                results.append(row);print(row,flush=True)
    paths=(Path(__file__),Path(continuous.__file__),Path(native.__file__),native.LIBRARY,
           HERE/'engineering_surfaces.json',HERE/'engineering_surfaces.csv',HERE/'study_engineering_surfaces.py')
    (HERE/'continuous_state_comparison.json').write_text(json.dumps(dict(rows=results,
        source_sha256={str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths},
        note='Long-region model comparison; not a finite-distance proof.'),indent=2)+'\n')


if __name__=='__main__': main()
