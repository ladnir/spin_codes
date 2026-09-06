"""Add two exact signed weight identities to the small OA21 model."""
import json
from pathlib import Path
from export_scaled_rational_lp import export_scaled
from run_higher_endpoint_preflight import write_new,sha
ROOT=Path(__file__).resolve().parents[1]
FOLDER=ROOT/'generated/oa21_mod4_probe'
FOLDER.mkdir(exist_ok=False)
source=ROOT/'generated/oa21_closure_probe/model.json'
model=json.loads(source.read_text())
evidence=ROOT/'generated/bch256_quadratic_weight_sums.json'
sums=json.loads(evidence.read_text())['results']
for prefix,label in (('q','Q'),('h','H')):
    model['constraints'].append(dict(name=f'exact_mod4_{prefix}',
        coeffs={f'{prefix}_{w}':str((1 if w==128 else 2)*(-1)**(w//2)) for w in range(0,129,2)},
        sense='eq',rhs=sums[label]['signed_sum']))
model['metadata']['signed_weight_source_sha256']={str(p.relative_to(ROOT)):sha(p) for p in (Path(__file__),source,evidence)}
write_new(FOLDER/'model.json',model)
export_scaled(FOLDER/'model.json','h_38',FOLDER/'h_38.lp')
print('Prepared 404-row exact OA21 plus modulo-four LP')
