"""Selected low-arity screen without the added OUTER pair layer.

The existing global paired route and inner remain unchanged. This is a
different construction from the crossed-layer proposal; no measured cost
or certificate transfers. It tests whether changed column wiring alone
can support the cheaper first mixers under the same regional measure.
"""
from hashlib import sha256
import json
from math import log2
from pathlib import Path
from time import perf_counter
import numpy as np

from hd_pair_scheduled_control import load_measure, prepare_data
from hd_pair_identity_kernel import mask_log_matrices, branch_log_matrices
from hd_pair_state_tail import word_tail
from hd_pair_singleton_word import whole_word_roots

DIRECTORY=Path(__file__).resolve().parent
CASES=(('Hbin','input-xor'),('Hbin','input-inverse'),
       ('lambda','static'),('lambda','input-xor'),('lambda','input-inverse'))


def run(*,max_seconds=120,progress=None):
    started=perf_counter()
    def check_time():
        if perf_counter()-started>max_seconds:
            raise TimeoutError('Prescribed no-pair screen exceeded its budget')
    saved,measure,prefix=load_measure()
    table=mask_log_matrices(measure,18)
    cases=[]
    for mixer,schedule in CASES:
        arities=[]
        for arity in (1,2):
            data=prepare_data(arity,schedule,mixer)
            blocks=[]
            for b in range(8):
                check_time()
                matrix=branch_log_matrices(data['first'][b],data['second'][b],table)
                blocks.append(matrix.reshape(len(data['sources']),data['branch_count'],2,2))
            norms=word_tail(blocks,data['log_counts'],data['input_weights'],18,9*arity*18,
                            check_time=check_time)/18
            methods={'aggregate_tail':norms.copy()}
            if arity==1:
                whole=whole_word_roots(blocks,18,check_time=check_time)
                methods['whole_word']=whole
                norms=np.minimum(norms,whole)
            roots=norms-8*arity*log2(255)
            row=dict(arity=arity,labelled_rowsets=len(data['sources']),
                method_bounds_log2={method:prefix+18*float(np.logaddexp2.reduce(v-8*arity*log2(255)))
                                   for method,v in methods.items()},
                bound_log2=prefix+18*float(np.logaddexp2.reduce(roots)),
                largest_fixed_A=list(data['sources'][int(np.argmax(roots))]))
            arities.append(row)
            if progress:
                progress(dict(mixer=mixer,schedule=schedule,**row,elapsed_seconds=perf_counter()-started))
        mixed=prefix+18*float(np.logaddexp2.reduce([(row['bound_log2']-prefix)/18 for row in arities]))
        cases.append(dict(mixer=mixer,schedule=schedule,arities=arities,mixed_arities_1_2_log2=mixed))
    names=('hd_scheduled_no_pair_control.py','hd_pair_scheduled_control.py','hd_pair_scheduled_outer.py',
        'hd_lambda_outer.py','hd_pair_identity_kernel.py','hd_pair_state_kernel.py',
        'hd_pair_state_tail.py','hd_pair_singleton_word.py','hd_pair_scalar_control.py')
    return dict(outer_pair_layer=False,global_paired_route_unchanged=True,
        lambda_value=2,groups=18,source_class='Fmax<=1',mark=saved['mark'],
        source_count=saved['source_count'],common_measure=measure,common_prefix_log2=prefix,
        controls=cases,directed=False,whole_code_certificate=False,
        native_implementation_changed=False,elapsed_seconds=perf_counter()-started,
        source_pins={name:sha256((DIRECTORY/name).read_bytes()).hexdigest() for name in names})


if __name__=='__main__':
    import argparse
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args()
    if args.output.exists():parser.error('Refuse to overwrite evidence')
    result=run(progress=lambda row:print(json.dumps(row),flush=True))
    with args.output.open('x',encoding='utf-8') as handle:
        json.dump(result,handle,indent=2,allow_nan=False);handle.write('\n')
    print(json.dumps(dict(elapsed_seconds=result['elapsed_seconds'],
        mixed=[(r['mixer'],r['schedule'],r['mixed_arities_1_2_log2']) for r in result['controls']])),flush=True)
