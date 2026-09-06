"""Extend the labeled frontier with a full certificate and a Q1 model point."""
import argparse
from fractions import Fraction as F
import math
from pathlib import Path
import numpy as np
import bridge as base
import dual_track_q1 as model


def modeled_point(m):
    _,_,spectrum,_=model.maps.load('t64_s20');rows=(1<<m)//128
    grid=np.arange(-40,81)/10
    coarse=model.values('refresh',spectrum,rows,grid);chosen=np.argmin(coarse,axis=0)
    fine=np.unique(np.concatenate([grid[chosen[w]]+np.arange(-30,31)/500 for w in (38,40,42,44,46,48,50)]))
    coefficients=np.concatenate((coarse,model.values('refresh',spectrum,rows,fine))).min(axis=0)
    counts=model.even_binomial();value=model.aggregate(counts,coefficients)
    return dict(margin_bits=value['margin_bits'],dominant_weight=value['dominant_weight'],
        status='BINARY64_MODELED_Q1_ONLY',modeled_full_margin_bits=None)


def run(output,previous,certificate):
    assert not output.exists();old=base.read(previous);full=base.read(certificate)
    assert full['full_distance_proved'] is True
    for name,digest in full['input_sha256'].items():assert base.sha(base.ROOT/name)==digest
    m=full['message_exponent'];bound=base.decode(full['failure_upper'])
    margin=math.log2(bound.denominator)-math.log2(bound.numerator)
    assert margin==full['margin_bits'] and bound<F(1,1<<40)
    point=modeled_point(m);rows={r['message_exponent']:r for r in old['rows']}
    assert m not in rows,'Do not overwrite a frontier point'
    rows[m]=dict(message_exponent=m,message_bits=1<<m,full_certificate_margin_bits=margin,
        partial_certificate=None,q1_constraint_diagnostic_bits=None,modeled_q1_bits=point['margin_bits'],
        modeled_full_margin_bits=None)
    certified=sorted((k,r['full_certificate_margin_bits']) for k,r in rows.items() if r['full_certificate_margin_bits'] is not None)
    slopes=[dict(from_exponent=a,to_exponent=b,bits_lost_per_doubling=(x-y)/(b-a))
            for (a,x),(b,y) in zip(certified,certified[1:])]
    paths=[previous.resolve(),certificate.resolve(),Path(__file__),Path(model.__file__),Path(model.refresh.__file__),
           Path(model.maps.__file__),base.HERE/'generated/larger_state_inputs_v1/t64_s20_selection.json',
           base.HERE/'MIGRATION_MANIFEST.json']
    base.write_new(output,dict(status='SEPARATED_SIZE_MARGIN_FRONTIER',configuration='t64_s20',
        rows=[rows[k] for k in sorted(rows)],adjacent_certified_slopes=slopes,
        added_q1_model=point,limitations=['Full certificates and modeled Q1 values are different quantities.',
            'No modeled full-margin curve or statistical confidence interval is claimed.',
            'Finite-difference slopes are descriptive, not an extrapolation theorem.',
            'K20 uses its original stronger cutoff; larger certificates use floor(N/10).',
            'Rounded higher-occupancy certificates can be looser than the modeled Q1 contribution.'],
        source_sha256={p.relative_to(base.ROOT).as_posix():base.sha(p) for p in paths}))
    print('K',m,'full margin',margin,'modeled Q1',point['margin_bits'],'slopes',slopes,flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True)
    p.add_argument('--previous',type=Path,required=True);p.add_argument('--certificate',type=Path,required=True)
    a=p.parse_args();run(a.output,a.previous,a.certificate)
