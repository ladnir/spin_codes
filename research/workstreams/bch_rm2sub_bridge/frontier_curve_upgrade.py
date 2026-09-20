"""Add a full certificate to a previously diagnostic-only frontier point."""
import argparse
from fractions import Fraction as F
import math
from pathlib import Path
import bridge as base
import frontier_ledger as ledger
import frontier_curve_extend as previous_model


def upgrade(rows,m,margin,model_bits):
    result={r['message_exponent']:dict(r) for r in rows}
    assert m in result and result[m]['full_certificate_margin_bits'] is None
    prior=dict(result[m])
    result[m]=dict(message_exponent=m,message_bits=1<<m,full_certificate_margin_bits=margin,
        partial_certificate=None,q1_constraint_diagnostic_bits=None,modeled_q1_bits=model_bits,
        modeled_full_margin_bits=None,prior_uncertified_point=prior)
    return [result[k] for k in sorted(result)]


def run(output,previous,certificate):
    output=output.resolve();previous=previous.resolve();certificate=certificate.resolve()
    assert not output.exists();old=base.read(previous);full=base.read(certificate)
    assert full['full_distance_proved'] is True
    m=full['message_exponent']
    assert ledger.build(m,[base.ROOT/r['file'] for r in full['certificates']])==full
    bound=base.decode(full['failure_upper']);assert bound<F(1,1<<40)
    margin=math.log2(bound.denominator)-math.log2(bound.numerator);assert margin==full['margin_bits']
    point=previous_model.modeled_point(m)
    rows=upgrade(old['rows'],m,margin,point['margin_bits'])
    certified=[(r['message_exponent'],r['full_certificate_margin_bits']) for r in rows if r['full_certificate_margin_bits'] is not None]
    slopes=[dict(from_exponent=a,to_exponent=b,bits_lost_per_doubling=(x-y)/(b-a))
            for (a,x),(b,y) in zip(certified,certified[1:])]
    model=previous_model.model
    paths=[previous,certificate,Path(__file__),Path(ledger.__file__),Path(previous_model.__file__),
        Path(model.__file__),Path(model.refresh.__file__),Path(model.maps.__file__),
        base.HERE/'generated/larger_state_inputs_v1/t64_s20_selection.json',base.HERE/'MIGRATION_MANIFEST.json']
    base.write_new(output,dict(status='SEPARATED_SIZE_MARGIN_FRONTIER',configuration='t64_s20',rows=rows,
        adjacent_certified_slopes=slopes,added_q1_model=point,limitations=old['limitations'],
        source_sha256={p.relative_to(base.ROOT).as_posix():base.sha(p) for p in paths}))
    print('K',m,'full margin',margin,'modeled Q1 only',point['margin_bits'],'slopes',slopes,flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True)
    p.add_argument('--previous',type=Path,required=True);p.add_argument('--certificate',type=Path,required=True)
    a=p.parse_args();run(a.output,a.previous,a.certificate)
