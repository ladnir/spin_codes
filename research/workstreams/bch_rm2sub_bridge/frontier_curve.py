"""Separate full certificates, partial certificates, and modeled Q1 values."""
import argparse
import math
from pathlib import Path
import bridge as base


def run(output):
    assert not output.exists()
    names=['refresh_q1_outward_v1.json','refresh_q1_outward_v1_replay.json',
           'frontier_k24_full_v1.json','dual_track_q1_v1.json','k30_partial_coverage_v1.json']
    paths=[base.HERE/'generated'/n for n in names]
    old,replay,new,diagnostic,partial=map(base.read,paths)
    assert replay['producer_sha256']==base.sha(paths[0])
    assert new['full_distance_proved'] is True
    for name,digest in new['input_sha256'].items():assert base.sha(base.ROOT/name)==digest
    for name,digest in old['source_sha256'].items():assert base.sha(base.ROOT/name)==digest
    assert partial['full_distance_proved'] is False
    rows={r['message_exponent']:dict(message_exponent=r['message_exponent'],message_bits=1<<r['message_exponent'],
        full_certificate_margin_bits=None,partial_certificate=None,
        q1_constraint_diagnostic_bits=r['refresh_caps_margin_bits'],modeled_q1_bits=r['modeled']['margin_bits'],
        modeled_full_margin_bits=None) for r in diagnostic['rows']}
    value=base.decode(old['combined_full_upper'])
    rows[20]['full_certificate_margin_bits']=math.log2(value.denominator)-math.log2(value.numerator)
    rows[24]['full_certificate_margin_bits']=new['margin_bits']
    rows[30]=dict(message_exponent=30,message_bits=1<<30,full_certificate_margin_bits=None,
        partial_certificate=dict(occupancy_range=partial['covered_occupancy_range'],margin_bits=partial['partial_margin_bits']),
        q1_constraint_diagnostic_bits=None,modeled_q1_bits=None,modeled_full_margin_bits=None)
    slope=(rows[20]['full_certificate_margin_bits']-rows[24]['full_certificate_margin_bits'])/4
    base.write_new(output,dict(status='SEPARATED_SIZE_MARGIN_FRONTIER',configuration='t64_s20',
        rows=[rows[m] for m in sorted(rows)],observed_certified_bits_lost_per_doubling_20_to_24=slope,
        limitations=['Only K20 and K24 have full certificates in this table.',
            'Modeled values are Q1 diagnostics, not a complete failure curve or confidence interval.',
            'The two-point certified slope is descriptive, not an extrapolation theorem.',
            'K20 full certificate uses cutoff 209716; its Q1 diagnostic uses floor(N/10)=209715.'],
        input_sha256={p.relative_to(base.ROOT).as_posix():base.sha(p) for p in paths},
        source_sha256=base.sha(Path(__file__))))
    print('Full certified bits lost per doubling, K20 to K24:',slope,flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True)
    run(p.parse_args().output)
