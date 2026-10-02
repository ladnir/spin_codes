"""Retarget completed pair-mixer POINT diagnostics, not a certificate.

The source receipt's outward endpoints are retained evidence from the local
screen, not freshly replayed proofs. The only new operation is a rigorous
cutoff adjustment exp(lambda*(D_new-D_old)). Optional fixed-inner ablation
regenerates the inner and checks the exact stored shell majorant, but does
not authenticate the stored outer CDF. No whole-code claim is supported.
"""
import argparse
import hashlib
import json
from fractions import Fraction as Q
from pathlib import Path
from flint import arb, ctx


def retarget(point, old, new):
    m,e = point['upper']
    lam = Q(point['witness']['parameters'][0])
    if (any(type(v) is not int for v in (m,e,old,new)) or m <= 0 or lam <= 0
            or not 0 <= old < 1<<21 or not 0 <= new < 1<<21):
        raise ValueError('positive dyadic bound and tilt, valid integer cutoffs required')
    value = arb(m)*arb(2)**e*(arb(lam.numerator)*(new-old)/lam.denominator).exp()
    return dict(threshold=new, upper=[int(v) for v in value.upper().man_exp()],
        log2_upper=str(value.log()/arb(2).log()))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('source', type=Path)
    parser.add_argument('--distances', nargs='+', default=['.095','.09'])
    parser.add_argument('--updates', type=int, choices=(1,2,3,4))
    parser.add_argument('--mean', default='.096')
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists(): parser.error('new output path required')
    raw = args.source.read_bytes(); source = json.loads(raw)
    if source['schema'] != 'shared-gf16-pair-mixer-diagnostic-1' or not source['probes']:
        parser.error('completed pair-mixer point diagnostics required')
    ctx.prec = 256
    old = int(Q(source['distance'])*(1<<21))
    result = dict(schema='shared-gf16-pair-mixer-retarget-1',
        source_sha256=hashlib.sha256(raw).hexdigest(), precision=256,
        note=__doc__, points=[])
    for point in source['probes']:
        if point.get('empty'): continue
        row = dict(mean=point['mean'], source_log2=point['log2_upper'], retargets=[])
        for distance in map(Q,args.distances):
            if not 0 < distance < Q(1,2): parser.error('distances in (0,1/2) required')
            value = dict(distance=str(distance), **retarget(point,old,int(distance*(1<<21))))
            row['retargets'].append(value)
            print('PAIR RETARGET', point['mean'], str(distance), value['log2_upper'], flush=True)
        result['points'].append(row)
    def save():
        args.output.parent.mkdir(parents=True,exist_ok=True)
        args.output.write_text(json.dumps(result,indent=2)+'\n')
    save()
    if args.updates is not None:
        import shared_mixture
        import birth_classes
        import scalar_cover as sc
        mixture = [(Q(v['mass']),Q(v['activity'])) for v in source['mixture']]
        shared_mixture.verify(list(map(Q,source['caps'])),mixture)
        data = birth_classes.actual(args.updates)
        ctx.prec = 256
        model = sc.Model(shared_mixture.as_components(mixture),data,old,33,Q(3,16),
            inner=birth_classes,variance_shuffle=True,variance_bins=16,regional_count=True)
        mean = Q(args.mean)
        point = next(v for v in source['probes'] if Q(v['mean']) == mean)
        upper = model.outward((mean,mean),point['witness'])
        result['fixed_inner_ablation'] = dict(mean=str(mean),updates=args.updates,
            distance=source['distance'],witness=point['witness'],
            upper=[int(v) for v in upper.upper().man_exp()],
            log2_upper=str(upper.log()/arb(2).log()))
        print('PAIR FIXED INNER', args.updates, str(mean), result['fixed_inner_ablation']['log2_upper'],flush=True)
        save()


if __name__ == '__main__': main()
