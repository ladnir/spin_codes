"""Fresh pointwise regional bounds with jointly constrained count laws.

The normalized law is the base-tilted Poisson-binomial count. The factor
tau**(-J) belongs in each objective, never in the normalization constraint.
Unselected variance parts retain their original outward bounds. Thus each
completed point still includes the entire variance partition, not just the
parts whose LP bounds improved. A point is not a full dense certificate.
"""
import argparse
import copy
from fractions import Fraction as Q
import hashlib
import json
from pathlib import Path
import sys

from flint import arb, arb_mat, ctx
import cell_search

SCHEMA = 'packed-gl32-joint-count-point-1'


def validate_scope(scope):
    expected = dict(schema='packed-canonical-gl32-hill-dense-context-1',
        ensemble='canonical-gl32-width8-shared4-r2', K=1 << 20, N=1 << 21,
        updates=2, block_width=8, minimum_groups=33, maximum_groups=2048,
        comparison='direct-expected-shell-majorant', last_lp=104, refined=True,
        base_tilt='3/16', regional_count=True)
    if (any(scope.get(k) != value for k, value in expected.items())
            or type(scope.get('variance_bins')) is not int
            or not 1 <= scope['variance_bins'] <= 64
            or scope.get('threshold') != cell_search.dense.cutoff(scope.get('distance', 0))):
        raise ValueError('matching canonical GL32 R2 refined-count context required')


def fresh_model(scope, precision):
    validate_scope(scope)
    from outer_hill_incidence import authenticated_bch_cdf
    from monotone import transport_shells
    from local_models import full_block
    caps, premises = authenticated_bch_cdf()
    shells = transport_shells(premises['canonical_cdf'], full_block(8))
    if (premises != scope['outer_premises']
            or cell_search.dense.fingerprint(caps) != scope['expected_cdf_sha256']
            or cell_search.dense.fingerprint(shells) != scope['comparison_caps_sha256']):
        raise ValueError('fresh count premises do not match the point source')
    parent = Path(__file__).resolve().parents[1]
    sys.path.insert(0, str(parent))
    sys.path.insert(0, str(parent/'gf16_packets'))
    mixture, _ = cell_search.dense.exact_mixture(shells, scope['mixture'])
    import shared_mixture
    import scalar_cover as sc
    import birth_classes
    data = birth_classes.actual(2)
    ctx.prec = precision
    model = sc.Model(shared_mixture.as_components(mixture), data, scope['threshold'],
        scope['minimum_groups'], Q(3, 16), inner=birth_classes,
        variance_shuffle=True, variance_bins=scope['variance_bins'], regional_count=True)
    if tuple(map(Q, scope['root'])) != model.root:
        raise ValueError('fresh comparison domain differs from the source')
    return model


def count_moments(groups, cell, interval):
    """Mean and centered-second-moment inequalities for the same count law."""
    lo, hi = map(Q, cell)
    v_lo, v_hi = map(Q, interval)
    if (type(groups) is not int or groups < 1 or not 0 < lo <= hi < 1
            or not 0 <= v_lo <= v_hi <= Q(1, 4)):
        raise ValueError('valid mean/variance intervals and group count required')
    center = groups*(lo+hi)/2
    counts = tuple(map(Q, range(groups+1)))
    radius2 = max((groups*lo-center)**2, (groups*hi-center)**2)
    return [(counts, groups*hi),
        (tuple(-j for j in counts), -groups*lo),
        (tuple((j-center)**2 for j in counts), groups*v_hi+radius2),
        (tuple(-(j-center)**2 for j in counts), -groups*v_lo)]


def selected_parts(selection, count):
    selected = list(range(count)) if selection is None else list(selection)
    if (len(set(selected)) != len(selected)
            or any(type(i) is not int or not 0 <= i < count for i in selected)):
        raise ValueError('distinct valid variance-part indices required')
    return set(selected)


def outward(model, cell, witness, *, selection=None, moments=False, replay=None, progress=None):
    import regional_count as rc
    import variance_partition as variance
    import scalar_cover as sc
    import joint_count_dual as joint
    if witness.get('regional_direct_counts') is not True:
        raise ValueError('direct-count regional witness required')
    if type(moments) is not bool or Q(witness['parameters'][0]) <= 0:
        raise ValueError('boolean moment option and positive output tilt required')
    parts, checked = rc.prepare_witness(model, cell, witness)
    selected = selected_parts(selection, len(parts))
    if replay is not None and (len(replay) != len(parts)
            or any((row is not None) != (i in selected) for i, row in enumerate(replay))):
        raise ValueError('one matching joint-dual block per selected variance part required')
    precision = ctx.prec
    region = rc.placement(rc.local_operators(model, witness))
    if len(region) != sc.G+1:
        raise ArithmeticError('regional count geometry mismatch')
    size = region[0].nrows()
    # The entry objectives are shared across every variance part.
    count_weights = [rc.aq(model.tilt)**(-j) for j in range(sc.G+1)]
    values = [[[rc.up(weight*matrix[a, b])
                for weight, matrix in zip(count_weights, region)]
                for b in range(size)] for a in range(size)]
    cs, _, _ = model.family(model.tilt)
    lam = Q(witness['parameters'][0])
    total, old_total = arb(0), arb(0)
    rows, dual_blocks = [], []
    for index, ((interval, dual), part) in enumerate(zip(parts, checked['regional_count_parts'])):
        cap = variance.factor(model, cell, interval[0], witness['variance_dual'])
        atom_caps = rc.count_mass_caps(model.features, model.active, cell, interval,
            model.q_min, sc.G, rc.aq(cap), part['mgf_witnesses'])
        old_matrix = arb_mat([[rc.up(sum((c*v for c, v in zip(atom_caps, values[a][b])), arb(0)))
            for b in range(size)] for a in range(size)])
        matrix, block = old_matrix, None
        if index in selected:
            # mgf_upper returns the raw log MGF. The atom correction inside
            # count_mass_caps is NOT a valid multiplier on this expectation.
            raw_mgfs = [(Q(w['tilt']), rc.mgf_upper(model.features, model.active,
                cell, interval, model.q_min, sc.G, w)) for w in part['mgf_witnesses']]
            prepared = joint.prepare(atom_caps, raw_mgfs,
                moment_bounds=count_moments(sc.G, cell, interval) if moments else ())
            block = [[None]*size for _ in range(size)]
            entries = [[arb(0)]*size for _ in range(size)]
            if replay is not None and (len(replay[index]) != size
                    or any(len(row) != size for row in replay[index])):
                raise ValueError('square matrix of joint duals required')
            for a in range(size):
                for b in range(size):
                    if replay is None:
                        upper, entry_witness = prepared.bound(values[a][b])
                    else:
                        entry_witness = replay[index][a][b]
                        upper = prepared.verify(values[a][b], entry_witness)
                    if not upper.is_finite() or upper < 0:
                        raise ArithmeticError('nonnegative finite joint bound required')
                    entries[a][b] = min(rc.up(upper), old_matrix[a, b])
                    block[a][b] = entry_witness
            matrix = arb_mat(entries)
        power, old_power = matrix**sc.REGIONS, old_matrix**sc.REGIONS
        moment = sum((power[0, j] for j in range(size)), arb(0))
        old_moment = sum((old_power[0, j] for j in range(size)), arb(0))
        eta, mu, gamma = dual
        count = sum((rc.aq(c)*rc.aq(eta*f+mu*a+gamma*f*(1-f)).exp()
            for c, f, a in zip(cs, model.features, model.active)), arb(0))
        exponent = sc.G*(min(eta*x for x in cell)+min(gamma*v for v in interval))+mu*model.q_min
        factor = (sc.G*count.log()-rc.aq(exponent)+rc.aq(lam)*model.threshold).exp()
        contribution, old_contribution = rc.up(factor*moment), rc.up(factor*old_moment)
        total, old_total = rc.up(total+contribution), rc.up(old_total+old_contribution)
        rows.append(dict(index=index, interval=list(map(str, interval)), joint=index in selected,
            log2_upper=str(contribution.log()/arb(2).log()) if contribution else '-inf',
            baseline_log2=str(old_contribution.log()/arb(2).log()) if old_contribution else '-inf'))
        dual_blocks.append(block)
        if progress:
            progress(rows[-1], dual_blocks)
    if ctx.prec != precision:
        raise ArithmeticError('precision changed during joint-count replay')
    return total, dict(base_witness=checked, selected_parts=sorted(selected), moments=moments,
        joint_duals=dual_blocks, parts=rows,
        baseline_upper=[int(v) for v in old_total.upper().man_exp()])


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('source', type=Path)
    parser.add_argument('--point', type=int, default=0)
    parser.add_argument('--parts', type=int, nargs='*')
    parser.add_argument('--moments', action='store_true')
    parser.add_argument('--precision', type=int, default=256)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists() or args.precision < 128 or args.point < 0:
        parser.error('fresh output, nonnegative point index and precision >=128 required')
    raw = args.source.read_bytes()
    source = json.loads(raw)
    if source.get('schema') != 'packed-canonical-gl32-dense-hill-points-1':
        parser.error('hill point receipt required')
    scope = source['scope']
    point = source['points'][args.point]
    cell = tuple(map(Q, point['cell']))
    witness = copy.deepcopy(point['checked']['witness'])
    model = fresh_model(scope, args.precision)
    if not model.root[0] <= cell[0] <= cell[1] <= model.root[1]:
        parser.error('point outside comparison domain')
    record = dict(schema=SCHEMA, scope=scope, cell=list(map(str, cell)),
        precision=args.precision, source=dict(path=str(args.source.resolve()),
            sha256=hashlib.sha256(raw).hexdigest(), point=args.point),
        partial_only=True, complete_point=False, whole_code_certificate=False,
        selected_parts=args.parts, moments=args.moments, parts=[])

    def save():
        if args.source.read_bytes() != raw:
            raise ValueError('point source changed during verification')
        args.output.parent.mkdir(parents=True, exist_ok=True)
        temporary = args.output.with_name(args.output.name+'.tmp')
        temporary.write_text(json.dumps(record, indent=2)+'\n')
        temporary.replace(args.output)

    def progress(row, blocks):
        record['parts'].append(row)
        record['joint_duals'] = blocks
        save()
        print('JOINT COUNT PART', row, flush=True)

    save()
    upper, checked = outward(model, cell, witness, selection=args.parts,
        moments=args.moments, progress=progress)
    record.update(complete_point=True, upper=[int(v) for v in upper.upper().man_exp()],
        log2_upper=str(upper.log()/arb(2).log()), checked=checked)
    save()
    print('JOINT COUNT POINT', record['log2_upper'], flush=True)


if __name__ == '__main__':
    main()
