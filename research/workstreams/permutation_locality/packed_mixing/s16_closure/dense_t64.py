"""Fresh dense bounds for two physical t64/s16 steps per t128 macro.

Only outer comparison geometry and exact count machinery are retained.
The inner map, physical-step law and every numerical bound are rebuilt.
Point searches do not certify a domain. A search cover still needs replay.
"""
import argparse
import copy
from fractions import Fraction as Q
import hashlib
import json
from pathlib import Path
from time import monotonic
from types import FunctionType

from flint import arb, ctx
import dense_model as prior
import kernel_t64
import regional_count as regional
import scalar_cover as sc
import shared_mixture
import variance_partition
import search_strategy as geometry

DEFAULT_SOURCE = Path('tmp/packed-hill/dense-r4-ekr-d10-complete-p256.json')
SCHEMA = 'packed-gl32-t64-s16-dense-search-1'


def local_operators(model, witness, *, _proposal_scratch=None):
    kernel_t64.authenticate(model.data)
    lam = Q(witness['parameters'][0])
    activity = Q(witness.get('t64_birth_density_activity', '1/2'))
    if lam <= 0 or not 0 <= activity <= 1:
        raise ValueError('positive output tilt and valid row-selection activity required')
    key = (model.data['map_sha256'], ctx.prec, lam, activity,
           model.data.get('birth_density', 'capped'))
    saved = None if _proposal_scratch is None else _proposal_scratch.get('t64_macro_local')
    if saved is None or saved[0] != key:
        saved = key, kernel_t64.local_operators(model.data, lam, activity)
        if _proposal_scratch is not None:
            _proposal_scratch['t64_macro_local'] = saved
    return saved[1]


_bindings = dict(regional.__dict__, local_operators=local_operators)
_regional_propose = FunctionType(regional.propose.__code__, _bindings,
    't64_regional_propose', regional.propose.__defaults__)
_regional_outward = FunctionType(regional.outward.__code__, _bindings,
    't64_regional_outward', regional.outward.__defaults__)


class Model(sc.Model):
    def __init__(self, mixture, data, threshold=209715, variance_bins=16):
        super().__init__(shared_mixture.as_components(mixture), data, threshold,
            33, Q(3, 16), inner=kernel_t64, variance_shuffle=True,
            variance_bins=variance_bins, regional_count=False)
        self.regional_count = True
        self.proposal_stop_bits = 54

    def outward(self, cell, witness):
        if 'regional_count_parts' in witness:
            return _regional_outward(self, cell, witness)[0]
        return super().outward(cell, witness)

    def candidates(self, cell, scales=(Q(1), Q(9,10), Q(4,5), Q(39,40)),
                   callback=None, stop=False):
        best = self.propose_with(cell, self.tilt)
        if callback:
            callback(dict(method='iid', score=float(best[0])), best)
        if stop and best[0] < -self.proposal_stop_bits:
            return best
        if not 0 < cell[0] <= cell[1] < 1:
            return best
        variance = variance_partition.propose(self, cell, best)
        if variance[0] < best[0]:
            best = variance
        if callback:
            callback(dict(method='variance', score=float(variance[0])), best)
        if stop and best[0] < -self.proposal_stop_bits or cell[1]-cell[0] > Q(1,1024):
            return best
        witness = copy.deepcopy(variance[1])
        weights, _ = self.weights(cell, self.tilt)
        witness.update(regional_tilted_atom=True, regional_fine_tilts=True,
            regional_tilted_variance=True, regional_direct_counts=True,
            t64_birth_density_activity=str(weights[1]/sum(weights)))
        lam = Q(witness['parameters'][0])
        for scale in scales:
            if Q(scale) <= 0:
                raise ValueError('positive output-tilt scale required')
            trial = copy.deepcopy(witness)
            trial['parameters'][0] = str(lam*Q(scale))
            score, checked = _regional_propose(self, cell, trial)
            if score < best[0]:
                best = score, checked
            if callback:
                callback(dict(method='regional', scale=str(scale),
                    output_tilt=trial['parameters'][0], score=float(score)), best)
            if stop and best[0] < -self.proposal_stop_bits:
                break
        return best

    def proposal(self, cell):
        return self.candidates(cell, stop=True)


def fresh_model(source=DEFAULT_SOURCE, *, precision=256, variance_bins=16,
                distance=Q(1,10), birth_density='capped'):
    if (type(precision) is not int or precision < 128
            or type(variance_bins) is not int or not 1 <= variance_bins <= 64
            or not 0 < Q(distance) < Q(1,2)
            or birth_density not in ('classes','capped')):
        raise ValueError('valid exact claim, precision, variance and density option required')
    source = Path(source).resolve(); raw = source.read_bytes()
    old = json.loads(raw)['scope']
    if (old.get('K') != 1 << 20 or old.get('N') != 1 << 21
            or old.get('block_width') != 8 or old.get('minimum_groups') != 33
            or old.get('maximum_groups') != 2048
            or old.get('comparison') != 'direct-expected-shell-majorant'):
        raise ValueError('matching canonical GL32 outer-mixture proposal required')
    print('T64 freshly authenticating outer count premises', flush=True)
    caps, premises = prior.authenticated_bch_cdf()
    comparison = prior.transport_shells(premises['canonical_cdf'], prior.full_block(8))
    fingerprint = prior.cell_search.dense.fingerprint
    if (premises != old['outer_premises'] or fingerprint(caps) != old['expected_cdf_sha256']
            or fingerprint(comparison) != old['comparison_caps_sha256']):
        raise ValueError('fresh outer premises differ from proposed mixture scope')
    mixture, verification = prior.cell_search.dense.exact_mixture(comparison, old['mixture'])
    print('T64 outer comparison checked; preparing physical maps and macro kernel', flush=True)
    data, maps = kernel_t64.prepare(birth_density=birth_density)
    ctx.prec = precision
    threshold = int(Q(distance)*(1 << 21))
    model = Model(mixture, data, threshold, variance_bins)
    scope = dict(schema='packed-gl32-t64-s16-dense-context-1',
        K=1 << 20, N=1 << 21, physical_t=64, state_bits=16,
        physical_steps=32768, macro_t=128, macro_windows=32, macro_steps=16384,
        physical_steps_per_macro=2, regions=256, macros_per_region=64,
        outer='fixed BCH[256,128]', block_rows=4, block_width=8, groups=2048,
        mixing='independent uniform GL32 on each canonical four-row/eight-column block',
        routing='shared uniform column shuffle per four-row group; independent uniform packet shuffle per region',
        inner_recurrence='y=x+Aq; next_q=Mq+Cx; independent uniform GL16 each physical step; zero initial state; no flush',
        maps=json.loads(json.dumps(maps)), maps_sha256=prior.fingerprint(maps),
        birth_density=birth_density, distance=str(Q(distance)), threshold=threshold,
        minimum_groups=33, maximum_groups=2048, root=list(map(str,model.root)),
        base_tilt=str(model.tilt), variance_bins=variance_bins,
        comparison='direct-expected-shell-majorant', outer_premises=premises,
        expected_cdf_sha256=fingerprint(caps), comparison_caps_sha256=fingerprint(comparison),
        mixture=[dict(mass=str(c), activity=str(p)) for c,p in mixture],
        mixture_verification=verification, outer_authenticated=True)
    if source.read_bytes() != raw:
        raise ValueError('outer proposal source changed during preparation')
    return model, scope, dict(path=str(source), sha256=hashlib.sha256(raw).hexdigest())


def endpoint(value):
    if not value.is_finite() or not value > 0:
        raise ArithmeticError('positive finite outward endpoint required')
    return list(map(int, value.upper().man_exp()))


def read_record(path):
    raw = Path(path).read_bytes()
    return json.loads(raw), dict(path=str(Path(path).resolve()),
        sha256=hashlib.sha256(raw).hexdigest())


def unchanged(source):
    if hashlib.sha256(Path(source['path']).read_bytes()).hexdigest() != source['sha256']:
        raise ValueError('source changed during computation')


def validate_cover_record(record, complete=False):
    if record.get('schema') != SCHEMA or record['scope'].get('schema') != 'packed-gl32-t64-s16-dense-context-1':
        raise ValueError('explicit t64/s16 dense search scope required')
    cover = record['cover']
    cells = geometry.partition(record['scope']['root'], cover['leaves'], cover['unresolved'])
    if complete and cover['unresolved']:
        raise ValueError('full dense domain must be covered before replay')
    if any(not isinstance(row.get('witness'),dict) for row in cover['leaves'].values()):
        raise ValueError('rational witnesses required for accepted cells')
    return cells


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('mode', choices=('points','search','replay'))
    parser.add_argument('--source', type=Path, default=DEFAULT_SOURCE)
    parser.add_argument('--input', type=Path, help='complete t64 search to replay')
    parser.add_argument('--geometry', type=Path, help='reuse only a complete exact partition; not old witnesses/bounds')
    parser.add_argument('--means', nargs='+', default=['.032','.104','.114'])
    parser.add_argument('--scales', nargs='+', default=['1','9/10','4/5'])
    parser.add_argument('--precision', type=int, default=256)
    parser.add_argument('--variance-bins', type=int, default=16)
    parser.add_argument('--distance', type=Q, default=Q(1,10))
    parser.add_argument('--birth-density', choices=('classes','capped'), default='capped')
    parser.add_argument('--outward', action='store_true')
    parser.add_argument('--cell-target-bits', type=int, default=52)
    parser.add_argument('--max-cells', type=int, default=800)
    parser.add_argument('--max-depth', type=int, default=24)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists() or not 20 <= args.cell_target_bits <= 256:
        parser.error('fresh output and valid cell target required')
    if (args.mode == 'replay') != (args.input is not None):
        parser.error('replay requires --input; other modes reject it')
    if args.geometry and args.mode != 'search':
        parser.error('--geometry is search-only')
    start = monotonic()
    model, scope, source = fresh_model(args.source, precision=args.precision,
        variance_bins=args.variance_bins, distance=args.distance, birth_density=args.birth_density)
    model.proposal_stop_bits = args.cell_target_bits+2
    record = dict(schema=SCHEMA, scope=scope, source=source, precision=args.precision,
        cell_target_bits=args.cell_target_bits, whole_code_certificate=False)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    def save():
        unchanged(source)
        record['elapsed_seconds'] = monotonic()-start
        args.output.write_text(json.dumps(record,indent=2)+'\n')
    if args.mode == 'points':
        record.update(schema='packed-gl32-t64-s16-dense-points-1', points=[],
            proof_status='Selected points only; no dense-domain or whole-code certificate.')
        save()
        for mean in map(Q,args.means):
            if not model.root[0] <= mean < model.root[1]:
                raise ValueError('point outside interior mean domain')
            point = dict(mean=str(mean),trials=[]); record['points'].append(point)
            def checkpoint(row, best):
                point['trials'].append(row); point['best_proposal'] = float(best[0]); save()
                print('T64 TRIAL',str(mean),row,flush=True)
            ctx.prec = args.precision
            score, witness = model.candidates((mean,mean),tuple(map(Q,args.scales)),checkpoint)
            point.update(witness=witness,best_proposal=float(score))
            if args.outward:
                upper = model.outward((mean,mean),witness)
                if ctx.prec != args.precision:
                    raise ArithmeticError('precision changed during point replay')
                point.update(upper=endpoint(upper),log2_upper=str(upper.log()/arb(2).log()))
                print('T64 CHECKED',str(mean),point['log2_upper'],flush=True)
            save()
    elif args.mode == 'search':
        assigned = ['']
        if args.geometry:
            old, original = read_record(args.geometry)
            old_root = old.get('scope',old).get('root')
            if tuple(map(Q,old_root)) != model.root:
                raise ValueError('geometry must cover the same global mean domain')
            cover = old['cover']
            cells = geometry.partition(model.root,cover['leaves'],cover['unresolved'])
            assigned = sorted(cells)
            record['geometry_source'] = original
        def checkpoint(cover):
            geometry.partition(model.root,cover['leaves'],cover['unresolved'])
            record.update(cover=cover,complete_search_partition=not cover['unresolved'],
                final_replay_required=True)
            save()
            print('T64 COVER',cover['visited'],'accepted',len(cover['leaves']),
                  'pending',len(cover['unresolved']),flush=True)
        prior.cell_search.search(model,assigned,precision=args.precision,
            target_bits=args.cell_target_bits,max_cells=args.max_cells,max_depth=args.max_depth,
            checkpoint_every=5,checkpoint=checkpoint)
        if args.geometry:
            unchanged(record['geometry_source'])
    else:
        old, original = read_record(args.input)
        cells = validate_cover_record(old,complete=True)
        if old['scope'] != scope:
            raise ValueError('saved t64 scope differs from freshly authenticated construction/counts')
        record.update(schema='packed-gl32-t64-s16-dense-replay-1',input_source=original,
            cells=[],complete_dense=False,proof_status='Fresh dense replay; sparse proof remains separate.')
        total = arb(0); save()
        for path in sorted(cells):
            ctx.prec = args.precision
            upper = model.outward(cells[path],old['cover']['leaves'][path]['witness'])
            if ctx.prec != args.precision:
                raise ArithmeticError('precision changed during dense replay')
            row = dict(path=path,cell=list(map(str,cells[path])),upper=endpoint(upper))
            total = arb((total+upper).upper())
            record['cells'].append(row); record['aggregate_upper'] = endpoint(total)
            record['complete_dense'] = len(record['cells']) == len(cells)
            unchanged(original); save()
            print('T64 REPLAY',len(record['cells']),'/',len(cells),'log2',upper.log()/arb(2).log(),flush=True)
        record['aggregate_below_2_minus_40'] = bool(total < arb(2)**-40)
        record['aggregate_margin'] = str(-total.log()/arb(2).log()); save()


if __name__ == '__main__':
    main()
