"""Disjoint coverage accounting for unsigned (r,a,b) support profiles.

Geometry, accepted probability mass, and unresolved work are separate.
No number of passing subdomains implies a whole-code certificate while
an unresolved box still contains a feasible profile.
"""
from dataclasses import dataclass
from copy import deepcopy
from fractions import Fraction as F
from math import prod
import heapq

from flint import arb,ctx

from hd_shifted_alias_certificate import rational
from hd_shifted_coverage import rational_ceiling
from hd_shifted_cover_backend import validate_box,profile_count,CachedBound,CONSTRUCTION
from hd_randomized_outer_proof import Parameters,ball_log2,dense_class_log2


@dataclass(frozen=True)
class Box:
    ranges: tuple

    def __post_init__(self):
        ranges = tuple(tuple(row) for row in self.ranges)
        if len(ranges)!=3 or any(len(row)!=2 or any(not isinstance(n,int) for n in row)
                                or row[0]>row[1] for row in ranges):
            raise ValueError('Three nonempty integer intervals required')
        object.__setattr__(self,'ranges',ranges)

    @property
    def volume(self):
        return prod(hi-lo+1 for lo,hi in self.ranges)

    def intersection(self,other):
        ranges = tuple((max(a,c),min(b,d)) for (a,b),(c,d) in zip(self.ranges,other.ranges))
        return None if any(lo>hi for lo,hi in ranges) else Box(ranges)

    def subtract(self,other):
        """At most six disjoint slabs partition self outside other."""
        overlap = self.intersection(other)
        if overlap is None:
            return [self]
        middle,slabs = list(self.ranges),[]
        for axis,(lo,hi) in enumerate(overlap.ranges):
            a,b = middle[axis]
            if a<lo:
                piece = middle.copy()
                piece[axis] = (a,lo-1)
                slabs.append(Box(tuple(piece)))
            if hi<b:
                piece = middle.copy()
                piece[axis] = (hi+1,b)
                slabs.append(Box(tuple(piece)))
            middle[axis] = (lo,hi)
        return slabs

    def split(self):
        widths = [hi-lo for lo,hi in self.ranges]
        axis = max(range(3),key=widths.__getitem__)
        if not widths[axis]:
            raise ValueError('Cannot split a singleton')
        lo,hi = self.ranges[axis]
        mid = (lo+hi)//2
        left,right = list(self.ranges),list(self.ranges)
        left[axis],right[axis] = (lo,mid),(mid+1,hi)
        return Box(tuple(left)),Box(tuple(right))


class Ledger:
    construction = CONSTRUCTION

    def __init__(self,root=Box(((1,8192),(0,4096),(0,4096))),*,pairs=4096):
        validate_box(root.ranges,pairs)
        self.root,self.pairs = root,pairs
        self.unresolved = [root]
        self.accepted,self.empty = [],[]

    def _remove(self,box):
        # A proposed leaf can cross earlier subtraction slabs. Their disjoint
        # intersections must cover it exactly; overlap with already accepted
        # or discarded geometry is rejected before mutating the partition.
        overlaps = [parent.intersection(box) for parent in self.unresolved]
        if sum(piece.volume for piece in overlaps if piece is not None)!=box.volume:
            raise ValueError('New box must lie entirely in unresolved geometry')
        following = []
        for parent,piece in zip(self.unresolved,overlaps):
            following.extend(parent.subtract(box) if piece is not None else [parent])
        self.unresolved = following

    def accept(self,receipt):
        if receipt.get('directed') is not True or receipt.get('covers_all_even_odd_splits') is not True:
            raise ValueError('A directed all-lane-balance receipt is required')
        if (receipt.get('construction')!=self.construction or receipt.get('field_order')!=256
                or receipt.get('forbidden_weight')!=734003):
            raise ValueError('Receipt construction or target does not match this ledger')
        box = Box(tuple(receipt['ranges']))
        F(receipt['log2_upper_rational'])
        self._remove(box)
        self.accepted.append(deepcopy(receipt))

    def discard_empty(self,box):
        if profile_count(box.ranges,self.pairs):
            raise ValueError('Cannot discard a feasible support class')
        self._remove(box)
        self.empty.append(box)

    def split(self,box):
        index = self.unresolved.index(box)
        self.unresolved[index:index+1] = box.split()

    def summary(self,precision=256):
        # Empty pending boxes may remain after geometric subtraction. Report
        # them separately; they never manufacture a full-code certificate.
        feasible = [box for box in self.unresolved if profile_count(box.ranges,self.pairs)]
        geometric = sum(Box(tuple(item['ranges'])).volume for item in self.accepted)
        geometric += sum(box.volume for box in self.empty+self.unresolved)
        if geometric!=self.root.volume:
            raise ArithmeticError('Coverage volume invariant failed')
        with ctx.workprec(precision):
            total = sum((arb(2)**rational(F(item['log2_upper_rational'])) for item in self.accepted),arb(0))
            ceiling = str(rational_ceiling(total.log()/arb(2).log())) if self.accepted else None
            below_target = bool(total<arb(2)**-40)
        complete = not feasible
        global_root = self.root==Box(((1,2*self.pairs),(0,self.pairs),(0,self.pairs)))
        accepted_count = sum(profile_count(tuple(tuple(row) for row in item['ranges']),self.pairs)
                             for item in self.accepted)
        unresolved_count = sum(profile_count(box.ranges,self.pairs) for box in feasible)
        if global_root:
            m = self.pairs//2
            expected = 4*m*(m+1)*(2*m+1)//6+(m+1)**2*(1 if self.pairs%2==0 else 3)-self.pairs
        else:
            expected = profile_count(self.root.ranges,self.pairs)
        if accepted_count+unresolved_count!=expected:
            raise ArithmeticError('Support-profile count invariant failed')
        highest = max((min(box.ranges[0][1],self.pairs+box.ranges[1][1],self.pairs+box.ranges[2][1])
                       for box in feasible),default=self.root.ranges[0][0]-1)
        top = (highest+1,self.root.ranges[0][1]) if highest<self.root.ranges[0][1] else None
        return dict(construction=self.construction,root=self.root.ranges,accepted_boxes=len(self.accepted),
                    unresolved_boxes=len(self.unresolved),unresolved_feasible_boxes=len(feasible),
                    unresolved_profile_triples=unresolved_count,
                    accepted_profile_triples=accepted_count,total_profile_triples=expected,
                    accepted_union_log2_upper=ceiling,complete_local_cover=complete,
                    complete_top_weight_range=top,
                    whole_code_certificate=bool(global_root and self.pairs==4096 and complete and below_target),
                    volume_accounting=True,profile_accounting=True)

    def band_summary(self,low,high,precision=256):
        band = Box(((low,high),(0,self.pairs),(0,self.pairs)))
        if self.root.intersection(band)!=band:
            raise ValueError('Weight band must lie within this ledger root')
        unresolved = sum(profile_count(piece.ranges,self.pairs)
                         for box in self.unresolved if (piece:=box.intersection(band)) is not None)
        accepted,masses = 0,[]
        for receipt in self.accepted:
            piece = Box(tuple(receipt['ranges'])).intersection(band)
            if piece is not None:
                count = profile_count(piece.ranges,self.pairs)
                accepted += count
                if count:
                    # Restricting an event only decreases its probability;
                    # retain the full leaf bound rather than prorating it.
                    masses.append(F(receipt['log2_upper_rational']))
        expected = profile_count(band.ranges,self.pairs)
        if accepted+unresolved!=expected:
            raise ArithmeticError('Weight-band profile accounting failed')
        with ctx.workprec(precision):
            total = sum((arb(2)**rational(value) for value in masses),arb(0))
            ceiling = str(rational_ceiling(total.log()/arb(2).log())) if masses else None
        return dict(weight_range=(low,high),total_profile_triples=expected,
                    accepted_profile_triples=accepted,unresolved_profile_triples=unresolved,
                    complete_local_cover=unresolved==0,accepted_union_log2_upper=ceiling,
                    whole_code_certificate=False)

    def import_band(self,source,low,high,*,label=None):
        """Import a recomputed source ledger, retaining whole leaf bounds.

        The assigned band must be completely unclaimed. Source gaps remain
        gaps; restricting accepted events never prorates their probability.
        """
        band = Box(((low,high),(0,self.pairs),(0,self.pairs)))
        if (source.construction!=self.construction or source.pairs!=self.pairs or self.root.intersection(band)!=band
                or source.root.intersection(band)!=band):
            raise ValueError('Band must belong to both compatible ledger roots')
        available = sum(piece.volume for parent in self.unresolved
                        if (piece:=parent.intersection(band)) is not None)
        if available!=band.volume:
            raise ValueError('Imported band overlaps previously claimed geometry')
        source.summary()  # Check the source partition before importing it.
        staged = deepcopy(self)
        count = 0
        for receipt in source.accepted:
            original = Box(receipt['ranges'])
            piece = original.intersection(band)
            if piece is None or not profile_count(piece.ranges,self.pairs):
                continue
            restricted = deepcopy(receipt)
            restricted['ranges'] = piece.ranges
            restricted['integer_profile_triples'] = profile_count(piece.ranges,self.pairs)
            restricted['restricted_from_ranges'] = original.ranges
            restricted['source_label'] = label
            # No change to log2_upper_rational: this is a subset event.
            staged.accept(restricted)
            count += 1
        staged.summary()
        self.unresolved,self.accepted,self.empty = staged.unresolved,staged.accepted,staged.empty
        return dict(source_label=label,imported_boxes=count,**self.band_summary(low,high))


def terminal_receipt(precision=256):
    p = Parameters(field=256,state=16,limit=734003)
    with ctx.workprec(precision):
        ball = ball_log2(p,mark=F(734003,347602995))
        value = dense_class_log2(7649,p,ball=ball)
        ceiling = rational_ceiling(value)
    return dict(ranges=((7649,8192),(0,4096),(0,4096)),
                log2_upper_rational=str(ceiling),directed=True,covers_all_even_odd_splits=True,
                method='inverse-image-terminal-range',precision=precision,whole_code_certificate=False,
                construction=CONSTRUCTION,field_order=256,forbidden_weight=734003)


def seeded_ledger(precision=256,progress=None):
    ledger = Ledger()
    ledger.accept(terminal_receipt(precision))
    witnesses = (
        (((4032,4160),(1536,1664),(1536,1664)),F('7722805725266083/144115188075855872'),
         (F('244777405/305135501'),F('1189188654/621493969'))),
        (((6016,6272),(2372,2628),(2372,2628)),F('7151016493358955/576460752303423488'),
         (F('973542133/268519977'),F('9096903675/859563932'))),
    )
    alpha = F('4097437497366089/36028797018963968')
    for ranges,z,row in witnesses:
        backend = CachedBound(z,alpha,(row,row),precision=precision,progress=progress)
        receipt = backend.receipt(ranges)
        receipt['witness'] = backend.witness()
        ledger.accept(receipt)
    return ledger


def refine(ledger,backends,*,threshold=-80,max_nodes=4096,minimum_r=1,progress=None):
    """Bound largest pending boxes first, retaining all unfinished work.

    All accepted bounds are recomputed by a directed backend. A node
    limit controls work only; it never turns pending boxes into coverage.
    """
    if not backends or max_nodes<1 or not 1<=minimum_r<=2*ledger.pairs:
        raise ValueError('Invalid refinement request')
    # Cut the scheduling boundary exactly; the lower part stays unresolved.
    for box in list(ledger.unresolved):
        lo,hi = box.ranges[0]
        if lo<minimum_r<=hi:
            low = Box(((lo,minimum_r-1),)+box.ranges[1:])
            high = Box(((minimum_r,hi),)+box.ranges[1:])
            index = ledger.unresolved.index(box)
            ledger.unresolved[index:index+1] = [low,high]
    queue,counter = [],0
    def enqueue(box):
        nonlocal counter
        if box.ranges[0][0]>=minimum_r:
            count = profile_count(box.ranges,ledger.pairs)
            heapq.heappush(queue,(-count,counter,box))
            counter += 1
    for box in ledger.unresolved:
        enqueue(box)
    visited,evaluations,accepted,empty = 0,0,0,0
    stalled = []
    while queue and visited<max_nodes:
        negative_count,_,box = heapq.heappop(queue)
        visited += 1
        if not negative_count:
            ledger.discard_empty(box)
            empty += 1
            continue
        best = None
        for index,backend in enumerate(backends):
            receipt = backend.receipt(box.ranges)
            evaluations += 1
            value = F(receipt['log2_upper_rational'])
            if best is None or value<best[0]:
                best = value,index,receipt
            if value<=threshold:
                receipt['witness_index'] = index
                ledger.accept(receipt)
                accepted += 1
                break
        else:
            if box.volume==1:
                stalled.append(dict(ranges=box.ranges,log2_upper_rational=str(best[0])))
            else:
                children = box.split()
                ledger.split(box)
                for child in children:
                    enqueue(child)
        if progress and visited%128==0:
            progress(dict(stage='directed-cover-refinement',visited=visited,evaluations=evaluations,
                          accepted=accepted,queued=len(queue),stalled=len(stalled)))
    return dict(visited=visited,evaluations=evaluations,accepted=accepted,empty=empty,
                queued=len(queue),stalled=stalled,minimum_r=minimum_r,
                witnesses=[backend.witness() for backend in backends],summary=ledger.summary())


def replay_plan(plan,*,precision=256,threshold=-80,progress=None,
                backend_type=CachedBound,ledger_type=Ledger):
    """Treat a FLOAT plan as proposals; rebuild every accepted bound."""
    if plan.get('profile_coordinates')!=['r','a','b'] or plan.get('all_even_odd_splits') is not True:
        raise ValueError('Expected an unsigned support proposal plan')
    if (backend_type.construction!=ledger_type.construction
            or plan.get('construction',ledger_type.construction)!=ledger_type.construction):
        raise ValueError('Proposal, backend, and ledger constructions must agree')
    rlo,rhi = plan['requested_weight_range']
    ledger = ledger_type(Box(((rlo,rhi),(0,4096),(0,4096))))
    grouped = {}
    for item in plan['accepted_boxes']:
        grouped.setdefault(item['witness_id'],[]).append(item)
    replays,failed = [],[]
    # Grouping equal z/alpha makes the unmarked coefficient cache useful.
    order = sorted(grouped,key=lambda i:(plan['witnesses'][i]['mark'],plan['witnesses'][i]['alpha'],i))
    for completed,index in enumerate(order,1):
        proposal = plan['witnesses'][index]
        markers = tuple(tuple(F(float(2**v)).limit_denominator(10**9) for v in row)
                        for row in proposal['marker_log_odds'])
        grid = proposal.get('grid')
        method = proposal.get('coefficient_method','positive' if grid is None else 'positive-alias')
        if method!=('positive' if grid is None else 'positive-alias'):
            raise ValueError('Inconsistent proposed coefficient method and grid')
        backend = backend_type(proposal['mark'],proposal['alpha'],markers,precision=precision,grid=grid,
                              uniform_aliases=proposal.get('uniform_aliases',False),
                              run_markers=proposal.get('run_markers'))
        for item in grouped[index]:
            ranges = tuple(tuple(row) for row in item['ranges'])
            receipt = backend.receipt(ranges)
            receipt['witness_index'] = index
            if F(receipt['log2_upper_rational'])<=threshold:
                ledger.accept(receipt)
            else:
                failed.append(receipt)
        replays.append(dict(index=index,**backend.witness()))
        if progress and (completed%10==0 or completed==len(order)):
            progress(dict(stage='directed-plan-replay',witnesses=completed,total_witnesses=len(order),
                          accepted=len(ledger.accepted),failed=len(failed)))
    bands = [ledger.band_summary(low,min(low+1023,rhi),precision)
             for low in range(rlo,rhi+1,1024)]
    return ledger,dict(summary=ledger.summary(precision),bands=bands,failed=failed,witnesses=replays,
                       all_proposed_leaves_replayed=len(ledger.accepted)+len(failed)==len(plan['accepted_boxes']))


def replay_cover(plan_bands,*,precision=256,threshold=-80,progress=None):
    """Join freshly replayed plans on disjoint bands with the dense tail.

    Each entry is (label, plan, low, high). There is no trusted-JSON receipt
    import: every bound is recomputed by replay_plan before composition.
    """
    plan_bands = tuple(plan_bands)
    bands = [Box(((low,high),(0,4096),(0,4096))) for _,_,low,high in plan_bands]
    tail = Box(((7649,8192),(0,4096),(0,4096)))
    for index,band in enumerate(bands):
        validate_box(band.ranges)
        if band.intersection(tail) or any(band.intersection(prior) for prior in bands[:index]):
            raise ValueError('Assigned plan bands and the terminal range must be disjoint')
    ledger = Ledger()
    ledger.accept(terminal_receipt(precision))
    sources = []
    for label,plan,low,high in plan_bands:
        source,result = replay_plan(plan,precision=precision,threshold=threshold,progress=progress)
        imported = ledger.import_band(source,low,high,label=label)
        sources.append(dict(label=label,replay=result,imported_band=imported))
        if progress:
            progress(dict(stage='directed-cover-import',label=label,**imported))
    summary = ledger.summary(precision)
    top = summary['complete_top_weight_range']
    return ledger,dict(summary=summary,sources=sources,
                       complete_top_band=ledger.band_summary(*top,precision) if top else None)


if __name__=='__main__':
    import argparse,json
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--precision',type=int,default=256)
    parser.add_argument('--details',action='store_true')
    selection = parser.add_mutually_exclusive_group()
    selection.add_argument('--plan',help='Replay accepted leaves of a saved FLOAT proposal plan')
    selection.add_argument('--band-plan',nargs=3,action='append',metavar=('PATH','LOW','HIGH'),
                           help='Replay each plan and import its assigned disjoint weight band, plus the dense tail')
    args = parser.parse_args()
    progress = lambda item:print(json.dumps(item),flush=True)
    if args.band_plan:
        from hashlib import sha256
        inputs,hashes = [],{}
        for path,low,high in args.band_plan:
            with open(path,'rb') as source:
                plan_bytes = source.read()
            inputs.append((path,json.loads(plan_bytes),int(low),int(high)))
            hashes[path] = sha256(plan_bytes).hexdigest()
        ledger,result = replay_cover(inputs,precision=args.precision,progress=progress)
        result['proposal_sha256'] = hashes
        if not args.details:
            for source in result['sources']:
                source['replay']['witness_count'] = len(source['replay'].pop('witnesses'))
    elif args.plan:
        from hashlib import sha256
        with open(args.plan,'rb') as source:
            plan_bytes = source.read()
        plan = json.loads(plan_bytes)
        ledger,result = replay_plan(plan,precision=args.precision,progress=progress)
        result['proposal_sha256'] = sha256(plan_bytes).hexdigest()
        if not args.details:
            result['witness_count'] = len(result.pop('witnesses'))
    else:
        ledger = seeded_ledger(args.precision,progress)
        result = ledger.summary(args.precision)
    if args.details:
        result['accepted'] = ledger.accepted
        result['unresolved'] = [box.ranges for box in ledger.unresolved]
    print(json.dumps(result),flush=True)
