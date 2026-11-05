"""Explicit coverage and inheritance for the independently swapped route.

Fixed-lane receipts are never silently relabeled. Only unsigned-envelope
receipts and the route-independent terminal bound use the proved transfer.
New swap proposals are replayed by the swap backend before acceptance.
"""
from copy import deepcopy
from fractions import Fraction as F

from hd_shifted_cover_backend import CONSTRUCTION, validate_box
from hd_shifted_cover_ledger import Box, Ledger, replay_plan, terminal_receipt
from hd_shifted_endpoint import swap_endpoint_receipt
from hd_shifted_swap_certificate import SWAP_CONSTRUCTION, SwapCachedBound


UNSIGNED_METHODS = frozenset((
    'matching-profile-envelope', 'whole-box-cyclic-runs',
    'tilted-cyclic-run-sum', 'joint-profile-count',
    'projected-normalization-with-cyclic-run-tilt',
))


class SwapLedger(Ledger):
    construction = SWAP_CONSTRUCTION

    def accept_inherited(self, receipt):
        """Apply the unsigned domination lemma to a directed old receipt.

        The caller supplies a freshly recomputed receipt, not an untrusted
        saved JSON bound. This helper records the exact transfer used.
        """
        if (receipt.get('construction') != CONSTRUCTION
                or receipt.get('covers_all_even_odd_splits') is not True
                or receipt.get('directed') is not True):
            raise ValueError('A directed original-route unsigned receipt is required')
        terminal = receipt.get('method') == 'inverse-image-terminal-range'
        if not terminal and receipt.get('support_count_method') not in UNSIGNED_METHODS:
            raise ValueError('This receipt is not an established unsigned-envelope bound')
        transferred = deepcopy(receipt)
        transferred.update(
            construction=self.construction,
            inherited_from=CONSTRUCTION,
            inheritance=('route-independent-inverse-image-bound' if terminal
                         else 'independent-swap-unsigned-coefficient-domination'),
            whole_code_certificate=False,
        )
        self.accept(transferred)


def inherit_unsigned_ledger(source):
    """Transfer a recomputed old ledger, preserving gaps and whole leaf mass."""
    if source.construction != CONSTRUCTION:
        raise ValueError('Expected the original fixed-lane unsigned ledger')
    source.summary()
    target = SwapLedger(source.root, pairs=source.pairs)
    for receipt in source.accepted:
        target.accept_inherited(receipt)
    for box in source.empty:
        target.discard_empty(box)
    target.summary()
    return target


def endpoint_receipt(precision=256):
    """The two alternating supports comprise the whole unsigned (4096,0,0) class."""
    receipt = swap_endpoint_receipt(precision=precision)
    receipt.update(
        construction=SWAP_CONSTRUCTION,
        ranges=((4096,4096),(0,0),(0,0)),
        integer_profile_triples=1,
        covers_all_even_odd_splits=True,
        method='exact-alternating-support-swap-average',
    )
    return receipt


def seed_ledger(precision=256):
    """Create only the proven endpoint and route-independent terminal range."""
    ledger = SwapLedger()
    ledger.accept_inherited(terminal_receipt(precision))
    ledger.accept(endpoint_receipt(precision))
    return ledger


def replay_swap_plan(plan, *, precision=256, threshold=-80, progress=None):
    """Replay a FLOAT swap proposal with explicit ensemble identity."""
    if plan.get('construction') != SWAP_CONSTRUCTION:
        raise ValueError('An explicit swapped-route proposal identity is required')
    for witness in plan.get('witnesses', []):
        if (witness.get('pair_swaps') is not True
                or witness.get('orientation_aggregation') != 'independent-fair-swaps'):
            raise ValueError('A proposal witness must select independent fair swaps')
    return replay_plan(plan, precision=precision, threshold=threshold, progress=progress,
                       backend_type=SwapCachedBound, ledger_type=SwapLedger)


def replay_sparse_family(last=64, *, precision=256, progress=None):
    """Recompute the whole sparse band, including every matching profile."""
    from hd_shifted_swap_sparse import family_receipt
    receipt = family_receipt(last,precision=precision,progress=progress)
    ledger = SwapLedger(Box(((1,last),(0,4096),(0,4096))))
    ledger.accept(receipt)
    return ledger,dict(summary=ledger.summary(precision),receipt=receipt)


def replay_sparse_batch(first,last,z,alpha,maximum, *, precision=256,progress=None):
    """Rebuild a stated band using one rational coefficient witness."""
    from hd_shifted_swap_sparse import batch_receipt
    receipt = batch_receipt(first,last,F(z),F(alpha),maximum=maximum,
                            precision=precision,progress=progress)
    ledger = SwapLedger(Box(((first,last),(0,4096),(0,4096))))
    ledger.accept(receipt)
    return ledger,dict(summary=ledger.summary(precision),receipt=receipt)


def replay_swap_cover(unsigned_bands=(), swap_bands=(), *, sparse_last=0,
                      sparse_batches=(), potential_bands=(), potential_max_weights=33,projected=False,
                      precision=256, threshold=-80, progress=None):
    """Join only freshly recomputed bounds on disjoint assigned domains.

    Plans contain proposals, not trusted receipts. Sparse coefficients and
    projected checkpoints are recomputed too. All uncovered profiles stay
    in the returned ledger; accepted probability mass is never prorated.
    """
    unsigned_bands,swap_bands = tuple(unsigned_bands),tuple(swap_bands)
    sparse_batches = tuple(sparse_batches)
    potential_bands = tuple(potential_bands)
    if not isinstance(potential_max_weights,int) or not 1<=potential_max_weights<=8192:
        raise ValueError('Potential replay requires an explicit positive weight budget')
    if not isinstance(sparse_last,int) or not 0<=sparse_last<=64:
        raise ValueError('Sparse range must be omitted or end between 1 and 64')
    assigned = [Box(((7649,8192),(0,4096),(0,4096)))]
    if sparse_last:
        assigned.append(Box(((1,sparse_last),(0,4096),(0,4096))))
    for first,last,z,alpha,maximum in sparse_batches:
        if (not all(isinstance(n,int) for n in (first,last,maximum))
                or not 1<=first<=last<=maximum<=256 or not 0<F(z)<1 or not 0<F(alpha)<=1):
            raise ValueError('Invalid explicit sparse batch witness')
        assigned.append(Box(((first,last),(0,4096),(0,4096))))
    if projected:
        assigned.extend((Box(((500,524),(188,212),(188,212))),
                         Box(((1012,1036),(388,412),(388,412)))))
    for _,_,low,high in unsigned_bands+swap_bands:
        assigned.append(Box(((low,high),(0,4096),(0,4096))))
    for _,proposal in potential_bands:
        low,high = proposal['requested_weight_range']
        if high-low+1>potential_max_weights:
            raise ValueError('Assigned potential band exceeds the replay weight budget')
        assigned.append(Box(((low,high),(0,4096),(0,4096))))
    for i,box in enumerate(assigned):
        validate_box(box.ranges)
        if any(box.intersection(prior) for prior in assigned[:i]):
            raise ValueError('Assigned proof domains must be disjoint')
    ledger,sources = SwapLedger(),[]
    ledger.accept_inherited(terminal_receipt(precision))
    if sparse_last:
        source,result = replay_sparse_family(sparse_last,precision=precision,progress=progress)
        imported = ledger.import_band(source,1,sparse_last,label='exact-sparse-coefficients')
        sources.append(dict(label='exact-sparse-coefficients',summary=result['summary'],imported=imported))
    for first,last,z,alpha,maximum in sparse_batches:
        source,result = replay_sparse_batch(first,last,z,alpha,maximum,precision=precision,progress=progress)
        label = f'exact-sparse-batch-{first}-{last}'
        imported = ledger.import_band(source,first,last,label=label)
        sources.append(dict(label=label,summary=result['summary'],imported=imported,
                            witness=result['receipt']['witness']))
    if projected:
        from hd_shifted_projected_bound import checkpoint_receipts
        for receipt in checkpoint_receipts(precision):
            if F(receipt['log2_upper_rational'])>threshold:
                raise ArithmeticError('Fixed projected checkpoint no longer passes its threshold')
            ledger.accept(receipt)
    for label,proposal in potential_bands:
        from hd_shifted_potential_band import replay_saved_band
        receipt = replay_saved_band(proposal,precision=precision,progress=progress,
                                     max_weights=potential_max_weights)
        if F(receipt['log2_upper_rational'])>threshold:
            raise ArithmeticError('Replayed potential band misses its acceptance threshold')
        ledger.accept(receipt)
        sources.append(dict(label=label,receipt=receipt))
    for inherited,inputs in ((True,unsigned_bands),(False,swap_bands)):
        for label,plan,low,high in inputs:
            if inherited:
                source,result = replay_plan(plan,precision=precision,threshold=threshold,progress=progress)
                source = inherit_unsigned_ledger(source)
            else:
                source,result = replay_swap_plan(plan,precision=precision,threshold=threshold,progress=progress)
            imported = ledger.import_band(source,low,high,label=label)
            sources.append(dict(label=label,replay=result,imported=imported))
            if progress:
                progress(dict(stage='directed-swap-cover-import',label=label,**imported))
    endpoint = Box(((4096,4096),(0,0),(0,0)))
    if any(box.intersection(endpoint)==endpoint for box in ledger.unresolved):
        ledger.accept(endpoint_receipt(precision))
    return ledger,dict(summary=ledger.summary(precision),sources=sources)


if __name__ == '__main__':
    import argparse
    import json
    from hashlib import sha256

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--plan', help='Recompute one swap FLOAT proposal plan')
    parser.add_argument('--unsigned-band',nargs=3,action='append',metavar=('PATH','LOW','HIGH'))
    parser.add_argument('--swap-band',nargs=3,action='append',metavar=('PATH','LOW','HIGH'))
    parser.add_argument('--sparse',type=int,default=0,metavar='LAST')
    parser.add_argument('--sparse-batch',nargs=5,action='append',
                        metavar=('FIRST','LAST','Z','ALPHA','MAXIMUM'),
                        help='Fresh directed replay of one exact-coefficient witness on a disjoint band')
    parser.add_argument('--potential-band',action='append',metavar='PATH',
                        help='Replay every weight from one exact fixed-atlas proposal')
    parser.add_argument('--potential-max-weights',type=int,default=33,
                        help='Explicit per-band replay resource budget; proposals cannot raise it')
    parser.add_argument('--projected',action='store_true')
    parser.add_argument('--precision', type=int, default=256)
    args = parser.parse_args()
    combined = bool(args.unsigned_band or args.swap_band or args.sparse or args.sparse_batch
                    or args.potential_band or args.projected)
    if args.plan and combined:
        parser.error('--plan cannot be combined with cover options')
    if combined:
        hashes = {}
        def load_bands(items):
            inputs = []
            for path,low,high in items or []:
                with open(path,'rb') as source:
                    data = source.read()
                hashes[path] = sha256(data).hexdigest()
                inputs.append((path,json.loads(data),int(low),int(high)))
            return inputs
        unsigned,swapped = load_bands(args.unsigned_band),load_bands(args.swap_band)
        batches = [(int(first),int(last),F(z),F(alpha),int(maximum))
                   for first,last,z,alpha,maximum in args.sparse_batch or []]
        potential_bands = []
        for path in args.potential_band or []:
            with open(path,'rb') as source:
                data = source.read()
            potential_bands.append((path,json.loads(data)))
            hashes[path] = sha256(data).hexdigest()
        ledger,result = replay_swap_cover(unsigned,swapped,sparse_last=args.sparse,
                                          sparse_batches=batches,potential_bands=potential_bands,
                                          potential_max_weights=args.potential_max_weights,
                                          projected=args.projected,precision=args.precision,
                                          progress=lambda item:print(json.dumps(item),flush=True))
        result['proposal_sha256'] = hashes
        for source in result['sources']:
            if 'replay' in source:
                source['replay']['witness_count'] = len(source['replay'].pop('witnesses'))
    elif args.plan:
        with open(args.plan, 'rb') as source:
            data = source.read()
        ledger, result = replay_swap_plan(json.loads(data), precision=args.precision,
                                         progress=lambda item: print(json.dumps(item), flush=True))
        result['proposal_sha256'] = sha256(data).hexdigest()
    else:
        ledger = seed_ledger(args.precision)
        result = dict(summary=ledger.summary(args.precision), accepted=ledger.accepted)
    print(json.dumps(result), flush=True)
