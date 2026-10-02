"""Replay retained joint-count witnesses, then propagate exact moments.

All input cap arrays are regenerated from the BCH premises. Saved numerical
LP objectives are ignored: only exact dual multipliers are replayed. Later
primal/dual and containment steps can propagate each accepted improvement.
"""
import argparse
import hashlib
import json
from fractions import Fraction as Q
from math import log2
from pathlib import Path

import shared_support
from count_refinements import refine


def actual(paths,iterations=8):
    if type(iterations) is not int or not 0<=iterations<=30:
        raise ValueError('zero through thirty exact propagation iterations required')
    paths=[Path(p) for p in paths]
    provenance=[dict(path=str(p),sha256=hashlib.sha256(p.read_bytes()).hexdigest()) for p in paths]
    counts,ranks=shared_support.shared_counts(True,True,count_witnesses=paths)
    stages=[dict(name='retained-exact-duals',shell_caps=counts,rank_cdfs=ranks)]
    if iterations:
        ranks=refine(ranks,shared_support.authenticated_caps(),iterations)
        counts=[sum(row[u] for row in ranks) for u in range(257)]
        stages.append(dict(name='propagated-exact-duals',shell_caps=counts,rank_cdfs=ranks))
    assert counts[0]==0 and counts[-1]==(1<<512)-1
    assert all(a<=b for a,b in zip(counts,counts[1:]))
    if any(hashlib.sha256(p.read_bytes()).hexdigest()!=row['sha256']
           for p,row in zip(paths,provenance)):
        raise RuntimeError('a count-witness file changed during exact replay')
    return counts,dict(witnesses=provenance,iterations=iterations,stages=stages)


def build_components(paths,*,iterations=0,step=8,zero_bits=64,cost_tilt=Q(1,4),
                     saved_mixture=None,starting_mixture=None):
    """Regenerate counts and verify a shared-row positive comparison.

    Setup/search uses the existing deterministic envelope and exact-repaired
    pruning. Replay should supply the saved rational mixture instead: checking
    every shell against freshly regenerated caps avoids reliance on numerical
    LP reproducibility. A changed file or insufficient majorant is rejected.

    A starting mixture can preserve a previously useful comparison geometry.
    Pruning only reduces its component masses. It is a search input, not an
    additional premise: both the start and result are checked against fresh
    counts. Saved-mixture replay does not need the original search input.

    Return (components, caps, mixture, metadata). The metadata authenticates
    the caps and witness files, but supplies no inner or distance certificate.
    """
    import shared_mixture
    import positive_prune
    if type(step) is not int or not 1<=step<=64 or not 0<Q(cost_tilt)<=1:
        raise ValueError('valid support step and positive bounded cost tilt required')
    if zero_bits is not None and (type(zero_bits) is not int or not 0<=zero_bits<=512):
        raise ValueError('valid optional empty-component budget required')
    if saved_mixture is not None and starting_mixture is not None:
        raise ValueError('saved replay and a pruning start are mutually exclusive')
    caps,source=actual(paths,iterations)
    centers=sorted({Q(u,256) for u in range(38,257,step)}|{Q(1)})

    def decode(rows):
        result=[(Q(row['mass']),Q(row['activity'])) if isinstance(row,dict)
                else (Q(row[0]),Q(row[1])) for row in rows]
        if (not result or len({p for _,p in result})!=len(result)
                or any(p not in centers for _,p in result)):
            raise ValueError('distinct saved activities on the declared center grid required')
        shared_mixture.verify(caps,result)
        return result

    if saved_mixture is None:
        baseline=(shared_mixture.envelope(caps,centers,zero_bits,Q(cost_tilt))
                  if starting_mixture is None else decode(starting_mixture))
        mixture,pruning=positive_prune.prune(caps,baseline,Q(cost_tilt))
        old_masses={p:c for c,p in baseline}
        if any(p not in old_masses or c>old_masses[p] for c,p in mixture):
            raise ArithmeticError('componentwise pruning increased a comparison mass')
        mode=('searched-and-exact-repaired' if starting_mixture is None
              else 'saved-start-componentwise-pruned')
    else:
        mixture=decode(saved_mixture)
        pruning=None
        mode='saved-rationals-freshly-verified'
    shared_mixture.verify(caps,mixture)
    if zero_bits is not None and any(c*(1-p)**256>Q(2)**-zero_bits for c,p in mixture):
        raise ArithmeticError('a saved component exceeds the declared empty-input budget')
    cap_hash=hashlib.sha256(json.dumps([str(c) for c in caps]).encode()).hexdigest()
    metadata=dict(schema='shared-relaxed-authenticated-mixture-1',cap_sha256=cap_hash,
        count_witnesses=source['witnesses'],count_refinement_iterations=iterations,
        step=step,zero_bits=zero_bits,cost_tilt=str(Q(cost_tilt)),pruned=True,
        mixture_mode=mode,pruning=pruning,
        scope='shared four-row coordinate shuffle and independent GF16* labels; outer comparison only')
    return shared_mixture.as_components(mixture),caps,mixture,metadata


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--count-witnesses',nargs='+',type=Path,required=True)
    parser.add_argument('--iterations',type=int,default=8)
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args()
    if args.output.exists():parser.error('refusing to overwrite a receipt')
    caps,provenance=actual(args.count_witnesses,args.iterations)
    for stage in provenance['stages']:
        print('EXACT COUNT STAGE',stage['name'],[(u,log2(stage['shell_caps'][u]))
             for u in (96,112,114,115,120,128,144,160,192,209,210)],flush=True)
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(dict(schema='shared-relaxed-propagated-counts-1',
        shell_caps=caps,provenance=provenance,proof_status='exact shared-group CDF caps; not a distance certificate'),indent=2)+'\n')


if __name__=='__main__':main()
