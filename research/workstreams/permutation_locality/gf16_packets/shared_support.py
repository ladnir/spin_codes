"""Occupancy-one screen for shared row shuffles and independent GF16 labels.

All four BCH rows in a group use the same coordinate permutation. This is
NOT the independent-row ensemble. No full certificate is imported or claimed.
"""
import argparse
import json
from fractions import Fraction as Q
from pathlib import Path
from flint import arb, ctx
import single_group
import sparse_cover
import occupancy_rank
from bch_joint_support import authenticated_caps, support_caps
from basis_lattice import improve_caps
from shortened_bound import dimension_caps


def shared_counts(refined=False, coupled=False, joint=False, wide=False, count_witnesses=()):
    if (any(type(v) is not bool for v in (refined, coupled, joint, wide))
            or (coupled and not refined) or (joint and not coupled) or (wide and not joint)
            or (count_witnesses and not coupled)):
        raise ValueError('boolean count-refinement option required')
    spectrum = authenticated_caps()
    dimensions = dimension_caps()
    if refined:
        from shortening_polynomial import improve_dimensions
        from dual_shortening import improve_dimensions as dual_dimensions
        dimensions = dual_dimensions(improve_dimensions(dimensions))
    ranks = improve_caps(support_caps(spectrum, g=4, dimensions=dimensions), spectrum)
    if refined:
        from shortening_moments import improve
        from dual_moments import refine_bch
        ranks, _ = improve(ranks, dimensions)
        ranks = refine_bch(ranks, dimensions)
    if coupled:
        if joint or count_witnesses:
            from count_refinements import refine_pair
            from constraint_counts import refine, Constraints
            from dual_moments import dual_shell_caps
            from bch_joint_support import rank_total
            ranks, dual, dims, ddims = refine_pair(ranks, spectrum)
            checked = refine(ranks, dual, dims, ddims, spectrum) if joint else [row[:] for row in ranks]
            if wide:
                wider = refine(ranks, dual, dims, ddims, spectrum, last=192, supports=range(96, 137, 2))
                checked = [[min(a, b) for a, b in zip(left, right)] for left, right in zip(checked, wider)]
            for path in count_witnesses:
                record = json.loads(Path(path).read_text())
                symmetry = record['complement_symmetry']
                if record['schema'] != 'shared-gf16-constraint-counts-1' or type(symmetry) is not bool:
                    raise ValueError('valid joint-count witness record required')
                system = Constraints(ranks, dual, dims, ddims, 128, record['last'],
                    spectra=(spectrum, dual_shell_caps()), ones=(symmetry, symmetry))
                accepted = 0
                for row in record['rows']:
                    if 'dual_witness' not in row:
                        continue
                    witness = row['dual_witness']
                    cap = system.verify_witness(witness)
                    h, u = witness['rank'], witness['support']
                    checked[h-1][u] = min(checked[h-1][u], cap*rank_total(h, 4, h))
                    accepted += 1
                if not accepted:
                    raise ValueError('count witness file has no exact multipliers to replay')
                print('EXACT COUNT WITNESSES replayed', accepted, 'from', path, flush=True)
            for row in checked:
                for u in range(255, -1, -1):
                    row[u] = min(row[u], row[u+1])
            ranks = checked
        else:
            from count_refinements import refine
            ranks = refine(ranks, spectrum)
    counts = [sum(row[u] for row in ranks) for u in range(257)]
    assert counts[0] == 0 and counts[-1] == (1 << 512)-1
    assert all(a <= b for a, b in zip(counts, counts[1:]))
    return counts, ranks


def run(updates, thresholds, tilts, precision=256, output=None):
    if (not updates or any(type(r) is not int or r not in (2,3,4) for r in updates)
            or not thresholds or any(type(t) is not int or not 0 <= t < 1 << 21 for t in thresholds)
            or not tilts or any(Q(t) <= 0 for t in tilts) or precision < 128):
        raise ValueError('updates 2..4, valid cutoffs, positive tilts, precision >=128 required')
    ctx.prec = precision
    counts, ranks = shared_counts()
    results = []
    for r in updates:
        args = sparse_cover.build_args(1, tilts, precision, 0, 40, None, r)
        operators = occupancy_rank.build_operators(args)
        moments = {}
        for tilt in tilts:
            exact = operators[tilt,'1'][0]
            moments[tilt] = single_group.support_moments(exact[0], exact[1])
        for threshold in thresholds:
            best = [arb(1)]*257
            choices = [None]*257
            for tilt, weights in moments.items():
                multiplier = (occupancy_rank.aq(Q(tilt))*threshold).exp()
                for u, value in enumerate(weights):
                    value = occupancy_rank.up(value*multiplier)
                    if value < best[u]:
                        best[u], choices[u] = value, tilt
            rank_uppers = [occupancy_rank.up(2048*single_group.fold_cdf(row,best)) for row in ranks]
            upper = occupancy_rank.up(2048*single_group.fold_cdf(counts,best))
            margin = -upper.log()/arb(2).log()
            print('SHARED GF16 ONE GROUP updates',r,'cutoff',threshold,'margin',margin,flush=True)
            for h, value in enumerate(rank_uppers, 1):
                print('  rank',h,'log2 upper',value.log()/arb(2).log(),flush=True)
            pair = lambda x: [int(v) for v in x.upper().man_exp()]
            results.append(dict(updates=r, threshold=threshold, upper=pair(upper),
                                rank_uppers=list(map(pair,rank_uppers)),
                                support_probability_uppers=list(map(pair,best)), choices=choices))
            if output:
                record = dict(schema='shared-gf16-single-group-screen-1', precision=precision,
                              tilts=tilts, results=results,
                              note='Shared coordinate shuffle per four rows. Occupancy one only, all supports. Not a full-code certificate.')
                output.parent.mkdir(parents=True,exist_ok=True)
                output.write_text(json.dumps(record,indent=2)+'\n')
    return results


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--updates',type=int,nargs='+',default=[2,4])
    parser.add_argument('--thresholds',type=int,nargs='+',default=[188743,193986,199229,209715])
    parser.add_argument('--tilts',nargs='+',default=['.00016','.00024','.00028','.00032','.0004','.00064','.001','.0032','.008','.016','.032','.064','.096'])
    parser.add_argument('--precision',type=int,default=256)
    parser.add_argument('--output',type=Path)
    args = parser.parse_args()
    run(args.updates,args.thresholds,args.tilts,args.precision,args.output)
