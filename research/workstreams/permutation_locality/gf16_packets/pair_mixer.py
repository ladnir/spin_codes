"""Exact support diagnostics for a proposed local GF16 pair mixer.

This is a NEW ensemble, not a certificate for the existing shared4 code.
Each four-row group uses a uniform matching of its 256 packet positions.
On each matched pair, sample a,b,e independently from GF16* and output
(a*x+b*y, e*(a*x+alpha*b*y)), for fixed alpha not in {0,1}. A final
independent uniform column shuffle restores exchangeability. Different
pairs and groups use independent setup randomness.

The extra postmultiplier e is necessary for this proof interface: conditional
on output support, all nonzero packet labels are independent uniform GF16*.
Without e, this script's support law still holds but the inner proof cannot
be reused through its iid-label interface.

All kernel and CDF transport calculations are rational. Optional dense
point probes are outward diagnostics, never a whole-code certificate.
"""
import argparse
import json
from collections import Counter
from fractions import Fraction as Q
from math import comb
from pathlib import Path


def multiply(a, b):
    """Polynomial-basis GF16 multiplication modulo x^4+x+1."""
    out = 0
    for _ in range(4):
        if b & 1:
            out ^= a
        b >>= 1
        a = (a << 1) ^ (0x13 if a & 8 else 0)
    return out


def pair_law(x, y, alpha=2, postmultiply=True):
    if any(type(v) is not int or not 0 <= v < 16 for v in (x, y, alpha)) or alpha < 2:
        raise ValueError('GF16 input values and alpha outside {0,1} required')
    result = Counter()
    scalars = range(1, 16)
    for a in scalars:
        for b in scalars:
            u, v = multiply(a, x), multiply(b, y)
            for e in scalars if postmultiply else (1,):
                result[u ^ v, multiply(e, u ^ multiply(alpha, v))] += 1
    total = sum(result.values())
    return {pair: Q(count, total) for pair, count in result.items()}


def support_kernel(n=256):
    """P[output support=v | input support=u] for one uniformly matched layer."""
    if type(n) is not int or n < 2 or n % 2:
        raise ValueError('positive even packet count required')
    m = n // 2
    result = []
    powers2 = [2**j for j in range(m+1)]
    powers13 = [13**j for j in range(m+1)]
    powers15 = [15**j for j in range(m+1)]
    for u in range(n+1):
        # A uniformly shuffled support occupies j full pairs and u-2j
        # singleton pairs. In each full pair one output vanishes with
        # probability 2/15, independently between pairs.
        denominator = comb(n, u) * powers15[u//2]
        numerators = [0] * (n+1)
        for j in range(max(0, u-m), u//2+1):
            single = u-2*j
            matching = comb(m, j)*comb(m-j, single)*powers2[single]
            matching *= powers15[u//2-j]
            for loss in range(j+1):
                v = 2*u-2*j-loss
                numerators[v] += matching*comb(j, loss)*powers2[loss]*powers13[j-loss]
        assert sum(numerators) == denominator
        result.append([Q(x, denominator) for x in numerators])
    return result


def mean_support(u, n=256):
    return 2*u-Q(16*u*(u-1), 15*(n-1))


def gf256_support_kernel(n=256):
    """Alternative: uniform GF256* multiplication on each matched pair.

    Every nonzero eight-bit pair becomes uniform over 255 nonzero pairs.
    The resulting two GF16 labels are independent uniform conditional on
    their support. Unlike the MDS mixer, its support kernel is monotone:
    under nested input supports, each newly occupied pair adds an
    independent positive output-support contribution.
    """
    if type(n) is not int or n < 2 or n % 2:
        raise ValueError('positive even packet count required')
    m = n//2
    result = []
    for u in range(n+1):
        denominator = comb(n,u)*17**min(u,m)
        numerators = [0]*(n+1)
        for j in range(max(0,u-m),u//2+1):
            single, occupied = u-2*j, u-j
            matching = comb(m,j)*comb(m-j,single)*2**single
            matching *= 17**(min(u,m)-occupied)
            for loss in range(occupied+1):
                numerators[2*occupied-loss] += matching*comb(occupied,loss)*2**loss*15**(occupied-loss)
        assert sum(numerators) == denominator
        result.append([Q(v,denominator) for v in numerators])
    return result


def transport_cdf(cdf, kernel):
    """Upper expected output CDF from an upper input CDF.

    Let B be an upper CDF with the correct total mass. Its increments form
    a measure stochastically BELOW the unknown input measure. For each
    output threshold v, replace K[u, <=v] by its suffix maximum in u.
    This is a decreasing function, so integrating against increments of
    B safely bounds the unknown measure. No stochastic monotonicity of
    the actual pair kernel is assumed (it is false near full support).
    """
    n = len(cdf)-1
    if (n < 2 or len(kernel) != n+1 or cdf[0] != 0
            or any(Q(v) < 0 for v in cdf)
            or any(a > b for a, b in zip(cdf, cdf[1:]))
            or any(len(row) != n+1 or sum(row) != 1 or min(row) < 0 for row in kernel)):
        raise ValueError('monotone nonnegative zero-free CDF and matching stochastic kernel required')
    increments = [Q(cdf[0])] + [Q(b)-Q(a) for a, b in zip(cdf, cdf[1:])]
    kernel_cdf = [Q(0)]*(n+1)
    output = []
    for v in range(n+1):
        for u in range(n+1):
            kernel_cdf[u] += kernel[u][v]
        suffix = Q(0)
        upper = Q(0)
        for u in range(n, -1, -1):
            suffix = max(suffix, kernel_cdf[u])
            upper += increments[u]*suffix
        output.append(upper)
    assert output[0] == 0 and output[-1] == cdf[-1]
    assert all(a <= b for a, b in zip(output, output[1:]))
    return output


def ceil_cdf(cdf):
    return [(v.numerator+v.denominator-1)//v.denominator for v in map(Q, cdf)]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--layers', type=int, choices=(1, 2), default=1)
    parser.add_argument('--kernel', choices=('mds', 'gf256'), default='mds')
    parser.add_argument('--updates', type=int, choices=(1,2,3,4), default=2)
    parser.add_argument('--ceil-cdf', action='store_true',
        help='Conservative integer rounding used by the first retained point screen')
    parser.add_argument('--means', nargs='*', default=[])
    parser.add_argument('--distance', default='.1')
    parser.add_argument('--precision', type=int, default=256)
    parser.add_argument('--count-witnesses', nargs='*', type=Path, default=[])
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    if args.output and args.output.exists():
        parser.error('new diagnostic path required; existing results are preserved')
    if args.precision < 128 or not 0 < Q(args.distance) < Q(1, 2):
        parser.error('precision >=128 and distance in (0,1/2) required')
    if any(not 0 <= Q(x) <= 1 for x in args.means):
        parser.error('mean fractions in [0,1] required')
    from math import log2
    import shared_support
    import shared_mixture
    from flint import arb, ctx
    kernel = (support_kernel if args.kernel == 'mds' else gf256_support_kernel)()
    caps, _ = shared_support.shared_counts(True, True, count_witnesses=args.count_witnesses)
    original = caps
    for _ in range(args.layers):
        caps = transport_cdf(caps, kernel)
    # These are EXPECTED counts over the fresh mixer setup. Preserve tiny
    # fractional tails; rounding them to one would invent low-support mass.
    if args.ceil_cdf:
        caps = ceil_cdf(caps)
    print('PAIR MIXER support means', args.kernel, [(u, float(sum(v*p for v,p in enumerate(kernel[u])))) for u in (38, 76, 114, 120, 128, 192, 256)], flush=True)
    print('PAIR MIXER CDF log2 caps', [(u, log2(original[u]), log2(caps[u])) for u in (96, 114, 128, 144, 160, 176, 192, 208, 224)], flush=True)
    centers = sorted({Q(u, 256) for u in range(1, 257, 4)} | {Q(1)})
    mixture = shared_mixture.envelope(caps, centers, zero_bits=64, cost_tilt=Q(1,4))
    from positive_prune import prune
    mixture, pruning = prune(caps, mixture, Q(1,4))
    shared_mixture.verify(caps, mixture)
    components = shared_mixture.as_components(mixture)
    print('EXACT PAIR MIXER ENVELOPE', len(mixture), 'components', flush=True)
    record = dict(schema='shared-gf16-pair-mixer-diagnostic-1', layers=args.layers,
        kernel=args.kernel, updates=args.updates,
        cdf_rounding=('ceil' if args.ceil_cdf else 'none'),
        postmultiplier=(args.kernel == 'mds'), distance=args.distance, precision=args.precision,
        caps=list(map(str, caps)), original_caps=list(map(str, original)),
        mixture=[dict(mass=str(c), activity=str(p)) for c,p in mixture], pruning=pruning,
        note='New ensemble. Exact support transport and shell majorant; any point bounds cover only those points. Not a distance certificate.', probes=[])
    def save():
        if args.output:
            args.output.parent.mkdir(parents=True, exist_ok=True)
            args.output.write_text(json.dumps(record, indent=2)+'\n')
    save()
    if args.means:
        import birth_classes
        import scalar_cover as sc
        data = birth_classes.actual(args.updates)
        ctx.prec = args.precision
        model = sc.Model(components, data, int(Q(args.distance)*sc.N), 33, Q(3,16),
            inner=birth_classes, variance_shuffle=True, variance_bins=16, regional_count=True)
        model.proposal_stop_bits = 42
        for mean in map(Q, args.means):
            cell = (mean, mean)
            if model.empty(cell):
                record['probes'].append(dict(mean=str(mean), empty=True))
                continue
            score, witness = model.proposal(cell)
            ctx.prec = args.precision
            upper = model.outward(cell, witness)
            row = dict(mean=str(mean), proposal=score, witness=witness,
                upper=[int(x) for x in upper.upper().man_exp()],
                log2_upper=str(upper.log()/arb(2).log()))
            record['probes'].append(row)
            print('PAIR MIXER POINT', args.layers, args.distance, mean, row['log2_upper'], flush=True)
            save()
    save()


if __name__ == '__main__':
    main()
