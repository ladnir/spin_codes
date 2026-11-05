"""Directed support-neighborhood receipts, not a whole-code cover.

The exact support-count envelope, scalar chord envelope, and interval FFT
coefficient envelopes each cover every feasible profile in the given box.
"""
from fractions import Fraction as F

from flint import arb,ctx

from hd_shifted_general_certificate import general_coefficients
from hd_shifted_general_count import ProfileBox,box_count_upper
from hd_shifted_general_sweep import CHECKPOINTS
from hd_shifted_good_pair_certificate import good_coefficients,regional_numerator
from hd_shifted_box_normalizer import scalar_envelope,histogram_ranges
from hd_sumrank_paired_certificate import scalar_power


DEFAULT_COUNT_MARKERS = (F('78107/24499979'),F('59183/54662624'),
                         F('79848470422/92864511'),F('17740896315/66165169'))


def receipt(box,z,alpha,markers,count_markers=DEFAULT_COUNT_MARKERS,*,grid=64,
            precision=256,bits=60,threshold=-40,progress=None):
    if precision<128:
        raise ValueError('At least128 bits are required')
    z,alpha = F(z),F(alpha)
    if not 0<z<1 or not 0<alpha<=1:
        raise ValueError('Invalid mark or fractional exponent')
    markers = tuple(tuple(map(F,row)) for row in markers)
    if len(markers)!=2 or any(len(row)!=3 or min(row)<=0 for row in markers):
        raise ValueError('Two positive marker triples are required')
    # A zero-message profile must never be accidentally counted as evidence
    # for a positive-distance claim. Split it away before invoking this API.
    if box.even[0]+box.odd[0]==0:
        raise ValueError('Box must exclude the zero message')
    count = box_count_upper(4096,box,count_markers)
    if count<=0:
        raise ValueError('Empty support box')
    histograms = histogram_ranges(box)
    with ctx.workprec(precision):
        if progress:
            progress(dict(stage='unmarked-single-coefficients'))
        single = general_coefficients(z,alpha,(F(1),)*3,bits=bits,progress=progress)
        if progress:
            progress(dict(stage='ordered-good-good-coefficients'))
        paired = good_coefficients(z,alpha,bits)
        regions,diagnostics = [],[]
        for parity,(ranges,row) in enumerate(zip(histograms,markers)):
            if progress:
                progress(dict(stage='regional-neighborhood',parity=parity))
            target = tuple(hi for lo,hi in ranges[1:])
            region,diagnostic = regional_numerator(single,paired,row,target,grid,progress,
                                                    target_ranges=ranges[1:])
            regions.append(region)
            diagnostics.append(diagnostic)
        first,second = regions
        paired_region = tuple(first[2*r]*second[s]+first[2*r+1]*second[2+s]
                              for r in range(2) for s in range(2))
        product = scalar_power(paired_region,128)
        scalar = scalar_envelope(4096,box,z,alpha,markers,precision=precision)
        log2 = lambda value:value.log()/arb(2).log()
        bound = log2(arb(count))+scalar+log2(product[0]+product[1])
        ceiling = F(int((bound.upper()*10**12).upper().ceil().fmpz()),10**12)
        return dict(box=box.ranges,mark=str(z),alpha=str(alpha),
                    markers=[list(map(str,row)) for row in markers],
                    count_markers=list(map(str,count_markers)),
                    support_count_upper=str(count),grid=grid,precision=precision,
                    packet_mantissa_bits=bits,log2_upper_enclosure=str(bound),
                    log2_upper_rational=str(ceiling),accepted=bool(bound<threshold),
                    threshold=threshold,regional_diagnostics=diagnostics,
                    directed=True,whole_code_certificate=False)


if __name__=='__main__':
    import argparse,json
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--radius',type=int,default=2)
    parser.add_argument('--grid',type=int,default=64)
    parser.add_argument('--precision',type=int,default=256)
    args = parser.parse_args()
    if args.radius<0:
        raise ValueError('Radius must be nonnegative')
    profile,z,alpha,odds = CHECKPOINTS['asymmetric']
    markers = tuple(tuple(F(float(2**v)).limit_denominator(10**9) for v in row) for row in odds)
    box = ProfileBox(*((max(0,n-args.radius),min(4096,n+args.radius)) for n in profile))
    print(json.dumps(receipt(box,F(z),F(alpha),markers,grid=args.grid,precision=args.precision,
        progress=lambda item:print(json.dumps(item),flush=True))),flush=True)
