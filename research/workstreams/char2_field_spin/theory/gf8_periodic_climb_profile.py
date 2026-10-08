"""Identity of a public periodic retained-output pattern, not a certificate."""
from fractions import Fraction as F
from gf8_rate_climb_profile import metadata as fixed_metadata


def metadata(dimension=84,pattern=(6,6,6,7),k_log=20):
    pattern=tuple(pattern)
    if not pattern or any(not isinstance(x,int) or not 1<=x<=8 for x in pattern):
        raise ValueError('nonempty fixed packet output pattern required')
    g=fixed_metadata(dimension,max(pattern),k_log)
    if g['rows']%(8*len(pattern)):
        raise ValueError('period must divide every regional packet count')
    average=F(sum(pattern),len(pattern))
    n=g['input_n']//(8*len(pattern))*sum(pattern)
    if g['padded_k']>=n:
        raise ValueError('rate below one required')
    del g['outputs'],g['kept_outputs']
    pattern_name='_'.join(map(str,pattern))
    g.update(construction=f'gf256-full16-cauchy-w8-s8-keep-period-{pattern_name}-grs128-{dimension}-regional-nonzero-diagonals',
        packet_period=len(pattern),output_pattern=list(pattern),
        average_outputs=str(average),kept_outputs_per_phase=[list(range(x)) for x in pattern],
        pattern_phase_origin='packet zero; repeated in chronological order without resetting state',
        output_n=n,output_threshold=(n+9)//10-1,
        nominal_rate=str(F(g['padded_k'],n)),requested_rate=str(F(2**k_log,n)))
    return g
