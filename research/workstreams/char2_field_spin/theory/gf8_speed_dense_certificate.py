"""Directed dense replay for the length64 speed profile.

Reuses the exact-MDS positive-excess/Maclaurin membership theorem, with
64 regional conditioning factors and the omitted-output completion union.
"""
from fractions import Fraction as F
from hashlib import sha256

from flint import ctx
from gf8_exact_spectrum_dense import interval_margin,propose_excess_tilt,support_counts
from gf8_outer_rate_dense import lower_integer
from gf8_speed_profile import metadata


def certify(dimension=40,k_log=20,first_density=F(4,5),cell_bits=128,
            initial_grid=200,precision=192,include_witnesses=False):
    if (not isinstance(first_density,F) or not 0<first_density<=1
            or cell_bits<40 or initial_grid<1 or precision<80):
        raise ValueError('invalid directed dense replay parameters')
    g=metadata(dimension,k_log)
    if not F(10,9)**100*2295**10<2**127:
        raise ArithmeticError('retained-output ball comparison failed')
    n=g['input_n']
    first=(n*first_density.numerator+first_density.denominator-1)//first_density.denominator
    previous=ctx.prec
    try:
        ctx.prec=precision
        span=(n+initial_grid-1)//initial_grid
        pending=[(lo,min(n,lo+span-1)) for lo in range(first,n+1,span)][::-1]
        entries=[]
        while pending:
            lo,hi=pending.pop()
            tilt=propose_excess_tilt(F(lo,n),F(hi,n),*g['outer'])
            bits=lower_integer(interval_margin(g,lo,hi,tilt))
            if bits>=cell_bits:
                entries.append((lo,hi,tilt,bits))
            elif lo==hi:
                raise ArithmeticError(f'dense weight{lo} has only{bits}bits')
            else:
                mid=(lo+hi)//2
                pending.extend(((mid+1,hi),(lo,mid)))
        if (entries[0][0]!=first or entries[-1][1]!=n
                or any(a[1]+1!=b[0] for a,b in zip(entries,entries[1:]))):
            raise ArithmeticError('dense interval coverage gap')
        bits=min(entry[3] for entry in entries)
        digest=sha256(''.join(f'{lo}:{hi}:{tilt}:{margin}\n' for lo,hi,tilt,margin in entries).encode('ascii')).hexdigest()
        spectrum=sha256(','.join(map(str,support_counts(*g['outer']))).encode('ascii')).hexdigest()
        result=dict(**g,weights=[first,n],covered_weights=n-first+1,first_density=str(first_density),
            checked_cells=len(entries),precision_bits=precision,initial_grid=initial_grid,target_cell_bits=cell_bits,
            ball_exponent='127/100',conditioning='(rows+1)^64',omitted_symbols=n-g['output_n'],
            completion_factor_log2=8*(n-g['output_n']),margin_bits_lower=bits,
            failure_upper_dyadic_exponent=bits,failure_upper_dyadic_hex=['0x1',hex(1<<bits)],
            witness_sha256=digest,exact_support_spectrum_sha256=spectrum,
            dense_method='exact-MDS-positive-excess-Maclaurin',
            dense_inner_hypothesis='invertible full causal inner; independent nonzero input diagonals',
            complete=True,rigorous=True,whole_code=False,
            scope='complete dense inverse-input interval including omitted-output completions')
        if include_witnesses:
            result['witnesses']=[[lo,hi,str(tilt),margin] for lo,hi,tilt,margin in entries]
        return result
    finally:
        ctx.prec=previous


if __name__=='__main__':
    import json
    result=certify()
    del result['failure_upper_dyadic_hex']
    print(json.dumps(result),flush=True)
