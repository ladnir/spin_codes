"""Bounded pre-H Toeplitz2 comparison; no full proof or native changes.

The iid-majorant depends only on the input support.  Artificial independent
input diagonals may therefore be inserted before it, retaining the known
RS counting bound although the proposed native map has no such diagonals.
The exact-support mode is ONLY a uniform-input diagnostic: it does not
establish that measure bridge for actual RS words.
"""
from fractions import Fraction as F
from functools import lru_cache
from math import comb, log2

from hd_two_sum_input_weight import exact_support_count
from hd_two_sum_toeplitz_pairs import checked_matching, checked_q, pair_kernel, xor_matching, column_contraction


def pre_pair_support_masses(active, *, q=256, input_marker=F(1), mode='iid-majorant'):
    """Unnormalized support masses after the comparison on uniform F^A.

    The weight v^wt(input) is attached before the pair map. Conditional on
    each output support, amplitudes are iid-NZ under this comparison.
    Exact mode has that property only after averaging uniform-F^A input.
    """
    checked_q(q)
    if type(active) is not int or not 0 <= active <= 3:
        raise ValueError('A two-coordinate input mask is required')
    marker = F(input_marker)
    if marker < 1 or mode not in ('exact-support', 'iid-majorant', 'diagonal', 'uniform-cap'):
        raise ValueError('A lower-weight marker>=1 and declared pair law are required')
    kernel = pair_kernel(q, mode)
    return tuple(sum(((F(q-1)*marker)**incoming.bit_count()*kernel[incoming][outgoing]
                      for incoming in range(4) if incoming & ~active == 0), F(0))
                 for outgoing in range(4))


def expanded_support_polynomial(active, matching, *, group=16, input_marker=F(1), mode='iid-majorant'):
    """Positive degree-b mass before H; symmetric scalar output marks only."""
    matching = checked_matching(matching, group)
    if group not in (8, 16) or type(active) is not int or not 0 <= active < 1 << group:
        raise ValueError('A fixed g8/g16 row set is required')
    coefficients = [F(1)]
    for left, right in matching:
        local = ((active >> left) & 1) | (((active >> right) & 1) << 1)
        support = pre_pair_support_masses(local, input_marker=input_marker, mode=mode)
        factor = support[0], support[1]+support[2], support[3]
        result = [F(0)]*(len(coefficients)+2)
        for i, value in enumerate(coefficients):
            for j, weight in enumerate(factor):
                result[i+j] += value*weight
        coefficients = result
    return tuple(coefficients)


@lru_cache(None)
def h_nz_scalar_moment(arity, marker, group=16):
    """Exact E[z^wt(Hx)] for iid NZ values on a fixed arity-set."""
    marker = F(marker)
    if (group not in (8, 16) or type(arity) is not int or not 0 <= arity <= group or marker < 0):
        raise ValueError('Valid arity and nonnegative scalar output marker required')
    numerator = F(0)
    for weight in range(group+1):
        count = sum(comb(arity, overlap)*comb(group-arity, weight-overlap)
                    *exact_support_count(arity, overlap, weight, group)
                    for overlap in range(max(0, arity+weight-group), min(arity, weight)+1))
        numerator += count*marker**weight
    return numerator/255**arity


def pre_column_scalar(active, matching, marker, *, group=16, input_marker=F(1), mode='iid-majorant'):
    """Raw positive column mass; whole-word255^(-8a) is NOT applied here."""
    coefficients = expanded_support_polynomial(active, matching, group=group, input_marker=input_marker, mode=mode)
    return sum((coefficient*h_nz_scalar_moment(arity, F(marker), group)
                for arity, coefficient in enumerate(coefficients)), F(0))


def uniform_space_mixture(active, matching, *, group=16):
    """Exact normalized v=1 mixture by dimension of uniform coordinate space.

    A singleton pair has weights(q-1)/q^2 at0,1/q^2 at1,(q-1)/q at2.
    This is a uniform-input diagnostic, not the native RS measure bridge.
    """
    matching = checked_matching(matching, group)
    if group not in (8, 16) or type(active) is not int or not 0 <= active < 1 << group:
        raise ValueError('A fixed g8/g16 input space is required')
    values = [F(1)]
    for left, right in matching:
        occupancy = (active >> left & 1)+(active >> right & 1)
        factor = ((F(1), F(0), F(0)) if occupancy == 0 else
                  (F(0), F(0), F(1)) if occupancy == 2 else
                  (F(255, 256**2), F(1, 256**2), F(255, 256)))
        updated = [F(0)]*(len(values)+2)
        for i, value in enumerate(values):
            for j, probability in enumerate(factor):
                updated[i+j] += value*probability
        values = updated
    return tuple(values)


def diagnostics():
    """Exactly the five inherited scalar-mark controls; no search/provider."""
    rows = []
    matching, z = xor_matching(8), F(1, 2)
    for active in (1, 3, 15, 255, 65535):
        pre = {mode: pre_column_scalar(active, matching, z, mode=mode)
               for mode in ('exact-support', 'iid-majorant', 'uniform-cap')}
        post = {mode: column_contraction(active, matching, (z,)*16, mode=mode)
                for mode in ('exact-support', 'iid-majorant', 'diagonal')}
        a = active.bit_count()
        rows.append({'active': active, 'arity': a, 'matching_direction': 8,
                     'output_marker': str(z), 'input_marker': '1',
                     'log2_uniform_input_moment_pre': {mode: log2(value)-8*a for mode, value in pre.items()},
                     'log2_uniform_input_moment_post': {mode: log2(value)-8*a for mode, value in post.items()},
                     'safe_pre_improvement_over_diagonal_bits': log2(post['diagonal'])-log2(pre['iid-majorant'])})
    return {'method': 'pre-H-Toeplitz2-positive-support-mixture',
            'exact_mode_is_uniform_input_diagnostic_only': True,
            'iid_majorant_has_RS_measure_bridge': True,
            'bridge': 'pair majorant invariant under independent NZ input diagonals; insert artificial diagonals',
            'same_native_input_diagonals_retained': False, 'native_output_diagonals_retained': True,
            'regional_banks_built': 0, 'witness_searches': 0,
            'directed_distance_receipt': False, 'whole_code_certificate': False, 'results': rows}


if __name__ == '__main__':
    import argparse
    from hashlib import sha256
    import json
    from pathlib import Path
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        parser.error('Output must be a new artifact')
    result = diagnostics()
    result['source_sha256'] = sha256(Path(__file__).read_bytes()).hexdigest()
    with args.output.open('x', encoding='utf-8', newline='\n') as destination:
        destination.write(json.dumps(result, indent=2, allow_nan=False)+'\n')
    print(json.dumps({'output': str(args.output), 'sha256': sha256(args.output.read_bytes()).hexdigest(),
                      'results': result}, allow_nan=False))
