"""All-region pair parity conditioned on any one complete region pair type.

Both single-message weights in the conditioned type must be even. The
bound is uniform over fixed occupied supports and compatible bit patterns.
"""
import math
from fractions import Fraction as F
from pathlib import Path
import bridge as base
from certify_tail_parity_mixing import kraw_values


def run():
    family_path = base.HERE/'generated/weight80_family.json'
    family = base.read(family_path)
    for name, digest in family['local_sha256'].items():
        assert base.sha(base.HERE/name) == digest
    close = base.decode(family['close_pair_probability_upper'])
    single = max(F(abs(value), math.comb(255, j))
                 for w in [79, 80] for j, value in enumerate(kraw_values(255, w)) if 0 < j < 255)
    far_tables = [kraw_values(255, c) for c in range(82, 161, 2)]
    far = max(F(abs(values[j]), math.comb(255, j))
              for values in far_tables for j in range(1, 255))
    assert single == F(97, 255) and far == F(91, 255)
    conditioned_seconds = [far+(1-far)*close/marginal**2 for marginal in [F(11, 16), F(5, 16)]]
    paired = max(conditioned_seconds)
    assert single*single <= paired < F(3, 8)
    characters = (1 << 255)-2
    epsilon = characters*single**2620+F(characters**2, 4)*paired**2620
    assert epsilon < F(1, 1 << 3200)
    result = dict(status='EXACT_ONE_PAIR_TYPE_CONDITIONED_GLOBAL_PARITY', occupation=2620,
        conditioned_region_type='arbitrary compatible type with both marginal weights even',
        uniform_over_occupied_supports=True, reference_probability=base.encode(F(1, 1 << 508)),
        single_character_upper=base.encode(single), conditioned_character_seconds=[base.encode(v) for v in conditioned_seconds],
        joint_character_upper=base.encode(paired), relative_error_upper=base.encode(epsilon),
        full_second_moment_certified=False,
        local_sha256={str(p.relative_to(base.HERE)):base.sha(p) for p in
                      [Path(__file__), family_path, base.HERE/'certify_tail_parity_mixing.py']})
    path = base.HERE/'generated/pair_type_conditioned_parity.json'
    if path.exists():
        assert result == base.read(path)
    else:
        base.write_new(path, result)
    print('Conditional pair-type parity: joint gap <=', float(paired),
          'relative error <2^-3200 around2^-508', flush=True)


if __name__ == '__main__':
    run()
