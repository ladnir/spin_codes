"""An admissible setup with a large bad-message cluster, not its prevalence.

This excludes a useful sup-count bound combined with the retained mean lower.
It does not refute the target probability for a random setup.
"""
import functools
import math
import operator
from pathlib import Path
from fractions import Fraction as F
import bridge as base


def run():
    base.load_map('t128_s15')
    selection_path = base.HERE/'inputs/t128_s15_selection.json'
    columns = [int(x, 16) for x in base.read(selection_path)['selected']['B_columns_hex']]
    found = None
    for a in range(1, 128):
        for b in range(a+1, 128):
            plane = [0, a, b, a ^ b]
            if functools.reduce(operator.xor, (columns[j] for j in plane)) == 0:
                cosets, unused = [], set(range(128))
                while unused:
                    shift = min(unused)
                    coset = sorted(shift ^ j for j in plane)
                    cosets.append(coset)
                    unused.difference_update(coset)
                if all(functools.reduce(operator.xor, (columns[j] for j in coset)) == 0
                       for coset in cosets):
                    found = (plane, cosets)
                    break
        if found:
            break
    assert found is not None
    plane, cosets = found
    assert len(cosets) == 32 and sorted(j for c in cosets for j in c) == list(range(128))
    quartets = [[128*block+j for j in coset] for block in range(64) for coset in cosets][:655]
    assert len(set(j for quartet in quartets for j in quartet)) == 2620
    family_path = base.HERE/'generated/weight80_family.json'
    family = base.read(family_path)
    for name, digest in family['local_sha256'].items():
        assert base.sha(base.HERE/name) == digest
    cluster = family['shell_size_lower']**655
    mean = base.decode(family['first_moment_lower'])
    assert cluster > 1 << 64000 and mean/cluster < F(1, 1 << 44000)
    result = dict(status='ADMISSIBLE_ALIGNED_SETUP_LARGE_CLUSTER',
                  configuration='t128_s15', row_permutations='identity', region_permutations='identity',
                  multipliers='arbitrary nonzero, irrelevant on zero-state branch',
                  plane=plane, kernel_quartet_partition=cosets,
                  selected_row_quartets=quartets, occupation=2620, output_weight=209600,
                  actual_bad_message_count_lower=cluster,
                  retained_mean_lower_over_cluster_lower=base.encode(mean/cluster),
                  random_setup_failure_above_2_to_minus40_proved=False,
                  local_sha256={str(p.relative_to(base.HERE)):base.sha(p)
                                for p in [Path(__file__), selection_path, family_path]})
    path = base.HERE/'generated/aligned_setup_cluster.json'
    if path.exists():
        assert result == base.read(path)
    else:
        base.write_new(path, result)
    print('Aligned admissible setup:', plane, 'partitioned into32 kernel quartets;',
          'bad cluster log2 lower', math.log2(cluster),
          'retained mean / cluster <2^-44000', flush=True)


if __name__ == '__main__':
    run()
