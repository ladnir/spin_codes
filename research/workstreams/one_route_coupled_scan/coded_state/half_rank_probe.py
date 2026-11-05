"""Bounded exact diagnostic of AB half-block ranks; not a certificate."""
from collections import defaultdict
import itertools
import time
from check_maps import Field, Maps


def main():
    started = time.monotonic()
    field = Field(4,0x13)
    maps = Maps(field,4)
    matrix = [[field.xor(field.table[a][b] for a,b in zip(row,column))
               for column in zip(*maps.b_matrix)] for row in maps.a_matrix]
    profiles = defaultdict(dict)
    for h in range(1,4):
        for holes in itertools.combinations(range(8),h):
            hl,hr = set(holes)&set(range(4)),set(holes)&set(range(4,8))
            for mask in range(1,256):
                support = [j for j in range(8) if mask>>j&1]
                wl,wr = set(support)&set(range(4)),set(support)&set(range(4,8))
                actual = field.rank([[matrix[i][j] for j in support] for i in holes])
                profile = len(hl),len(hr),len(wl),len(wr),len(hl&wl),len(hr&wr)
                profiles[profile].setdefault(actual,(holes,tuple(support)))
        if time.monotonic()-started>30:
            raise RuntimeError('30-second cap')
    ambiguous = [(key,examples) for key,examples in profiles.items() if len(examples)>1]
    print(f'{len(profiles)} half/intersection profiles; '
          f'{len(ambiguous)} have more than one actual rank.',flush=True)
    for key,examples in ambiguous[:3]:
        print('Profile',key,'rank examples',examples,flush=True)
    print(f'Elapsed {time.monotonic()-started:.3f}s. No large-field rank claim.',flush=True)


if __name__=='__main__':
    main()
