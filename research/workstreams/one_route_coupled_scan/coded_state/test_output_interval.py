"""Output enumeration and state-constraint checks. Serial and file-free."""
from collections import Counter
import itertools
import math
import random
import time
from check_maps import Field, Maps, weight
from output_interval import touched_regions, dense_mixed_certificate


def main():
    started = time.monotonic()
    def cap():
        if time.monotonic()-started>50:
            raise RuntimeError('50-second test cap')
    region_cases = 0
    for region in range(1,9):
        for length in range(1,65):
            actual = max((offset+length-1)//region+1 for offset in range(region))
            assert actual==touched_regions(length,region)
            region_cases +=1

    field = Field(2,0x7)
    maps = Maps(field,1)
    counts = Counter()
    for flat in itertools.product(range(4),repeat=8):
        if not any(flat):
            continue
        state,intervals,quiet,begin = 0,[],[],None
        for t in range(4):
            y = flat[2*t:2*t+2]
            scales = [1+(t%3),1+((t+1)%3)]
            next_state = state ^ maps.bmap([field.table[r][z] for r,z in zip(scales,y)])[0]
            if not state:
                if next_state:
                    begin = t
                elif any(y):
                    quiet.append(t)
            elif not next_state:
                intervals.append((begin,t+1))
                begin = None
            state = next_state
        if state:
            intervals.append((begin,4))
        counts[tuple(intervals),tuple(quiet),weight(flat)] +=1
    for (intervals,quiet,w),amount in counts.items():
        v = sum(b-a for a,b in intervals)
        j = len(quiet)
        assert w>=2*j
        upper = math.comb(2*(v+j),w)*4**(w-j)
        assert amount<=upper,(intervals,quiet,w,amount,upper)
    print(f'{region_cases} interval-region checks;65535 GF4 outputs in '
          f'{len(counts)} exact trajectory/weight groups: PASS',flush=True)
    cap()

    # Exact ranks of A_H B_W; column scaling by R does not change these ranks.
    field = Field(4,0x13)
    maps = Maps(field,4)
    matrix = [[field.xor(field.table[a][b] for a,b in zip(row,column))
               for column in zip(*maps.b_matrix)] for row in maps.a_matrix]
    local,strict = 0,0
    selected = Counter()
    for h in range(4):
        for holes in itertools.combinations(range(8),h):
            for mask in range(256):
                support = [j for j in range(8) if mask>>j&1]
                w = len(support)
                actual = field.rank([[matrix[i][j] for j in support] for i in holes]) if h and w else 0
                lower = max(0,h+min(4,w)-4)
                assert actual>=lower
                local +=1
                strict +=actual>lower
                if h==3 and w==3:
                    selected[actual] +=1
        cap()
    print(f'{local} exact A_H B_W rank checks PASS; {strict} strict; '
          f'h=w=3 rank histogram={dict(sorted(selected.items()))}',flush=True)

    # Fresh output columns certify a lower-triangular rank contribution.
    maps = Maps(field,2)
    rng = random.Random(3917)
    global_cases = 0
    for iteration in range(512):
        outputs = [[rng.randrange(16) if iteration%4 or t!=1 else 0 for _ in range(4)]
                   for t in range(3)]
        state = [0]*2
        forms = [[0]*12 for _ in range(2)]
        rows,lower = [],0
        for t,y in enumerate(outputs):
            scales = [rng.randrange(1,16) for _ in range(4)]
            old = state[:]
            state = [a^b for a,b in zip(state,maps.bmap([field.table[r][z] for r,z in zip(scales,y)]))]
            for i in range(2):
                for j in range(4):
                    forms[i][4*t+j] ^= field.table[maps.b_matrix[i][j]][scales[j]]
            if not any(y):
                assert state==old
            if t==2:
                continue  # The terminal state is free: impose no constraints.
            if not any(state):
                rows.extend(row[:] for row in forms)
                lower +=min(2,weight(y))
            else:
                aq = maps.a(state)
                holes = [j for j,z in enumerate(aq) if not z]
                lower +=max(0,len(holes)+min(2,weight(y))-2)
                rows.extend([field.xor(field.table[a][forms[i][k]] for i,a in enumerate(maps.a_matrix[j]))
                             for k in range(12)] for j in holes)
        support = [4*t+j for t,y in enumerate(outputs) for j,z in enumerate(y) if z]
        actual = field.rank([[row[j] for j in support] for row in rows]) if rows and support else 0
        assert actual>=lower
        global_cases +=1
    print(f'{global_cases} joint state-constraint ranks and zero-output persistence checks PASS',flush=True)
    dense_mixed_certificate()
    print(f'Elapsed {time.monotonic()-started:.3f}s; no full certificate claimed.',flush=True)


if __name__=='__main__':
    main()
