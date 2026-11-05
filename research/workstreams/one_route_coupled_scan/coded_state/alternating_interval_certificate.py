"""Complete restricted trajectory classes, not a whole-code certificate."""
import time
from half_update_bound import pair_weight_bound


def certificate(end,bonus,numerator,denominator,expected_cells,expected_worst):
    n,N,M = 16384,1<<21,8192
    threshold = 491*N//1000
    worst,where,cells = None,None,0
    for v in range(1,end+1):
        live = v-1
        touched = (live+126)//64 if live else 0
        dimension = M*max(0,touched-128)
        pairs = live//2
        qualifying = (numerator*pairs+denominator-1)//denominator
        for outside_zero in range(33):
            quiet = n-v-outside_zero
            if quiet<0 or 65*quiet>threshold:
                continue
            log_bound = 128*(v+quiet)+128*(threshold-64*quiet+dimension
                                          -65*live-bonus*qualifying)
            if worst is None or log_bound>worst:
                worst,where = log_bound,(v,outside_zero,quiet,qualifying)
            cells +=1
    assert cells==expected_cells and worst==expected_worst,(cells,worst,where)
    # At most n^2 intervals, 33 outside-zero counts and n^32 position tuples.
    # The count contains all v,z choices: no extra cell multiplier is needed.
    assert n*n*33*n**32 < 2**482
    assert worst+482 < -40
    print(f'v<= {end}, bonus={bonus}, qualifying >=ceil({numerator}/{denominator} *pairs): '
          f'{cells} cells; max log2={worst}, at(v,z,j,g)={where}; '
          f'full class failure <=2^-{-worst-482}',flush=True)


if __name__=='__main__':
    start = time.monotonic()
    assert pair_weight_bound(64,1)[0]-130==20
    certificate(8700,20,9,10,269742,-706688)
    certificate(8800,20,1,1,273042,-770688)
    certificate(10000,63,3,4,312642,-73216)
    # A diagnostic of this relaxation, not a counterexample to the code.
    live = 16384-1
    feedback_floor = 65*live+63*(live//2)
    remaining = 491*(1<<21)//1000+(1<<20)-feedback_floor
    assert feedback_floor==1580928 and remaining==497349
    print(f'Full-length containing-interval relaxation remains vacuous: '
          f'{remaining} positive field-symbol exponent before support counting.')
    print(f'Elapsed {time.monotonic()-start:.3f}s. Only these restricted classes are certified.')
