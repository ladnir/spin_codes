"""Exact toy checks; serial, <45 seconds, no generated files."""
from collections import Counter
from fractions import Fraction
import itertools
import random
import time
from check_maps import Field
from compact_reset_obstruction import occupancy_coefficients,gate,intervals
from random_feedback_bridge import output_exponent


def main():
    started = time.monotonic()
    f = Field(2,0x7)
    mul = f.table
    # Uniform 2x1 feedback over GF4, for all fixed inputs and nonzero states.
    for Q in range(1,4):
        for X in itertools.product(range(4),repeat=2):
            outputs = Counter()
            states = Counter()
            for T in itertools.product(range(4),repeat=2):
                Y = tuple(x^mul[t][Q] for x,t in zip(X,T))
                outputs[Y] += 1
                states[Q^Y[0]^Y[1]] += 1
            assert len(outputs)==16 and set(outputs.values())=={1}
            assert set(states.values())=={4}
    print('48 exhaustive fixed-input uniform feedback/state laws PASS.',flush=True)

    # Enumerate every output for three b=2,s=1 batches. Count trajectories
    # by their episodes, quiet batches and total symbol weight.
    counts = Counter()
    for y in itertools.product(range(4),repeat=6):
        Q,v,c,j = 0,0,0,0
        for t in range(3):
            a,b = y[2*t:2*t+2]
            new = Q^a^b
            if Q:
                v += 1
            elif new:
                c += 1
                v += 1
            elif a or b:
                j += 1
            Q = new
        counts[v,c,j,sum(x!=0 for x in y)] += 1
    comparisons = 0
    for v,c,j,_ in counts:
        if not c:
            continue
        for W in range(7):
            actual = sum(counts[v,c,j,w] for w in range(W+1))
            for beta in (Fraction(0),Fraction(1,2),Fraction(3,4),Fraction(1)):
                exponent = 2*(v+j)+2*output_exponent(v,c,j,2,1,W,beta)
                if exponent>=0:
                    assert actual**exponent.denominator<=2**exponent.numerator
                else:
                    assert actual==0
                comparisons += 1
    print(f'4096 exact output trajectories; {comparisons} weighted-count comparisons PASS.',flush=True)

    # Exhaust routes of four [4,2] blocks with two selected labels per
    # region. This checks the expectation without block independence.
    choices = list(itertools.combinations(range(4),2))
    values = []
    for selected in itertools.product(choices,repeat=4):
        t = [sum(i in region for region in selected) for i in range(4)]
        values.append(sum(max(0,x-2) for x in t))
    coeff,denom = occupancy_coefficients([2]*4,4)
    mean = Fraction(4*sum(max(0,j-2)*a for j,a in enumerate(coeff)),denom)
    assert mean==Fraction(sum(values),len(values))
    lower = (mean-1)/(4-1)
    assert Fraction(sum(x>1 for x in values),len(values))>=lower
    print('1296 exact route restrictions: dimension expectation/failure lower bound PASS.',flush=True)

    # Direct low-weight witnesses: two RS[4,2] blocks over GF4, state1,
    # b=1, U=the first six routed coordinates, one imposed reset.
    # d(C_U)=2 > one reset equation for every one of the 16 routes.
    rng = random.Random(20261022)
    for swaps in itertools.product((0,1),repeat=4):
        for _ in range(8):
            T = [rng.randrange(4) for _ in range(8)]
            witnesses = 0
            for msg in itertools.product(range(4),repeat=4):
                if not any(msg):
                    continue
                words = [[msg[2*i]^mul[msg[2*i+1]][x] for x in range(4)] for i in range(2)]
                X = [words[i^swaps[r]][r] for r in range(4) for i in range(2)]
                if any(X[6:]):
                    continue
                Q,Y = 0,[]
                for x,t in zip(X,T):
                    y = x^mul[t][Q]
                    Y.append(y)
                    Q ^= y
                if not any(Y[6:]):
                    assert any(Y)
                    witnesses += 1
            assert witnesses>=3
    print('128 exact route/inner fixtures contain reset-supported nonzero witnesses PASS.',flush=True)
    for state in (64,128,138):
        size,mean,rank,lower = gate(state)
        assert size<=491*2**21//1000 and mean>rank and lower>Fraction(1,2**40)
    assert len(intervals())==128 and sum(b-a for a,b in intervals())*128==gate(64)[0]
    assert time.monotonic()-started<45
    print(f'K20 exact obstruction decisions PASS; elapsed {time.monotonic()-started:.3f}s.')


if __name__=='__main__':
    main()
