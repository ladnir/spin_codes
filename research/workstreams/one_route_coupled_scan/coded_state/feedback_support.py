"""Exact, small regional-moment and output-rank interfaces."""
from fractions import Fraction
from math import comb


def one_block_moment(counts, region_size, redundancy, tau):
    """Poisson-binomial bound base, BEFORE raising to number of outer blocks."""
    assert tau>=1 and region_size>0 and 0<=redundancy<=len(counts)
    probabilities = [Fraction(1)]
    for selected in counts:
        assert 0<=selected<=region_size
        p = Fraction(selected,region_size)
        row = [Fraction(0)]*(len(probabilities)+1)
        for j,mass in enumerate(probabilities):
            row[j] += mass*(1-p)
            row[j+1] += mass*p
        probabilities = row
    return sum(mass*tau**max(0,j-redundancy) for j,mass in enumerate(probabilities))


def binomial_moment(length, redundancy, numerator, denominator, tau):
    """Exact positive numerator; preserves even tiny excesses above one."""
    assert 0<=numerator<=denominator and denominator>0 and tau>=1
    assert 0<=redundancy<=length
    absent = denominator-numerator
    total = denominator**length
    for j in range(redundancy+1,length+1):
        total += comb(length,j)*numerator**j*absent**(length-j)*(tau**(j-redundancy)-1)
    return Fraction(total,denominator**length)


def fresh_dimension(state, holes, weight):
    assert 0<=holes<state and 0<=weight<=2*state
    return weight-max(0,holes+min(state,weight)-state)
