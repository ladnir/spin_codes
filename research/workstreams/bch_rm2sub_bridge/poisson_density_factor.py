"""Sharper integer density factor; not used by frozen K24/K26 producers."""
import math


def density_factor(length):
    assert type(length) is int and length>=1
    root=math.isqrt(length)
    ceiling=root+int(root*root<length)
    return min(length+1,4*ceiling)
