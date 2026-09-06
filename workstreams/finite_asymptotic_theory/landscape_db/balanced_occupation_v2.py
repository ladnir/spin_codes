"""Batched balanced witnesses and reusable power-of-two region coefficients."""
import math

import numpy as np

import activation_occupation as general
import balanced_occupation as original
from composition_occupation import terminal_logs


def region_ladder(epoch, t, lengths, maximum):
    lengths = sorted(set(lengths))
    if (not lengths or maximum < 1 or len(epoch) != min(t,maximum)+1
            or any(L < t or L % t or (L//t)&(L//t-1) for L in lengths)):
        raise ValueError('power-of-two native lengths and complete coefficients required')
    power = epoch+np.array([math.log(math.comb(t,j)) for j in range(len(epoch))])[:,None,None]
    result = {}; length = t
    while length <= lengths[-1]:
        if length in lengths:
            limit = min(maximum,length)
            result[length] = power[:limit+1]-np.array([math.log(math.comb(length,j)) for j in range(limit+1)])[:,None,None]
        length *= 2
        if length <= lengths[-1]:
            power = general.polynomial_product(power,power,maximum)
    return result


class PreparedModel:
    def __init__(self, counts, block, scales, band_count=0):
        self.block = block; self.scales = list(scales)
        if not self.scales: raise ValueError('at least one scale required')
        self.envelopes = []
        for scale in scales:
            bands, roots, lp, ln = original.density_roots(counts,block,scale,band_count)
            self.envelopes.append(original.Envelope(roots,lp,ln))
        self.log_bands = math.log(len(bands))
        self.log_mass = math.log(sum(counts.values()))

    def bounds(self, regions, length, cutoff, lam):
        maximum = len(regions)-1
        starts = np.empty((len(self.scales),maximum,3,3))
        for index,envelope in enumerate(self.envelopes):
            current = regions
            for q in range(maximum):
                current = envelope.apply(current[:-1],current[1:])
                starts[index,q] = current[0]
        moments = terminal_logs(starts.reshape(-1,3,3),self.block).reshape(len(self.scales),maximum)
        occupations = np.arange(1,maximum+1)
        choose = np.array([math.log(math.comb(length,int(q))) for q in occupations])
        values = moments+choose+occupations*self.log_bands+cutoff*lam
        np.minimum(values,choose+occupations*self.log_mass,out=values)
        selected = np.argmin(values,axis=0)
        return values[selected,np.arange(maximum)],np.array(self.scales)[selected]
