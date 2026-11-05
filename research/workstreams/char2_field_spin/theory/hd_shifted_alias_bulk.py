"""FLOAT box witnesses using one uniform regional coefficient envelope.

A complete positive alias table bounds every exact coefficient: each
coefficient belongs to one nonnegative residue sum. Taking the entrywise
maximum over all residues removes target dependence from the numerator.
This module proposes witnesses; CachedBound performs directed replay.
"""
from dataclasses import dataclass
from math import log2

import numpy as np

from hd_shifted_bulk_search import PositiveWitness,LOG_DENSITY
from hd_grouped_packet_screen import LIMIT
from hd_randomized_outer_state_profile import multiply,power,total
from hd_shifted_unoriented_screen import unsigned_coefficients,low_degree_frequency
from hd_shifted_good_pair_screen import good_step_coefficients,good_pair_coefficients,paired_frequency


def uniform_regional(z,alpha,odds,grid=128,cushion=1e-12):
    """Return log2 of a target-independent FLOAT regional matrix envelope."""
    if (not 0<z<1 or not 0<alpha<=1 or len(odds)!=2
            or grid<64 or grid&(grid-1) or cushion<0):
        raise ValueError('Valid marks and a power-of-two grid at least64 required')
    dirty = unsigned_coefficients(z,alpha).copy()
    dirty[:,:2,:] = 0
    dirty *= np.exp2(np.arange(33)[None,:,None]*odds[0]+np.arange(33)[None,None,:]*odds[1])
    dirty = np.fft.fftn(dirty,s=(grid,grid),axes=(1,2))
    # No degree truncation: the good-good polynomial has degree at most64.
    good = low_degree_frequency(good_step_coefficients(z,alpha),odds,grid,(64,64))
    composed = low_degree_frequency(good_pair_coefficients(z,alpha),odds,grid,(64,64))
    values = paired_frequency(dirty,good,composed)
    a,b,c,d = values[:,0,0].real
    eigenvalue = (a+d+np.sqrt((a-d)**2+4*b*c))/2
    vector = np.array([b,eigenvalue-a] if eigenvalue-a>=eigenvalue-d else [eigenvalue-d,c])
    if not eigenvalue>0 or not np.all(vector>0) or not np.all(np.isfinite(vector)):
        raise ArithmeticError('Degenerate similarity scaling')
    vector /= max(vector)
    values *= (vector[None,:]/vector[:,None]/eigenvalue).ravel()[:,None,None]
    for _ in range(6):
        a,b,c,d = values
        trace,bc = a+d,b*c
        values = np.stack((a*a+bc,b*trace,c*trace,d*d+bc))
    aliases = np.fft.ifftn(values,axes=(1,2))
    maxima = np.max(aliases.real,axis=(1,2))
    if np.max(abs(aliases.imag))>1e-10 or np.min(maxima)<100*max(1e-12,cushion):
        raise ArithmeticError('FLOAT alias cancellation diagnostic failed')
    logs = np.log2(maxima+cushion).reshape((2,2))
    logs += 64*log2(eigenvalue)+np.log2(vector[:,None]/vector[None,:])
    return logs


@dataclass(frozen=True)
class AliasWitness(PositiveWitness):
    grid: int = 128

    @classmethod
    def create(cls,z,alpha,odds,centre=(0,0,0),*,grid=128):
        odds = tuple(tuple(map(float,row)) for row in odds)
        cached = {row:uniform_regional(z,alpha,row,grid) for row in set(odds)}
        regions = [tuple(map(tuple,cached[row])) for row in odds]
        run = power(multiply(*regions),128)
        numerator = float(total(run[0]))
        intercept = numerator-alpha*(log2(255)+LIMIT*log2(z))
        (s0,f0),(s1,f1) = odds
        slopes = (alpha*LOG_DENSITY+1-128*(s0+s1),-2+256*s0-128*f0,256*s1-128*f1)
        return cls(float(z),float(alpha),odds,tuple(centre),intercept,slopes,numerator,grid)

    def serializable(self):
        return dict(super().serializable(),grid=self.grid,uniform_aliases=True,
                    coefficient_method='positive-alias')
