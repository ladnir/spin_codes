"""FLOAT proposals from normalized regional alias tables for the swap route.

This is the screening analogue of ProjectedBound, not a certificate. The
regional matrix entries use their individual residue sums before histogram
normalization. Every accepted proposal requires a fresh directed replay.
"""
from fractions import Fraction as F
from functools import lru_cache
from math import log2

import numpy as np

from hd_grouped_packet_screen import LIMIT
from hd_randomized_outer_state_profile import multiply,power,total
from hd_shifted_bulk_search import LOG_DENSITY,LOG_FACTORIAL
from hd_shifted_cover_backend import validate_box,profile_count
from hd_shifted_good_pair_screen import (good_step_coefficients,
    good_pair_coefficients,paired_frequency)
from hd_shifted_swap_screen import CONSTRUCTION,swap_coefficients
from hd_shifted_unoriented_screen import low_degree_frequency


@lru_cache(16)
def alias_table(z,alpha,odds,grid=128):
    """Return four full unnormalized log2 alias tables; no target truncation."""
    if not 0<z<1 or not 0<alpha<=1 or len(odds)!=2 or grid<64 or grid&(grid-1):
        raise ValueError('Invalid swap alias proposal')
    dirty = swap_coefficients(z,alpha).copy()
    dirty[:,:2,:] = 0
    dirty *= np.exp2(np.arange(33)[None,:,None]*odds[0]
                     +np.arange(33)[None,None,:]*odds[1])
    dirty = np.fft.fftn(dirty,s=(grid,grid),axes=(1,2))
    good = low_degree_frequency(good_step_coefficients(z,alpha),odds,grid,(64,64))
    composed = low_degree_frequency(good_pair_coefficients(z,alpha),odds,grid,(64,64))
    values = paired_frequency(dirty,good,composed)
    peak = float(np.max(values[:,0,0].real))
    if not 0<peak<float('inf') or not np.all(np.isfinite(values)):
        raise ArithmeticError('Invalid positive matrix evaluation')
    values /= peak
    a,b,c,d = values[:,0,0].real
    rho = (a+d+np.sqrt((a-d)**2+4*b*c))/2
    vector = np.array([b,rho-a] if rho-a>=rho-d else [rho-d,c])
    if not rho>0 or not np.all(vector>0) or not np.all(np.isfinite(vector)):
        raise ArithmeticError('Invalid positive similarity scaling')
    vector /= max(vector)
    values *= (vector[None,:]/vector[:,None]/rho).ravel()[:,None,None]
    for _ in range(6):
        a,b,c,d = values
        trace,bc = a+d,b*c
        values = np.stack((a*a+bc,b*trace,c*trace,d*d+bc))
    aliases = np.fft.ifftn(values,axes=(1,2))
    if np.max(abs(aliases.imag))>1e-10 or np.min(aliases.real)<-1e-10:
        raise ArithmeticError('FLOAT alias diagnostic failed')
    # The cushion is a screening safeguard, not a rigorous error bound.
    logs = np.log2(np.maximum(aliases.real,0)+1e-12)
    logs += 64*(log2(peak)+log2(rho))
    logs += np.log2(vector[:,None]/vector[None,:]).ravel()[:,None,None]
    logs.setflags(write=False)
    return logs


class ProjectedScreen:
    def __init__(self,z,alpha,odds,count_markers,*,grid=128):
        self.z,self.alpha,self.grid = float(z),float(alpha),grid
        self.odds = tuple(tuple(map(float,row)) for row in odds)
        self.count_markers = tuple(map(F,count_markers))
        if len(self.odds)!=2 or len(self.count_markers)!=2 or min(self.count_markers)<=0:
            raise ValueError('Two regional marker pairs and two positive count markers required')
        self.tables = tuple(alias_table(self.z,self.alpha,row,grid) for row in self.odds)
        self.log_x,self.log_y = map(log2,self.count_markers)
        trace = power(((0.,self.log_x+self.log_y),(0.,self.log_x)),8192)
        self.log_count = total((trace[0][0],trace[1][1]))
        self.constant = -self.alpha*(log2(255)+LIMIT*log2(self.z))

    def regional(self,parity,r,f):
        """Vectorized normalized/tilted regional entries, with entries first."""
        r,f = np.broadcast_arrays(np.asarray(r),np.asarray(f))
        if (parity not in (0,1) or r.dtype.kind not in 'iu' or f.dtype.kind not in 'iu'):
            raise ValueError('Parity and histogram counts must be integers')
        r,f = r.astype(np.int64),f.astype(np.int64)
        empty,singleton = 4096-r+f,r-2*f
        if np.any(empty<0) or np.any(singleton<0) or np.any(f<0):
            raise ValueError('Invalid regional histogram')
        normalizer = (LOG_FACTORIAL[4096]-LOG_FACTORIAL[empty]
                      -LOG_FACTORIAL[singleton]-LOG_FACTORIAL[f]
                      +singleton*self.odds[parity][0]+f*self.odds[parity][1])
        tilt = ((self.alpha*LOG_DENSITY-self.log_x-self.log_y)*r/256
                +self.log_y*f/128-normalizer)
        return self.tables[parity][:,singleton%self.grid,f%self.grid]+tilt

    def evaluate(self,ranges):
        ranges = validate_box(ranges)
        if not profile_count(ranges):
            raise ValueError('Empty support box')
        regions,counts = [],[]
        for parity in (0,1):
            r,f = np.meshgrid(np.arange(ranges[0][0],ranges[0][1]+1),
                              np.arange(ranges[parity+1][0],ranges[parity+1][1]+1),indexing='ij')
            valid = (4096-r+f>=0)&(r-2*f>=0)
            r,f = r[valid],f[valid]
            entries = self.regional(parity,r,f)
            regions.append(tuple(map(tuple,np.max(entries,axis=1).reshape(2,2))))
            counts.append(len(r))
        run = power(multiply(*regions),128)
        return dict(ranges=ranges,log2_upper_screen=float(self.log_count+total(run[0])+self.constant),
                    projected_histogram_counts=counts,directed=False,whole_code_certificate=False,
                    construction=CONSTRUCTION)

    def serializable(self):
        return dict(mark=self.z,alpha=self.alpha,marker_log_odds=self.odds,
                    count_markers=list(map(str,self.count_markers)),grid=self.grid,
                    coefficient_method='positive-alias',projected_normalization=True,
                    construction=CONSTRUCTION,pair_swaps=True,
                    orientation_aggregation='independent-fair-swaps')
