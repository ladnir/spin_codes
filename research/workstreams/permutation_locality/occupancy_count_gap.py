"""Isolate outer-count slack; OA results require the dual-distance premise.

Exact integer count comparisons only. This is not a full distance verifier.
"""
from math import log2

from bch_joint_support import authenticated_caps, support_caps, rank_total
from shortened_bound import dimension_caps
from basis_lattice import improve_caps
from joint_oa_caps import shell_caps, self_test


def oa_refinement(caps):
    """Upper shells and CDFs for four-tuples, conditional on OA strength 29."""
    shells=[]
    cumulative=[]
    for h,row in enumerate(caps,1):
        # Each h-space has |GL(h,2)| ordered bases, and this many
        # spanning ordered four-tuples. Dependent h-tuples in the OA
        # population only enlarge its upper bound.
        numerator=rank_total(h,4,h)
        denominator=rank_total(h,h,h)
        bounds=shell_caps(256,128,h)
        shell=[min(row[u],numerator*b//denominator) for u,b in enumerate(bounds)]
        cdf=[]
        running=0
        for u,b in enumerate(shell):
            running+=b
            cdf.append(min(row[u],running))
        for u in range(255,-1,-1):
            cdf[u]=min(cdf[u],cdf[u+1])
        assert all(a<=b for a,b in zip(cdf,row))
        assert all(a<=b for a,b in zip(cdf,cdf[1:]))
        shells.append(shell)
        cumulative.append(cdf)
    return shells,cumulative


def main():
    self_test()
    spectrum=authenticated_caps()
    caps=improve_caps(support_caps(spectrum,g=4,dimensions=dimension_caps()),spectrum)
    shells,cdfs=oa_refinement(caps)
    print('u old_CDF_log2 OA_CDF_log2 OA_shell_log2 saved_bits_16_groups',flush=True)
    for u in (38,57,72,76,80,83,84,96,128,176,216,256):
        old=sum(row[u] for row in caps)
        cdf=sum(row[u] for row in cdfs)
        shell=sum(row[u] for row in shells)
        print(u,log2(old),log2(cdf),log2(shell),16*(log2(old)-log2(shell)),flush=True)
    print('OA count comparisons are conditional on dual distance >=30; no inner probability is changed.')


if __name__=='__main__':
    main()
