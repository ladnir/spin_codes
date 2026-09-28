"""Outward parameter diagnostics at selected row-mixture endpoints.

These points do not constitute a full occupancy cover. In particular,
changed update counts cannot inherit the existing two-update sparse proof.
"""
import argparse
from fractions import Fraction as Q

from flint import arb,ctx

import model,iid_kernel
from affine_mixture import AffineModel
from row_mixture import envelope,pair_components
from bch_joint_support import authenticated_caps


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--minimum-groups',type=int,default=401)
    parser.add_argument('--theta',nargs='+',default=['2/5'])
    parser.add_argument('--updates',type=int,nargs='+',default=[2,3,4,6])
    parser.add_argument('--threshold',type=int,default=model.THRESHOLD)
    parser.add_argument('--central-bits',type=int,default=130)
    parser.add_argument('--component',default='0,3')
    parser.add_argument('--clip-domain',action='store_true')
    args=parser.parse_args()
    if (not 1<=args.minimum_groups<=4096 or not 0<=args.threshold<2097152
            or any(not 1<=r<=32 for r in args.updates) or args.central_bits<1):
        parser.error('bounded group, update, and weight parameters required')
    if any(not 0<Q(t)<Q(1,2) for t in args.theta):parser.error('bias must lie strictly between zero and one half')
    caps=authenticated_caps();ctx.prec=192
    images,columns,_=model.maps();data=iid_kernel.prepare(images,columns,19)
    for raw in args.theta:
        theta=Q(raw);components=pair_components(envelope(caps,1<<args.central_bits,theta))
        selected=[i for i,row in enumerate(components) if row[0]==args.component and row[-1]]
        if len(selected)!=1:parser.error('unknown active component')
        for updates in args.updates:
            instance=AffineModel(components,{**data,'updates':updates},args.minimum_groups,
                                 [1,Q(1,4),Q(1,4)],threshold=args.threshold)
            instance.clip_moments=args.clip_domain
            mean=[Q(args.minimum_groups,4096)*x for x in instance.features[selected[0]]]
            cell=tuple(x for coordinate in mean for x in (coordinate,coordinate))
            if instance.empty(cell):raise ArithmeticError('constructed endpoint was incorrectly pruned')
            _,witness=instance.proposal(cell);bound=instance.outward(cell,witness)
            print('SELECTED ENDPOINT q/component',args.minimum_groups,args.component,'theta',theta,
                  'updates',updates,'cutoff',args.threshold,'log2 upper',bound.log()/arb(2).log(),flush=True)
    print('Selected endpoints only; changed updates require new sparse bounds. No full-code certificate.',flush=True)


if __name__=='__main__':main()
