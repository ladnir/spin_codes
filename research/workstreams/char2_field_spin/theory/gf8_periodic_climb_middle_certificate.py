"""Directed middle certificates for chronological periodic puncturing.

Each packet retains its stated prefix of the eight full-Cauchy outputs.
The state persists across phases. Exact positive matrix products preserve
the order of the puncture pattern; no averaged packet matrix is used.
"""
from __future__ import annotations

import argparse
from fractions import Fraction as F
import json

from gf8_periodic_climb_profile import metadata
from gf8_rate_climb_middle_certificate import certify_profile
from return_cover_certificate import Witness


def witnesses(entries):
    return tuple((F(anchor),Witness(*(F(value) for value in values)))
                 for anchor,values in entries)


WITNESSES_84_6667=witnesses((
    ('3/100',('2977/9730','3928/5977','41649/6467','41611/9661')),
    ('1/20',('5493/9287','3551/7434','60999/8633','18585/7466')),
    ('1/10',('6644/3833','1871/9608','21001/2121','10281/9635')),
    ('7/50',('22569/6937','814/8829','114301/8128','6356/9977')),
    ('1/5',('11939/1712','163/5106','82055/3238','193/606')),
    ('3/10',('26527/1290','23/3631','576454/7711','1013/9572')),
    ('2/5',('541361/9044','10/7787','1990765/7807','21/625')),
    ('1/2',('624592/9925','1/1036','255','122/5423')),
    ('3/5',('529743/8156','4/5169','255','146/9689')),
    ('7/10',('366167/5515','6/9301','255','80/8209')),
    ('4/5',('274401/4067','3/5423','255','55/9572')),
    ('17/20',('308252/4539','2/3873','255','39/9428'))))

WITNESSES_74_56=witnesses((
    ('3/100',('1758/5663','3584/5771','55585/7877','15009/2663')),
    ('1/20',('5139/8596','3290/7587','49527/6302','23377/7200')),
    ('1/10',('14442/8489','1603/9859','114589/9996','5689/4113')),
    ('7/50',('22757/7440','607/8133','166112/9999','532/643')),
    ('1/5',('1229/203','122/4701','32520/1081','1664/3947')),
    ('3/10',('75707/5013','41/7347','330537/3982','254/1667')),
    ('2/5',('145/4','7/5520','255','251/4569')),
    ('1/2',('166997/4408','2/2099','255','5/136')),
    ('3/5',('366977/9421','7/9164','255','11/447')),
    ('7/10',('347394/8749','6/9407','255','34/2139')),
    ('4/5',('214487/5326','1/1828','255','85/9104')),
    ('9/10',('114350/2809','1/2088','255','5/1176'))))


def certify(dimension=84,pattern=(6,6,6,7),first=F(3,100),last=None,k_log=20,
            initial_grid=1000,cell_bits=128,rational_witnesses=None):
    pattern=tuple(pattern)
    defaults={(84,(6,6,6,7)):(WITNESSES_84_6667,F(17,20)),
              (74,(5,6)):(WITNESSES_74_56,F(9,10))}
    if rational_witnesses is None:
        if (dimension,pattern) not in defaults:
            raise ValueError('supply rational witnesses for this nondefault periodic profile')
        rational_witnesses,default_last=defaults[dimension,pattern]
        if last is None:
            last=default_last
    if last is None:
        raise ValueError('specify the exclusive final density')
    return certify_profile(metadata(dimension,pattern,k_log),first,last,rational_witnesses,
                           initial_grid,cell_bits)


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--dimension',type=int,default=84)
    parser.add_argument('--pattern',type=int,nargs='+',default=[6,6,6,7])
    parser.add_argument('--first',type=F,default=F(3,100))
    parser.add_argument('--last',type=F)
    args=parser.parse_args()
    print(json.dumps(certify(args.dimension,tuple(args.pattern),args.first,args.last)),flush=True)
