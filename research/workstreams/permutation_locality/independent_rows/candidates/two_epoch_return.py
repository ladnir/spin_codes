"""Exact zero-to-zero moment for two consecutive single-packet epochs.

Local building block only: this does not bound the other output states,
arbitrary occupancies, or the complete code's distance.
"""
import argparse
from collections import Counter
from math import comb, log2

from flint import arb, ctx

import mass_density_screen  # Establish the existing research import paths.
from group_moment import maps


def census(images, columns, width=4):
    if (type(width) is not int or not 1<=width<=4 or not columns
        or len(columns)%width or len(images)<2 or len(images)&(len(images)-1)):
        raise ValueError('valid packet width and power-of-two state space required')
    if images[0]!=0:
        raise ValueError('zero state must have zero expansion')
    packets={b:[] for b in range(1,width+1)}
    seen=set()
    for start in range(0,len(columns),width):
        for mask in range(1,1<<width):
            syndrome=0
            for bit in range(width):
                if mask>>bit&1:syndrome^=columns[start+bit]
            if not 0<syndrome<len(images) or syndrome in seen:
                raise ValueError('single-packet feedback must be nonzero and globally injective')
            seen.add(syndrome)
            packets[mask.bit_count()].append((mask<<start,syndrome))
    result={}
    for a,left in packets.items():
        assert len(left)==len(columns)//width*comb(width,a)
        for b,right in packets.items():
            all_pairs=Counter();equal=Counter()
            for _,s in left:
                for word,t in right:
                    weight=a+(images[s]^word).bit_count()
                    all_pairs[weight]+=1
                    if s==t:equal[weight]+=1
            denominator=len(left)*len(right)
            assert sum(all_pairs.values())==denominator
            assert sum(equal.values())==(len(left) if a==b else 0)
            assert all(n<=all_pairs[w] for w,n in equal.items())
            result[a,b]={'all':all_pairs,'equal':equal,'denominator':denominator}
    return result


def moments(data,states,tilt,penalty='1',rounds=2,width=4):
    if (type(states) is not int or states<1 or type(rounds) is not int or rounds<1
        or arb(tilt)<0 or not 0<arb(penalty)<=1):
        raise ValueError('nonnegative tilt, positive state/update counts, penalty in (0,1] required')
    alpha=arb(2)**(-rounds)
    result={}
    for (a,b),record in data.items():
        def moment(hist):
            return sum((n*(-arb(tilt)*w).exp() for w,n in hist.items()),arb(0))
        scale=arb(penalty)**(int(a==width)+int(b==width))
        result[a,b]=scale*(alpha*moment(record['equal'])
                          +(1-alpha)*moment(record['all'])/states)/record['denominator']
    return result


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--tilt',default='.056')
    parser.add_argument('--penalty',default='.9')
    parser.add_argument('--rounds',type=int,default=2)
    args=parser.parse_args();ctx.prec=192
    images,columns,_=maps()
    data=census(images,columns)
    values=moments(data,len(images)-1,args.tilt,args.penalty,args.rounds)
    for (a,b),record in data.items():
        print('LOCAL TWO-EPOCH RETURN a/b',a,b,
              'equal/total',sum(record['equal'].values()),record['denominator'],
              'moment',values[a,b],flush=True)
    # Separate activation and return maxima, solely as a local control.
    activation={a:(-arb(args.tilt)*a).exp()*arb(args.penalty)**int(a==4)
                for a in range(1,5)}
    coupled=max(float(v.upper()) for v in values.values())
    separate=max(float(v.upper()) for v in activation.values())*max(
        float((v/activation[a]).upper()) for (a,b),v in values.items())
    print('LOCAL CONTROL joint maximum',coupled,'separate maxima',separate,
          'log2 ratio',log2(separate/coupled),flush=True)
    print('Exact integer counts and Arb local moments; no global integration or distance certificate.')


if __name__=='__main__':main()
