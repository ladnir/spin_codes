"""Exact centered second overlap moments from forced-coordinate characters.

This independent candidate is not used by an encoder or a verifier.
It uses unbounded integers for polynomial divisions and signed sums.
"""
from collections import Counter
from fractions import Fraction as Q
from itertools import combinations
from math import comb,factorial,prod
from numbers import Integral
from pathlib import Path
import sys

import numpy as np

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'overlap_mean'))
from mean import kraw,polynomial,shape_tuples


def expansion_characters(expansion,length):
    return [sum(((image>>i)&1)<<b for b,image in enumerate(expansion)) for i in range(length)]


def small_character_data(expansion,feedback,maximum,width):
    """Small exact reference polynomial table; production uses compression."""
    if len(expansion)>8:
        raise ValueError('large maps require the checked compressed character table')
    windows=len(feedback)//width;shapes=shape_tuples(maximum,width)
    table=[]
    for a in range(1<<len(expansion)):
        histogram=Counter(sum((a&feedback[width*w+b]).bit_count()%2 for b in range(width))
                          for w in range(windows))
        coefficients=polynomial(tuple(histogram[h] for h in range(width+1)),width,maximum)
        table.append([coefficients[s] for s in shapes])
    return shapes,np.array(table,dtype=object).T,np.arange(len(table),dtype=np.int64)


def moments(expansion,feedback,maximum,*,width=4,character_data=None):
    """shape -> (D, E[(2H_0-W)^2], upper for every translated target)."""
    if (type(width) is not int or width<1 or not feedback or len(feedback)%width
            or type(maximum) is not int or not 1<=maximum<=len(feedback)//width
            or not expansion or any(type(x) is not int or x<0 for x in (*expansion,*feedback))
            or any(x>=1<<len(feedback) for x in expansion)
            or any(x>=1<<len(expansion) for x in feedback)):
        raise ValueError('matching maps and a complete packet-window geometry required')
    windows=len(feedback)//width
    expected=shape_tuples(maximum,width)
    shapes,table,inverse=(small_character_data(expansion,feedback,maximum,width)
                          if character_data is None else character_data)
    if (shapes!=expected or table.ndim!=2 or table.shape[0]!=len(shapes)
            or table.dtype.kind not in 'iO' or inverse.dtype.kind not in 'iu'
            or (table.dtype.kind=='O' and any(not isinstance(v,Integral) or isinstance(v,bool) for v in table.flat))
            or len(inverse)!=1<<len(expansion) or np.any(inverse<0)
            or np.any(inverse>=table.shape[1])):
        raise ValueError('complete compatible character coefficients required')
    index={s:i for i,s in enumerate(shapes)}
    degree=[sum(s) for s in shapes]
    limits={m:sum(d<=m for d in degree) for m in range(maximum+1)}
    previous=[[(index[tuple(c-int(i==b) for i,c in enumerate(s))],b)
               for b in range(width) if s[b]] for s in shapes]
    pairs=list(combinations(range(width),2))+[(b,b) for b in range(width)]
    previous_two=[]
    for s in shapes:
        row=[]
        for g,(b,c) in enumerate(pairs):
            if s[b] and s[c] and (b!=c or s[b]>=2):
                child=tuple(n-int(i==b)-int(i==c) for i,n in enumerate(s))
                row.append((index[child],g))
        previous_two.append(row)
    denominators=[comb(windows,sum(s))*factorial(sum(s))//prod(factorial(n) for n in s)
                  *prod(comb(width,b)**s[b-1] for b in range(1,width+1)) for s in shapes]
    weights=[sum((b+1)*n for b,n in enumerate(s)) for s in shapes]
    assert all(int(table[i,int(inverse[0])])==D for i,D in enumerate(denominators))
    assert all(int(table[i].min())>=-D and int(table[i].max())<=D for i,D in enumerate(denominators))
    characters=expansion_characters(expansion,len(feedback))
    signed=[0]*len(shapes);absolute=[0]*len(shapes)
    cache={};checked=0

    def quotient(column,removed):
        key=column,removed
        if key in cache:return cache[key]
        limit=limits[max(0,maximum-len(removed))]
        if len(removed)==1:
            source=table[:limit,column]
        else:
            source=quotient(column,removed[:-1])
        factors=[kraw(width,removed[-1],b) for b in range(1,width+1)]
        result=[0]*limit
        for i in range(limit):
            result[i]=int(source[i])-sum(result[p]*factors[b] for p,b in previous[i])
        cache[key]=result
        return result

    for i in range(len(feedback)):
        for l in range(i):
            character=characters[i]^characters[l]
            signs=[(character&column).bit_count()%2 for column in feedback]
            wi,wl=i//width,l//width
            hi=sum(signs[width*wi:width*(wi+1)])
            hl=sum(signs[width*wl:width*(wl+1)])
            if wi==wl:
                poly=quotient(int(inverse[character]),(hi,))
                factors=[(-1)**(signs[i]+signs[l])*kraw(width-2,hi-signs[i]-signs[l],b-2)
                         for b in range(1,width+1)]
                parents=previous
            else:
                poly=quotient(int(inverse[character]),tuple(sorted((hi,hl))))
                first=[(-1)**signs[i]*kraw(width-1,hi-signs[i],b-1) for b in range(1,width+1)]
                second=[(-1)**signs[l]*kraw(width-1,hl-signs[l],b-1) for b in range(1,width+1)]
                factors=[first[b]*second[c]+(first[c]*second[b] if b!=c else 0) for b,c in pairs]
                parents=previous_two
            for j,row in enumerate(parents):
                value=sum(poly[p]*factors[b] for p,b in row)
                assert abs(value)<=denominators[j]
                signed[j]+=value;absolute[j]+=abs(value)
            checked+=1
            if checked%1024==0:
                print('EXACT OVERLAP SECOND coordinate pairs',checked,'quotients',len(cache),flush=True)
    result={}
    for s,W,D,S,A in zip(shapes,weights,denominators,signed,absolute):
        if not W:continue
        exact=Q(W)+Q(2*S,D)
        upper=min(Q(W*W),Q(W)+Q(2*A,D))
        assert 0<=exact<=upper<=W*W
        shape=tuple(b+1 for b,n in enumerate(s) for _ in range(n))
        result[shape]=D,exact,upper
    print('EXACT OVERLAP SECOND complete pairs',checked,'shapes',len(result),flush=True)
    return result


def check_all_full(result,expansion,feedback,maximum,width=4):
    """Independent univariate check when every active packet is all one."""
    windows=len(feedback)//width
    columns=[]
    for w in range(windows):
        column=0
        for b in range(width):column^=feedback[width*w+b]
        columns.append(column)
    characters=expansion_characters(expansion,len(feedback))
    signed=[0]*(maximum+1);absolute=signed.copy()
    for i in range(len(feedback)):
        for l in range(i):
            a=characters[i]^characters[l]
            signs=[(a&column).bit_count()%2 for column in columns]
            forced={i//width,l//width};removed=sum(signs[w] for w in forced)
            for j in range(len(forced),maximum+1):
                value=(-1)**removed*kraw(windows-len(forced),sum(signs)-removed,j-len(forced))
                signed[j]+=value;absolute[j]+=abs(value)
    for j in range(1,maximum+1):
        D=comb(windows,j);W=width*j
        expected=D,Q(W)+Q(2*signed[j],D),min(Q(W*W),Q(W)+Q(2*absolute[j],D))
        assert result[(width,)*j]==expected
    print('SECOND independent univariate all-full checks',maximum,flush=True)


def build(maximum=10):
    from group_moment import maps
    from feedback_character_census import character_polynomials
    images,feedback,_=maps()
    shapes,coefficients,denominators,_,inverse=character_polynomials(maximum)
    expansion=[images[1<<b] for b in range(19)]
    result=moments(expansion,feedback,maximum,character_data=(shapes,coefficients,inverse))
    for shape,(D,exact,upper) in result.items():
        counts=tuple(shape.count(b) for b in range(1,5))
        assert D==denominators[counts]
    check_all_full(result,expansion,feedback,maximum)
    for j in range(1,maximum+1):
        print('SECOND all-four',j,'exact',result[(4,)*j][1],'universal',result[(4,)*j][2],flush=True)
    return result


if __name__=='__main__':build()
