"""Independent literal-map, basis, sampler, and generated-SIMD checks for s15.

SIMD lanes below hold symbolic binary supports. Executing the actual C++
intrinsic expressions therefore checks every input bit at once, without
requiring AVX512 hardware. Compiled validation remains a separate check.
"""
from collections import Counter
from itertools import product
from pathlib import Path
import random
import re
import sys
from types import SimpleNamespace
import unittest

HERE=Path(__file__).resolve().parent
REPO=HERE.parents[3]
KERNEL=REPO/'spin/experiments/k16_codesign_100us/kernel'
sys.path.insert(0,str(HERE.parent/'proof'))
import screen_paired_restriction as restriction

PERM=(0,1,2,7,4,5,6,3,8,9,10,11,12,13,14,15)


def vx(a,b):return tuple(x^y for x,y in zip(a,b))
def shuffle(a,b,imm):
    return (a[imm&3],a[(imm>>2)&3],b[(imm>>4)&3],b[(imm>>6)&3])
def mask(src,bits,value,width=2):
    selectors=[(bits>>(width*j))&((1<<width)-1) for j in range(4)]
    assert all(s in (0,(1<<width)-1) for s in selectors)
    return tuple(v if s else x for s,x,v in zip(selectors,src,value))
def permute(a,indices,b):
    assert len(indices)==8
    assert all(indices[2*j]%2==0 and indices[2*j+1]==indices[2*j]+1 for j in range(4))
    return tuple((a+b)[indices[2*j]//2] for j in range(4))
def ternary(a,b,c,imm):
    assert imm==0x96
    return vx(vx(a,b),c)
def transpose4(a,b,c,d):
    columns=list(zip(a,b,c,d))
    return SimpleNamespace(a=columns[0],b=columns[1],c=columns[2],d=columns[3])


OPS={
    '_mm512_xor_si512':vx,
    '_mm512_shuffle_i32x4':shuffle,
    '_mm512_mask_xor_epi64':lambda src,k,a,b:mask(src,k,vx(a,b)),
    '_mm512_mask_shuffle_i32x4':lambda src,k,a,b,i:mask(src,k,shuffle(a,b,i),4),
    '_mm512_ternarylogic_epi64':ternary,
    '_mm512_permutex2var_epi64':permute,
    '_mm512_maskz_permutex2var_epi64':lambda k,a,i,b:mask((0,)*4,k,permute(a,i,b)),
    '_mm512_setr_epi64':lambda *indices:indices,
    '_mm512_loadu_si512':lambda value:value,
    'transpose4':transpose4,
    'totals4':lambda a,b,c,d:tuple(v[0]^v[1]^v[2]^v[3] for v in (a,b,c,d)),
}


class Raw:
    def __init__(self,values):self.values=values
    def __add__(self,offset):return tuple(self.values[offset:offset+4])


def body(source,name):
    source=re.sub(r'//[^\n]*','',source)
    start=source.index('{',source.index(name+'('))+1
    depth=1
    for i in range(start,len(source)):
        depth+=(source[i]=='{')-(source[i]=='}')
        if not depth:return source[start:i]
    raise ValueError('unterminated function body')


def evaluate(expression,env):
    expression=expression.strip()
    if expression.startswith('RawFeedback?'):
        yes,no=expression[len('RawFeedback?'):].split(':',1)
        return evaluate(yes if env['RawFeedback'] else no,env)
    return eval(expression,{'__builtins__':{},**OPS},env)


def execute(source,env):
    """Execute the small straight-line intrinsic subset in the two functions."""
    for statement in source.split(';'):
        statement=statement.strip().lstrip('{} \n\r').strip()
        if not statement or statement.startswith('__m512i '):continue
        statement=re.sub(r'^(?:const auto|auto)\s+','',statement)
        if statement.startswith('emit('):
            key,value=statement[len('emit('):-1].split(',',1)
            env['emitted'][evaluate(key,env)]=evaluate(value,env)
        else:
            left,right=statement.split('=',1)
            left=left.strip()
            value=evaluate(right,env)
            indexed=re.fullmatch(r'(high|feedback\.v)\[(\d+)\]',left)
            if indexed:
                obj=env['high'] if indexed[1]=='high' else env['feedback'].v
                obj[int(indexed[2])]=value
            else:env[left]=value


def rank(rows,bits):
    rows=list(rows);pivot=0
    for column in range(bits):
        selected=next((i for i in range(pivot,len(rows)) if rows[i]>>column&1),None)
        if selected is None:continue
        rows[pivot],rows[selected]=rows[selected],rows[pivot]
        for row in range(pivot+1,len(rows)):
            if rows[row]>>column&1:rows[row]^=rows[pivot]
        pivot+=1
    return pivot


def apply(rows,value):
    return sum(((row&value).bit_count()&1)<<j for j,row in enumerate(rows))


def basis(value):
    return sum(((value>>PERM[j])&1)<<j for j in range(16))


class Paired15Interface(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.wrapper,cls.record=restriction.prepare((10,))
        cls.rows=[int(x,16) for x in cls.record['expansion_rows_hex']]
        cls.columns=cls.record['feedback_columns']
        cls.header=(KERNEL/'T64Paired15Map.h').read_text()
        cls.implementation=(KERNEL/'K16Paired15.cpp').read_text()
        cls.expansion=body(cls.header,'streamHigh')
        cls.feedback=body(cls.implementation,'finishWide')

    def test_literal_deletion_and_complete_state_census(self):
        _,parent=restriction.paired.prepare()
        expected=parent['expansion_rows_hex'][:10]+parent['expansion_rows_hex'][11:]
        self.assertEqual(self.record['expansion_rows_hex'],expected)
        literal=re.search(r'columns\[64\]=\{([^}]+)\}',self.header).group(1)
        self.assertEqual([int(x,16) for x in literal.split(',')],self.columns)
        self.assertEqual(rank(self.rows,64),15)
        self.assertTrue(all((a&b).bit_count()%2==0 for a in self.rows for b in self.rows))
        weights=Counter(sum((state&c).bit_count()%2 for c in self.columns) for state in range(1<<15))
        self.assertEqual({str(k):v for k,v in sorted(weights.items())},self.record['spectrum'])
        self.assertEqual(weights[0],1)
        self.assertTrue(all(c<1<<15 for c in self.columns))

    def test_actual_generated_expansion_and_dummy_independence(self):
        # Coordinate15 is symbolic/nonzero here: expansion must still ignore it.
        state=[1<<PERM[j] for j in range(16)]
        for raw_feedback in (False,True):
            env=dict(state=SimpleNamespace(v=[tuple(state[i:i+4]) for i in range(0,16,4)]),
                     raw=Raw([0]*64),packetBase=0,RawFeedback=raw_feedback,
                     high={},emitted={})
            execute(self.expansion,env)
            emitted=sum((env['emitted'][j] for j in range(16)),())
            self.assertEqual(list(emitted),self.columns)

    def test_actual_wide_feedback_network(self):
        # Every raw input position has a distinct symbolic bit. The generated
        # stream must form exactly C*raw in the permuted state basis.
        for raw_feedback in (False,True):
            env=dict(state=SimpleNamespace(v=[(0,)*4]*4),
                     raw=Raw([1<<p for p in range(64)]),packetBase=0,
                     RawFeedback=raw_feedback,high={},emitted={})
            execute(self.expansion,env)
            env['feedback']=SimpleNamespace(v={})
            execute(self.feedback,env)
            actual=sum((env['feedback'].v[j] for j in range(4)),())
            padded=self.rows+[0]
            self.assertEqual(actual,tuple(padded[j] for j in PERM))
            self.assertEqual(actual[-1],0)

    def test_basis_conjugation_and_dummy_coordinate(self):
        rng=random.Random(913)
        for _ in range(24):
            while True:
                rows=[rng.getrandbits(15) for _ in range(15)]
                if rank(rows,15)==15:break
            rows.append(1<<15)
            fast=[basis(rows[PERM[j]]) for j in range(16)]
            self.assertEqual(fast[15],1<<15)
            self.assertTrue(all(not(r>>15&1) for r in fast[:15]))
            for column in range(16):
                literal=1<<column
                self.assertEqual(apply(fast,basis(literal)),basis(apply(rows,literal)))
            for out in range(2):
                for inp in range(2):
                    # Exact little-byte row packing used by preparePairedBasis.
                    word=sum(((fast[8*out+j]>>(8*inp))&255)<<(8*j) for j in range(8))
                    self.assertEqual([(word>>(8*j))&255 for j in range(8)],
                                     [(fast[8*out+j]>>(8*inp))&255 for j in range(8)])

    def test_uniform_rejection_sampling_small_exact_analogue(self):
        # All equiprobable binary matrices are rejected exactly by singularity;
        # conditioning yields uniform GL. Enumerate the same rank test at n=3.
        matrices=[rows for rows in product(range(8),repeat=3) if rank(rows,3)==3]
        self.assertEqual(len(matrices),(8-1)*(8-2)*(8-4))
        for value in range(1,8):
            self.assertEqual(Counter(apply(rows,value) for rows in matrices),
                             Counter({image:24 for image in range(1,8)}))
        source=(KERNEL/'K16Paired15Setup.cpp').read_text()
        self.assertIn('words()&0x7fffU',source)
        self.assertIn('rows[15]=0x8000U',source)
        self.assertIn('while(!fullRank15(rows))',source)


if __name__=='__main__':
    unittest.main()
