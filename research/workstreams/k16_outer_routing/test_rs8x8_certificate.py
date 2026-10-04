"""Independent integer/physical-step audit of the composed RS8x8 proof."""
from fractions import Fraction as F
import hashlib
import json
from math import comb
from pathlib import Path
import sys
import unittest

HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parent/'k16_codesign_100us'/'proof'))
import screen_paired_restriction as maps
from flint import arb,ctx


def rational(pair):
    m,e=pair
    return F(m<<max(e,0),1<<max(-e,0))


def normalize(m,e):
    while not m&1:m>>=1;e+=1
    return m,e


def shells():
    q=2**32
    symbols=[comb(8,h)*sum((-1)**j*comb(h,j)*q**max(0,h-j-4)
        for j in range(h+1)) for h in range(1,9)]
    packet=[0]+[comb(8,k)*15**k for k in range(1,9)]
    power=[1];counts=[F(0)]*65
    for h,number in enumerate(symbols,1):
        nxt=[0]*(len(power)+8)
        for i,a in enumerate(power):
            for j,b in enumerate(packet):nxt[i+j]+=a*b
        power=nxt
        for v,num in enumerate(power):counts[v]+=F(number*num,(q-1)**h)
    assert sum(counts)==2**128-1
    return counts


class CertificateTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        ctx.prec=256
        cls.record=json.loads((HERE/'rs8x8-s15-complete-p256.json').read_text())
        cls.old=json.loads(Path(cls.record['old_receipt']['path']).read_text())
        cls.wrapper,cls.map_record=maps.prepare((10,))

    def test_sources_and_scope(self):
        r=self.record;ref=r['old_receipt']
        self.assertEqual(hashlib.sha256(Path(ref['path']).read_bytes()).hexdigest(),ref['sha256'])
        for path,digest in {**self.old['source_sha256'],**r['new_source_sha256']}.items():
            self.assertEqual(hashlib.sha256(Path(path).read_bytes()).hexdigest(),digest,path)
        self.assertEqual(set(r['occupancy_uppers']),set(map(str,range(1,513))))
        self.assertFalse(r['fresh_replay'])
        self.assertTrue(r['q1_fresh'] and r['q2_fresh'] and r['whole_code_certificate'])
        self.assertEqual(r['map_record'],json.loads(json.dumps(self.map_record)))

    def test_each_exact_rational_tail_transfer(self):
        ratio=F(65535,65537)**4
        self.assertEqual(ratio,F(self.record['tail_transfer']['ratio']))
        self.assertEqual(ratio,F(self.record['tail_transfer']['beta_new'])/F(self.record['tail_transfer']['beta_old']))
        for q in range(3,513):
            expected=rational(self.old['occupancy_uppers'][str(q)])*ratio**q
            m,e=self.record['occupancy_uppers'][str(q)]
            self.assertLessEqual(expected,rational((m,e)))
            self.assertLess(rational((m-1,e)),expected)
            self.assertEqual(m.bit_length(),256)

    def test_exact_union_and_target(self):
        r=self.record;values=list(r['occupancy_uppers'].values());base=min(e for _,e in values)
        total=sum(m<<(e-base) for m,e in values)
        self.assertEqual(normalize(total,base),(int(r['union_exact_mantissa_hex'],16),r['union_exact_exponent']))
        shift=max(0,total.bit_length()-256);ceil=(total+(1<<shift)-1)>>shift
        self.assertEqual(normalize(ceil,base+shift),normalize(*r['union_upper']))
        self.assertLess(rational(r['union_upper']),F(1,2**40))

    def test_fresh_q1_direct_physical_steps(self):
        q1=maps.screen.q1;kernel=q1.kernel_t64;physical=kernel.authenticate(self.wrapper)
        counts=shells();best=[arb(1)]*65
        for tilt in self.record['q1_tilts']:
            z=(-kernel.aq(F(tilt))).exp()
            local=kernel.sparse_kernel.outward_at_z(physical,z)
            local=kernel.kernel_birth_density.refine_local(physical,local,z,F(1,2))
            regional=q1.placement(local,epochs=32,windows=16,maximum_groups=1,rounding=q1.rounded)
            factor=(kernel.aq(F(tilt))*13107).exp()
            moments=q1.support_moments(regional[0],regional[1],64)
            best=[min(a,kernel.up(factor*b)) for a,b in zip(best,moments)]
        fresh=kernel.up(512*sum((kernel.aq(c)*m for c,m in zip(counts,best)),arb(0)))
        saved=kernel.aq(rational(self.record['occupancy_uppers']['1']))
        self.assertLess(abs(fresh/saved-1),arb(2)**-200)
        folded=512*sum((c*rational(p) for c,p in zip(counts,self.record['q1_support_uppers'])),F())
        self.assertLessEqual(folded,rational(self.record['occupancy_uppers']['1']))


if __name__=='__main__':unittest.main()
