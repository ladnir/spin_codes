"""Portable identities behind the exact native-field outer variants."""
from pathlib import Path
import re
import unittest

HERE=Path(__file__).resolve().parent


def unpack(a,b,high):
    # AVX512 unpacklo/hi_epi64 acts separately in every 128-bit lane.
    return sum((a[16*i+8*high:16*i+8*high+8]+b[16*i+8*high:16*i+8*high+8]
                for i in range(4)),[])


def affine(value,matrix):
    return sum(((value & (matrix >> (8*(7-i)) & 255)).bit_count() & 1) << i
               for i in range(8))


class OuterIdentities(unittest.TestCase):
    def test_pack_byte_sources(self):
        a,b=list(range(64)),list(range(64,128))
        index=[8*(2*((7-c)%4)+(7-c)//4)+byte for byte in range(8) for c in range(8)]
        for high in (0,1):
            expected=[(a+b)[16*(7-c)+8*high+byte] for byte in range(8) for c in range(8)]
            gathered=unpack(a,b,high)
            self.assertEqual([gathered[i] for i in index],expected)

    def test_unpack_byte_destinations(self):
        a,b=list(range(64)),list(range(64,128))
        index=[8*byte+lane//2+4*(lane%2) for lane in range(8) for byte in range(8)]
        x,y=[a[i] for i in index],[b[i] for i in index]
        for high in (0,1):
            expected=sum(([8*byte+coordinate for byte in range(8)]+
                          [64+8*byte+coordinate for byte in range(8)]
                          for coordinate in range(4*high,4*high+4)),[])
            self.assertEqual(unpack(x,y,high),expected)

    def test_shared_gfni_constants_all_byte_inputs(self):
        source=(HERE/'OuterVariants.cpp').read_text()
        constants={int(k):int(v,16) for k,v in
                   re.findall(r'm(\d+)=_mm512_set1_epi64\((0x[0-9a-f]+)ULL\)',source)}
        for x in range(256):
            t2,t4,t8=[affine(x,constants[k]) for k in (2,4,8)]
            for c,result in ((5,t4^x),(11,t8^t2^x),(12,t8^t4),(13,t8^t4^x)):
                self.assertEqual(affine(x,constants[c]),result)

    def test_non_temporal_variant_keeps_arithmetic_and_aligned_stores(self):
        original=(HERE/'OuterVariants.cpp').read_text().split('template<unsigned Plane,bool Shared>')[0]
        nt=(HERE/'OuterNonTemporal.cpp').read_text()
        copied=nt.split('// Retain the exact shared-loop arithmetic and scheduling.')[0]
        self.assertEqual(copied.replace('OuterNonTemporal.h','OuterVariants.h'),original)
        self.assertIn('reinterpret_cast<std::uintptr_t>(out)&63U',nt)
        self.assertIn('fieldLoopSharedParity(in,out,p,t);',nt)
        self.assertIn('_mm_sfence();',nt)
        for group in range(3):
            for plane in (0,2):
                for symbol in range(8):
                    for half in (0,1):
                        address=16*(128*group+32*(plane+half)+4*symbol)
                        self.assertEqual(address%64,0)


if __name__=='__main__':
    unittest.main()
