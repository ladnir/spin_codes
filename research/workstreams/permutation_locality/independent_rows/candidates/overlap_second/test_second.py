"""Exact direct enumeration for the centered overlap second moment."""
from fractions import Fraction as Q
import unittest
from unittest.mock import patch

from second import moments


def linear(columns,x):
    value=0
    for i,c in enumerate(columns):
        if x>>i&1:value^=c
    return value


class OverlapSecond(unittest.TestCase):
    def test_all_small_shapes_and_translated_states(self):
        cases=[((0b101101,0b011110,0b110001),(1,2,3,4,5,6),2),
               ((0b10101010,0b11001100,0b11110000),(1,3,5,7,2,4,6,1),2),
               ((0b00001111,0b11001100),(1,2,3,1,2,3,2,1),4),
               ((0,0),(0,0,0,0),2)]
        for expansion,feedback,width in cases:
            maximum=min(3,len(feedback)//width)
            with patch('builtins.print'):
                result=moments(expansion,feedback,maximum,width=width)
            words={s:[] for s in result}
            for x in range(1,1<<len(feedback)):
                shape=tuple(sorted(b for b in ((x>>(width*w)&((1<<width)-1)).bit_count()
                                               for w in range(len(feedback)//width)) if b))
                if shape in words:words[shape].append(x)
            for shape,group in words.items():
                D,exact,upper=result[shape];W=sum(shape)
                self.assertEqual(D,len(group))
                for target in range(1<<len(expansion)):
                    moment=Q(sum((2*(x&linear(expansion,target^linear(feedback,x))).bit_count()-W)**2
                                 for x in group),D)
                    if target==0:self.assertEqual(moment,exact)
                    self.assertLessEqual(moment,upper)
                self.assertLessEqual(upper,W*W)

    def test_bad_geometry(self):
        for width,maximum in ((0,1),(2,0),(2,3)):
            with self.assertRaises(ValueError):moments((1,2),(1,2,1,2),maximum,width=width)


if __name__=='__main__':unittest.main()
