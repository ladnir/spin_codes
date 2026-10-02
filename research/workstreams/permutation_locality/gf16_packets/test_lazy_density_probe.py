from fractions import Fraction as Q
from math import comb
import unittest
from unittest.mock import patch
from flint import arb,ctx
import lazy_density_probe as probe
import lazy_density
import occupancy_birth_classes as birth
import single_packet
import return_moment
from test_scalar_cover import endpoint


class LazyDensityTests(unittest.TestCase):
    def test_actual_refinement_regenerates_census_and_checks_maps(self):
        ctx.prec=192;images=list(range(8));columns=[1,2,4,3,5,7,6,1]
        data=birth.prepare(images,columns,3,3);z=arb(3)/4
        local=birth.outward_at_z(data,z)
        with patch('group_moment.maps',return_value=(images,columns,None)), \
             patch.object(single_packet,'census',wraps=single_packet.census) as census:
            first=lazy_density.refine_actual(data,local,z,2)
            second=lazy_density.refine_actual(data,local,z,2)
            self.assertEqual(census.call_count,2)
            self.assertEqual(first,second)
        with patch('group_moment.maps',return_value=(images,[0]*8,None)):
            with self.assertRaises(ArithmeticError):lazy_density.refine_actual(data,local,z,2)

    def test_caps_and_mixed_products_against_exact_transitions(self):
        ctx.prec=192;z=Q(3,4)
        images=[sum(((a>>i)&1)*v for i,v in enumerate((1,6,120))) for a in range(8)]
        for columns in ([1,2,4,3,5,7,6,1],[0]*8,[1]*8):
            single=single_packet.census(images,columns,arb(3)/4)
            census=return_moment.census(images,columns,3,2)
            for updates in (2,3,4):
                data=birth.prepare(images,columns,3,updates)
                local=birth.outward_at_z(data,arb(3)/4)
                caps=probe.density_caps(data,arb(3)/4,single)
                plain=probe.density_caps(data,arb(3)/4)
                density=[[Q(0)]*8 for _ in range(3)]
                exact=[[[Q(0)]*8 for _ in range(8)] for _ in range(3)]
                for x in range(256):
                    j=int(bool(x&15))+int(bool(x>>4));chance=Q(1,comb(2,j)*15**j)
                    feedback=0
                    for bit,c in enumerate(columns):
                        if x>>bit&1:feedback^=c
                    for state in range(8):
                        weight=chance*z**(images[state]^x).bit_count()
                        if not state:exact[j][state][feedback]+=weight
                        else:
                            density[j][state^feedback]+=weight
                            exact[j][state][state^feedback]+=weight/Q(2**updates)
                            for refreshed in range(1,8):
                                exact[j][state][refreshed^feedback]+=weight*(1-Q(1,2**updates))/7
                for j,cap in enumerate(caps):
                    self.assertGreaterEqual(endpoint(cap),max(density[j][1:]))
                    self.assertLessEqual(endpoint(cap),endpoint(plain[j]))
                refined=return_moment.refine_class_returns(data,local,census,arb(3)/4)
                alternatives=[probe.candidate(data,source,arb(3)/4,caps,through)
                              for source in (local,refined) for through in range(3)]
                for matrices in alternatives:
                    bounds=[[[endpoint(m[i,k]) for k in range(m.ncols())] for i in range(m.nrows())] for m in matrices]
                    for sequence in ((0,0,0),(1,1,1),(2,2,2),(0,1,2),(2,0,1),(1,0,0)):
                        values=[Q(1)]*8;upper=[Q(1)]*len(bounds[0])
                        for j in sequence:
                            values=[sum(a*v for a,v in zip(row,values)) for row in exact[j]]
                            upper=[sum(a*v for a,v in zip(row,upper)) for row in bounds[j]]
                            self.assertGreaterEqual(upper[0],values[0])
                            self.assertGreaterEqual(upper[1],max(values[1:]))
                            self.assertGreaterEqual(upper[2],sum(values[1:])/7)
                            for k,level in enumerate(data['birth_class_levels'],3):
                                self.assertGreaterEqual(upper[k],max(values[a] for a in range(1,8) if images[a].bit_count()==int(level)))

    def test_invalid_geometry_and_output_weight(self):
        data={'windows':2,'bits':3,'updates':2}
        for through in (-1,3):
            with self.assertRaises(ValueError):probe.candidate(data,[],arb(1),[],through)
        for z in (arb(0),arb(2)):
            with self.assertRaises(ValueError):probe.density_caps(data,z)


if __name__=='__main__':unittest.main()
