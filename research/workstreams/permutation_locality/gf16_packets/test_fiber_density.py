from fractions import Fraction as Q
from math import comb
from itertools import product
import unittest
from flint import ctx
import birth_classes
import fiber_density
from test_scalar_cover import endpoint


class FiberDensityTests(unittest.TestCase):
    def test_class_allocation_composes_for_positive_terminal_costs(self):
        ctx.prec=192;z=Q(3,4)
        images=[0x93*(a&1)^0x65*((a>>1)&1)^0x3c*((a>>2)&1) for a in range(8)]
        for columns in ([1,2,4,3,6,5,7,1],[0]*8):
            for updates in (2,3,4):
                data=fiber_density.attach(birth_classes.prepare(images,columns,3,updates))
                local=__import__('occupancy_birth_classes').outward_at_z(data,fiber_density.aq(z))
                actual=fiber_density.candidate(data,local,fiber_density.aq(z),1,2,allocation='classes',include_uniform=True)
                levels=list(data['birth_class_levels']);alpha=Q(1,2**updates)
                exact=[[[Q(0)]*8 for _ in range(8)] for _ in range(3)]
                for x in range(256):
                    j=int(bool(x&15))+int(bool(x>>4));chance=Q(1,comb(2,j)*15**j)
                    feedback=0
                    for i,c in enumerate(columns):
                        if x>>i&1:feedback^=c
                    for a in range(8):
                        weight=chance*z**(images[a]^x).bit_count()
                        if not a:exact[j][a][feedback]+=weight
                        else:
                            exact[j][a][a^feedback]+=alpha*weight
                            for b in range(1,8):exact[j][a][b^feedback]+=(1-alpha)*weight/7
                bounds=[[[endpoint(m[i,k]) for k in range(m.ncols())] for i in range(m.nrows())] for m in actual]
                for sequence in product(range(3),repeat=3):
                    values=[1+Q(a,7) for a in range(8)]
                    upper=[values[0],max(values[1:]),sum(values[1:])/7,
                        *(max(values[a] for a in range(1,8) if images[a].bit_count()==level) for level in levels)]
                    for j in sequence:
                        values=[sum(c*v for c,v in zip(row,values)) for row in exact[j]]
                        upper=[sum(c*v for c,v in zip(row,upper)) for row in bounds[j]]
                        self.assertGreaterEqual(upper[0],values[0])
                        self.assertGreaterEqual(upper[1],max(values[1:]))
                        self.assertGreaterEqual(upper[2],sum(values[1:])/7)
                        for i,level in enumerate(levels,3):
                            self.assertGreaterEqual(upper[i],max(values[a] for a in range(1,8) if images[a].bit_count()==level))

    def test_every_feedback_target_including_zero(self):
        ctx.prec=192;bits=3;z=Q(2,3);width=8
        images=[0x93*(a&1)^0x65*((a>>1)&1)^0x3c*((a>>2)&1) for a in range(8)]
        for columns in ([1,2,4,3,6,5,7,1],[0]*8,[1]*8):
            data=fiber_density.attach(birth_classes.prepare(images,columns,bits,3))
            caps=fiber_density.profile_caps(data,fiber_density.aq(z))
            profiles=[tuple(map(int,h)) for h in data['histograms']]
            tables=[[[Q(0)]*8 for _ in range(3)] for _ in range(8)]
            for a in range(1,8):
                for x in range(1<<width):
                    feedback=0
                    for i,c in enumerate(columns):
                        if x>>i&1:feedback^=c
                    j=int(bool(x&15))+int(bool(x>>4))
                    tables[a][j][feedback]+=z**(images[a]^x).bit_count()/(comb(2,j)*15**j)
                histogram=tuple(sum(((images[a]>>(4*w))&15).bit_count()==v for w in range(2)) for v in range(5))
                profile=profiles.index(histogram)
                for j in range(3):
                    for target in range(8):
                        self.assertGreaterEqual(endpoint(caps[j][profile]),tables[a][j][target])

    def test_uniform_class_allocation_preserves_already_allocated_branch(self):
        from flint import arb
        ctx.prec=192;z=fiber_density.aq(Q(3,4))
        images=[0x93*(a&1)^0x65*((a>>1)&1)^0x3c*((a>>2)&1) for a in range(8)]
        data=fiber_density.attach(birth_classes.prepare(images,[1,2,4,3,6,5,7,1],3,3))
        local=__import__('occupancy_birth_classes').outward_at_z(data,z)
        # Model an earlier valid allocation of the U lazy branch.
        local[1][2,2]=fiber_density.up(local[1][2,2]+arb(1)/8);local[1][2,1]=arb(0)
        actual=fiber_density.candidate(data,local,z,1,2,allocation='classes',include_uniform=True)
        for k in range(local[1].ncols()):self.assertEqual(actual[1][2,k],local[1][2,k])
        for value,allocation in (('true','classes'),(True,'density')):
            with self.assertRaises(ValueError):fiber_density.candidate(data,local,z,1,2,allocation=allocation,include_uniform=value)

    def test_uniform_replacement_selects_whole_rows_and_preserves_other_coordinates(self):
        from flint import arb
        ctx.prec=192;z=fiber_density.aq(Q(3,4))
        images=[0x93*(a&1)^0x65*((a>>1)&1)^0x3c*((a>>2)&1) for a in range(8)]
        data=fiber_density.attach(birth_classes.prepare(images,[1,2,4,3,6,5,7,1],3,3))
        source=__import__('occupancy_birth_classes').outward_at_z(data,z)
        local=fiber_density.candidate(data,source,z,1,2,allocation='classes')
        local[1][2,2]=fiber_density.up(local[1][2,2]+arb(1)/8);local[1][2,1]=arb(0)
        expected=fiber_density.candidate(data,source,z,1,1,allocation='classes',include_uniform=True)
        actual=fiber_density.replace_uniform_classes(data,local,source,z,1,1)
        for j,matrix in enumerate(actual):
            for i in range(matrix.nrows()):
                for k in range(matrix.ncols()):
                    truth=expected[j][i,k] if j==1 and i==2 else local[j][i,k]
                    self.assertEqual(matrix[i,k],truth)
        with self.assertRaises(ValueError):fiber_density.replace_uniform_classes(data,local[:-1],source,z,1,1)

    def test_reallocated_lazy_branch_dominates_each_outgoing_state(self):
        ctx.prec=192;z=Q(3,4);images=[0x93*(a&1)^0x65*((a>>1)&1)^0x3c*((a>>2)&1) for a in range(8)]
        for columns in ([1,2,4,3,6,5,7,1],[0]*8):
            for updates in (2,3,4):
                data=fiber_density.attach(birth_classes.prepare(images,columns,3,updates))
                local=__import__('occupancy_birth_classes').outward_at_z(data,fiber_density.aq(z))
                actual=fiber_density.candidate(data,local,fiber_density.aq(z),1,2)
                alpha=Q(1,2**updates);L=7
                for j in (1,2):
                    self.assertEqual(actual[j][1,1],0)
                    for a in range(1,8):
                        masses=[Q(0)]*8
                        for x in range(256):
                            if int(bool(x&15))+int(bool(x>>4))!=j:continue
                            feedback=0
                            for i,c in enumerate(columns):
                                if x>>i&1:feedback^=c
                            masses[a^feedback]+=alpha*z**(images[a]^x).bit_count()/(comb(2,j)*15**j)
                        level=images[a].bit_count();source=3+list(data['birth_class_levels']).index(level)
                        for row in (1,source):
                            added=actual[j][row,2]-local[j][row,2]
                            for target in range(1,8):
                                self.assertGreaterEqual(endpoint(added/L),masses[target])
                    self.assertEqual(actual[j][2,1],local[j][2,1])
                self.assertEqual(actual[0],local[0])

    def test_combined_constraints_bound_each_outgoing_class_mass(self):
        ctx.prec=192;z=Q(3,4)
        for basis in ((0x93,0x65,0x3c),(0x11,0x22,0x44)):
            images=[basis[0]*(a&1)^basis[1]*((a>>1)&1)^basis[2]*((a>>2)&1) for a in range(8)]
            for columns in ([1,2,4,3,6,5,7,1],[0]*8):
                for updates in (2,3,4):
                    data=fiber_density.attach(birth_classes.prepare(images,columns,3,updates))
                    local=__import__('occupancy_birth_classes').outward_at_z(data,fiber_density.aq(z))
                    actual=fiber_density.candidate(data,local,fiber_density.aq(z),1,2,allocation='classes')
                    levels=list(data['birth_class_levels']);alpha=Q(1,2**updates)
                    for j in (1,2):
                        for a in range(1,8):
                            masses=[Q(0)]*len(levels)
                            for x in range(256):
                                if int(bool(x&15))+int(bool(x>>4))!=j:continue
                                feedback=0
                                for i,c in enumerate(columns):
                                    if x>>i&1:feedback^=c
                                b=a^feedback
                                if b:masses[levels.index(images[b].bit_count())]+=alpha*z**(images[a]^x).bit_count()/(comb(2,j)*15**j)
                            source=3+levels.index(images[a].bit_count())
                            for row in (1,source):
                                self.assertEqual(actual[j][row,1],0)
                                self.assertEqual(actual[j][row,2],local[j][row,2])
                                for k,mass in enumerate(masses,3):
                                    added=actual[j][row,k]-local[j][row,k]
                                    self.assertGreaterEqual(endpoint(added),mass)
                    self.assertEqual(actual[0],local[0])


if __name__=='__main__':unittest.main()
