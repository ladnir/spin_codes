from itertools import combinations
from random import Random
import unittest

import gf8_packet_maps as maps


class PacketMapTests(unittest.TestCase):
    def test_all_1976_width8_square_minors(self):
        count=0
        for phase in range(4):
            a,_=maps.matrices(8,4,phase)
            for size in range(1,5):
                for rows in combinations(range(8),size):
                    for columns in combinations(range(4),size):
                        direct=maps.determinant([[a[i][j] for j in columns] for i in rows])
                        self.assertNotEqual(direct,0)
                        self.assertEqual(direct,maps.cauchy_minor(8,4,phase,rows,columns))
                        count+=1
        self.assertEqual(count,1976)

    def test_zero_feedback_and_extended_invertibility(self):
        for width,state in ((8,4),(16,8),(32,8),(32,16)):
            for phase in range(4):
                a,c=maps.matrices(width,state,phase)
                self.assertEqual(c,tuple(zip(*a)))
                self.assertTrue(all(value==0 for row in maps.matmul(c,a) for value in row))
                block=[tuple(int(i==j) for j in range(width))+tuple(maps.MUL[2][v] for v in a[i])
                       for i in range(width)]
                block += [tuple(c[i])+tuple(int(i==j) for j in range(state)) for i in range(state)]
                self.assertEqual(maps.rank(block),width+state)
                returned=maps.returned_matrix(width,state,phase)
                self.assertEqual(maps.matmul(returned,returned),
                                 tuple(tuple(int(i==j) for j in range(width)) for i in range(width)))

    def test_additive_power_sums(self):
        for width,state in ((8,4),(16,8),(32,8),(32,16)):
            for phase in range(4):
                labels,_,_=maps.geometry(width,state,phase)
                powers=[1]*width
                for degree in range(2*state-1):
                    self.assertEqual(maps.xor_sum(powers),0)
                    powers=[maps.MUL[a][b] for a,b in zip(powers,labels)]

    def test_circulant_gram_and_singleton_returns(self):
        for width,state in ((8,4),(16,8),(32,8),(32,16)):
            for phase in range(4):
                a,c=maps.matrices(width,state,phase)
                for first in range(0,width,state):
                    kernel=a[first]
                    self.assertEqual(maps.xor_sum(kernel),1)
                    for i in range(state):
                        for j in range(state):self.assertEqual(a[first+i][j],kernel[i^j])
                gram=maps.matmul(a,c)
                for i in range(width):
                    for j in range(width):
                        if i//state==j//state:self.assertEqual(gram[i][j],int(i==j))
                        else:self.assertNotEqual(gram[i][j],0)
                returned=maps.returned_matrix(width,state,phase)
                for column in zip(*returned):self.assertEqual(sum(value!=0 for value in column),width-state+1)
                self.assertEqual(returned,maps.returned_matrix(width,state,0))

    def test_width8_all_shortened_return_rank_bounds(self):
        # A nonzero vector on I cannot have output supported in allowed,
        # because the remaining rows restrict injectively to I.
        returned=maps.returned_matrix(8,4)
        minima=maps.proved_return_minima(8,4)
        self.assertEqual(minima,(5,5,4,2,1,1,1,1))
        count=0
        for weight,minimum in enumerate(minima,1):
            for inputs in combinations(range(8),weight):
                for allowed in combinations(range(8),minimum-1):
                    forbidden=[i for i in range(8) if i not in allowed]
                    self.assertEqual(maps.rank([[returned[i][j] for j in inputs] for i in forbidden]),weight)
                    count+=1
        self.assertEqual(count,6309)

    def test_wide_minor_formula_samples(self):
        rng=Random(83216)
        for width,state in ((16,8),(32,8),(32,16)):
            for phase in range(4):
                a,_=maps.matrices(width,state,phase)
                for size in range(1,state+1):
                    for _ in range(4):
                        rows=sorted(rng.sample(range(width),size))
                        columns=sorted(rng.sample(range(state),size))
                        self.assertEqual(maps.determinant([[a[i][j] for j in columns] for i in rows]),
                                         maps.cauchy_minor(width,state,phase,rows,columns))
                        self.assertNotEqual(maps.cauchy_minor(width,state,phase,rows,columns),0)

    def test_packet_adjoint(self):
        rng=Random(804)
        dot=lambda x,y:maps.xor_sum(maps.MUL[a][b] for a,b in zip(x,y))
        for width,state in ((8,4),(16,8),(32,8),(32,16)):
            for phase in range(4):
                a,c=maps.matrices(width,state,phase)
                for _ in range(8):
                    x=tuple(rng.randrange(256) for _ in range(width))
                    s=tuple(rng.randrange(256) for _ in range(state))
                    v=tuple(rng.randrange(256) for _ in range(width))
                    h=tuple(rng.randrange(256) for _ in range(state))
                    diagonal=tuple(rng.randrange(1,256) for _ in range(state))
                    r=tuple(maps.MUL[d][value] for d,value in zip(diagonal,s))
                    y=tuple(value^maps.MUL[2][e] for value,e in zip(x,maps.matvec(a,r)))
                    successor=tuple(value^e for value,e in zip(r,maps.matvec(c,x)))
                    dual_input=tuple(value^e for value,e in zip(v,maps.matvec(a,h)))
                    dual_state=tuple(maps.MUL[d][value^maps.MUL[2][e]]
                                     for d,value,e in zip(diagonal,h,maps.matvec(c,v)))
                    self.assertEqual(dot(y,v)^dot(successor,h),dot(x,dual_input)^dot(s,dual_state))

    def test_costs_and_validation(self):
        self.assertEqual([maps.packed_cost(w,s)['packet_gfni512'] for w,s in
                          ((8,4),(16,8),(32,8),(32,16))],[19,70,138,268])
        self.assertEqual(maps.packed_cost(8,4)['packet_xor_or_xor3'],11)
        self.assertEqual(maps.packed_cost(8,4)['reused_record_broadcasts'],12)
        self.assertEqual(maps.packed_cost(16,8)['packet_xor_or_xor3'],38)
        with self.assertRaises(ValueError):maps.matrices(8,8)
        with self.assertRaises(ValueError):maps.matrices(8,4,4)
        with self.assertRaises(ValueError):maps.returned_matrix(8,4,eta=1)

    def test_cubic_actual_subfields_and_old_geometry(self):
        from cauchy_diagonal_packet_reference import LABELS,NODES
        self.assertEqual(set(maps.GF4),{v for v in range(256) if maps.power(v,4)==v})
        self.assertEqual(NODES,maps.CUBIC_NODES)
        union=set()
        for phase in range(8):
            labels,nodes=maps.cubic_geometry(phase)
            self.assertEqual(labels,LABELS[phase//2][8*(phase%2):8*(phase%2+1)])
            self.assertTrue(union.isdisjoint(labels))
            union.update(labels)
            self.assertTrue(set(labels).isdisjoint(nodes))
        self.assertEqual(len(union),64)

    def test_cubic_all_7904_square_minors(self):
        count=0
        old_bad_counts=[]
        for phase in range(8):
            a,c=maps.cubic_matrices(phase)
            _,old=maps.cubic_matrices(phase,injection_shift=0)
            old=tuple(zip(*old))
            bad=0
            for size in range(1,5):
                for rows in combinations(range(8),size):
                    for columns in combinations(range(4),size):
                        self.assertNotEqual(maps.determinant([[a[i][j] for j in columns] for i in rows]),0)
                        self.assertNotEqual(maps.determinant([[c[j][i] for j in columns] for i in rows]),0)
                        bad+=maps.determinant([[old[i][j] for j in columns] for i in rows])==0
                        count+=2
            old_bad_counts.append(bad)
        self.assertEqual(count,7904)
        self.assertEqual(old_bad_counts,[4,4,0,0,4,4,0,0])
        _,old=maps.cubic_matrices(0,injection_shift=0)
        self.assertEqual(maps.determinant([[old[j][i] for j in (1,2)] for i in (0,6)]),0)

    def test_cubic_feedback_pairing_and_invertibility(self):
        expected=maps.cubic_returned_matrix()
        for phase in range(8):
            a,c=maps.cubic_matrices(phase)
            labels,nodes=maps.cubic_geometry(phase)
            self.assertTrue(all(v==0 for row in maps.matmul(c,a) for v in row))
            ac=maps.matmul(a,c)
            for i,u in enumerate(labels):
                for j,v in enumerate(labels):self.assertEqual(ac[i][j],maps.power(u^v,3)^3)
                denominator=maps.power(u,4)^u^maps.power(16,4)^16
                for j,v in enumerate(nodes):self.assertEqual(maps.MUL[a[i][j]][u^v],denominator)
            block=[tuple(int(i==j) for j in range(8))+tuple(maps.MUL[2][v] for v in a[i]) for i in range(8)]
            block += [tuple(c[i])+tuple(int(i==j) for j in range(4)) for i in range(4)]
            self.assertEqual(maps.rank(block),12)
            returned=maps.cubic_returned_matrix(phase)
            self.assertEqual(returned,expected)
            self.assertEqual(maps.matmul(returned,returned),tuple(tuple(int(i==j) for j in range(8)) for i in range(8)))

    def test_cubic_return_minima_exhaustive(self):
        returned=maps.cubic_returned_matrix()
        self.assertEqual(maps.exact_return_minima(returned),maps.CUBIC_RETURN_MINIMA)
        count=0
        # Independent simpler lower-bound test: all vectors supported on I,
        # including proper subsets, are excluded from every smaller output J.
        for weight,minimum in enumerate(maps.CUBIC_RETURN_MINIMA,1):
            for inputs in combinations(range(8),weight):
                for allowed in combinations(range(8),minimum-1):
                    self.assertEqual(maps.rank([[returned[i][j] for j in inputs]
                                               for i in range(8) if i not in allowed]),weight)
                    count+=1
        self.assertEqual(count,11329)

    def test_cubic_packed_maps_basis_scalars(self):
        for phase in range(8):
            a,c=maps.cubic_matrices(phase)
            for adjoint in (False,True):
                emission,injection=(tuple(zip(*c)),tuple(zip(*a))) if adjoint else (a,c)
                for column in range(12):
                    for scalar in range(256):
                        values=tuple(scalar if j==column else 0 for j in range(8))
                        state=tuple(scalar if j+8==column else 0 for j in range(4))
                        emitted,injected=maps.cubic_packed_maps(values,state,phase,adjoint)
                        self.assertEqual(emitted,maps.matvec(emission,state))
                        self.assertEqual(injected,maps.matvec(injection,values))

    def test_cubic_full_packet_and_adjoint(self):
        rng=Random(843)
        dot=lambda x,y:maps.xor_sum(maps.MUL[a][b] for a,b in zip(x,y))
        for phase in range(8):
            a,c=maps.cubic_matrices(phase)
            for _ in range(32):
                values=tuple(rng.randrange(256) for _ in range(8))
                state=tuple(rng.randrange(256) for _ in range(4))
                outputs=tuple(rng.randrange(256) for _ in range(8))
                terminal=tuple(rng.randrange(256) for _ in range(4))
                dx=tuple(rng.randrange(1,256) for _ in range(8))
                ds=tuple(rng.randrange(1,256) for _ in range(4))
                y,t=maps.cubic_step(values,state,dx,ds,phase)
                x=tuple(maps.MUL[d][v] for d,v in zip(dx,values))
                r=tuple(maps.MUL[d][v] for d,v in zip(ds,state))
                self.assertEqual(y,tuple(v^maps.MUL[2][e] for v,e in zip(x,maps.matvec(a,r))))
                self.assertEqual(t,tuple(v^e for v,e in zip(r,maps.matvec(c,x))))
                xd,sd=maps.cubic_step(outputs,terminal,dx,ds,phase,adjoint=True)
                self.assertEqual(dot(y,outputs)^dot(t,terminal),dot(values,xd)^dot(state,sd))
        self.assertEqual(maps.cubic_packed_cost()['packet_gfni512'],9)
        self.assertEqual(maps.cubic_packed_cost()['packet_xor_or_xor3'],15)


if __name__=='__main__':unittest.main()
