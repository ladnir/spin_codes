"""Exact toy distinction between one-region invariance and global dependence."""
import itertools
import random
import unittest


class RegionTests(unittest.TestCase):
    def test_one_region_fixed_count_does_not_imply_global_fixed_count(self):
        # T is coordinate-transitive but is not invariant under every permutation.
        tail=(0b0011,0b1100);length=4;rows=4;q=2
        messages=[]
        for support in itertools.combinations(range(rows),q):
            for values in itertools.product(tail,repeat=q):
                message=[0]*rows
                for i,value in zip(support,values):message[i]=value
                messages.append(message)
        rng=random.Random(403901);single=set();joint=set()
        for _ in range(100):
            row_perms=[rng.sample(range(length),length) for _ in range(rows)]
            region_perms=[rng.sample(range(rows),rows) for _ in range(length)]
            one=all_regions=0
            for message in messages:
                events=[]
                for region in range(length):
                    bits=[(message[i]>>row_perms[i][region])&1 for i in region_perms[region]]
                    events.append(bits[0]==bits[1] and bits[2]==bits[3])
                one+=events[0];all_regions+=all(events)
            single.add(one);joint.add(all_regions)
        self.assertEqual(single,{8})
        self.assertGreater(len(joint),1)
        print('Exact toy: single-region count',single,'; all-region counts',sorted(joint))


if __name__=='__main__':unittest.main()
