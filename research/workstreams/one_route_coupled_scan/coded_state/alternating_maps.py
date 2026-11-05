"""Inner-only candidate: alternate two disjoint evaluation halves."""
from check_maps import Maps


class AlternatingMaps(Maps):
    def a_at(self,q,phase):
        assert phase in (0,1)
        coefficients = self.f.fft(q,self.b,True)
        return self.f.fft(coefficients+[0]*self.s,phase*self.b)


def run_hole_bound(state_size,length):
    """Fixed nonzero state across length feedback batches, alternating A."""
    return (state_size-1)*((length+1)//2)
