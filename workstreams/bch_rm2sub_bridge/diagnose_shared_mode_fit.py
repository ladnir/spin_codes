"""Inspect pointwise slack of the saved shared-rate discovery fit."""
import math
import numpy as np
from scipy.special import logsumexp
import bridge as base
import region_log_cache as cache


def run():
    screen=base.read(base.HERE/'generated/shared_mode_slopes_m10_g32_screen.json')
    grid=screen['grid'];region=cache.get(screen['tilt_tenth']).reshape(-1,9)
    modes=np.full((2*grid+1,9),-np.inf);modes[0]=region[0]
    for segment in screen['segments']:
        lo,hi,i=segment['lower'],segment['upper'],segment['label']
        coefficient=np.max(region[lo:hi+1]-np.arange(lo,hi+1)[:,None]*math.log(i/grid),axis=0)
        modes[i]=np.maximum(modes[i],coefficient)
    for j in [128,256,384,512,640,768,800,896,1024,1200,1536,2048,4096,6144,8191]:
        fitted=logsumexp(modes[1:]+j*np.log(np.arange(1,2*grid+1)/grid)[:,None],axis=0)
        print('j',j,'log slack Z->Z',round(float(fitted[0]-region[j,0]),3),
            'max finite-entry slack',round(float(np.max(fitted-region[j])),3))


if __name__=='__main__':run()
