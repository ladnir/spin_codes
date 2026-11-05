"""Bounded floating diagnostic. No certificate, no files, no encoder changes."""
import math
import time

N, K, M = 1 << 21, 1 << 20, 8192
U, W, MASK_N, HMAX = 8193*128, 491*N//1000, N-128, 8192*63


def logchoose(n, k):
    if not 0 <= k <= n:
        return -math.inf
    return (math.lgamma(n+1)-math.lgamma(k+1)-math.lgamma(n-k+1))/math.log(2)


def route_log(h, bits):
    """Complemented moment after equalizing all 256 regional fractions."""
    p = (h+128)/N
    terms = [logchoose(256,j)+j*math.log2(p)+(256-j)*math.log2(1-p)
             + bits*(j-128)+math.log2(-math.expm1(-bits*(j-128)*math.log(2)))
             for j in range(129,257)]
    peak = max(terms)
    excess = peak+math.log2(sum(2**(a-peak) for a in terms))
    log_moment = (excess+math.log1p(2**(-excess))/math.log(2) if excess>0
                  else math.log1p(2**excess)/math.log(2))
    return M*log_moment


def point(h):
    v = min(W,U-h-(64*(U-W)+126)//127)
    # All free-output support sizes: binomial volume with a dimension ceiling.
    output = v+(math.log2(U+1)+logchoose(U,max(v,U//2)))/128
    mask = min(logchoose(MASK_N,h),local_mask_log(h))
    return min((mask+math.log2(tau)*(output-K)+route_log(h,math.log2(tau)),tau)
               for tau in (2,4,8,9,10,11,12,16,32,64,128,256,1024,4096,65536))


def local_mask_log(h):
    if h==0:
        return 0.0
    if h==HMAX:
        return 8191*logchoose(256,63)+logchoose(128,63)
    def tilted(bits):
        ans = -h*bits
        for length,multiplicity in ((256,8191),(128,1)):
            a = [logchoose(length,j)+j*bits for j in range(64)]
            peak = max(a)
            ans += multiplicity*(peak+math.log2(sum(2**(x-peak) for x in a)))
        return ans
    return min(tilted(bits) for bits in (-16,-12,-8,-6,-4,-3,-2,-1,0,1,2,3,4,6,8,12,16))


if __name__ == '__main__':
    started = time.monotonic()
    for h in (0,1000,10000,U-W,30000,60000,100000,200000,300000,400000,480000,HMAX):
        print(f'H={h}: {point(h)}',flush=True)
        if time.monotonic()-started>20:
            raise RuntimeError('20-second diagnostic cap')
    print('FLOAT diagnostics only; no complete coverage or directed arithmetic.')
