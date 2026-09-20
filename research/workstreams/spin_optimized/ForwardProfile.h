#pragma once
#include <chrono>
#include <cstdio>
namespace forward_profile {
using Clock=std::chrono::steady_clock;
struct State {
    unsigned calls=0;double ms[3]{};
    void add(unsigned i,Clock::time_point start) {
        if(calls>3) ms[i]+=std::chrono::duration<double,std::milli>(Clock::now()-start).count();
    }
    ~State() {
        if(calls>3) std::fprintf(stderr,"forward profile: calls=%u outer_ms=%.6f route_ms=%.6f inner_gather_ms=%.6f\n",
                               calls-3,ms[0]/(calls-3),ms[1]/(calls-3),ms[2]/(calls-3));
    }
};
inline State& state() {static thread_local State s;return s;}
}
