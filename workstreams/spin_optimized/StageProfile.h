#pragma once
// Diagnostic-only Linux counters. Never included by a normal build.
#include <array>
#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <cstdint>
#include <stdexcept>
#include <linux/perf_event.h>
#include <sys/ioctl.h>
#include <sys/syscall.h>
#include <unistd.h>
namespace spin_profile {
using Clock=std::chrono::steady_clock;
struct Reading {std::uint64_t nr,enabled,running,value[4];};
struct Token {Reading counters;Clock::time_point time;};
struct State {
    bool active=std::getenv("SPIN_PROFILE")!=nullptr;
    int fd[4]={-1,-1,-1,-1};
    unsigned calls=0;
    double ms[3]={};
    double counts[3][4]={};
    double enabled[3]={},running[3]={};
    State() {
        if(!active) return;
        constexpr std::uint64_t events[]={PERF_COUNT_HW_CPU_CYCLES,PERF_COUNT_HW_INSTRUCTIONS,
            PERF_COUNT_HW_CACHE_REFERENCES,PERF_COUNT_HW_CACHE_MISSES};
        for(unsigned i=0;i<4;++i) {
            perf_event_attr a{};a.size=sizeof(a);a.type=PERF_TYPE_HARDWARE;a.config=events[i];
            a.disabled=i==0;a.exclude_kernel=1;a.exclude_hv=1;
            a.read_format=PERF_FORMAT_GROUP|PERF_FORMAT_TOTAL_TIME_ENABLED|PERF_FORMAT_TOTAL_TIME_RUNNING;
            fd[i]=int(syscall(SYS_perf_event_open,&a,0,-1,i?fd[0]:-1,0));
            if(fd[i]<0) throw std::runtime_error("stage perf_event_open failed");
        }
        if(ioctl(fd[0],PERF_EVENT_IOC_ENABLE,PERF_IOC_FLAG_GROUP)) throw std::runtime_error("enable perf failed");
    }
    Reading readCounters() {
        Reading r{};
        if(read(fd[0],&r,sizeof(r))!=sizeof(r) || r.nr!=4) throw std::runtime_error("read perf group failed");
        return r;
    }
    Token begin() {
        if(!active) return {};
        auto r=readCounters();return {r,Clock::now()};
    }
    void end(unsigned stage,const Token& t) {
        if(!active) return;
        const auto stop=Clock::now();const auto r=readCounters();
        if(calls<=3) return; // The harness has exactly three warmup calls.
        ms[stage]+=std::chrono::duration<double,std::milli>(stop-t.time).count();
        for(unsigned j=0;j<4;++j) counts[stage][j]+=r.value[j]-t.counters.value[j];
        enabled[stage]+=r.enabled-t.counters.enabled;running[stage]+=r.running-t.counters.running;
    }
    ~State() {
        if(!active) return;
        if(calls>3) {
            const unsigned n=calls-3;const char* names[]={"inner_route","scatter","bch"};
            std::fprintf(stderr,"{\"profile_calls\":%u,\"stages\":[",n);
            for(unsigned i=0;i<3;++i) std::fprintf(stderr,
                "%s{\"name\":\"%s\",\"mean_ms\":%.9f,\"cycles\":%.3f,\"instructions\":%.3f,"
                "\"cache_references\":%.3f,\"cache_misses\":%.3f,\"counter_running_fraction\":%.9f}",
                i?",":"",names[i],ms[i]/n,counts[i][0]/n,counts[i][1]/n,
                counts[i][2]/n,counts[i][3]/n,enabled[i]?running[i]/enabled[i]:0);
            std::fprintf(stderr,"]}\n");
        }
        for(int f:fd) if(f>=0) close(f);
    }
};
inline State& state() {static thread_local State s;return s;}
}
