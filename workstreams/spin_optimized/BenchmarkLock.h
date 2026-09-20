#pragma once
#include <algorithm>
#include <cctype>
#include <filesystem>
#include <stdexcept>
#include <string>
#include <fcntl.h>
#include <sched.h>
#include <sys/file.h>
#include <unistd.h>
static void lockBenchmarks() {
    for(const char* path:{"/tmp/prindal-addition-encoder-benchmark.lock","/tmp/bare-spin-benchmark.lock"}) {
        const int fd=open(path,O_CREAT|O_RDWR,0600);
        if(fd<0 || flock(fd,LOCK_EX|LOCK_NB)) throw std::runtime_error("benchmark lock unavailable");
    }
    for(const auto& e:std::filesystem::directory_iterator("/proc")) {
        const auto id=e.path().filename().string();
        if(id.empty() || !std::all_of(id.begin(),id.end(),[](unsigned char c){return std::isdigit(c);}) ||
           std::stoul(id)==static_cast<unsigned long>(getpid())) continue;
        std::error_code error;
        auto name=std::filesystem::read_symlink(e.path()/"exe",error).filename().string();
        if(error) continue;
        std::transform(name.begin(),name.end(),name.begin(),[](unsigned char c){return std::tolower(c);});
        if(name.find("bench")!=std::string::npos) throw std::runtime_error("another benchmark is active: "+name);
    }
    cpu_set_t cpus;CPU_ZERO(&cpus);CPU_SET(15,&cpus);
    if(sched_setaffinity(0,sizeof(cpus),&cpus)) throw std::runtime_error("cannot pin CPU 15");
}
