#include "perf_gate.hpp"
#include <algorithm>
#include <fstream>
#include <iostream>
#include <numeric>
#include <random>
#include <sched.h>
#include <vector>

int main() {
    cpu_set_t set; CPU_ZERO(&set); CPU_SET(24,&set);
    if(sched_setaffinity(0,sizeof(set),&set)) throw std::runtime_error("affinity");
    PerfGate pmu;
    std::vector<uint64_t> memory((64ULL<<20)/8);
    for(size_t i=0;i<memory.size();++i) memory[i]=i*7+1;
    std::vector<size_t> order(memory.size()/8);
    std::iota(order.begin(),order.end(),0);
    std::mt19937_64 rng(20261022);std::shuffle(order.begin(),order.end(),rng);
    std::cout<<"condition,repeat,loads,cycles,l2_miss,dram_io_fills,checksum\n";
    pmu.start();
    for(int mode=0;mode<2;++mode) for(int repeat=0;repeat<5;++repeat) {
        size_t n=mode?order.size():256;
        if(mode) {
            for(size_t i=0;i<memory.size();i+=8) _mm_clflush(&memory[i]);
        } else {
            for(size_t i=0;i<256;++i) {
                auto v=((volatile uint64_t*)memory.data())[i*8];
                asm volatile(""::"r"(v):"memory");
            }
        }
        _mm_mfence();
        uint64_t sum=0;
        auto p0=pmu.snapshot();
        if(mode) for(auto i:order) sum+=((volatile uint64_t*)memory.data())[i*8];
        else for(int pass=0;pass<1024;++pass) for(size_t i=0;i<256;++i)
            sum+=((volatile uint64_t*)memory.data())[i*8];
        auto p1=pmu.snapshot();
        std::cout<<(mode?"cold_clflush_64MiB":"hot_16KiB")<<','<<repeat<<','
            <<(mode?n:n*1024)<<','<<p1[0]-p0[0]<<','<<p1[1]-p0[1]<<','
            <<p1[2]-p0[2]<<','<<sum<<'\n';
    }
    pmu.stop();
    std::cerr<<"PMU enabled="<<pmu.enabled<<" running="<<pmu.running<<'\n';
}
