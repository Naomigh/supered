#include "../tests/support.hpp"
#include "strong/myers_v2.hpp"
#include "strong/banded.hpp"
#include "initial_band.hpp"
#include <iostream>
#include <chrono>
#include <sched.h>
using namespace strong;
struct alignas(64) Input {PackedDNA128 q[4],t[4];};
template<int Algorithm>
__attribute__((noinline)) uint64_t kernel(const Input* inputs,uint64_t rounds) {
    uint64_t sum=0;
    for(uint64_t r=0;r<rounds;++r) for(int i=0;i<32;++i) {
        asm volatile("":::"memory");const auto& p=inputs[i];
        Quad q{p.q,p.q+1,p.q+2,p.q+3},t{p.t,p.t+1,p.t+2,p.t+3};
        Distances d;
        if constexpr(Algorithm==0) d=v2::myers_2x64x4<50>(q,t);
        if constexpr(Algorithm==1) d=initial_band::myers_banded64x4<50>(q,t);
        if constexpr(Algorithm==2) d=myers_banded64x4<50,1>(q,t);
        if constexpr(Algorithm==3) d=myers_banded64x4<50,2>(q,t);
        if constexpr(Algorithm==4) d=myers_banded64x4<50,4>(q,t);
        if constexpr(Algorithm==5) d=myers_banded64x4<50,8>(q,t);
        for(int j=0;j<4;++j)sum+=d[j]+1;
    }
    return sum;
}
int main(int argc,char**argv) {
    cpu_set_t set;CPU_ZERO(&set);CPU_SET(24,&set);sched_setaffinity(0,sizeof(set),&set);
    std::mt19937_64 rng(0x20261022);Input inputs[32];
    for(auto& x:inputs)for(int j=0;j<4;++j){auto q=test_support::random_dna(rng);x.q[j]=encode(q);x.t[j]=encode(test_support::mutate(q,j*15,rng,true));}
    uint64_t rounds=std::stoull(argv[2]);int a=std::stoi(argv[1]);
    auto start=std::chrono::steady_clock::now();uint64_t sum=0;
    switch(a) {
        case 0:sum=kernel<0>(inputs,rounds);break;
        case 1:sum=kernel<1>(inputs,rounds);break;
        case 2:sum=kernel<2>(inputs,rounds);break;
        case 3:sum=kernel<3>(inputs,rounds);break;
        case 4:sum=kernel<4>(inputs,rounds);break;
        case 5:sum=kernel<5>(inputs,rounds);break;
    }
    auto ns=std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::steady_clock::now()-start).count();
    std::cout<<"algorithm="<<a<<" pairs="<<rounds*128<<" ns="<<ns<<" checksum="<<sum<<" ns/pair="<<double(ns)/(rounds*128)<<'\n';
}
