#include "../tests/support.hpp"
#include "strong/myers_v1.hpp"
#include "exploration.hpp"
#include "diagnostic.hpp"
#include <iostream>
#include <chrono>
#include <sched.h>
using namespace strong;
struct alignas(64) Input { PackedDNA128 q[4],t[4]; };
template<int Algorithm>
__attribute__((noinline)) uint64_t profile_kernel(const Input* inputs, uint64_t rounds) {
    uint64_t sum=0;
    for(uint64_t r=0;r<rounds;++r) for(int i=0;i<32;++i) {
        asm volatile("":::"memory");
        const auto& p=inputs[i];
        Quad q{p.q,p.q+1,p.q+2,p.q+3},t{p.t,p.t+1,p.t+2,p.t+3};
        Distances d;
        if constexpr(Algorithm==0) d=v1::myers_2x64x4<128>(q,t);
        if constexpr(Algorithm==1) d=myers_128x4<128>(q,t);
        if constexpr(Algorithm==2) d=exploration::myers_2x64x4_bits<128,false>(q,t);
        if constexpr(Algorithm==3) d=exploration::myers_2x64x4_bits<128,true>(q,t);
        if constexpr(Algorithm==4) d=exploration::myers_2x64x4_skew<128>(q,t);
        if constexpr(Algorithm==5) d=exploration::myers_2x64x4_peeled<128,false>(q,t);
        if constexpr(Algorithm==6) d=exploration::myers_2x64x4_peeled<128,true>(q,t);
        if constexpr(Algorithm==7) d=exploration::myers_2x64x4_raw<128,false>(q,t);
        if constexpr(Algorithm==8) d=exploration::myers_2x64x4_raw<128,true>(q,t);
        if constexpr(Algorithm==9) d=exploration::myers_2x64x4_raw<128,true,2>(q,t);
        if constexpr(Algorithm==10) d=exploration::myers_2x64x4_raw<128,true,4>(q,t);
        if constexpr(Algorithm==11) d=exploration::myers_2x64x4_raw<128,true,8>(q,t);
        if constexpr(Algorithm==12) d=exploration::myers_2x64x4_raw<128,true,4,true>(q,t);
        if constexpr(Algorithm==13) d=exploration::myers_2x64x4_raw<128,true,1,true>(q,t);
        if constexpr(Algorithm==14) d=exploration::myers_2x64x4_raw<128,true,2,true>(q,t);
        if constexpr(Algorithm==15) d=exploration::myers_2x64x4_raw<128,true,8,true>(q,t);
        if constexpr(Algorithm==17) d=diagnostic::myers_128x4_diagnostic<128,true,false>(q,t);
        if constexpr(Algorithm==18) d=diagnostic::myers_128x4_diagnostic<128,false,true>(q,t);
        if constexpr(Algorithm==19) d=diagnostic::myers_128x4_diagnostic<128,true,true>(q,t);
        if constexpr(Algorithm==16) d=myers_2x64x4<128>(q,t);
        for(int j=0;j<4;++j) sum+=d[j];
    }
    return sum;
}
int main(int argc,char**argv) {
    cpu_set_t set;CPU_ZERO(&set);CPU_SET(24,&set);sched_setaffinity(0,sizeof(set),&set);
    std::mt19937_64 rng(0x20261022);Input inputs[32];
    for(auto& x:inputs) for(int j=0;j<4;++j) {
        auto q=test_support::random_dna(rng);
        x.q[j]=encode(q);x.t[j]=encode(test_support::mutate(q,j*3,rng,true));
    }
    int a=std::stoi(argv[1]);uint64_t rounds=std::stoull(argv[2]);
    auto start=std::chrono::steady_clock::now();
    uint64_t sum=0;
    switch(a) {
        case 0:sum=profile_kernel<0>(inputs,rounds);break;
        case 1:sum=profile_kernel<1>(inputs,rounds);break;
        case 2:sum=profile_kernel<2>(inputs,rounds);break;
        case 3:sum=profile_kernel<3>(inputs,rounds);break;
        case 4:sum=profile_kernel<4>(inputs,rounds);break;
        case 5:sum=profile_kernel<5>(inputs,rounds);break;
        case 6:sum=profile_kernel<6>(inputs,rounds);break;
        case 7:sum=profile_kernel<7>(inputs,rounds);break;
        case 8:sum=profile_kernel<8>(inputs,rounds);break;
        case 9:sum=profile_kernel<9>(inputs,rounds);break;
        case 10:sum=profile_kernel<10>(inputs,rounds);break;
        case 11:sum=profile_kernel<11>(inputs,rounds);break;
        case 12:sum=profile_kernel<12>(inputs,rounds);break;
        case 13:sum=profile_kernel<13>(inputs,rounds);break;
        case 14:sum=profile_kernel<14>(inputs,rounds);break;
        case 15:sum=profile_kernel<15>(inputs,rounds);break;
        case 17:sum=profile_kernel<17>(inputs,rounds);break;
        case 18:sum=profile_kernel<18>(inputs,rounds);break;
        case 19:sum=profile_kernel<19>(inputs,rounds);break;
        case 16:sum=profile_kernel<16>(inputs,rounds);break;
    }
    auto ns=std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::steady_clock::now()-start).count();
    std::cout<<"algorithm="<<a<<" pairs="<<rounds*128<<" ns="<<ns<<" checksum="<<sum<<" ns/pair="<<double(ns)/(rounds*128)<<'\n';
}
