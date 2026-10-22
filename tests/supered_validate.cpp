#include "supered/supered.hpp"
#include "../src/references.hpp"
#include <iostream>
#include <tuple>

using namespace test_support;

template<int K>
void check(supered::Context<K>& context, const supered::PackedDNA128* q,
           const supered::PackedDNA128* t, const supered::Distances& exact,
           const supered::Distances& shared_exact) {
    supered::Quad qp{q,q+1,q+2,q+3}, tp{t,t+1,t+2,t+3};
    auto four = context.distance4(qp,tp);
    auto shared = context.distance4(q[0],tp);
    int pairs[4], many[4];
    context.distance_pairs(q,t,4,pairs);
    context.distance_many(q[0],t,4,many);
    for(int lane=0; lane<4; ++lane) {
        int expected=exact[lane]<=K?exact[lane]:-1;
        int shared_expected=shared_exact[lane]<=K?shared_exact[lane]:-1;
        if(four[lane]!=expected || pairs[lane]!=expected || shared[lane]!=shared_expected
           || many[lane]!=shared_expected || context.distance(q[lane],t[lane])!=expected)
            throw std::runtime_error("supered API/original Edlib mismatch at K="+std::to_string(K));
    }
}

int main() {
    constexpr int count=131072;
    std::tuple<supered::Context<0>,supered::Context<1>,supered::Context<2>,
        supered::Context<3>,supered::Context<4>,supered::Context<5>,
        supered::Context<50>,supered::Context<63>,supered::Context<64>,
        supered::Context<128>> contexts;
    std::mt19937_64 rng(0x20261022eda);
    uint64_t checksum=0,dp_pairs=0;
    for(int group=0;group<count/4;++group) {
        Sequence queries[4],targets[4];
        supered::PackedDNA128 q[4],t[4];
        supered::Distances exact,shared_exact;
        for(int lane=0;lane<4;++lane) {
            int id=4*group+lane;
            auto a=random_dna(rng),b=a;
            switch(id%6) {
                case 0:b=mutate(a,id%8,rng,false);break;
                case 1:b=mutate(a,id%129,rng,false);break;
                case 2:{
                    int shift=1+(id%64);
                    std::rotate(b.begin(),b.begin()+shift,b.end());
                    for(int j=128-shift;j<128;++j)b[j]="ACGT"[rng()%4];
                    break;
                }
                case 3:b=random_dna(rng);break;
                case 4:a.fill('A');b.fill('T');break;
                default:break;
            }
            queries[lane]=Sequence(a);targets[lane]=Sequence(b);
            q[lane]=queries[lane].packed;t[lane]=targets[lane].packed;
            exact[lane]=reference::api(queries[lane],targets[lane]);
            checksum+=exact[lane];
            if(id%127==0) {
                if(scalar_dp(a,b)!=exact[lane])throw std::runtime_error("independent DP mismatch");
                ++dp_pairs;
            }
        }
        for(int lane=0;lane<4;++lane)
            shared_exact[lane]=reference::api(queries[0],targets[lane]);
        std::apply([&](auto&... ctx){(check(ctx,q,t,exact,shared_exact),...);},contexts);
    }
    std::cout<<"PASS supered independent_pairs="<<count<<" shared_query_pairs="<<count
        <<" thresholds=0,1,2,3,4,5,50,63,64,128 all_public_methods=verified"
        <<" independent_DP_pairs="<<dp_pairs<<" checksum="<<checksum<<'\n';
}
