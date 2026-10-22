#include "../src/references.hpp"
#include "strong/banded.hpp"
#include <iostream>
using namespace strong;
using namespace test_support;
template<int K>
void check(const Quad& q,const Quad& t,const Distances& exact,uint64_t group) {
    auto a=myers_banded64x4<K,1>(q,t),b=myers_banded64x4<K,2>(q,t);
    auto c=myers_banded64x4<K,4>(q,t),d=myers_banded64x4<K,8>(q,t);
    for(int i=0;i<4;++i) {
        int e=exact[i]<=K?exact[i]:-1;
        if(a[i]!=e||b[i]!=e||c[i]!=e||d[i]!=e) {
            std::cerr<<"group="<<group<<" lane="<<i<<" K="<<K<<" expected="<<e
                <<" variants="<<a[i]<<','<<b[i]<<','<<c[i]<<','<<d[i]<<'\n';
            throw std::runtime_error("banded oracle mismatch");
        }
    }
}
int main() {
    constexpr uint64_t exhaustive=3*65536,total=exhaustive+32768;
    std::mt19937_64 rng(0x20261022b64);uint64_t checksum=0,dp_count=0;
    for(uint64_t g=0;g<total/4;++g) {
        Sequence q[4],t[4];Quad qp,tp;Distances exact;
        for(int lane=0;lane<4;++lane) {
            uint64_t id=4*g+lane;auto a=random_dna(rng),b=a;
            if(id<exhaustive) {
                int pos=(id/65536)*60;unsigned x=(id%65536)/256,y=id%256;
                for(int j=0;j<8;++j) {a[pos+j]="AC"[(x>>j)&1];b[pos+j]="AC"[(y>>j)&1];}
            } else {
                unsigned h=(id-exhaustive)%65,lo=rng()%(129-h);
                unsigned hi=lo+h+rng()%(129-lo-h);
                if(id%2) {
                    std::rotate(b.begin()+lo,b.begin()+lo+h,b.begin()+hi);
                    for(unsigned j=hi-h;j<hi;++j)b[j]="ACGT"[rng()%4];
                } else {
                    std::rotate(b.begin()+lo,b.begin()+hi-h,b.begin()+hi);
                    for(unsigned j=lo;j<lo+h;++j)b[j]="ACGT"[rng()%4];
                }
            }
            q[lane]=Sequence(a);t[lane]=Sequence(b);qp[lane]=&q[lane].packed;tp[lane]=&t[lane].packed;
            exact[lane]=reference::api(q[lane],t[lane]);checksum+=exact[lane];
            if(id>=exhaustive || id%257==0) {
                if(scalar_dp(a,b)!=exact[lane])throw std::runtime_error("independent DP mismatch");
                ++dp_count;
            }
        }
        check<0>(qp,tp,exact,g);check<2>(qp,tp,exact,g);check<5>(qp,tp,exact,g);
        check<50>(qp,tp,exact,g);check<62>(qp,tp,exact,g);check<63>(qp,tp,exact,g);
        check<64>(qp,tp,exact,g);check<128>(qp,tp,exact,g);
    }
    std::cout<<"PASS banded_pairs="<<total<<" exhaustive_binary8_pairs="<<exhaustive
        <<" long_indel_pairs=32768 max_indel_excursion=64 thresholds=0,2,5,50,62,63,64,128"
        <<" unroll_variants=1,2,4,8 independent_DP_pairs="<<dp_count<<" checksum="<<checksum<<'\n';
}
