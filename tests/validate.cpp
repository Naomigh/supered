#include "../src/references.hpp"
#include "strong/banded.hpp"
#include <fstream>
#include <iostream>
#include <map>
#include <stdexcept>

using namespace strong;
using namespace test_support;

template<int K>
void check_group(const Sequence* q, const Sequence* t, const std::array<int,4>& exact,
                 SIMD_ED& leap, uint64_t group) {
    Quad qp{&q[0].packed,&q[1].packed,&q[2].packed,&q[3].packed};
    Quad tp{&t[0].packed,&t[1].packed,&t[2].packed,&t[3].packed};
    auto a = myers_2x64x4<K>(qp,tp), b = myers_128x4<K>(qp,tp);
    auto c=myers_banded64x4<K>(qp,tp);
    for (int lane = 0; lane < 4; ++lane) {
        int expected = exact[lane] <= K ? exact[lane] : -1;
        int original = reference::api(q[lane], t[lane], K);
        int hot = reference::hot<K,false>(q[lane],t[lane]);
        int stack = reference::hot<K,true>(q[lane],t[lane]);
        int l = K <= 5 ? reference::leap(leap,q[lane].packed,t[lane].packed) : expected;
        if (a[lane] != expected || b[lane] != expected || c[lane]!=expected || original != expected
            || hot != expected || stack != expected || l != expected) {
            std::cerr << "FAIL group="<<group<<" lane="<<lane<<" K="<<K
                <<" exact="<<exact[lane]<<" original="<<original<<" hot="<<hot
                <<" stack="<<stack<<" two64="<<a[lane]<<" wide="<<b[lane]<<" LEAP="<<l<<'\n'
                <<std::string(q[lane].dna.begin(),q[lane].dna.end())<<'\n'
                <<std::string(t[lane].dna.begin(),t[lane].dna.end())<<'\n';
            throw std::runtime_error("distance mismatch");
        }
    }
}
int main(int argc, char** argv) {
    const uint64_t groups = argc > 1 ? std::stoull(argv[1]) / 4 : 32768;
    std::mt19937_64 rng(0x1ed11b512);
    for(int iteration=0;iteration<10000;++iteration) {
        alignas(64) uint64_t a[8],b[8],sum[8],shifted[8];
        for(int j=0;j<8;++j) { a[j]=rng();b[j]=rng(); }
        _mm512_store_si512(sum,detail::add128x4(_mm512_load_si512(a),_mm512_load_si512(b)));
        _mm512_store_si512(shifted,detail::shl128x4(_mm512_load_si512(a)));
        for(int j=0;j<8;j+=2) {
            unsigned __int128 x=(static_cast<unsigned __int128>(a[j+1])<<64)|a[j];
            unsigned __int128 y=(static_cast<unsigned __int128>(b[j+1])<<64)|b[j];
            auto z=x+y;
            if(sum[j]!=uint64_t(z)||sum[j+1]!=uint64_t(z>>64)
                ||shifted[j]!=uint64_t(x<<1)||shifted[j+1]!=uint64_t((x<<1)>>64))
                throw std::runtime_error("128-bit primitive mismatch");
        }
    }
    SIMD_ED l0,l2,l5,unused;
    l0.init_levenshtein(0,ED_GLOBAL,false);
    l2.init_levenshtein(2,ED_GLOBAL,false);
    l5.init_levenshtein(5,ED_GLOBAL,false);
    const int budgets[] = {0,1,2,3,4,5,6,10,49,50,51,80,128};
    std::map<int,uint64_t> histogram;
    uint64_t checksum = 0, dp_checked = 0;
    for (uint64_t g = 0; g < groups; ++g) {
        Sequence q[4],t[4];
        std::array<int,4> exact;
        for (int lane=0;lane<4;++lane) {
            uint64_t id=4*g+lane;
            DNA a=random_dna(rng),b;
            if (id < 16) {
                a.fill("ACGT"[id/4]);b.fill("ACGT"[id%4]);
            } else if (id < 512) {
                a.fill("ACGT"[(id/128)%4]); b=a;
                b[id%128]="ACGT"[(base_code(b[id%128])+1)%4];
            } else if (id < 1024) {
                for (int p=0;p<128;++p) a[p]="ACGT"[(p+id)%4];
                b=mutate(a,2+(id%5),rng,true);
            } else if (id%7 == 0) b=random_dna(rng);
            else b=mutate(a,budgets[id%13],rng,id%2);
            q[lane]=Sequence(a);t[lane]=Sequence(b);
            if (decode(q[lane].packed)!=a || decode(t[lane].packed)!=b)
                throw std::runtime_error("encoding mismatch");
            exact[lane]=reference::api(q[lane],t[lane]);
            if (id < 2048 || id%257 == 0) {
                if (scalar_dp(a,b)!=exact[lane]) throw std::runtime_error("independent DP mismatch");
                ++dp_checked;
            }
            ++histogram[exact[lane]]; checksum+=exact[lane];
        }
        check_group<0>(q,t,exact,l0,g);
        check_group<2>(q,t,exact,l2,g);
        check_group<5>(q,t,exact,l5,g);
        check_group<50>(q,t,exact,unused,g);
        check_group<128>(q,t,exact,unused,g);
    }
    std::cout<<"PASS pairs="<<4*groups<<" thresholds=0,2,5,50,128 oracle=original_edlib_API"
        <<" independent_DP_pairs="<<dp_checked<<" primitive_pair_tests=40000 checksum="<<checksum<<'\n';
    std::cout<<"distance,count\n";
    for (auto [d,n]:histogram) std::cout<<d<<','<<n<<'\n';
}
