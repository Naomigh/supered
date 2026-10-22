#include "../src/references.hpp"
#include "strong/banded.hpp"
#include <iostream>
#include <stdexcept>

// Target the positive K=2 indel path, with deletion/insertion across word 63/64.
// This independent program is not linked into the benchmark or production core.
int main() {
    using namespace strong;
    using namespace test_support;
    std::mt19937_64 rng(0x128b0172);
    SIMD_ED leap;
    leap.init_levenshtein(2,ED_GLOBAL,false);
    uint64_t count=0,checksum=0;
    for(int lo=0;lo<64;++lo) for(int hi=64;hi<128;++hi) {
        Sequence q[4],t[4];
        Quad qp,tp;
        for(int lane=0;lane<4;++lane) {
            auto a=random_dna(rng),b=a;
            if(lane%2) {
                std::rotate(b.begin()+lo,b.begin()+hi,b.begin()+hi+1);
                b[lo]="ACGT"[rng()%4];
            } else {
                std::rotate(b.begin()+lo,b.begin()+lo+1,b.begin()+hi+1);
                b[hi]="ACGT"[rng()%4];
            }
            q[lane]=Sequence(a);t[lane]=Sequence(b);
            qp[lane]=&q[lane].packed;tp[lane]=&t[lane].packed;
        }
        auto a=myers_banded64x4<2>(qp,tp),b=myers_128x4<2>(qp,tp);
        for(int lane=0;lane<4;++lane) {
            int d=reference::api(q[lane],t[lane]);
            int dp=scalar_dp(q[lane].dna,t[lane].dna);
            int l=reference::leap(leap,q[lane].packed,t[lane].packed);
            if(d>2||d!=dp||a[lane]!=d||b[lane]!=d||l!=d)
                throw std::runtime_error("positive K=2 cross-word indel mismatch");
            ++count;checksum+=d;
        }
    }
    std::cout<<"PASS targeted_indel_pairs="<<count
        <<" K=2 independent_DP_pairs="<<count<<" checksum="<<checksum<<'\n';
}
