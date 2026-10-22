#pragma once
#include "myers.hpp"
namespace strong {
// Equal-length global alignment: a path drifting d positions must pay >=2*d
// indels before returning to (128,128). For K<64, width 64 around the diagonal
// contains every qualifying path. Larger K uses the exact dense kernel.
template<int K,int Unroll=4>
inline Distances myers_banded64x4(const Quad& q,const Quad& t) {
    if constexpr(K>=64) return myers_2x64x4<K>(q,t);
    else {
        using namespace detail;
        const __m256i zero=_mm256_setzero_si256(),one=_mm256_set1_epi64x(1);
        const __m256i low=_mm256_set1_epi64x(0xffffffffULL);
        auto qp=transpose4(q),tp=transpose4(t);
        // Window origin is query row -32; artificial negative rows have cost |i|.
        __m256i pv=_mm256_set1_epi64x(0xffffffff00000000ULL),mv=low;
        __m256i q0=_mm256_slli_epi64(qp.q0l,32),q1=_mm256_slli_epi64(qp.q1l,32);
        __m256i n0=_mm256_shldi_epi64(qp.q0h,qp.q0l,32);
        __m256i n1=_mm256_shldi_epi64(qp.q1h,qp.q1l,32);
        __m256i drops=_mm256_setzero_si256();
        for(int half=0;half<2;++half) {
            auto t0=half?tp.q0h:tp.q0l,t1=half?tp.q1h:tp.q1l;
            __m256i hits=zero;
            for(int j=0;j<64/Unroll;++j) for(int u=0;u<Unroll;++u) {
                auto b0=_mm256_sub_epi64(zero,_mm256_and_si256(t0,one));
                auto b1=_mm256_sub_epi64(zero,_mm256_and_si256(t1,one));
                auto eq=equality4(q0,q1,b0,b1);
                auto xv=_mm256_or_si256(eq,mv);
                auto sum=_mm256_add_epi64(_mm256_and_si256(eq,pv),pv);
                auto xh=_mm256_ternarylogic_epi64(sum,pv,eq,0xbe);
                auto ph=_mm256_ternarylogic_epi64(mv,xh,pv,0xf1);
                auto mh=_mm256_and_si256(pv,xh);
                // Moving origin advances one row and one column. Its score
                // increases by 1-(Xv&1); collect the subtracted bits in-register.
                hits=_mm256_shrdi_epi64(hits,xv,1);
                xv=_mm256_srli_epi64(xv,1);
                // PV[63]=1 and MV[63]=0 imply Ph[63]=0. After the window
                // shift Xv[63]=0, so the new top PV bit is already one.
                pv=_mm256_ternarylogic_epi64(mh,xv,ph,0xf1);
                mv=_mm256_and_si256(ph,xv);
                q0=_mm256_shrdi_epi64(q0,n0,1);q1=_mm256_shrdi_epi64(q1,n1,1);
                n0=_mm256_srli_epi64(n0,1);n1=_mm256_srli_epi64(n1,1);
                t0=_mm256_srli_epi64(t0,1);t1=_mm256_srli_epi64(t1,1);
            }
            drops=_mm256_add_epi64(drops,_mm256_popcnt_epi64(hits));
            n0=_mm256_srli_epi64(qp.q0h,32);n1=_mm256_srli_epi64(qp.q1h,32);
        }
        // Final origin is row 96; the first 32 differences lead to row 128.
        auto score=_mm256_sub_epi64(_mm256_set1_epi64x(160),drops);
        score=_mm256_add_epi64(score,_mm256_popcnt_epi64(_mm256_and_si256(pv,low)));
        score=_mm256_sub_epi64(score,_mm256_popcnt_epi64(_mm256_and_si256(mv,low)));
        return finish4<K>(score);
    }
}
}
