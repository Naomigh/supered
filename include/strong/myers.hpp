#pragma once
#include "packed.hpp"
#include <immintrin.h>

namespace strong {
namespace detail {
inline __m256i limb4(const Quad& s, int plane, int half) {
    return _mm256_setr_epi64x(
        (plane ? s[0]->bit1 : s[0]->bit0)[half],
        (plane ? s[1]->bit1 : s[1]->bit0)[half],
        (plane ? s[2]->bit1 : s[2]->bit0)[half],
        (plane ? s[3]->bit1 : s[3]->bit0)[half]);
}
inline __m512i plane4(const Quad& s, int plane) {
    auto ptr = [plane](const PackedDNA128* p) { return plane ? p->bit1 : p->bit0; };
    __m512i v = _mm512_castsi128_si512(_mm_loadu_si128((const __m128i*)ptr(s[0])));
    v = _mm512_inserti64x2(v, _mm_loadu_si128((const __m128i*)ptr(s[1])), 1);
    v = _mm512_inserti64x2(v, _mm_loadu_si128((const __m128i*)ptr(s[2])), 2);
    return _mm512_inserti64x2(v, _mm_loadu_si128((const __m128i*)ptr(s[3])), 3);
}
inline __m512i add128x4(__m512i a, __m512i b) {
    __m512i sum = _mm512_add_epi64(a, b);
    __mmask8 carry = _mm512_mask_cmp_epu64_mask(0x55, sum, a, _MM_CMPINT_LT);
    return _mm512_mask_add_epi64(sum, _kshiftli_mask8(carry, 1), sum,
                                _mm512_set1_epi64(1));
}
inline __m512i shl128x4(__m512i v) {
    __m512i neighbor = _mm512_maskz_permutexvar_epi64(
        0xaa, _mm512_setr_epi64(0, 0, 2, 2, 4, 4, 6, 6), v);
    return _mm512_shldi_epi64(v, neighbor, 1);
}
template<int K> inline Distances finish4(__m256i score) {
    if constexpr (K < 128) {
        __mmask8 over = _mm256_cmp_epi64_mask(score, _mm256_set1_epi64x(K), _MM_CMPINT_GT);
        score = _mm256_mask_mov_epi64(score, over, _mm256_set1_epi64x(-1));
    }
    Distances d;
    _mm_storeu_si128((__m128i*)d.data(), _mm256_cvtepi64_epi32(score));
    return d;
}
}

namespace detail {
struct Edge4 {__m256i pos,neg;};
inline __m256i equality4(__m256i q0,__m256i q1,__m256i b0,__m256i b1) {
    return _mm256_ternarylogic_epi64(_mm256_xor_si256(q0,b0),q1,b1,0x09);
}
struct Planes4 { __m256i q0l,q0h,q1l,q1h; };
inline Planes4 transpose4(const Quad& p) {
    auto a=_mm256_loadu_si256((const __m256i*)p[0]);
    auto b=_mm256_loadu_si256((const __m256i*)p[1]);
    auto c=_mm256_loadu_si256((const __m256i*)p[2]);
    auto d=_mm256_loadu_si256((const __m256i*)p[3]);
    auto ab0=_mm256_unpacklo_epi64(a,b),ab1=_mm256_unpackhi_epi64(a,b);
    auto cd0=_mm256_unpacklo_epi64(c,d),cd1=_mm256_unpackhi_epi64(c,d);
    return {_mm256_permute2x128_si256(ab0,cd0,0x20),_mm256_permute2x128_si256(ab1,cd1,0x20),
            _mm256_permute2x128_si256(ab0,cd0,0x31),_mm256_permute2x128_si256(ab1,cd1,0x31)};
}
template<bool First>
inline Edge4 block_raw4(__m256i& pv,__m256i& mv,__m256i eq,Edge4 in) {
    __m256i xv=_mm256_or_si256(eq,mv);
    __m256i e=First?eq:_mm256_or_si256(eq,_mm256_srli_epi64(in.neg,63));
    __m256i sum=_mm256_add_epi64(_mm256_and_si256(e,pv),pv);
    __m256i xh=_mm256_ternarylogic_epi64(sum,pv,e,0xbe);
    __m256i ph=_mm256_ternarylogic_epi64(mv,xh,pv,0xf1),mh=_mm256_and_si256(pv,xh);
    Edge4 out{ph,mh};
    if constexpr(First) {
        ph=_mm256_or_si256(_mm256_slli_epi64(ph,1),_mm256_set1_epi64x(1));
        mh=_mm256_slli_epi64(mh,1);
    } else {
        ph=_mm256_shldi_epi64(ph,in.pos,1);
        mh=_mm256_shldi_epi64(mh,in.neg,1);
    }
    pv=_mm256_ternarylogic_epi64(mh,xv,ph,0xf1);mv=_mm256_and_si256(ph,xv);
    return out;
}
}
template<int K=128,bool CountBits=false,int Unroll=1,bool FinalScore=false,bool Transpose=false>
inline Distances myers_2x64x4_raw(const Quad& q,const Quad& t) {
    using namespace detail;
    const __m256i zero=_mm256_setzero_si256(),one=_mm256_set1_epi64x(1);
    Planes4 qp,tp;
    if constexpr(Transpose) {qp=transpose4(q);tp=transpose4(t);}
    __m256i q0l=Transpose?qp.q0l:limb4(q,0,0),q1l=Transpose?qp.q1l:limb4(q,1,0);
    __m256i q0h=Transpose?qp.q0h:limb4(q,0,1),q1h=Transpose?qp.q1h:limb4(q,1,1);
    __m256i pl=_mm256_set1_epi64x(-1),ph=pl,ml=zero,mh=zero,score=_mm256_set1_epi64x(128);
    __m256i t0=Transpose?tp.q0l:limb4(t,0,0),t1=Transpose?tp.q1l:limb4(t,1,0),positive=zero,negative=zero;
    auto expand=[&](__m256i v) {return _mm256_sub_epi64(zero,_mm256_and_si256(v,one));};
    auto b0=expand(t0),b1=expand(t1);
    Edge4 edge=block_raw4<true>(pl,ml,equality4(q0l,q1l,b0,b1),{zero,zero});
    __m256i previous=equality4(q0h,q1h,b0,b1);
    t0=_mm256_srli_epi64(t0,1);t1=_mm256_srli_epi64(t1,1);
    auto count=[&](Edge4 out) {
        if constexpr(FinalScore) {
            // PV/MV encode vertical +/-1 differences; recover the final column
            // score once after the complete global recurrence.
        } else if constexpr(CountBits) {
            positive=_mm256_shldi_epi64(positive,out.pos,1);
            negative=_mm256_shldi_epi64(negative,out.neg,1);
        } else {
            score=_mm256_sub_epi64(_mm256_add_epi64(score,_mm256_srli_epi64(out.pos,63)),_mm256_srli_epi64(out.neg,63));
        }
    };
    auto fold=[&] {
        if constexpr(CountBits && !FinalScore) {
            score=_mm256_sub_epi64(_mm256_add_epi64(score,_mm256_popcnt_epi64(positive)),_mm256_popcnt_epi64(negative));
            positive=zero;negative=zero;
        }
    };
    auto step=[&] {
        auto x0=expand(t0),x1=expand(t1);
        auto out=block_raw4<false>(ph,mh,previous,edge);
        edge=block_raw4<true>(pl,ml,equality4(q0l,q1l,x0,x1),{zero,zero});
        previous=equality4(q0h,q1h,x0,x1);
        count(out);
        t0=_mm256_srli_epi64(t0,1);t1=_mm256_srli_epi64(t1,1);
    };
    for(int j=0;j<63/Unroll;++j) {
        for(int u=0;u<Unroll;++u) step();
    }
    for(int j=0;j<63%Unroll;++j) step();
    fold();
    t0=Transpose?tp.q0h:limb4(t,0,1);t1=Transpose?tp.q1h:limb4(t,1,1);
    for(int j=0;j<64/Unroll;++j) {
        for(int u=0;u<Unroll;++u) step();
    }
    for(int j=0;j<64%Unroll;++j) step();
    fold();
    auto out=block_raw4<false>(ph,mh,previous,edge);
    count(out);fold();
    if constexpr(FinalScore) {
        score=_mm256_add_epi64(score,_mm256_add_epi64(_mm256_popcnt_epi64(pl),_mm256_popcnt_epi64(ph)));
        score=_mm256_sub_epi64(score,_mm256_add_epi64(_mm256_popcnt_epi64(ml),_mm256_popcnt_epi64(mh)));
    }
    return finish4<K>(score);
}
template<int K=128>
inline Distances myers_2x64x4(const Quad& q,const Quad& t) {
    return myers_2x64x4_raw<K,true,2,true,true>(q,t);
}

// ZMM qwords [pair0.lo, pair0.hi, ..., pair3.lo, pair3.hi].
template<int K = 128>
inline Distances myers_128x4(const Quad& q, const Quad& t) {
    using namespace detail;
    const __m512i zero = _mm512_setzero_si512(), one = _mm512_set1_epi64(1);
    const __m512i ones = _mm512_set1_epi64(-1);
    const __m512i first = _mm512_setr_epi64(1, 0, 1, 0, 1, 0, 1, 0);
    __m512i q0 = plane4(q, 0), q1 = plane4(q, 1);
    __m512i tp0 = plane4(t, 0), tp1 = plane4(t, 1);
    __m512i pv = ones, mv = zero, score = _mm512_set1_epi64(128);
    for (int half = 0; half < 2; ++half) {
        __m512i idx = _mm512_setr_epi64(half, half, half+2, half+2,
                                       half+4, half+4, half+6, half+6);
        __m512i t0 = _mm512_permutexvar_epi64(idx, tp0);
        __m512i t1 = _mm512_permutexvar_epi64(idx, tp1);
        for (int j = 0; j < 64; ++j) {
            __m512i b0 = _mm512_sub_epi64(zero, _mm512_and_si512(t0, one));
            __m512i b1 = _mm512_sub_epi64(zero, _mm512_and_si512(t1, one));
            __m512i eq = _mm512_xor_si512(_mm512_or_si512(
                _mm512_xor_si512(q0, b0), _mm512_xor_si512(q1, b1)), ones);
            __m512i xv = _mm512_or_si512(eq, mv);
            __m512i sum = add128x4(_mm512_and_si512(eq, pv), pv);
            __m512i xh = _mm512_ternarylogic_epi64(sum, pv, eq, 0xbe);
            __m512i ph = _mm512_ternarylogic_epi64(mv, xh, pv, 0xf1);
            __m512i mh = _mm512_and_si512(pv, xh);
            score = _mm512_add_epi64(score, _mm512_srli_epi64(ph, 63));
            score = _mm512_sub_epi64(score, _mm512_srli_epi64(mh, 63));
            ph = _mm512_or_si512(shl128x4(ph), first);
            mh = shl128x4(mh);
            pv = _mm512_ternarylogic_epi64(mh, xv, ph, 0xf1);
            mv = _mm512_and_si512(ph, xv);
            t0 = _mm512_srli_epi64(t0, 1);
            t1 = _mm512_srli_epi64(t1, 1);
        }
    }
    __m512i compact = _mm512_maskz_compress_epi64(0xaa, score);
    return finish4<K>(_mm512_castsi512_si256(compact));
}
}
