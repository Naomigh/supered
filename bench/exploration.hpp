#pragma once
#include "strong/packed.hpp"
#include <immintrin.h>

namespace strong::exploration {
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
inline __m256i block4(__m256i& pv, __m256i& mv, __m256i eq, __m256i hin) {
    const __m256i one = _mm256_set1_epi64x(1);
    __m256i neg = _mm256_srli_epi64(hin, 63);
    __m256i xv = _mm256_or_si256(eq, mv);
    eq = _mm256_or_si256(eq, neg);
    __m256i sum = _mm256_add_epi64(_mm256_and_si256(eq, pv), pv);
    __m256i xh = _mm256_ternarylogic_epi64(sum, pv, eq, 0xbe);
    __m256i ph = _mm256_ternarylogic_epi64(mv, xh, pv, 0xf1);
    __m256i mh = _mm256_and_si256(pv, xh);
    __m256i hout = _mm256_sub_epi64(_mm256_srli_epi64(ph, 63),
                                   _mm256_srli_epi64(mh, 63));
    ph = _mm256_or_si256(_mm256_slli_epi64(ph, 1),
                         _mm256_srli_epi64(_mm256_add_epi64(hin, one), 1));
    mh = _mm256_or_si256(_mm256_slli_epi64(mh, 1), neg);
    pv = _mm256_ternarylogic_epi64(mh, xv, ph, 0xf1);
    mv = _mm256_and_si256(ph, xv);
    return hout;
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

// Four independent pairs. AVX512VL: each YMM holds one 64-bit block of four pairs.
// No band/early rejection: always compute the exact distance, then apply K.
template<int K = 128>
inline Distances myers_2x64x4_signed(const Quad& q, const Quad& t) {
    using namespace detail;
    const __m256i zero = _mm256_setzero_si256(), one = _mm256_set1_epi64x(1);
    const __m256i ones = _mm256_set1_epi64x(-1);
    __m256i q0l = limb4(q, 0, 0), q1l = limb4(q, 1, 0);
    __m256i q0h = limb4(q, 0, 1), q1h = limb4(q, 1, 1);
    __m256i pl = ones, ph = ones, ml = zero, mh = zero;
    __m256i score = _mm256_set1_epi64x(128);
    for (int half = 0; half < 2; ++half) {
        __m256i t0 = limb4(t, 0, half), t1 = limb4(t, 1, half);
        for (int j = 0; j < 64; ++j) {
            __m256i b0 = _mm256_sub_epi64(zero, _mm256_and_si256(t0, one));
            __m256i b1 = _mm256_sub_epi64(zero, _mm256_and_si256(t1, one));
            __m256i e0 = _mm256_xor_si256(_mm256_or_si256(
                _mm256_xor_si256(q0l, b0), _mm256_xor_si256(q1l, b1)), ones);
            __m256i e1 = _mm256_xor_si256(_mm256_or_si256(
                _mm256_xor_si256(q0h, b0), _mm256_xor_si256(q1h, b1)), ones);
            __m256i hin = block4(pl, ml, e0, one);
            score = _mm256_add_epi64(score, block4(ph, mh, e1, hin));
            t0 = _mm256_srli_epi64(t0, 1);
            t1 = _mm256_srli_epi64(t1, 1);
        }
    }
    return finish4<K>(score);
}

namespace detail {
struct Edge4 { __m256i pos,neg; };
inline __m256i equality4(__m256i q0,__m256i q1,__m256i b0,__m256i b1) {
    return _mm256_ternarylogic_epi64(_mm256_xor_si256(q0,b0),q1,b1,0x09);
}
// Keep the two horizontal boundary bits instead of encoding/decoding a signed hin.
template<bool First>
inline Edge4 block_bits4(__m256i& pv,__m256i& mv,__m256i eq,Edge4 in) {
    __m256i xv=_mm256_or_si256(eq,mv);
    __m256i e=First?eq:_mm256_or_si256(eq,in.neg);
    __m256i sum=_mm256_add_epi64(_mm256_and_si256(e,pv),pv);
    __m256i xh=_mm256_ternarylogic_epi64(sum,pv,e,0xbe);
    __m256i ph=_mm256_ternarylogic_epi64(mv,xh,pv,0xf1);
    __m256i mh=_mm256_and_si256(pv,xh);
    Edge4 out{_mm256_srli_epi64(ph,63),_mm256_srli_epi64(mh,63)};
    ph=_mm256_or_si256(_mm256_slli_epi64(ph,1),in.pos);
    mh=First?_mm256_slli_epi64(mh,1):_mm256_or_si256(_mm256_slli_epi64(mh,1),in.neg);
    pv=_mm256_ternarylogic_epi64(mh,xv,ph,0xf1);
    mv=_mm256_and_si256(ph,xv);
    return out;
}
}

template<int K=128,bool Pipeline=false>
inline Distances myers_2x64x4_bits(const Quad& q,const Quad& t) {
    using namespace detail;
    const __m256i zero=_mm256_setzero_si256(),one=_mm256_set1_epi64x(1);
    __m256i q0l=limb4(q,0,0),q1l=limb4(q,1,0),q0h=limb4(q,0,1),q1h=limb4(q,1,1);
    __m256i pl=_mm256_set1_epi64x(-1),ph=pl,ml=zero,mh=zero,score=_mm256_set1_epi64x(128);
    Edge4 edge{one,zero}; __m256i previous=zero;
    for(int half=0;half<2;++half) {
        __m256i t0=limb4(t,0,half),t1=limb4(t,1,half);
        for(int j=0;j<64;++j) {
            __m256i b0=_mm256_sub_epi64(zero,_mm256_and_si256(t0,one));
            __m256i b1=_mm256_sub_epi64(zero,_mm256_and_si256(t1,one));
            __m256i el=equality4(q0l,q1l,b0,b1),eh=equality4(q0h,q1h,b0,b1);
            if constexpr(Pipeline) {
                // The previous column's high block and current low block are independent.
                // Column zero is peeled below by the compiler's constant first iteration.
                if(half!=0 || j!=0) {
                    auto out=block_bits4<false>(ph,mh,previous,edge);
                    score=_mm256_sub_epi64(_mm256_add_epi64(score,out.pos),out.neg);
                }
                edge=block_bits4<true>(pl,ml,el,{one,zero});
                previous=eh;
            } else {
                auto in=block_bits4<true>(pl,ml,el,{one,zero});
                auto out=block_bits4<false>(ph,mh,eh,in);
                score=_mm256_sub_epi64(_mm256_add_epi64(score,out.pos),out.neg);
            }
            t0=_mm256_srli_epi64(t0,1);t1=_mm256_srli_epi64(t1,1);
        }
    }
    if constexpr(Pipeline) {
        auto out=block_bits4<false>(ph,mh,previous,edge);
        score=_mm256_sub_epi64(_mm256_add_epi64(score,out.pos),out.neg);
    }
    return finish4<K>(score);
}

namespace detail {
inline __m256i reverse_bits4(__m256i v) {
    const __m256i bytes=_mm256_setr_epi8(7,6,5,4,3,2,1,0,15,14,13,12,11,10,9,8,
                                       7,6,5,4,3,2,1,0,15,14,13,12,11,10,9,8);
    v=_mm256_gf2p8affine_epi64_epi8(v,_mm256_set1_epi64x(0x8040201008040201ULL),0);
    return _mm256_shuffle_epi8(v,bytes);
}
}
template<int K=128,bool Reverse=false>
inline Distances myers_2x64x4_peeled(const Quad& q,const Quad& t) {
    using namespace detail;
    const __m256i zero=_mm256_setzero_si256(),one=_mm256_set1_epi64x(1);
    auto query=[&](int plane,int half) {
        if constexpr(Reverse) return reverse_bits4(limb4(q,plane,1-half));
        else return limb4(q,plane,half);
    };
    auto expand=[&](__m256i v) {
        if constexpr(Reverse) return _mm256_srai_epi64(v,63);
        else return _mm256_sub_epi64(zero,_mm256_and_si256(v,one));
    };
    auto advance=[](__m256i v) {
        if constexpr(Reverse) return _mm256_slli_epi64(v,1);
        else return _mm256_srli_epi64(v,1);
    };
    __m256i q0l=query(0,0),q1l=query(1,0),q0h=query(0,1),q1h=query(1,1);
    __m256i pl=_mm256_set1_epi64x(-1),ph=pl,ml=zero,mh=zero,score=_mm256_set1_epi64x(128);
    __m256i t0=limb4(t,0,Reverse?1:0),t1=limb4(t,1,Reverse?1:0);
    auto b0=expand(t0),b1=expand(t1);
    auto edge=block_bits4<true>(pl,ml,equality4(q0l,q1l,b0,b1),{one,zero});
    __m256i previous=equality4(q0h,q1h,b0,b1);
    t0=advance(t0);t1=advance(t1);
    auto step=[&] {
        auto x0=expand(t0),x1=expand(t1);
        auto out=block_bits4<false>(ph,mh,previous,edge);
        edge=block_bits4<true>(pl,ml,equality4(q0l,q1l,x0,x1),{one,zero});
        previous=equality4(q0h,q1h,x0,x1);
        score=_mm256_sub_epi64(_mm256_add_epi64(score,out.pos),out.neg);
        t0=advance(t0);t1=advance(t1);
    };
    for(int j=0;j<63;++j) step();
    t0=limb4(t,0,Reverse?0:1);t1=limb4(t,1,Reverse?0:1);
    for(int j=0;j<64;++j) step();
    auto out=block_bits4<false>(ph,mh,previous,edge);
    score=_mm256_sub_epi64(_mm256_add_epi64(score,out.pos),out.neg);
    return finish4<K>(score);
}

namespace detail {
// Raw Ph/Mh retain their boundary in bit 63, directly usable by VPSHLDQ.
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
    return myers_2x64x4_raw<K,true,4,true>(q,t);
}

// Two independent 64-bit block recurrences per pair, with high blocks one column
// behind low blocks. Layout: [lo0,lo1,lo2,lo3,hi0,hi1,hi2,hi3].
// This is not a 128-bit add/shift: every add and shift stays in its own qword.
template<int K=128>
inline Distances myers_2x64x4_skew(const Quad& q,const Quad& t) {
    using namespace detail;
    auto join=[](__m256i lo,__m256i hi) {
        return _mm512_inserti64x4(_mm512_castsi256_si512(lo),hi,1);
    };
    const __m256i z4=_mm256_setzero_si256(),o4=_mm256_set1_epi64x(1);
    const __m512i zero=_mm512_setzero_si512(),one=_mm512_set1_epi64(1);
    __m512i q0=join(limb4(q,0,0),limb4(q,0,1));
    __m512i q1=join(limb4(q,1,0),limb4(q,1,1));
    __m256i t0l=limb4(t,0,0),t1l=limb4(t,1,0);
    __m256i pl=_mm256_set1_epi64x(-1),ml=z4;
    auto e0=equality4(_mm512_castsi512_si256(q0),_mm512_castsi512_si256(q1),
        _mm256_sub_epi64(z4,_mm256_and_si256(t0l,o4)),
        _mm256_sub_epi64(z4,_mm256_and_si256(t1l,o4)));
    auto edge=block_bits4<true>(pl,ml,e0,{o4,z4}); // low block column zero
    __m512i pv=join(pl,_mm256_set1_epi64x(-1)),mv=join(ml,z4);
    __m512i pos=join(o4,edge.pos),neg=join(z4,edge.neg);
    __m512i score=_mm512_set1_epi64(128);
    // First segment: low columns 1..63, high columns 0..62.
    // Second segment: low columns 64..127 plus one unused, high columns 63..127.
    for(int half=0;half<2;++half) {
        __m256i a0=half?limb4(t,0,1):_mm256_srli_epi64(t0l,1);
        __m256i a1=half?limb4(t,1,1):_mm256_srli_epi64(t1l,1);
        __m256i b0=half?_mm256_or_si256(_mm256_slli_epi64(a0,1),_mm256_srli_epi64(t0l,63)):t0l;
        __m256i b1=half?_mm256_or_si256(_mm256_slli_epi64(a1,1),_mm256_srli_epi64(t1l,63)):t1l;
        __m512i t0=join(a0,b0),t1=join(a1,b1);
        for(int j=0;j<(half?65:63);++j) {
            // Low column 127's MSB becomes high column 127 in the final iteration.
            // The 65th high bit is injected separately below at the word boundary.
            if(half && j==64) {
                t0=_mm512_mask_set1_epi64(t0,0xf0,0);t1=_mm512_mask_set1_epi64(t1,0xf0,0);
                t0=_mm512_mask_mov_epi64(t0,0xf0,join(z4,_mm256_srli_epi64(a0,63)));
                t1=_mm512_mask_mov_epi64(t1,0xf0,join(z4,_mm256_srli_epi64(a1,63)));
            }
            __m512i x0=_mm512_sub_epi64(zero,_mm512_and_si512(t0,one));
            __m512i x1=_mm512_sub_epi64(zero,_mm512_and_si512(t1,one));
            __m512i eq=_mm512_ternarylogic_epi64(_mm512_xor_si512(q0,x0),q1,x1,0x09);
            __m512i xv=_mm512_or_si512(eq,mv);
            eq=_mm512_or_si512(eq,neg);
            __m512i sum=_mm512_add_epi64(_mm512_and_si512(eq,pv),pv);
            __m512i xh=_mm512_ternarylogic_epi64(sum,pv,eq,0xbe);
            __m512i ph=_mm512_ternarylogic_epi64(mv,xh,pv,0xf1);
            __m512i mh=_mm512_and_si512(pv,xh);
            __m512i op=_mm512_srli_epi64(ph,63),on=_mm512_srli_epi64(mh,63);
            score=_mm512_sub_epi64(_mm512_add_epi64(score,op),on);
            ph=_mm512_or_si512(_mm512_slli_epi64(ph,1),pos);
            mh=_mm512_or_si512(_mm512_slli_epi64(mh,1),neg);
            pv=_mm512_ternarylogic_epi64(mh,xv,ph,0xf1);mv=_mm512_and_si512(ph,xv);
            pos=_mm512_mask_mov_epi64(one,0xf0,_mm512_shuffle_i64x2(op,op,0x44));
            neg=_mm512_maskz_shuffle_i64x2(0xf0,on,on,0x44);
            t0=_mm512_srli_epi64(t0,1);t1=_mm512_srli_epi64(t1,1);
        }
    }
    return finish4<K>(_mm512_extracti64x4_epi64(score,1));
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
