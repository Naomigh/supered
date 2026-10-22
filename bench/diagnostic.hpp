#pragma once
#include "strong/myers.hpp"
namespace strong::diagnostic {
inline __m512i add_vector(__m512i a,__m512i b) {
    auto sum=_mm512_add_epi64(a,b);
    auto carry=_mm512_ternarylogic_epi64(a,b,sum,0xd4);
    carry=_mm512_srli_epi64(carry,63);
    auto high=_mm512_maskz_shuffle_epi32(0xcccc,carry,_MM_PERM_BABA);
    return _mm512_add_epi64(sum,high);
}
inline __m512i shift_local(__m512i v) {
    auto neighbor=_mm512_maskz_shuffle_epi32(0xcccc,v,_MM_PERM_BABA);
    return _mm512_shldi_epi64(v,neighbor,1);
}
// ZMM qwords [pair0.lo, pair0.hi, ..., pair3.lo, pair3.hi].
template<int K = 128,bool VectorCarry=false,bool LocalShift=false>
inline Distances myers_128x4_diagnostic(const Quad& q, const Quad& t) {
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
            __m512i sum = VectorCarry?add_vector(_mm512_and_si512(eq, pv),pv):add128x4(_mm512_and_si512(eq,pv),pv);
            __m512i xh = _mm512_ternarylogic_epi64(sum, pv, eq, 0xbe);
            __m512i ph = _mm512_ternarylogic_epi64(mv, xh, pv, 0xf1);
            __m512i mh = _mm512_and_si512(pv, xh);
            score = _mm512_add_epi64(score, _mm512_srli_epi64(ph, 63));
            score = _mm512_sub_epi64(score, _mm512_srli_epi64(mh, 63));
            ph = _mm512_or_si512((LocalShift?shift_local(ph):shl128x4(ph)), first);
            mh = (LocalShift?shift_local(mh):shl128x4(mh));
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
