#pragma once
#include "packed.hpp"
#include <immintrin.h>

namespace strong::v1 {
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
inline Distances myers_2x64x4(const Quad& q, const Quad& t) {
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
