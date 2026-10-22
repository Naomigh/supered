#pragma once
#include "../tests/support.hpp"
#include "../edlib_origin/edlib/src/edlib.cpp"
#include "../generated/original_nw_preallocated.hpp"
#include "SIMD_ED.h"

namespace reference {
using test_support::Sequence;
inline int api(const Sequence& q, const Sequence& t, int k = -1) {
    auto r = edlibAlign(q.dna.data(), 128, t.dna.data(), 128,
                       edlibNewAlignConfig(k, EDLIB_MODE_NW, EDLIB_TASK_DISTANCE, nullptr, 0));
    int d = r.editDistance;
    edlibFreeAlignResult(r);
    return d;
}
template<int K, bool Preallocated>
inline int hot(const Sequence& q, const Sequence& t) {
    int d, pos;
    AlignmentData* alignment = nullptr;
    if constexpr (Preallocated)
        original_nw_preallocated(q.peq, 0, 2, 128, t.codes.data(), 128,
                                 K, &d, &pos, false, &alignment, -1);
    else
        myersCalcEditDistanceNW(q.peq, 0, 2, 128, t.codes.data(), 128,
                                K, &d, &pos, false, &alignment, -1);
    return d;
}
// Unmodified corrected LEAP AVX2 backend. Zero-extension adapts exactly 16 bytes
// per plane to its existing 256-bit interface; no DNA re-encoding is done here.
inline int leap(SIMD_ED& engine, const strong::PackedDNA128& q,
                const strong::PackedDNA128& t) {
    engine.load_reads(
        _mm256_zextsi128_si256(_mm_loadu_si128((const __m128i*)q.bit0)),
        _mm256_zextsi128_si256(_mm_loadu_si128((const __m128i*)q.bit1)),
        _mm256_zextsi128_si256(_mm_loadu_si128((const __m128i*)t.bit0)),
        _mm256_zextsi128_si256(_mm_loadu_si128((const __m128i*)t.bit1)), 128);
    engine.calculate_masks();
    engine.reset();
    engine.run();
    return engine.check_pass() ? engine.get_ED() : -1;
}
}
