#pragma once
#include "strong/banded.hpp"
#include "SIMD_ED.h"
#include <cstddef>

namespace supered {
using PackedDNA128 = strong::PackedDNA128;
using DNA = strong::DNA;
using Quad = strong::Quad;
using Distances = strong::Distances;
using strong::encode;
using strong::decode;
using strong::broadcast;

namespace detail {
template<int K, bool Small = (K < 5)>
struct Backend {};

template<int K>
struct Backend<K, true> {
    SIMD_ED leap;
    Backend() { leap.init_levenshtein(K, ED_GLOBAL, false); }

    int distance(const PackedDNA128& q, const PackedDNA128& t) {
        leap.load_reads(
            _mm256_zextsi128_si256(_mm_loadu_si128((const __m128i*)q.bit0)),
            _mm256_zextsi128_si256(_mm_loadu_si128((const __m128i*)q.bit1)),
            _mm256_zextsi128_si256(_mm_loadu_si128((const __m128i*)t.bit0)),
            _mm256_zextsi128_si256(_mm_loadu_si128((const __m128i*)t.bit1)), 128);
        leap.calculate_masks();
        leap.reset();
        leap.run();
        return leap.check_pass() ? leap.get_ED() : -1;
    }
};
}

// Contract: valid prepacked 128 bp DNA; K in [0,128]. Each worker owns its Context.
// Distance is exact global unit-cost Levenshtein if <= K, otherwise -1.
template<int K = 128>
class Context {
    detail::Backend<K> backend_;
public:
    static constexpr int max_edit = K;
    static constexpr bool uses_leap = K < 5;
    Context() = default;
    Context(const Context&) = delete;
    Context& operator=(const Context&) = delete;
    Context(Context&&) = delete;
    Context& operator=(Context&&) = delete;

    Distances distance4(const Quad& queries, const Quad& targets) {
        if constexpr(K < 5) {
            Distances distances;
            for(int lane = 0; lane < 4; ++lane)
                distances[lane] = backend_.distance(*queries[lane], *targets[lane]);
            return distances;
        } else {
            return strong::myers_banded64x4<K>(queries, targets);
        }
    }

    Distances distance4(const PackedDNA128& query, const Quad& targets) {
        return distance4(broadcast(query), targets);
    }

    // Convenience for an isolated pair. Myers broadcasts to four lanes;
    // use distance4/distance_many/distance_pairs for measured batch throughput.
    int distance(const PackedDNA128& query, const PackedDNA128& target) {
        if constexpr(K < 5) return backend_.distance(query, target);
        else return distance4(broadcast(query), broadcast(target))[0];
    }

    // count is a multiple of four. Output storage is supplied by the caller.
    void distance_many(const PackedDNA128& query, const PackedDNA128* targets,
                       std::size_t count, int* output) {
        for(std::size_t i = 0; i < count; i += 4) {
            auto distances = distance4(query,
                Quad{targets+i, targets+i+1, targets+i+2, targets+i+3});
            for(int lane = 0; lane < 4; ++lane) output[i+lane] = distances[lane];
        }
    }

    void distance_pairs(const PackedDNA128* queries, const PackedDNA128* targets,
                        std::size_t count, int* output) {
        for(std::size_t i = 0; i < count; i += 4) {
            auto distances = distance4(
                Quad{queries+i, queries+i+1, queries+i+2, queries+i+3},
                Quad{targets+i, targets+i+1, targets+i+2, targets+i+3});
            for(int lane = 0; lane < 4; ++lane) output[i+lane] = distances[lane];
        }
    }
};
}
