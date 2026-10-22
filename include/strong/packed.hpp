#pragma once
#include <array>
#include <cstdint>

namespace strong {
// A=00, C=01, G=10, T=11; position i is bit i, little-endian words.
struct alignas(32) PackedDNA128 {
    uint64_t bit0[2];
    uint64_t bit1[2];
};
static_assert(sizeof(PackedDNA128) == 32);
using DNA = std::array<char, 128>;
using Quad = std::array<const PackedDNA128*, 4>;
using Distances = std::array<int, 4>;

inline unsigned base_code(char c) {
    return c == 'C' ? 1 : c == 'G' ? 2 : c == 'T' ? 3 : 0;
}
inline PackedDNA128 encode(const DNA& dna) {
    PackedDNA128 p{};
    for (unsigned i = 0; i < 128; ++i) {
        unsigned b = base_code(dna[i]);
        p.bit0[i / 64] |= uint64_t(b & 1) << (i % 64);
        p.bit1[i / 64] |= uint64_t(b >> 1) << (i % 64);
    }
    return p;
}
inline DNA decode(const PackedDNA128& p) {
    DNA dna;
    for (unsigned i = 0; i < 128; ++i)
        dna[i] = "ACGT"[((p.bit0[i / 64] >> (i % 64)) & 1)
                       | (((p.bit1[i / 64] >> (i % 64)) & 1) << 1)];
    return dna;
}
inline Quad broadcast(const PackedDNA128& p) { return {&p, &p, &p, &p}; }
}
