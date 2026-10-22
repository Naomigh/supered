#pragma once
#include "strong/packed.hpp"
#include "strong/myers.hpp"
#include <algorithm>
#include <array>
#include <random>
#include <string>
#include <vector>

namespace test_support {
using namespace strong;
struct alignas(64) Sequence {
    PackedDNA128 packed;
    std::array<unsigned char, 128> codes;
    DNA dna;
    uint64_t peq[8];
    explicit Sequence(const DNA& s = {}) : packed(encode(s)), dna(s) {
        for (int i = 0; i < 128; ++i) codes[i] = base_code(s[i]);
        for (int h = 0; h < 2; ++h) {
            uint64_t a = packed.bit0[h], b = packed.bit1[h];
            peq[0*2+h] = ~a & ~b;
            peq[1*2+h] =  a & ~b;
            peq[2*2+h] = ~a &  b;
            peq[3*2+h] =  a &  b;
        }
    }
};
inline DNA random_dna(std::mt19937_64& rng) {
    DNA s;
    for (auto& c : s) c = "ACGT"[rng() & 3];
    return s;
}
// A balanced deletion/insertion costs at most two, and preserves 128 bp.
inline DNA mutate(DNA s, int budget, std::mt19937_64& rng, bool indels) {
    int pairs = indels ? budget / 4 : 0;
    for (int i = 0; i < pairs; ++i) {
        int a = rng() % 127, b = a + 1 + rng() % (127 - a);
        std::rotate(s.begin() + a, s.begin() + a + 1, s.begin() + b + 1);
        s[b] = "ACGT"[rng() & 3];
    }
    std::array<int, 128> positions;
    for (int i = 0; i < 128; ++i) positions[i] = i;
    std::shuffle(positions.begin(), positions.end(), rng);
    for (int i = 0; i < budget - 2*pairs; ++i) {
        int p = positions[i];
        s[p] = "ACGT"[(base_code(s[p]) + 1 + rng()%3) & 3];
    }
    return s;
}
inline int scalar_dp(const DNA& a, const DNA& b) {
    int prev[129], next[129];
    for (int j = 0; j <= 128; ++j) prev[j] = j;
    for (int i = 1; i <= 128; ++i) {
        next[0] = i;
        for (int j = 1; j <= 128; ++j)
            next[j] = std::min({prev[j] + 1, next[j-1] + 1,
                               prev[j-1] + (a[i-1] != b[j-1])});
        std::copy(next, next+129, prev);
    }
    return prev[128];
}
}
