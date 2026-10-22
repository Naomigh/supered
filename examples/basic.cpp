#include "supered/supered.hpp"
#include <iostream>

int main() {
    supered::DNA dna;
    dna.fill('A');
    auto query = supered::encode(dna);
    supered::PackedDNA128 targets[4];
    for(int i = 0; i < 4; ++i) {
        auto target = dna;
        for(int j = 0; j < i; ++j) target[10 + 30*j] = 'C';
        targets[i] = supered::encode(target);
    }
    // Encoding and Context construction are outside the compute region.
    supered::Context<2> tiny;       // corrected LEAP, threshold 2
    supered::Context<5> bounded;    // four-pair moving-band Myers
    supered::Context<> complete;   // K=128, complete two-block Myers
    int distances[4];
    tiny.distance_many(query, targets, 4, distances);
    for(int d : distances) std::cout << d << ' ';
    std::cout << '\n';
    bounded.distance_many(query, targets, 4, distances);
    for(int d : distances) std::cout << d << ' ';
    std::cout << '\n';
    complete.distance_many(query, targets, 4, distances);
    for(int d : distances) std::cout << d << ' ';
    std::cout << '\n';
}
