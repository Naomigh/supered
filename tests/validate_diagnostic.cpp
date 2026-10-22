#include "../src/references.hpp"
#include "../bench/diagnostic.hpp"
#include <iostream>
int main(){using namespace strong;using namespace test_support;std::mt19937_64 rng(512);uint64_t count=0;
for(int g=0;g<4096;g++){Sequence q[4],t[4];Quad qp,tp;for(int l=0;l<4;l++){q[l]=Sequence(random_dna(rng));t[l]=Sequence(mutate(q[l].dna,(g+l)%129,rng,g%2));qp[l]=&q[l].packed;tp[l]=&t[l].packed;}
auto a=diagnostic::myers_128x4_diagnostic<128,true,false>(qp,tp),b=diagnostic::myers_128x4_diagnostic<128,false,true>(qp,tp),c=diagnostic::myers_128x4_diagnostic<128,true,true>(qp,tp);
for(int l=0;l<4;l++){int d=reference::api(q[l],t[l]);if(a[l]!=d||b[l]!=d||c[l]!=d)throw std::runtime_error("diagnostic mismatch");++count;}}
std::cout<<"PASS diagnostic_variants=3 pairs_each="<<count<<" oracle=original_Edlib_API\n";}
