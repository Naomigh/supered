# supered

**Fast, exact edit distance for 128 bp DNA.**

**English** | [简体中文](README.zh-CN.md)

supered is a C++ library for global, unit-cost Levenshtein distance on prepacked DNA sequences. It combines four-pair AVX-512 Myers kernels with a corrected LEAP backend for small edit thresholds, using the same **32-byte input format** throughout.

```cpp
// Inputs are already packed; reuse the context across calls.
supered::Context<80> engine;
auto distances = engine.distance4(
    read, {&anchor0, &anchor1, &anchor2, &anchor3});
// distances[i]: exact distance <= 80, or -1 if the distance exceeds 80.
```

In a cache-conditioned, single-thread experiment with **10 million reads × 20 anchors**, supered achieved **97.91 ns/pair**, or **10.21 million pairs/s**. See [Performance](#performance) for the workload, measurement boundaries, and comparison with original Edlib.

## Features

- **Exact distances:** global alignment with insertion, deletion, and substitution cost 1; matching cost 0.
- **Compact inputs:** two 128-bit planes, totaling 32 bytes per 128 bp sequence, shared by both backends.
- **Four-pair SIMD:** compute four independent distances per Myers call, including one read against four anchors.
- **Threshold specialization:** compile-time `Context<K>` selects corrected LEAP for `K < 5`, moving-band Myers for `5 <= K < 64`, and complete two-block Myers for larger thresholds.
- **Reusable workspace:** prepare inputs and construct a context once; reuse them in the compute loop.
- **CMake integration:** link the `supered::supered` target; required corrected LEAP sources are included.

The public API covers fixed-length **128 bp, uppercase A/C/G/T** inputs and returns distances. It provides global distance calculation without traceback, alignment coordinates, or CIGAR generation.

## Contents

- [Building](#building)
- [Using supered in your project](#using-supered-in-your-project)
- [Usage](#usage)
- [Input format](#input-format)
- [API](#api)
- [Performance](#performance)
- [Implementation](#implementation)
- [Validation](#validation)
- [Acknowledgements and source provenance](#acknowledgements-and-source-provenance)

## Building

Requirements: Linux/x86-64, CMake 3.20 or later, and a C++20 compiler. The tested toolchain is Clang 22.1.2 with LLVM lld. The current kernels require AVX2, BMI instructions, and AVX-512 extensions including **F, VL, VBMI2, and VPOPCNTDQ**.

```bash
git clone https://github.com/Naomigh/supered.git
cd supered

cmake -S . -B build \
  -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_CXX_COMPILER=clang++-22 \
  -DCMAKE_EXE_LINKER_FLAGS=-fuse-ld=lld
cmake --build build --target supered_example -j4
./build/supered_example
```

The included example prints:

```text
0 1 2 -1
0 1 2 3
0 1 2 3
```

These lines use `Context<2>`, `Context<5>`, and `Context<128>`, respectively. The last pair has distance 3 and exceeds the first context's threshold.

For a library-only build, configure with `-DSUPERED_BUILD_EXAMPLES=OFF` and build `--target supered`. The resulting static library is `build/libsupered.a`. A normal clone contains the required dependencies; submodule initialization is unnecessary.

The build uses `-march=native` and link-time optimization. Build and deploy on compatible CPUs, and use a compatible compiler/linker when linking the static library. CPU dispatch is selected at build time.

## Using supered in your project

Clone this repository into `third_party/supered`, then link its CMake target:

```cmake
cmake_minimum_required(VERSION 3.20)
project(my_alignment_tool LANGUAGES CXX)

set(CMAKE_INTERPROCEDURAL_OPTIMIZATION TRUE)
add_subdirectory(third_party/supered)

add_executable(my_alignment_tool main.cpp)
target_link_libraries(my_alignment_tool PRIVATE supered::supered)
```

Configure your application with a Release build and the same Clang/lld options shown above. The target propagates its include directories, C++20 requirement, and `-march=native`. Examples are disabled by default when supered is added as a subdirectory.

The public header is [`<supered/supered.hpp>`](include/supered/supered.hpp). Myers kernels are defined in the headers; the static library supplies the corrected LEAP backend. Keep both when distributing a built dependency.

## Usage

This complete example compares one read against four anchors. Encoding and context construction happen before the distance call; an application with existing packed inputs can pass them directly.

```cpp
#include <supered/supered.hpp>
#include <iostream>

int main() {
    supered::DNA dna;
    dna.fill('A');
    const auto read = supered::encode(dna);

    supered::PackedDNA128 anchors[4];
    for (int i = 0; i < 4; ++i) {
        auto target = dna;
        for (int j = 0; j < i; ++j) target[10 + 30 * j] = 'C';
        anchors[i] = supered::encode(target);
    }

    supered::Context<80> engine;
    const auto distances = engine.distance4(
        read, {&anchors[0], &anchors[1], &anchors[2], &anchors[3]});

    for (int d : distances) std::cout << d << ' ';
    std::cout << '\n';
}
```

Save as `example.cpp`, then compile and run from the repository root after building the library:

```bash
clang++-22 -O3 -std=c++20 -march=native -flto=thin -fuse-ld=lld \
  -Iinclude -Ileap_origin example.cpp build/libsupered.a -o example
./example
```

Output: `0 1 2 3`.

For four independent pairs, pass two pointer arrays:

```cpp
const auto distances = engine.distance4(
    {&query0, &query1, &query2, &query3},
    {&target0, &target1, &target2, &target3});
```

`distances[i]` corresponds to `query[i]` versus `target[i]`. Inputs may be noncontiguous and remain read-only during the call. To compare one read with 20 contiguous anchors, allocate `int output[20]` and call `engine.distance_many(read, anchors, 20, output)`; this makes five four-pair calls.

### Choosing a threshold

`K` is a compile-time integer from 0 through 128. Each result is the **exact distance when `d <= K`, or `-1` otherwise**. `Context<>` means `Context<128>` and returns the complete distance for every valid input pair.

Choose `Context<2>` to accept distances 0–2, `Context<50>` for 0–50, or `Context<>` when the distance is unknown. Backend selection follows the supplied threshold: `Context<50>` uses Myers even for a pair whose actual distance is 1.

Construct one context per worker and reuse it. Contexts cannot be copied or moved; pass them by reference. The library executes in the calling thread, so applications control parallel scheduling and output ownership.

## Input format

[`PackedDNA128`](include/strong/packed.hpp) has size 32 bytes and alignment 32 bytes:

```cpp
struct alignas(32) PackedDNA128 {
    uint64_t bit0[2];
    uint64_t bit1[2];
};
```

A/C/G/T map to `00/01/10/11`. `bit0` contains the low bit of each base, and `bit1` contains the high bit. Position `i` occupies bit `i % 64` of word `i / 64`; position 0 is the least significant bit. Memory order is `bit0[0], bit0[1], bit1[0], bit1[1]`.

Use `encode(DNA)` to prepare the two bit planes and `decode(packed)` to recover the 128 characters. `DNA` is `std::array<char, 128>`, without a terminating null character. Existing read/index data in this layout can be used directly as properly aligned, live `PackedDNA128` objects. Conventional packing with two adjacent bits per base requires conversion to the two-plane layout.

Callers provide valid 128 bp A/C/G/T inputs, valid pointers, aligned packed objects, and sufficient output storage. The hot path relies on these preconditions and omits runtime input validation. Array batch sizes must be multiples of four.

## API

The following declarations summarize the public interface:

```cpp
namespace supered {
using DNA = std::array<char, 128>;
using Quad = std::array<const PackedDNA128*, 4>;
using Distances = std::array<int, 4>;

PackedDNA128 encode(const DNA& dna);
DNA decode(const PackedDNA128& packed);
Quad broadcast(const PackedDNA128& packed);

template<int K = 128>
class Context {
public:
    static constexpr int max_edit = K;
    static constexpr bool uses_leap = K < 5;

    Distances distance4(const Quad& queries, const Quad& targets);
    Distances distance4(const PackedDNA128& query, const Quad& targets);
    int distance(const PackedDNA128& query, const PackedDNA128& target);

    void distance_many(const PackedDNA128& query,
                       const PackedDNA128* targets,
                       std::size_t count, int* output);
    void distance_pairs(const PackedDNA128* queries,
                        const PackedDNA128* targets,
                        std::size_t count, int* output);
};
}
```

| Method | Comparisons | Result |
| --- | --- | --- |
| `distance4(queries, targets)` | Four element-wise pairs | Four distances, in lane order |
| `distance4(query, targets)` | One query against four targets | Four distances, in target order |
| `distance(query, target)` | One pair | One distance |
| `distance_many(query, targets, count, output)` | One query against a contiguous target array | `output[i] = distance(query, targets[i])` |
| `distance_pairs(queries, targets, count, output)` | Two contiguous arrays, element by element | `output[i] = distance(queries[i], targets[i])` |

Every distance obeys the same threshold semantics. Both array APIs require `count` to be a multiple of four and `output` to hold at least `count` integers. `broadcast(packed)` creates four pointers to one existing packed object.

For `K >= 5`, `distance4()` processes four useful pairs concurrently. The convenience `distance()` broadcasts one pair across the Myers lanes and returns lane 0; use the batch APIs for throughput. For `K < 5`, `distance4()` makes four individual LEAP calls. Its small-threshold throughput is therefore a LEAP result, rather than four-pair AVX-512 execution.

## Performance

### Comparison with original Edlib

The following results use the **same independent input pairs and threshold** for both methods, measured on 2026-10-22. Units are **ns/pair**; `±` is the 95% confidence-interval half-width across five fitted slopes. The supered column is already normalized by all four useful pairs.

| `K` | Original Edlib | supered four-pair API | Speedup |
| ---: | ---: | ---: | ---: |
| 0 | 1344.37 ± 4.12 | 16.17 ± 0.12 | 83.13× |
| 2 | 1346.30 ± 2.59 | 66.87 ± 0.88 | 20.13× |
| 5 | 1354.27 ± 2.72 | 64.37 ± 0.16 | 21.04× |
| 10 | 1345.76 ± 6.27 | 64.42 ± 0.07 | 20.89× |
| 20 | 1374.90 ± 2.25 | 64.54 ± 0.10 | 21.30× |
| 30 | 1401.45 ± 0.79 | 64.34 ± 0.09 | 21.78× |
| 40 | 1417.34 ± 3.80 | 64.75 ± 0.82 | 21.89× |
| 50 | 1438.45 ± 1.67 | 64.35 ± 0.09 | 22.35× |
| 60 | 1446.94 ± 13.03 | 64.26 ± 0.12 | 22.52× |
| 70 | 1472.58 ± 1.72 | 99.50 ± 0.21 | 14.80× |
| 80 | 1477.02 ± 1.31 | 99.43 ± 0.20 | 14.86× |
| 90 | 1483.12 ± 1.23 | 99.55 ± 0.03 | 14.90× |

The shown workload contains 655,360 independent 128 bp pairs per threshold, all with true distance `<= K`. Queries are random DNA; targets contain substitutions and, where the edit budget permits, balanced insertions/deletions. `K=0` uses identical query/target sequences; `K=2` near inputs contain substitutions only. `K=0/2` use LEAP; the other listed thresholds use four-pair Myers.

Original Edlib means its **unmodified global NW distance kernel**, with prepared equality masks and target codes; its original block allocation remains inside the timed call. This comparison excludes the public `edlibAlign()` wrapper. supered is measured through its public `distance4()` API. The [complete results](docs/benchmarks/2026-10-22/independent_summary.csv) also include mixed passing/random inputs, single-pair calls, and a separately labeled Edlib reference with stack block storage.

Each repeat fits `T(N) = a + bN` at 131,072, 262,144, 393,216, 524,288, and 655,360 pairs; the slope `b` represents per-pair compute time. Inputs are prepared once, copied into cache-sized scratch batches, and touched before timing. Encoding, Edlib input preparation, context allocation, batch loading, and warming are excluded. Input loads, per-call SIMD setup, pair-dependent masks, state initialization, recurrence, thresholding, and result checksums are included.

### Read × anchor simulation

With **10,000,000 reads and 20 anchors**, all 128 bp, `Context<80>` computed **200,000,000 distances**. Every accepted read has true distance `<= 80` from **every** anchor, and every result was independently checked against original Edlib.

- **97.91 ± 0.25 ns/pair**, fitted across 40–200 million comparisons, three repeats.
- **391.65 ns per four-pair call**, equivalent to **10.21 million pairs/s**.
- **19.58 s** fitted compute time for 200 million pairs; the full measured compute passes took **19.88–19.89 s**.
- **20.77–20.78 s** full-trial wall time, including batch movement, checks, and cache-gate retries.

The 20 packed anchors occupy 640 bytes. Prepared reads enter scratch in batches of 16–80, outside the compute window; each read makes five four-pair calls against the resident anchors. Input generation, rejection sampling, and oracle validation are excluded from supered compute time. This is a separate workload from the Edlib comparison above.

### Machine and cache conditions

Both experiments used one compute thread pinned to logical CPU 24 in a Linux VMware guest reporting an AMD EPYC 9654, with Clang 22.1.2, C++20, `-O3`, `-march=native`, and LTO. These are single-thread measurements; the guest does not establish exclusive physical-core ownership.

A PMU gate accepted windows with zero observed DRAM/IO demand fills and retried the same inputs otherwise. L2/L3 hits were allowed, and L2 misses were retained. The results describe cache-conditioned computation and do not imply that every load hits L1. See the [benchmark notes, source identities, fit points, and summaries](docs/benchmarks/2026-10-22/README.md).

## Implementation

supered specializes Myers bit-vector dynamic programming for fixed 128 bp global alignment. Its packed input is distinct from the evolving Myers `Pv/Mv` state. Equality masks are derived from the two bit planes inside the distance calculation.

For `5 <= K < 64`, a 64-bit moving band uses the equal-length global-alignment bound `|i-j| <= floor(K/2)`. Larger thresholds use two 64-bit blocks with explicit cross-block carry and shift propagation. Four independent pairs occupy four 64-bit lanes; AVX-512 instructions, including 256-bit forms enabled by AVX-512VL, update these lanes concurrently. The implementation uses intrinsics and does not require each operation to use a 512-bit register.

For `K < 5`, the corrected AVX2 LEAP backend consumes the same packed inputs and reuses its context workspace. Threshold specialization selects the backend at compile time.

## Validation

The independent-pair benchmark validated **15,728,640 input pairs** against original Edlib and checked **98,304 pairs** with independent scalar dynamic programming. The read/anchor experiment checked all **200 million results** against original Edlib, with another **20,000 scalar-DP checks**. Timed windows also verified output checksums. These describe the archived experiments; they are not new tests run during documentation updates.

The distribution includes the runnable [`examples/basic.cpp`](examples/basic.cpp). The [published benchmark data](docs/benchmarks/2026-10-22/README.md) contain the source identities and measurement points underlying the reported figures. The full laboratory benchmark harness and per-window binary archives are maintained separately from this library distribution.

## Acknowledgements and source provenance

supered builds on [Myers' bit-vector algorithm](https://doi.org/10.1145/316542.316550), the engineering approach of [Edlib](https://github.com/Martinsos/edlib), and the [corrected LEAP branch](https://github.com/YourongZhou/LEAP/tree/codex/leap-correctness-fixes).

The included LEAP subset is pinned to commit [`6a711282`](https://github.com/YourongZhou/LEAP/commit/6a711282f3494b6b9886236e58c7dfba6efa885a). Original source notices are retained; [NOTICE](NOTICE) and [`leap_origin/UPSTREAM.json`](leap_origin/UPSTREAM.json) record its provenance. This documentation does not assign a new license to third-party source files.
