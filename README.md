# supered

**English** | [简体中文](README.zh-CN.md)

**supered is an ultrafast edit-distance computation library for fixed-length 128 bp DNA sequences.** It combines compact 32-byte sequence inputs, AVX-512 SIMD processing of four sequence pairs, specialized fixed-length kernels, and reusable working memory. A corrected AVX2-based LEAP backend handles very small edit thresholds through the same public API.

Get the code:

```bash
git clone https://github.com/Naomigh/supered.git
cd supered
```

This guide covers building, linking, and using the library. The public header is [`include/supered/supered.hpp`](include/supered/supered.hpp), and the CMake target is `supered::supered`.

## 1. Distance semantics and input contract

supered computes the global Levenshtein edit distance between two **128 bp DNA sequences**. Insertions, deletions, and substitutions each cost 1; a match costs 0. Input consists of uppercase `A/C/G/T` characters or data already packed in the format described below.

Choose the maximum accepted distance, `max_edit=K`, through `Context<K>`:

- If the true distance `d <= K`, return the exact distance `d`.
- If the true distance `d > K`, return `-1`.
- `Context<>` is equivalent to `Context<128>`. Any two 128 bp DNA sequences have distance at most 128, so this returns the full distance.

The output contains distances only, without a CIGAR, alignment coordinates, or an alignment path. The public API supports fixed-length, global, unit-cost alignment.

The caller guarantees valid lengths, characters, packed data, pointers, and output storage. The compute path does not validate these conditions at runtime. The contiguous-array batch APIs require a pair count that is a multiple of four.

## 2. Requirements and repository layout

### Build requirements

- Linux / x86-64.
- CMake 3.20 or later and a C++20 compiler.
- The build commands below use Clang 22.1.2 and LLVM lld.
- The general SIMD path requires AVX2 and extensions including AVX-512F, AVX-512VL, AVX-512VBMI2, and AVX-512VPOPCNTDQ. The LEAP path uses its AVX2 backend and BMI instructions.

Builds use `-march=native` to target the build machine. There is no runtime CPU dispatch; deploy the resulting binary on a machine that supports its build target. An AVX2-only machine does not meet the requirements of the current general four-pair kernel.

### Repository layout

```text
supered/
├── CMakeLists.txt
├── README.md                   English usage guide
├── README.zh-CN.md             Chinese usage guide
├── include/supered/supered.hpp  Public API
├── include/strong/             Core implementation
├── examples/basic.cpp          Runnable usage example
├── leap_origin/                Required corrected LEAP source and provenance
└── NOTICE                      Third-party source notice
```

The required LEAP source is included as ordinary files, so a normal clone is ready to build without submodule initialization. `leap_origin/UPSTREAM.json` records the pinned upstream commit and file checksums; original notices are retained. See [NOTICE](NOTICE) for provenance.

## 3. Build and run the included example

From the root of your cloned supered repository:

```bash
cmake -S . -B build_clang \
  -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_CXX_COMPILER=clang++-22 \
  -DCMAKE_EXE_LINKER_FLAGS=-fuse-ld=lld

cmake --build build_clang --target supered_example -j4
./build_clang/supered_example
```

Expected output:

```text
0 1 2 -1
0 1 2 3
0 1 2 3
```

The three lines use `Context<2>`, `Context<5>`, and `Context<128>` on the same four input pairs. The pair with true distance 3 returns `-1` when `K=2`.

To build only the static library:

```bash
cmake --build build_clang --target supered -j4
```

To configure a library-only build, add `-DSUPERED_BUILD_EXAMPLES=OFF` to the CMake configuration command. Standalone builds enable the example by default.

The output is `build_clang/libsupered.a`. The static library contains the compiled LEAP implementation; Myers kernels are inlined into the caller through the public headers. Provide both the headers and the static library when linking your application.

## 4. Quick start: compare four independent pairs in one call

This complete program prepares four `(query[i], target[i])` pairs and computes their distances with one `distance4()` call.

Save the following as `four_pairs.cpp`:

```cpp
#include <supered/supered.hpp>
#include <array>
#include <iostream>

int main() {
    constexpr char bases[] = "ACGT";
    std::array<supered::PackedDNA128, 4> queries;
    std::array<supered::PackedDNA128, 4> targets;

    // Prepare inputs. Production workflows can use existing packed data directly.
    for (int lane = 0; lane < 4; ++lane) {
        supered::DNA query;
        query.fill(bases[lane]);
        auto target = query;
        for (int j = 0; j < lane; ++j)
            target[10 + 30 * j] = bases[(lane + 1) % 4];
        queries[lane] = supered::encode(query);
        targets[lane] = supered::encode(target);
    }

    supered::Context<5> engine;
    supered::Quad query_ptrs{
        &queries[0], &queries[1], &queries[2], &queries[3]
    };
    supered::Quad target_ptrs{
        &targets[0], &targets[1], &targets[2], &targets[3]
    };

    // K=5 uses the AVX-512 four-pair kernel.
    const auto distances = engine.distance4(query_ptrs, target_ptrs);

    for (int d : distances) std::cout << d << ' ';
    std::cout << '\n';
}
```

Compile, link, and run from the supered directory:

```bash
clang++-22 -O3 -std=c++20 -march=native -flto=thin -fuse-ld=lld \
  -Iinclude -Ileap_origin \
  four_pairs.cpp build_clang/libsupered.a \
  -o four_pairs

./four_pairs
```

Expected output:

```text
0 1 2 3
```

Results correspond to `query[0]` versus `target[0]`, `query[1]` versus `target[1]`, and so on.

`distance4()` accepts four input pointers, so the queries and targets do not need to be contiguous. The pointed-to objects must remain valid throughout the call. The inputs are read-only.

## 5. Integrate with your CMake project

Assume your application contains `main.cpp`. Clone the dependency from your application's root directory:

```bash
git clone https://github.com/Naomigh/supered.git third_party/supered
```

Use the following CMake configuration:

```cmake
cmake_minimum_required(VERSION 3.20)
project(my_alignment_tool LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_INTERPROCEDURAL_OPTIMIZATION TRUE)

add_subdirectory(
    ${CMAKE_CURRENT_SOURCE_DIR}/third_party/supered
    ${CMAKE_BINARY_DIR}/supered-build
)

add_executable(my_alignment_tool main.cpp)
target_link_libraries(my_alignment_tool PRIVATE supered::supered)
```

Configure and build the application:

```bash
cmake -S . -B build \
  -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_CXX_COMPILER=clang++-22 \
  -DCMAKE_EXE_LINKER_FLAGS=-fuse-ld=lld

cmake --build build --target my_alignment_tool -j4
```

`supered::supered` propagates the required include directories, C++20 requirement, and `-march=native` option. The current build enables LTO; use a compatible Clang/lld toolchain when linking the existing LLVM ThinLTO static library.

When included through `add_subdirectory()`, the usage example is disabled by default. Build your application with `--target my_alignment_tool`; only the library dependency is added to that target.

## 6. Input format: 32-byte PackedDNA128

### Use existing preprocessed data directly

The public type has this layout:

```cpp
struct alignas(32) PackedDNA128 {
    uint64_t bit0[2];
    uint64_t bit1[2];
};
```

`sizeof(supered::PackedDNA128)` is 32, with 32-byte type alignment. Bases are encoded as A=`00`, C=`01`, G=`10`, and T=`11`:

- `bit0` stores the low bit of each base; `bit1` stores the high bit.
- Word `[0]` stores positions 0..63; word `[1]` stores positions 64..127.
- Base `i` occupies bit `i % 64` of word `i / 64`. The first base occupies the least significant bit.
- The four words appear in memory as `bit0[0], bit0[1], bit1[0], bit1[1]`.

If your read/index preprocessing already produces this layout, pass it directly as `PackedDNA128` objects. Allocate arrays using this type, or preserve the same type, alignment, and object lifetime within your own storage layout.

A conventional packed array with each base's two bits adjacent must first be converted into these two bit planes.

### Prepare inputs from DNA characters

`supered::DNA` is `std::array<char, 128>`. For example:

```cpp
#include <supered/supered.hpp>

supered::PackedDNA128 make_example_input() {
    supered::DNA dna;
    dna.fill('A');
    dna[10] = 'C';
    return supered::encode(dna);
}
```

If you already have a `std::string` containing exactly 128 characters:

```cpp
#include <supered/supered.hpp>
#include <algorithm>
#include <string>

// The caller provides exactly 128 characters, all from A/C/G/T.
supered::PackedDNA128 pack_read(const std::string& text) {
    supered::DNA dna;
    std::copy_n(text.begin(), 128, dna.begin());
    return supered::encode(dna);
}
```

To recover characters, call `supered::decode(packed)`, which returns `DNA`. This contains exactly 128 characters. You can construct a string with `std::string(dna.begin(), dna.end())`.

Encode each read or candidate sequence ahead of time and reuse the packed result. Already packed inputs can be passed directly to the distance API.

## 7. Choose max_edit and reuse a Context

`K` is a compile-time template argument, with valid integer values from 0 through 128:

```cpp
#include <supered/supered.hpp>

void context_examples() {
    supered::Context<0> exact_match;  // Accept distance 0 only.
    supered::Context<2> tiny;         // Accept distances 0..2.
    supered::Context<5> bounded;      // Accept distances 0..5.
    supered::Context<50> wider;       // Accept distances 0..50.
    supered::Context<> full;          // K=128: return the full distance.

    static_assert(decltype(tiny)::max_edit == 2);
    static_assert(decltype(tiny)::uses_leap);
    static_assert(!decltype(bounded)::uses_leap);
}
```

The backend is selected by K:

- `K=0..4`: corrected LEAP.
- `K=5..63`: AVX-512 four-pair banded Myers kernel.
- `K=64..128`: AVX-512 four-pair complete two-block Myers kernel.

Dispatch follows the supplied upper bound. A `Context<50>` uses its Myers backend even if a particular pair turns out to have distance 1.

For `K>=5`, `distance4()` uses four-pair SIMD. For `K<5`, the same API performs four individual LEAP calculations.

Use `Context<>` when you need the full distance for arbitrary inputs. For example, if the true distance is 8, `Context<5>` returns `-1`, while `Context<10>` and `Context<>` both return 8.

There is currently no overload that accepts an arbitrary runtime `int max_edit`. If your workflow uses several thresholds, construct the corresponding `Context<K>` instances beforehand and select the appropriate instance in your application.

Create a Context during worker initialization and reuse it. Context objects cannot be copied or moved; pass them between functions by reference or pointer.

## 8. Public API reference

The following declarations summarize the public API. All types are in the `supered` namespace:

```cpp
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
```

### `distance4(queries, targets)`: four independent pairs

Section 4 contains a complete example. `result[i]` is the distance between `queries[i]` and `targets[i]`, in lane order 0..3.

Each lane independently returns an exact distance or `-1`, for example `{0, 2, -1, 5}`. One call produces four pair results.

### `distance4(query, targets)`: one query against four candidates

Use this overload for a read with several candidate reference fragments. The following function can be used with already prepared data:

```cpp
#include <supered/supered.hpp>

supered::Distances compare_four_candidates(
    supered::Context<50>& engine,
    const supered::PackedDNA128& read,
    const supered::PackedDNA128& candidate0,
    const supered::PackedDNA128& candidate1,
    const supered::PackedDNA128& candidate2,
    const supered::PackedDNA128& candidate3) {
    return engine.distance4(read,
        {&candidate0, &candidate1, &candidate2, &candidate3});
}
```

The four returned values correspond to `read` versus `candidate0..3`. The candidates may reside at different addresses.

To obtain all 16 combinations of four queries and four targets, call this overload once for each query.

### `distance(query, target)`: one pair

```cpp
#include <supered/supered.hpp>

int compare_one_pair(supered::Context<>& engine,
                     const supered::PackedDNA128& query,
                     const supered::PackedDNA128& target) {
    return engine.distance(query, target);
}
```

This example uses the default K=128 and returns the full distance. With a bounded Context, the result is still the exact distance when it is <=K, or `-1` otherwise.

This is a convenience API for an isolated pair. The current Myers backend broadcasts that pair across four lanes and returns one result; the LEAP backend computes it once. For a stream of independent pairs, use the batch APIs to fill all four lanes.

### `distance_many(query, targets, count, output)`: one query against a contiguous array

`targets` points to a contiguous `PackedDNA128` array, `count` is the number of candidates, and `output` must hold at least `count` integers. `output[i]` is the distance between `query` and `targets[i]`. `count` must be a multiple of four.

```cpp
#include <supered/supered.hpp>
#include <vector>

// Candidates are already packed; their count is a multiple of four.
std::vector<int> compare_candidates(
    supered::Context<50>& engine,
    const supered::PackedDNA128& read,
    const std::vector<supered::PackedDNA128>& candidates) {
    std::vector<int> output(candidates.size());
    engine.distance_many(read, candidates.data(), candidates.size(), output.data());
    return output;
}
```

This wrapper allocates output for demonstration. In a high-frequency workflow, allocate output storage ahead of time and reuse it with direct `distance_many()` calls.

### `distance_pairs(queries, targets, count, output)`: pair two contiguous arrays element by element

`queries` and `targets` each contain at least `count` elements, and `output` must hold at least `count` integers. `output[i]` is the distance between `queries[i]` and `targets[i]`. `count` must be a multiple of four.

```cpp
#include <supered/supered.hpp>
#include <vector>

// Both arrays are packed, have equal lengths, and contain a multiple of four pairs.
std::vector<int> compare_pair_arrays(
    supered::Context<50>& engine,
    const std::vector<supered::PackedDNA128>& queries,
    const std::vector<supered::PackedDNA128>& targets) {
    std::vector<int> output(queries.size());
    engine.distance_pairs(queries.data(), targets.data(), queries.size(), output.data());
    return output;
}
```

This API automatically calls `distance4()` for each group of four pairs and writes results in the original array order.

### `broadcast(packed)`: construct four pointers to the same input

`broadcast(packed)` returns `{&packed, &packed, &packed, &packed}` for constructing a `Quad` with a shared input. It only constructs a pointer array; the original object must remain alive while those pointers are used.

For one query against four targets, you can usually call `distance4(query, targets)` directly.

## 9. Threads, memory, and high-frequency calls

- The library executes in the calling thread and does not create compute threads automatically.
- Give each compute thread its own Context, especially for the LEAP backend's mutable workspace.
- Workers may share read-only packed inputs; write their outputs to separate storage.
- Construct and reuse the Context outside the compute loop, and pass it by reference.
- Use `distance4()` for four independent pairs, `distance_many()` for one read against contiguous candidates, and `distance_pairs()` for two corresponding contiguous input arrays.
- Encode inputs and allocate output storage before the compute loop, then reuse them across calls.

The library handles the per-call state initialization and matching information needed for computation. Callers prepare the PackedDNA128 inputs and Context. Cache residency depends on the caller's working set and scheduling; the API does not read sequence files or lock an entire dataset into CPU cache.
