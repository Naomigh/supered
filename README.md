# supered

**supered is an ultrafast edit-distance computation library for fixed-length 128 bp DNA sequences.** It combines compact 32-byte sequence inputs, AVX-512 SIMD processing of four sequence pairs, specialized fixed-length kernels, and reusable working memory. A corrected AVX2-based LEAP backend handles very small edit thresholds through the same public API.

获取代码：

```bash
git clone https://github.com/Naomigh/supered.git
cd supered
```

本文介绍如何构建、链接和调用代码库。公共头文件是 [`include/supered/supered.hpp`](include/supered/supered.hpp)，CMake target 是 `supered::supered`。

## 1. 计算内容与使用约定

supered 计算两条 **128 bp DNA** 的全局 Levenshtein 编辑距离：插入、删除和替换各计 1，匹配计 0。输入为大写 `A/C/G/T`，或者已经按下文格式编码的数据。

调用者通过 `Context<K>` 指定最大允许距离 `max_edit=K`：

- 真实距离 `d <= K`：返回精确距离 `d`。
- 真实距离 `d > K`：返回 `-1`。
- `Context<>` 等价于 `Context<128>`；任意两条 128 bp DNA 的距离都不超过 128，因此它返回完整距离。

返回值只有距离，不包含 CIGAR、比对位置或比对路径。公开接口处理固定长度、全局、单位编辑代价的比对。

输入长度、字符、编码、指针和输出空间由调用者保证；计算路径不做这些运行时校验。连续数组批量接口的 pair 数必须是 4 的倍数。

## 2. 环境与目录

### 构建环境

- Linux / x86-64。
- CMake 3.20 或更新版本，以及支持 C++20 的编译器。
- Boost 头文件（LEAP 的 `mask.h` 使用 Boost.Preprocessor；Ubuntu/Debian 对应 `libboost-dev`）。
- 当前已经验证的工具链：Clang 22.1.2 与 LLVM lld。
- 常规 SIMD 路径需要 AVX2 以及 AVX-512F、AVX-512VL、AVX-512VBMI2、AVX-512VPOPCNTDQ 等扩展；LEAP 路径使用 AVX2 后端及 BMI 指令。

构建使用 `-march=native`，目标是构建机器的指令能力。库没有运行时 CPU 能力分派；部署时使用支持该构建目标的机器。仅支持 AVX2 的机器不满足当前常规四 pair 内核的要求。

### 当前目录布局

```text
supered/
├── CMakeLists.txt
├── include/supered/supered.hpp  公共接口
├── include/strong/             内部实现
├── src/                        参考实现适配
├── generated/                  已生成的参考代码
├── examples/basic.cpp          可运行的基本示例
├── tests/                      正确性验证程序
├── bench/                      性能实验程序
├── scripts/                    实验与分析脚本
├── edlib_origin/               原始 Edlib 参考源码与许可证
├── leap_origin/                修正版 LEAP 所需源码
├── THIRD_PARTY.md              上游来源及版本
└── results/README.md           结果保存说明
```

两个依赖目录以实际源码随仓库提供，克隆后即可构建，不需要初始化 submodule。
上游 commit 和文件校验值记录在各目录的 `UPSTREAM.json` 中。第三方来源见 [THIRD_PARTY.md](THIRD_PARTY.md)。

## 3. 构建与运行现有示例

在克隆后的 supered 根目录执行：

```bash
cmake -S . -B build_clang \
  -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_C_COMPILER=clang-22 \
  -DCMAKE_CXX_COMPILER=clang++-22 \
  -DCMAKE_EXE_LINKER_FLAGS=-fuse-ld=lld

cmake --build build_clang --target supered_example supered_validate -j4
./build_clang/supered_example
```

示例输出：

```text
0 1 2 -1
0 1 2 3
0 1 2 3
```

三行分别使用 `Context<2>`、`Context<5>` 和 `Context<128>` 比对相同的四组输入。真实距离为 3 的 pair 在 `K=2` 时返回 `-1`。

只构建静态库：

```bash
cmake --build build_clang --target supered -j4
```

输出为 `build_clang/libsupered.a`。库包含 LEAP 的编译实现；Myers 内核通过公共头文件内联到调用者中，因此使用时需要同时提供头文件和静态库。

## 4. 快速开始：一次比对四组独立 pair

下面是一个完整程序。它准备四组 `(query[i], target[i])`，通过一次 `distance4()` 调用计算四个距离。

将下面代码保存为 `four_pairs.cpp`：

```cpp
#include <supered/supered.hpp>
#include <array>
#include <iostream>

int main() {
    constexpr char bases[] = "ACGT";
    std::array<supered::PackedDNA128, 4> queries;
    std::array<supered::PackedDNA128, 4> targets;

    // 准备输入。实际工作流中可以直接使用已有的预编码数据。
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

    // K=5 使用 AVX-512 四 pair 内核。
    const auto distances = engine.distance4(query_ptrs, target_ptrs);

    for (int d : distances) std::cout << d << ' ';
    std::cout << '\n';
}
```

在 supered 目录下编译、链接和运行：

```bash
clang++-22 -O3 -std=c++20 -march=native -flto=thin -fuse-ld=lld \
  -Iinclude -Ileap_origin \
  four_pairs.cpp build_clang/libsupered.a \
  -o four_pairs

./four_pairs
```

输出：

```text
0 1 2 3
```

结果顺序为 `query[0]` 对 `target[0]`、`query[1]` 对 `target[1]`，依此类推。

`distance4()` 接收的是四个输入指针，所以四条 query 或 target 可以位于不同位置，不要求彼此连续。输入对象需要在本次调用期间保持有效。计算不会修改输入。

## 5. 接入自己的 CMake 工程

下面假设调用工程包含 `main.cpp`，先在该工程根目录克隆依赖：

```bash
git clone https://github.com/Naomigh/supered.git third_party/supered
```

```cmake
cmake_minimum_required(VERSION 3.20)
project(my_alignment_tool LANGUAGES C CXX)

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

配置并构建调用程序：

```bash
cmake -S . -B build \
  -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_C_COMPILER=clang-22 \
  -DCMAKE_CXX_COMPILER=clang++-22 \
  -DCMAKE_EXE_LINKER_FLAGS=-fuse-ld=lld

cmake --build build --target my_alignment_tool -j4
```

`supered::supered` 会传递所需的头文件目录、C++20 要求和 `-march=native` 编译选项。当前构建启用了 LTO；链接已有的 LLVM ThinLTO 静态库时使用兼容的 Clang/lld 工具链。

当前 CMake 项目也定义了验证和实验 target。指定 `--target my_alignment_tool` 可以只构建调用程序及其库依赖。

## 6. 输入数据：32-byte PackedDNA128

### 直接使用已有预处理结果

公开类型的布局为：

```cpp
struct alignas(32) PackedDNA128 {
    uint64_t bit0[2];
    uint64_t bit1[2];
};
```

`sizeof(supered::PackedDNA128)` 为 32，类型对齐为 32 bytes。每个碱基的编码为 A=`00`、C=`01`、G=`10`、T=`11`：

- `bit0` 保存每个碱基的低 bit，`bit1` 保存高 bit。
- `[0]` 保存位置 0..63，`[1]` 保存位置 64..127。
- 第 `i` 个碱基对应 word `i / 64` 中的 bit `i % 64`；第一个碱基位于最低 bit。
- 四个 word 的内存顺序为 `bit0[0], bit0[1], bit1[0], bit1[1]`。

已有的 read/index 预处理结果符合这一布局时，可以直接以 `PackedDNA128` 对象传入。使用该类型分配数组，或者在自有存储布局中保持相同的类型、对齐及对象生存期。

普通的“每个碱基两 bit 相邻存储”数组需要先转换为这里的两个位平面。

### 从 DNA 字符准备输入

`supered::DNA` 是 `std::array<char, 128>`。例如：

```cpp
#include <supered/supered.hpp>

supered::PackedDNA128 make_example_input() {
    supered::DNA dna;
    dna.fill('A');
    dna[10] = 'C';
    return supered::encode(dna);
}
```

如果已有长度恰好为 128 的 `std::string`：

```cpp
#include <supered/supered.hpp>
#include <algorithm>
#include <string>

// 调用者保证 text 长度为 128，且只含 A/C/G/T。
supered::PackedDNA128 pack_read(const std::string& text) {
    supered::DNA dna;
    std::copy_n(text.begin(), 128, dna.begin());
    return supered::encode(dna);
}
```

需要读取编码后的字符时，调用 `supered::decode(packed)`，返回 `DNA`。它包含恰好 128 个字符；构造字符串可以使用 `std::string(dna.begin(), dna.end())`。

对同一条 read 或候选序列预先编码并复用；已有合法编码数据时直接调用距离接口。

## 7. 选择 max_edit 与 Context

`K` 是编译期模板参数，合法范围为整数 0..128：

```cpp
#include <supered/supered.hpp>

void context_examples() {
    supered::Context<0> exact_match;  // 只接受距离 0
    supered::Context<2> tiny;         // 只接受距离 0..2
    supered::Context<5> bounded;      // 只接受距离 0..5
    supered::Context<50> wider;       // 只接受距离 0..50
    supered::Context<> full;          // K=128，返回完整距离

    static_assert(decltype(tiny)::max_edit == 2);
    static_assert(decltype(tiny)::uses_leap);
    static_assert(!decltype(bounded)::uses_leap);
}
```

后端选择由 K 决定：

- `K=0..4`：修正版 LEAP。
- `K=5..63`：AVX-512 四 pair Myers 带宽受限内核。
- `K=64..128`：AVX-512 四 pair Myers 完整双 block 内核。

路由依据是给定上限。使用 `Context<50>` 时，即使某个 pair 的实际距离为 1，也由该 Context 的 Myers 后端处理。

`distance4()` 在 `K>=5` 时使用四 pair SIMD；在 `K<5` 时以统一接口执行四次 LEAP 单 pair 计算。

需要任意输入的完整距离时使用 `Context<>`。例如真实距离为 8：`Context<5>` 返回 `-1`，`Context<10>` 和 `Context<>` 均返回 8。

当前没有接收任意运行时 `int max_edit` 的重载。工作流使用多个阈值时，可以预先构造对应的多个 `Context<K>`，由上层选择调用哪个实例。

在 worker 初始化阶段创建 Context，然后重复使用。Context 禁止复制和移动，函数间通过引用或指针传递。

## 8. 公共接口参考

以下声明概括实际公共接口；类型均位于 `supered` 命名空间：

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

### `distance4(queries, targets)`：四组独立 pair

第 4 节给出了完整示例。返回值 `result[i]` 对应 `queries[i]` 与 `targets[i]` 的距离，按 lane 0..3 排列。

每个 lane 独立返回精确距离或 `-1`，例如 `{0, 2, -1, 5}`。一次调用得到四个 pair 的结果。

### `distance4(query, targets)`：一个 query 对四条候选

适合一个 read 对多个候选片段。下面的函数可以直接接入已经准备好数据的工作流：

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

返回的四项依次对应 `read` 与 `candidate0..3`。这四条候选可以位于不同地址。

需要四条 query 与四条 target 的全部 16 个组合时，对每条 query 分别调用这一重载即可。

### `distance(query, target)`：单 pair

```cpp
#include <supered/supered.hpp>

int compare_one_pair(supered::Context<>& engine,
                     const supered::PackedDNA128& query,
                     const supered::PackedDNA128& target) {
    return engine.distance(query, target);
}
```

示例使用默认 K=128，所以返回完整距离。使用有界 Context 时，返回值仍遵循“距离 <=K 时返回距离，否则 -1”。

这是处理孤立 pair 的便利接口。当前 Myers 后端会将该 pair 广播到四个 lane，再取一个结果；LEAP 后端直接计算一次。
四 pair 的摊销吞吐与这里单 pair 调用的耗时需要分别测量。

### `distance_many(query, targets, count, output)`：一对多连续数组

`targets` 是连续的 `PackedDNA128` 数组，`count` 是候选数量，`output` 至少可写入 `count` 个 `int`。
`output[i]` 对应 `query` 与 `targets[i]`。`count` 必须是 4 的倍数。

```cpp
#include <supered/supered.hpp>
#include <vector>

// candidates 已编码，数量为 4 的倍数。
std::vector<int> compare_candidates(
    supered::Context<50>& engine,
    const supered::PackedDNA128& read,
    const std::vector<supered::PackedDNA128>& candidates) {
    std::vector<int> output(candidates.size());
    engine.distance_many(read, candidates.data(), candidates.size(), output.data());
    return output;
}
```

这个封装为了演示而在函数内分配输出。在高频工作流中，可以提前分配输出空间并复用，直接调用 `distance_many()`。

### `distance_pairs(queries, targets, count, output)`：两个连续数组逐项配对

`queries` 和 `targets` 各有至少 `count` 个元素，输出至少可写入 `count` 个 `int`。
`output[i]` 对应 `queries[i]` 与 `targets[i]`。`count` 必须是 4 的倍数。

```cpp
#include <supered/supered.hpp>
#include <vector>

// 两个数组均已编码、长度相同，且长度为 4 的倍数。
std::vector<int> compare_pair_arrays(
    supered::Context<50>& engine,
    const std::vector<supered::PackedDNA128>& queries,
    const std::vector<supered::PackedDNA128>& targets) {
    std::vector<int> output(queries.size());
    engine.distance_pairs(queries.data(), targets.data(), queries.size(), output.data());
    return output;
}
```

这个接口每四组 pair 自动调用一次 `distance4()`，结果按原数组顺序写入。

### `broadcast(packed)`：构造四个相同输入指针

`broadcast(packed)` 返回 `{&packed, &packed, &packed, &packed}`，用于手工构造共享输入的 `Quad`。
它只生成指针数组，指向的原始对象需要在使用这些指针期间保持有效。

通常，一个 query 对四条 target 直接使用 `distance4(query, targets)` 即可。

## 9. 线程、内存与高频调用

- 库自身在当前调用线程中执行，不自动创建计算线程。
- 每个计算线程持有自己的 Context，尤其是持有可变工作区的 LEAP Context。
- 只读的编码输入可以跨 worker 共享；输出分别写入各 worker 的独立空间。
- Context 在循环外构造并复用；函数间通过引用传递。
- 四个独立 pair 使用 `distance4()`；一条 read 对连续候选使用 `distance_many()`；两组连续输入逐项配对使用 `distance_pairs()`。
- 将编码、内存分配和新一批数据的准备放在计算计时之外；已经预处理的数据直接传给接口。

一次计算所需的状态初始化与匹配信息处理已由库完成，调用者只需要准备 PackedDNA128 输入和 Context。缓存驻留由调用方的工作集与调度方式决定；API 不负责从文件读取序列或将整个数据集锁定在 CPU cache 中。

## 10. 运行正确性验证

```bash
cmake --build build_clang --target supered_validate -j4
ctest --test-dir build_clang -R supered_api_correctness --output-on-failure
```

也可以直接运行：

```bash
./build_clang/supered_validate
```

该程序对照原始 Edlib，验证各个公开距离接口和多个 K，并用独立标量 DP 检查部分输入。成功时输出以 `PASS supered` 开头的验证摘要。

更多已有内核验证程序位于 [`tests/`](tests/)。

## 11. 性能实验工具与原始结果

使用库本身不需要性能计数器权限。`bench/` 中的实验程序另外依赖 Linux perf / RDPMC 和本实验机的 AMD 事件配置，并固定到 CPU 24。

构建基准程序：

```bash
cmake --build build_clang --target benchmark -j4
```

当前基础基准程序的命令行是：

```text
benchmark OUTPUT K near|cross TILES REPEATS [EXISTING_INPUTS]
```

其中 `OUTPUT` 为新的结果目录，`near` 为接近固定 query 的输入，`cross` 为固定 queries 对候选的全交叉，`EXISTING_INPUTS` 可复用同一程序保存的输入。
当前该程序的命令行 K 支持 0、2、5、50；库的 `Context<K>` 支持完整的 0..128 范围。

例如，比较原始 Edlib 预处理后的 NW 核心与 supered 四 pair 公共入口：

```bash
STRONG_ALGORITHMS=0,11 ./build_clang/benchmark \
  results/my_benchmark/k5_near 5 near 32768 5
```

结果目录需使用尚未保存过 trial 的新路径。上述命令在每组中使用五个计算量、重复五次；数据准备和 cache touch 位于计时外，计算窗口保留时间、PMU 和结果 checksum。
该基础基准中的 supered 项使用 `distance4()`；单 pair 接口的耗时应使用实际调用 `distance()` 的基准测量。

### 单 pair / 四 pair / 原始 Edlib 比较

新的独立入口覆盖 K=0、2、5、10、20、30、40、50、60、70、80、90：

```bash
cmake --build build_clang --target benchmark_single_pair cache_control -j4
python3 scripts/run_single_pair.py results/my_single_pair_run
```

该脚本顺序执行各组，保存源码和二进制快照、输入、逐窗计时及计数器，并自动运行分析。
输出目录使用新路径；脚本默认每组 4096 窗、5 次重复，可通过 `--windows` 与 `--repeats` 指定。
正式五轮拟合结果包括 `summary.csv`、`fits_per_repeat.csv`、`fit_points.csv`、输入及每窗原始记录。

已保存的数据也可重新分析：

```bash
python3 scripts/analyze_single_pair.py results/my_single_pair_run
```

原始开发环境中的大体积测速输入和历史结果没有纳入源码仓库；结果目录说明见 [results/README.md](results/README.md)。
历史实验脚本也保留在 `scripts/`，其中 `run_supered.py`、旧 profiling/report 脚本可能依赖原始工作目录及历史输入位置。

引用性能结果时，应注明 K、输入分布、单 pair 或四 pair 接口、CPU、缓存条件和计时范围。四 pair 结果中的 ns/pair 已按实际 pair 数归一化。
