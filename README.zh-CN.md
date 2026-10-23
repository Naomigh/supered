# supered

**面向 128 bp DNA 的快速、精确编辑距离计算库。**

[English](README.md) | **简体中文**

supered 是一个 C++ 全局单位代价 Levenshtein 距离计算库，直接接收预处理后的 DNA 输入。它将四 pair AVX-512 Myers 内核与适用于极小编辑距离阈值的修正版 LEAP 结合，统一使用每条序列 **32 字节**的输入格式。

```cpp
// 输入已完成编码；在多次调用之间复用 context。
supered::Context<80> engine;
auto distances = engine.distance4(
    read, {&anchor0, &anchor1, &anchor2, &anchor3});
// distances[i]：距离 <= 80 时返回精确距离，超过 80 时返回 -1。
```

在 **1,000 万条 read × 20 个 anchor** 的单线程、cache 内实验中，supered 达到 **97.91 ns/pair**，即 **每秒 1,021 万 pair**。[性能比较](#性能比较)说明了数据分布、计时边界，以及与原始 Edlib 的对比。

## 特性

- **精确距离：**全局比对；插入、删除、替换的代价均为 1，匹配为 0。
- **紧凑输入：**每条 128 bp 序列使用两个 128-bit 位平面，共 32 字节，两种后端共享输入。
- **四 pair SIMD：**Myers 每次并行计算四个独立 pair，支持一条 read 同时比较四个 anchor。
- **按阈值特化：**编译期 `Context<K>` 在 `K < 5` 时使用修正版 LEAP，`5 <= K < 64` 时使用移动带 Myers，更大阈值使用完整双 block Myers。
- **复用工作空间：**提前准备输入并创建 context，在计算循环中反复使用。
- **CMake 集成：**链接 `supered::supered`；仓库已包含所需修正版 LEAP 源码。

公开接口面向固定 **128 bp、大写 A/C/G/T** 输入，返回距离；提供全局距离计算，不输出 traceback、比对坐标或 CIGAR。

## 目录

- [构建](#构建)
- [在项目中使用](#在项目中使用)
- [使用示例](#使用示例)
- [输入格式](#输入格式)
- [API](#api)
- [性能比较](#性能比较)
- [实现原理](#实现原理)
- [正确性验证](#正确性验证)
- [致谢与源码来源](#致谢与源码来源)

## 构建

需要 Linux/x86-64、CMake 3.20 或更新版本，以及支持 C++20 的编译器。已验证的工具链是 Clang 22.1.2 与 LLVM lld。当前内核要求 AVX2、BMI 指令，以及包括 **F、VL、VBMI2、VPOPCNTDQ** 在内的 AVX-512 扩展。

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

内置示例的输出为：

```text
0 1 2 -1
0 1 2 3
0 1 2 3
```

三行依次使用 `Context<2>`、`Context<5>` 和 `Context<128>`。最后一对序列的距离为 3，超过第一个 context 的阈值。

仅构建库时，配置中加入 `-DSUPERED_BUILD_EXAMPLES=OFF`，然后构建 `--target supered`，生成 `build/libsupered.a`。普通 clone 即包含所需依赖，无需初始化 submodule。

构建使用 `-march=native` 与链接时优化。构建和部署机器需支持兼容的 CPU 指令集；链接静态库时使用兼容的编译器和链接器。CPU 指令路径在构建时确定。

## 在项目中使用

将仓库克隆至 `third_party/supered`，然后链接其 CMake target：

```cmake
cmake_minimum_required(VERSION 3.20)
project(my_alignment_tool LANGUAGES CXX)

set(CMAKE_INTERPROCEDURAL_OPTIMIZATION TRUE)
add_subdirectory(third_party/supered)

add_executable(my_alignment_tool main.cpp)
target_link_libraries(my_alignment_tool PRIVATE supered::supered)
```

应用程序使用 Release 配置，以及上文相同的 Clang/lld 选项。target 会传递头文件目录、C++20 要求和 `-march=native`。通过子目录引入时，默认不构建示例程序。

公开头文件为 [`<supered/supered.hpp>`](include/supered/supered.hpp)。Myers 内核定义在头文件中，静态库提供修正版 LEAP 后端；分发已构建的依赖时应同时保留两者。

## 使用示例

以下完整程序计算一条 read 与四个 anchor 的距离。编码和 context 构建位于距离调用之前；已有预处理输入的应用可以直接传入 packed 数据。

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

保存为 `example.cpp`，构建库后从仓库根目录编译运行：

```bash
clang++-22 -O3 -std=c++20 -march=native -flto=thin -fuse-ld=lld \
  -Iinclude -Ileap_origin example.cpp build/libsupered.a -o example
./example
```

输出为 `0 1 2 3`。

计算四个独立 pair 时，传入两组指针：

```cpp
const auto distances = engine.distance4(
    {&query0, &query1, &query2, &query3},
    {&target0, &target1, &target2, &target3});
```

`distances[i]` 对应 `query[i]` 与 `target[i]`。输入无需连续存储，调用期间保持只读。计算一条 read 与 20 个连续 anchor 的距离时，提前准备 `int output[20]`，调用 `engine.distance_many(read, anchors, 20, output)`；内部执行五次四 pair 调用。

### 选择阈值

`K` 是 0 至 128 的编译期整数。每个结果均满足：**真实距离 `d <= K` 时返回精确值 `d`，否则返回 `-1`**。`Context<>` 等价于 `Context<128>`，对任意合法输入 pair 返回完整距离。

`Context<2>` 接受距离 0–2，`Context<50>` 接受距离 0–50；距离未知时使用 `Context<>`。后端按给定阈值选择：即使某个 pair 的真实距离为 1，`Context<50>` 仍使用 Myers。

每个 worker 创建并复用自己的 context。Context 不可复制或移动，通过引用传递。库在调用线程内计算，由应用负责并行调度与输出存储。

## 输入格式

[`PackedDNA128`](include/strong/packed.hpp) 的大小为 32 字节，对齐要求为 32 字节：

```cpp
struct alignas(32) PackedDNA128 {
    uint64_t bit0[2];
    uint64_t bit1[2];
};
```

A/C/G/T 分别编码为 `00/01/10/11`。`bit0` 保存每个碱基的低位，`bit1` 保存高位。第 `i` 个碱基位于第 `i / 64` 个 word 的第 `i % 64` 位；第 0 个碱基位于最低位。内存排列顺序为 `bit0[0], bit0[1], bit1[0], bit1[1]`。

使用 `encode(DNA)` 生成两个位平面，`decode(packed)` 恢复 128 个字符。`DNA` 是 `std::array<char, 128>`，不含字符串终止符。read/index 预处理若已生成上述布局，可以直接使用满足对齐和对象生命周期要求的 `PackedDNA128`。每个碱基占相邻两位的常规压缩格式需要先转换为双位平面格式。

调用方保证输入为合法 128 bp A/C/G/T，指针有效，packed 对象正确对齐，输出空间充足。热路径依赖这些前提，不增加运行时输入校验。数组批量接口的数量必须为 4 的倍数。

## API

公开接口摘要如下：

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

| 方法 | 计算对象 | 输出 |
| --- | --- | --- |
| `distance4(queries, targets)` | 四个逐项对应的 pair | 按 lane 顺序返回四个距离 |
| `distance4(query, targets)` | 一条 query 与四条 target | 按 target 顺序返回四个距离 |
| `distance(query, target)` | 一个 pair | 一个距离 |
| `distance_many(query, targets, count, output)` | 一条 query 与连续 target 数组 | `output[i] = distance(query, targets[i])` |
| `distance_pairs(queries, targets, count, output)` | 两个连续数组逐项比较 | `output[i] = distance(queries[i], targets[i])` |

所有距离采用相同阈值语义。两个数组接口要求 `count` 是 4 的倍数，且 `output` 至少能容纳 `count` 个整数。`broadcast(packed)` 构造四个指向同一已有 packed 对象的指针。

`K >= 5` 时，`distance4()` 并行计算四个有效 pair。便捷接口 `distance()` 将一个 pair 广播到 Myers 各 lane，再返回 lane 0；追求吞吐时应使用批量接口。`K < 5` 时，`distance4()` 执行四次独立 LEAP 调用，因此小阈值数据反映的是 LEAP 性能，不属于四 pair AVX-512 并行结果。

## 性能比较

### 与原始 Edlib 对比

以下为 2026-10-22 的实验结果，两种方法使用**相同的独立输入 pair 和阈值**。单位为 **ns/pair**，`±` 表示五次拟合斜率的 95% 置信区间半宽。supered 一列已经按四个有效 pair 归一化。

| `K` | 原始 Edlib | supered 四 pair 接口 | 加速比 |
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

表中每个阈值使用 655,360 个独立 128 bp pair，真实距离均 `<= K`。query 为随机 DNA；target 包含替换，并在编辑预算允许时交替包含平衡的插入/删除。`K=0` 的 query 与 target 相同；`K=2` 的近似输入仅含替换。`K=0/2` 使用 LEAP，其余列出的阈值使用四 pair Myers。

原始 Edlib 指**未经修改的全局 NW 距离计算内核**，已提前准备 equality masks 和 target codes；原实现的 block 分配仍计入调用时间。该比较不包含公开 `edlibAlign()` 包装层。supered 通过公开 `distance4()` 接口测量。[完整数据](docs/benchmarks/2026-10-22/independent_summary.csv)还包括近似/随机混合输入、单 pair 调用，以及单独标注的采用栈上 block 存储的 Edlib 参考版本。

每轮对 131,072、262,144、393,216、524,288、655,360 个 pair 拟合 `T(N) = a + bN`，以斜率 `b` 表示每 pair 计算成本。输入提前准备，按适合 cache 的批次复制到 scratch 并预热。编码、Edlib 输入预处理、context 分配、批次加载和预热不计时；输入读取、每次调用的 SIMD 设置、依赖 pair 的 masks、状态初始化、递推、阈值判断和结果 checksum 均计时。

### Read × anchor 模拟

实验使用 **10,000,000 条 read 与 20 个 anchor**，长度均为 128 bp，通过 `Context<80>` 计算 **200,000,000 个距离**。每条保留的 read 与**全部** anchor 的真实距离均 `<= 80`；所有结果均独立与原始 Edlib 核对。

- **97.91 ± 0.25 ns/pair**：40–200 million 次比较、三轮斜率拟合。
- **391.65 ns/四 pair 调用**，相当于 **每秒 1,021 万 pair**。
- 2 亿次比较的拟合计算时间为 **19.58 s**；完整实测计算时间为 **19.88–19.89 s**。
- 包括批次搬运、检查、cache gate 重试的完整 trial 墙钟时间为 **20.77–20.78 s**。

20 个 packed anchor 共占 640 字节。预处理后的 read 以 16–80 条为一批，在计时外加载到 scratch；每条 read 对常驻 anchor 执行五次四 pair 调用。输入生成、拒绝采样、oracle 验证均不计入 supered 计算时间。该模拟与上面的 Edlib 对比使用不同工作负载。

### 机器与 cache 条件

两个实验均在 Linux VMware guest 中使用一个计算线程，绑定逻辑 CPU 24；guest 报告 CPU 为 AMD EPYC 9654。工具链为 Clang 22.1.2、C++20、`-O3`、`-march=native` 与 LTO。这些结果为单线程测量；虚拟机环境不能证明独占物理核心。

PMU gate 接受观测到 DRAM/IO demand fills 为零的窗口，否则重测同一输入。允许命中 L2/L3，并保留 L2 miss 记录。因此结果代表上述 cache 条件下的计算，不表示每次读取都命中 L1。[测速说明、源码身份、拟合点及汇总数据](docs/benchmarks/2026-10-22/README.md)提供了进一步信息。

## 实现原理

supered 为固定 128 bp 全局比对特化 Myers 位向量动态规划。双位平面输入与运行中的 Myers `Pv/Mv` 状态是不同概念；距离计算内部由两个位平面构造 equality masks。

`5 <= K < 64` 使用 64-bit 移动带，依据等长全局比对的 `|i-j| <= floor(K/2)` 边界裁剪。更大阈值使用两个 64-bit block，显式传播跨 block 的进位与移位。四个独立 pair 占据四个 64-bit lane；包括 AVX-512VL 支持的 256-bit 形式在内的 AVX-512 指令并行更新这些 lane。实现使用 intrinsic，不要求每个操作都使用 512-bit 寄存器。

`K < 5` 使用修正版 AVX2 LEAP，直接接收同一 packed 输入并复用 context 工作空间。后端由编译期阈值特化选择。

## 正确性验证

独立 pair 实验将 **15,728,640 个输入 pair** 与原始 Edlib 核对，另用独立标量动态规划验证 **98,304 个 pair**。read/anchor 实验将全部 **2 亿个结果**与原始 Edlib 核对，另有 **20,000 个标量 DP 验证**。每个计时窗口还检查结果 checksum。这些数字来自已保存的实验，并非文档更新期间新执行的测试。

发布仓库包含可运行的 [`examples/basic.cpp`](examples/basic.cpp)。[公开测速数据](docs/benchmarks/2026-10-22/README.md)提供支持本文数值的源码身份和测量点；完整实验 harness 与逐窗口二进制记录独立保存在实验归档中。

## 致谢与源码来源

supered 基于 [Myers 位向量算法](https://doi.org/10.1145/316542.316550)、[Edlib](https://github.com/Martinsos/edlib) 的工程方法，以及[修正版 LEAP 分支](https://github.com/YourongZhou/LEAP/tree/codex/leap-correctness-fixes)。

所包含的 LEAP 子集固定于提交 [`6a711282`](https://github.com/YourongZhou/LEAP/commit/6a711282f3494b6b9886236e58c7dfba6efa885a)。原始源码声明均已保留；[NOTICE](NOTICE) 与 [`leap_origin/UPSTREAM.json`](leap_origin/UPSTREAM.json)记录其来源。本文档不为第三方源码重新指定许可证。
