# Benchmark evidence — 2026-10-22

[English README](../../../README.md#performance) | [中文 README](../../../README.zh-CN.md#性能比较)

These files publish the measurement points and fitted summaries behind the README performance figures. They are copied from the completed experiment archives, without rerunning or changing the library. CSV/JSON numeric precision is retained; CSV line endings are normalized to LF. The full harness, large input matrices, executables, and per-window binary archives remain in the laboratory archives; this directory alone is not a standalone benchmark suite.

## Files

| File | Contents |
| --- | --- |
| [independent_summary.csv](independent_summary.csv) | All 96 workload/threshold/method summaries, with confidence intervals, throughput, passing fractions, and cache counters |
| [independent_fit_points.csv](independent_fit_points.csv) | 2,400 measured CPU-time points, fitted values, and residuals |
| [independent_fits.csv](independent_fits.csv) | 480 per-repeat fits |
| [anchor_summary.json](anchor_summary.json) | 200-million-pair experiment, validation results, timing boundaries, and input/output hashes |
| [anchor_fit_points.csv](anchor_fit_points.csv) | 15 measured CPU-time points and residuals |
| [anchor_fits.csv](anchor_fits.csv) | Three per-repeat fits |
| [provenance.json](provenance.json) | Compiler, guest CPU, upstream commits, production-header hashes, and evidence hashes |

## Independent-pair comparison

Each `near` or `mixed` group contains 655,360 independent 128 bp pairs. Queries are uniform random DNA. Near targets use edit budgets from 0 through K; substitution and balanced indel scripts are alternated where possible, with true distance checked by original Edlib. K=0 uses matching sequences; K=2 near inputs use substitutions only. Mixed groups contain half near targets and half independent uniform targets. At large K, some independent targets also pass: the CSV records actual passing fractions.

Methods are defined as follows:

- `edlib_original_hot`: original `myersCalcEditDistanceNW`, prepared Peq/target codes, fixed K, global distance only. Original block allocation is timed; the `edlibAlign()` wrapper is excluded.
- `edlib_stack_reference`: the same core recurrence with stack block storage. This is a modified reference and is distinct from original Edlib.
- `supered_batch4`: public `Context<K>::distance4()`, four independent query/target pairs and all four outputs consumed.
- `supered_single`: public `Context<K>::distance()`, called once per independent pair in a normal optimized loop. This is amortized call cost, not a timer around one isolated call.

Five repeats each fit `T=a+bN` across N=131,072/262,144/393,216/524,288/655,360. Each point comprises 4,096 windows of 32/64/96/128/160 pairs. Reported ns/pair is the mean slope; the confidence-interval half-width uses Student t with four degrees of freedom. Method/point order is randomized. The same logical pairs are used by every method.

Inputs are encoded and prepared before timing, copied to 24–120 KiB scratch batches, and touched. Context allocation is excluded. Timing includes input loads, SIMD setup, pair-dependent Eq/LEAP masks, fresh state/reset, recurrence, thresholding, results, and checksum work. In particular, preparation that must be repeated for each pair is not silently excluded.

All 15,728,640 input pairs matched the original Edlib full NW oracle; 98,304 additional scalar-DP checks passed. There were 9,830,400 accepted windows, 255,707 rejected attempts, and zero observed DRAM/IO demand fills in accepted windows. The minimum per-repeat R² was 0.99870235.

## Read/anchor experiment

Twenty independent uniform random anchors are retained in cache. Uniform random read candidates are accepted only if every distance to all 20 anchors is <=80. Of 10,044,397 candidates, 10,000,000 were accepted and 44,397 rejected. Original Edlib subsequently checked every one of the 200 million results; another 20,000 pairs were checked with scalar DP. Input and result SHA256 hashes are in `anchor_summary.json`.

The primary read-major layout calls `distance4(read, four_anchors)` five times per read. Packed anchors occupy 640 bytes. Read batches contain 16/32/48/64/80 sequences; fixed 125,000 windows produce 40/80/120/160/200 million comparisons per fit. Three repeats use the same complete accepted corpus. Preparation, read batch loading, and warm-up are outside timing; per-call setup, masks, recurrence, thresholds, integer result stores, and checksums are timed.

The mean fitted slope is 97.9131314617 ns/pair, with a t(2) 95% CI half-width of 0.2531141926 ns. The fitted 200-million-pair equivalent is 19.5826262923 seconds. Observed full compute trials are 19.883681589, 19.892829037, and 19.881434691 seconds; the positive fitted intercept explains why these exceed the slope-only equivalent. Full trial wall times are 20.782410387, 20.778723004, and 20.769829090 seconds.

A secondary anchor-major layout compares one anchor against four reads. Its one full compute trial took 19.898361893 seconds (99.491809465 ns/pair); it is not pooled into the primary three-repeat slope estimate.

Generation/filtering took 19.779052687 seconds and the full Edlib validation stage took 259.896208571 seconds. Both are outside the reported supered compute time. The validation-stage wall time has different boundaries and must not be used as an Edlib speed baseline.

## Hardware, cache gate, and interpretation

Both experiments ran on one compute thread pinned to logical CPU 24 in a VMware Linux guest reporting AMD EPYC 9654. They used Clang 22.1.2, C++20, `-O3`, `-march=native`, and LTO. Upstream identities and the production-header hashes are in `provenance.json`. The distributed production headers match the measured working-tree headers byte for byte.

An AMD PMU raw event (`0x4843`) monitors DRAM/IO demand fills. Windows with a nonzero observed count are recorded and retried with the same inputs; windows are not selected by elapsed time. PMU running/enabled fraction must be >=0.99. Hot/cold controls were retained in the original archives. L2 miss counts remain part of the audit; L2/L3 cache hits are permitted. The anchor experiment accepted 2,000,000 windows and recorded 40,988 rejected attempts.

These are cache-conditioned, single-thread compute measurements. They neither establish exclusive ownership of a physical core in the VM nor claim universal L1 residency. `supered_batch4` values are already divided by four useful pairs. K<5 executes individual LEAP calls; its speedup is not four-pair AVX-512 acceleration.

## 中文说明

这些文件发布 README 速度比较所依据的测量点与拟合结果；完整 harness、大型输入、可执行文件与逐窗口二进制记录保存在独立实验归档中。本目录不是独立测速工具。

`independent_*` 对应相同输入、相同阈值下的 Edlib/supered 对比，保留 near、mixed、单 pair、四 pair 与 stack reference 全部结果。`anchor_*` 对应一千万条 read 与二十个 anchor 的两亿次距离计算；每条 read 到全部 anchor 的真实距离均 <=80。

所有 ns/pair 均以实际 pair 数归一化，四 pair 数值不可再次除以四。编码、预处理、context 分配和批次入 cache 不计时；每 pair 必需的 masks、状态初始化与递推计时。斜率是 `T=a+bN` 中的 b，完整实测时间另列。gate 依据观测到的 DRAM/IO demand fills 拒绝并重测窗口，不按耗时筛选，允许 L2/L3 命中。Edlib 验证阶段的墙钟时间不用于速度比较。源码身份、CPU、工具链和证据文件哈希见 `provenance.json`。
