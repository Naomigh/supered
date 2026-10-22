# 第三轮：等长全局距离的 64-bit 移动带

接口为 `myers_banded64x4<K>(queries, targets)`，每条 DNA 仍为两个 128-bit 位平面，共 32 bytes。四个 qword 分别计算四个 pair。输入转置、query 窗口构造、Eq 与新鲜 DP 状态都在调用内计算，不增加预处理数据。

`K<64` 采用下面的移动带；`K>=64` 调用原有完整双 block 实现。`myers_2x64x4<K>` 本身保留原实现与接口。本轮没有把无界完整 Myers 的速度称为移动带的速度。

## 带宽为什么足够

两个长度都是 128 的序列，全局路径从 `(0,0)` 到 `(128,128)`。任何路径的 insertion 数与 deletion 数相等。若路径达到偏移 `d=i-j`，至少需要 `|d|` 次相应 indel 才能到达该偏移，随后还必须用 `|d|` 次反向 indel 回到偏移 0。因此总距离不超过 K 的路径满足：

```text
|i-j| <= floor(K/2)
```

K<=63 时最大偏移为 31。64-bit 窗口在列 j 使用 origin `r=j-32`，保存从行 r 到 r+64 的纵向差分；所有可能通过阈值的路径都在窗口边界以内。被省略的路径只可能产生超阈值结果。K=64 已可能需要偏移 32，因此选择完整双 block，不使用该移动带。

约束是 NW、全局、单位代价、两端长度均为 128。带宽推导不适用于任意不等长或半全局模式。本项目的输入契约已保证这些条件，生产代码不添加检查。

## 状态与递推

每列的 bit k 表示 `D(r+k+1,j)-D(r+k,j)`；PV 对应 +1，MV 对应 -1。初始 origin=-32，虚拟负行的列零代价为 `|i|`，所以：

```text
PV = ffffffff00000000
MV = 00000000ffffffff
S_origin = 32
```

query 位窗口以 32 个前导零开始；下一 word 和尾 word 在寄存器中拼接。每列用 VPSHRDQ 移动 query 窗口，必要的移位与拼装全部计时。

先对当前窗口计算普通 Myers 中间量：

```text
Xv = Eq | MV
Xh = (((Eq & PV) + PV) ^ PV) | Eq
Ph = MV | ~(Xh | PV)
Mh = PV & Xh
```

原本左移 Ph/Mh 后构造下一列 PV/MV，然后将窗口右移一行。两次相反方向的移位可以直接合并为：

```text
Xv_next = Xv >> 1
PV_next = Mh | ~(Xv_next | Ph)
MV_next = Xv_next & Ph
```

最高位不需要另加 OR：初始 `PV[63]=1, MV[63]=0`，因此 `Ph[63]=0`；`Xv_next[63]=0` 又保证新 PV 最高位仍是 1、新 MV 最高位仍是 0。该不变量让新窗口最下行自然由上一行付一次 deletion 延伸。生产路径没有每列边缘判断。

## origin 移动后的分数

窗口顶部的水平差分采用 +1。普通 Myers 更新在 bit 0 满足 `PV'=0, MV'=Xv[0]`，因此 origin 同时向右、向下移动时：

```text
S_next = S + 1 - (Xv & 1)
```

不需要每列进行整数 score 加减。用 VPSHRDQ 将 `Xv[0]` 收进 64-bit hits，每 64 列 popcount 一次。128 列后 origin 在 row 96；从该处到目标 row 128 还需 32 个纵向差分：

```text
D(128,128) = 160 - sum(hits)
             + popcount(PV & 0xffffffff)
             - popcount(MV & 0xffffffff)
```

最终结果才与 K 比较，返回精确距离或 -1。所有 pair 都执行完整 128 列，没有增加相等特判、提前拒绝或数据相关的循环次数。

## 不会因虚拟行产生假阳性

负行的初始代价是到 row 0 的距离。任何使用虚拟负行的路径，都可将其虚拟 diagonal 转为真实输入上的 insertion；移去对应初始代价后，转换后的代价不高于原路径代价。超出 query 末端的行不能逆向影响 row 128。边界上的 insertion/deletion 延伸也均对应合法付费路径。因此移动带不会产生低于真实 NW 距离的值；另一方面，所有 d<=K 的最优路径均完整保留，故阈值内距离与完整算法相等。

## 与 v2 的实现区别

v2 每列计算低、高两个 64-bit block，跨列交错调度，并在末尾对四个 PV/MV word popcount。新候选每列只计算一个 64-bit 带，不再传递两个 block 的 Ph/Mh 边界；代价是 query 每列需要滑动，并且必须用移动 origin 的公式恢复距离，不能套用完整列的 `128+popcount(PV)-popcount(MV)`。

初版显式 OR 最高 PV 位的代码留在 `bench/initial_band.hpp` 作对照。`profile_v3` 的 0 是 v2，1 是初版带，2 是消去 OR 的未展开版本，4 是四列展开版本。正式性能来自独立 cache 门控拟合，不使用这个热点循环代替正式结果。
