# Sign-free Hamiltonian catalog

**辅助场口径校正：新增21类中，仅R001–R005、R009共6类确认常规Gaussian-HS辅助场PfQMC正性；其余15类未确认。**
“verified”、21类总数及历史全部sign统计涵盖不同表示，不能当作21类辅助场PfQMC结果。
逐例标签和旧目录中的静态/非Gaussian换基/非厄米例外见 [确认状态](MODELS.md#auxiliary-field-status)。

[物理模型及简单数学](MODELS.md) · [原始笔记逐项对照](AUDIT.md) · [实测结果](RESULTS.md)

本目录把 `/home/han/notes/sign_free` 的剩余构造接入本仓库。排除单列的 Wei 两类及
interacting Kitaev chain；保留其余具体族及局域组合，并如实区分实辅助场、静态配置、
世界线和非厄米扩展。

持续搜索的新结果单独登记在 [MODELS.md §12](MODELS.md#12-持续搜索登记本轮新增构造2026-09-11)、
[搜索验算报告](SEARCH_RESULTS.md) 和 [计数登记](search_registry.json)。使用新目标
`sign_free_search`；本页下方 `sign_free_models` 的 27 个旧预设不计入新增结果阈值。
后续待验证线索见 [SEARCH_QUEUE.md](SEARCH_QUEUE.md)。

## Build and run

按仓库 [CONTEXT.md](../CONTEXT.md) 在 **WSL Ubuntu / han** 中激活 oneAPI 后：

```bash
cd /home/han/PfQMC
source /opt/intel/oneapi/setvars.sh
export PATH="$PWD/build/python-tools/cmake/data/bin:$PATH"
cmake -S . -B build/context-check -DCMAKE_CXX_COMPILER=mpiicpx \
  -DEIGEN3_INCLUDE_DIR="$PWD/build/eigen" \
  -DCMAKE_PREFIX_PATH=/home/han/miniconda3
cmake --build build/context-check -j 2
ctest --test-dir build/context-check --output-on-failure

build/context-check/sign_free_models --list
build/context-check/sign_free_models --model four_mode --samples 48 --sweeps 16 --seed 91011
build/context-check/sign_free_models --model sp_exterior_5_5 --samples 48 --sweeps 16 --seed 91011
build/context-check/sign_free_models --model checks --samples 48 --sweeps 16 --seed 91011
build/context-check/sign_free_models --model all --samples 48 --sweeps 16 --seed 91011 \
  --csv build/catalog_91011.csv
```

`--model MODEL` 接收 `--list` 的 27 个实 HS Hamiltonian 预设；
`checks` 运行静态 FK/BdG、五个非厄米扩展、pair hopping，以及正负/零权重控制；
`all` 两者都运行。`--samples N` 对每个普通模型产生 N 个独立 Gaussian 场配置，
场标准差轮流为 2、0.5、1，以覆盖不同幅度；这部分是压力验算，不是物理分布采样。
`--sweeps S` 则运行 S 次真实 Metropolis/PFQMC 左右交替扫描，记录初态和每次扫描后状态。
它是数值验算运行器，不把这些短链当作已热化的物理可观测量估计。
静态测试还各自穷举 3×3 的 512 个背景，并作 N×9 次单背景自旋翻转 Metropolis。
所有失败退出码非零，日志写 stderr，CSV 写 stdout 或指定文件。
`--samples>0`，`--sweeps>=2`。程序固定 MKL 单线程便于复现。

## 实现入口与使用方式

| 文件/API | 用途 |
| --- | --- |
| [catalog.h](../inc/catalog/catalog.h) | 公共 Hamiltonian、HS walker、静态及世界线 API |
| [catalog_models.cpp](../src/catalog_models.cpp) | 27 个预设、参数化模型工厂、非厄米生成元及静态条件 Hamiltonian |
| [catalog_core.cpp](../src/catalog_core.cpp) | Majorana 转换、精确 Fock oracle、带符号的 Gaussian 迹、通用连续场更新 |
| [catalog_worldline.cpp](../src/catalog_worldline.cpp) | 四次 pair hopping Hamiltonian、带占据投影的闭合路径 Pfaffian |
| [checks.cpp](checks.cpp) | 各模型独立权重对照、实际扫描、强脉冲和负/零权重测试 |
| [specials.cpp](specials.cpp) | 静态背景、非厄米三角结构、世界线求和与 Hamiltonian 核对 |

链接 CMake target `pfqmc_catalog` 即可使用；普通 HS 模型实际继承现有
`Spinless_tV`，传给已有 `PfQMC` 类，没有另造一个返回常数 sign 的模拟器：

```cpp
#include "catalog/catalog.h"
using namespace catalog;
rdGenerator rng(91011);
auto model = fiveOrbital(2, 2); // 或 adjointOdd(5), spinLattice(2,2,2), ...
Walker walker(model, 4, .9, &rng); // L_tau, beta
PfQMC simulation(&walker, 5);   // 每 5 个算符做稳定化
simulation.rightSweep();
simulation.leftSweep();
auto sampledFields = walker.fields();
auto signedTrace = traceProduct(model.pulses(4, .9, sampledFields));
// signedTrace 和 simulation.sign 的相位应一致。
```

`Hamiltonian` 存储 $H=-Q(q0)-\sum_C\lambda_C Q(a_C)^2/2+E_0$。
可以修改公开的 `q0`、`channels`、`scalarEnergy` 构造笔记允许的参数和图，但
`validate()` 只检查代数形式、维数、厄米性和非负 HS 耦合；**它不是自动证明任意新矩阵无符号的分类器**。
必须自己保持 [MODELS.md](MODELS.md) 的共同场和通道限制。
`exactHamiltonian()` 给出至多 10 复模的完整 Fock 矩阵，仅作小系统核对。
外幂构造使用显式稠密矩阵，集合枚举限辅助维数 d<=16，实际内存还取决于组合数。

静态、非厄米与世界线 API 的最小调用分别是：

```cpp
StaticModel fk; // 3x3，f 背景不沿时间独立变化
std::vector<int> s(9, 1);
auto wStatic = std::exp(fk.logScalar(s, .9))*traceProduct({.9*fk.generator(s)});

auto generator = makeGenerator("creation_only_spin");
std::vector<double> fields(generator.fieldCount(3), .4);
auto wTriangular = traceProduct(generator.pulses(3, .9, fields));

PairHopping pairs;
auto physicalH = pairs.exactHamiltonian();
auto wPath = pairs.pathPfaffian({3u,12u,192u,48u,3u}, .11);
```

## 符号如何计算，如何排除伪验证

`traceProduct` 使用单个 Gaussian 的归一化迹 $\eta=\operatorname{Tr}O/2^N$、
协方差 $G=2(I+B)^{-1}-I$，以及**有符号 Pfaffian**的乘积公式。
普通乘积用 UDT 重建协方差；遇到中间零迹等病态前缀时，退回整个时空反对称矩阵。
单脉冲零迹会显式报错并要求先拆分脉冲；厄米 HS 因子不遇到这个问题。
一般复脉冲另用矩阵指数，不能错误地送进厄米本征求解器。
这些公式见原始 [PfQMC 论文 Eqs. (4)–(6)](https://arxiv.org/html/2408.10311v1)。

`FieldOperator` 对标准正态场作独立正态提议，提议比正好消去 Gaussian 先验。
接受率使用实际 Pfaffian 权重比的绝对值，另外追踪其复相位；接受后解析更新 Green 矩阵。
单个厄米指数的 `getSignOfWeight()` 返回 1 是因为它的迹严格为正，**不等于整个乘积的 sign 返回 1**。
全乘积负号由 Pfaffian 计算。测试明确覆盖全乘积负权重、局部接受/拒绝决策以及两个扫描方向。

独立对照包括：小系统完整 Fock 指数；守粒子数模型的中心化实行列式（保留符号）；
伴随族的小矩阵字符；配对阵列的独立单胞 Fock 迹乘积；静态单个厄米 H 的热谱公式；
pair hopping 原始四次费米子矩阵与小系统全部离散路径求和。
没有用 $+\sqrt{\det(I+B)}$ 选择权重符号。

## 结果文件和适用范围

[results/](results/) 中保留两种子的完整逐配置 CSV、运行摘要和校验和。
CSV 的 `sign_real/sign_imag` 由 $Z/|Z|$ 得到；明确的零点记为 `nan`。
`negative` 和 `strong_commutator` 中有注明的负权重反例，不是模型验算失败。
一般复矩阵 oracle 的 `control` 行允许非实相位，不能混进平均 sign。
`attempts/accepted` 只用于真正的 Metropolis/PFQMC 行。

数值比较用 $|Z-Z_{\rm ref}|/\max(1,|Z_{\rm ref}|)$；由于 132 模的迹可达 $10^{43}$，
是否实数应该看相位误差，不应只看原始虚部的绝对值。
最大允差为 $2\times10^{-7}$，实际误差见报告。
这些是有限温度、短虚时间、有限模型尺寸的数值检查；通用 dense walker 主要为了验证模型，
没有承诺生产级大规模吞吐量或任意长虚时间稳定性。
时空 Pfaffian 分支限制矩阵维数不超过 8192，纯 Fock oracle 限 10 模。
零权重允许的模型不能把 sign“恒为 1”理解成零点上也有定义。


## 持续搜索中的专用分解

新增构造登记在 [MODELS.md §12](MODELS.md)，当前数量从 [search_registry.json](search_registry.json) 读取。
[SEARCH_RESULTS.md](SEARCH_RESULTS.md) 保留每轮独立数据；不要重跑第一轮历史汇总脚本覆盖后续注册表。
`sign_free_search` 运行 Gaussian 族及真实 PfQMC sweep。
`sign_free_projected` 通过 signed Pfaffian Fourier 和得到守恒扇区权重，
`sign_free_dynamic_ising` 计算完整 CT-INT／自旋路径配置权重；后两者没有生产级更新器。
源码 API 分别为 [search_models.h](../inc/catalog/search_models.h)、
[projected_models.h](../inc/catalog/projected_models.h)、[dynamic_ising.h](../inc/catalog/dynamic_ising.h)。

`sign_free_poisson` 实现局域高阶指数顶点的完整 CT 权重；
`sign_free_albert` 实现 27 轨道 Albert 单胞及邻胞 HS，并运行真实 PfQMC sweep。
API 为 [poisson_models.h](../inc/catalog/poisson_models.h) 与 [albert_models.h](../inc/catalog/albert_models.h)。

`sign_free_mobile_spin` 与 `sign_free_quantum_link` 验证动态pair和一维quantum-link路径，
API见 [mobile_spin.h](../inc/catalog/mobile_spin.h)、[quantum_link.h](../inc/catalog/quantum_link.h)。

第五轮专用CT验证器为 `sign_free_meron`、`sign_free_ring_exchange`、
`sign_free_color_projector`、`sign_free_frustrated_impurity`；API分别见
[meron_models.h](../inc/catalog/meron_models.h)、[ring_exchange.h](../inc/catalog/ring_exchange.h)、
[color_projector.h](../inc/catalog/color_projector.h)、[frustrated_impurity.h](../inc/catalog/frustrated_impurity.h)。

第六轮专用CT接口：

- `stabilizer`：[API](../inc/catalog/stabilizer_models.h)，验证器 `sign_free_stabilizer`。
- `jq`：[API](../inc/catalog/jq_models.h)，验证器 `sign_free_jq`。
- `bilayer`：[API](../inc/catalog/bilayer_models.h)，验证器 `sign_free_bilayer`。
- `zigzag`：[API](../inc/catalog/zigzag_models.h)，验证器 `sign_free_zigzag`。
- `tj_chain`：[API](../inc/catalog/tj_chain.h)，验证器 `sign_free_tj_chain`。
- `fredkin`：[API](../inc/catalog/fredkin_models.h)，验证器 `sign_free_fredkin`。


`sign_free_pxp`提供完整Fock的Rydberg条件翻转验证，[API](../inc/catalog/pxp_models.h)。
本轮最终21个新构造已严格去重，结果与范围见[最终汇总](SEARCH_RESULTS.md#最终汇总)。
在已激活oneAPI的WSL中，运行 `python sign_free_catalog/summarize_search_all.py` 可只读核对七轮统计，
不会覆盖注册表。原始CSV和逐次运行日志保留在本地，但已被 `.gitignore` 排除，不随 Git clone 提供。
汇总 JSON、CTest 记录、校验和及小型独立反例证据继续保留；复现各轮原始样本或从归档恢复后，
才能重新运行汇总脚本。各轮源码快照记录当时版本，不因后续文档或代码编辑而重写。
