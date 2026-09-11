# 持续搜索数值记录

**辅助场口径校正：新增21类中，仅R001–R005、R009共6类确认常规Gaussian-HS辅助场PfQMC正性；其余15类未确认。**
“verified”、21类总数及历史全部sign统计涵盖不同表示，不能当作21类辅助场PfQMC结果。
逐例标签和旧目录中的静态/非Gaussian换基/非厄米例外见 [确认状态](MODELS.md#auxiliary-field-status)。

原始逐样本 CSV 和运行日志是本地生成产物，已由 `.gitignore` 排除；下文相应链接仅在本地数据存在时可用。
Git 保留汇总、CTest 记录和校验和；新 clone 需按复现命令生成或恢复原始数据后才能重算统计。

本页对应 MODELS.md §12 的持续搜索，独立于前面的旧目录报告。
**最终计数：21 个独立 Hamiltonian 构造（R007/R010合并）。已超过20，停止新增并完成最终去重。**

## 第一轮验算（历史记录）

- WSL Ubuntu、Intel oneAPI/MKL、单线程；两个种子 9112026 / 9112027。
- 11 个单胞／阵列预设，每种子每预设 128 个独立实场配置和 16 次真实 PFQMC 左右交替扫描（含初态17个状态）。
- 共 3190 条应为正的配置全部通过 phase=1 检查；另有 4 条预期负号、2 条一般复相位、2 条精确零点控制。
- 全记录最大权重相对误差 `8.00169e-14`；正配置最大 `|sign-1|=1.19663e-15`。
- 完整 CTest 58/58 通过，包含之前 46 项及本轮 12 项。
- 数值样本不替代全参数证明；每个模型的证明、局域性及允许场限制直接写在 [MODELS.md](MODELS.md) §12。

| 验证预设 | 复模数 | 独立场配置 | 扫描状态 | 最大相对误差 |
| --- | ---: | ---: | ---: | ---: |
| `conjugate_tensor2` | 4 | 256 | 34 | 1.41e-14 |
| `conjugate_tensor3` | 9 | 256 | 34 | 3.05e-14 |
| `conjugate_tensor3_lattice` | 18 | 256 | 34 | 5.99e-14 |
| `lorentz4` | 4 | 256 | 34 | 9.46e-15 |
| `lorentz4_lattice` | 16 | 256 | 34 | 2.98e-14 |
| `lorentz5` | 5 | 256 | 34 | 1.04e-14 |
| `lorentz8` | 8 | 256 | 34 | 1.37e-14 |
| `symmetric_power3_2` | 6 | 256 | 34 | 7.1e-15 |
| `symmetric_power3_2_lattice` | 24 | 256 | 34 | 2.29e-14 |
| `symmetric_power3_3` | 10 | 256 | 34 | 2.25e-14 |
| `symmetric_power4_2` | 10 | 256 | 34 | 1.23e-14 |

## 实际检查覆盖

- 权重用现有 PFAPACK 的 signed Pfaffian；逐配置与中心化行列式和辅助小矩阵谱公式交叉检查。至多6复模的独立场记录还比较完整 Fock 指数。
- Walker 使用已有 PfQMC 稳定化和真实 Metropolis 更新；逐扫描核对 tracked sign、getSignRaw 和独立 Green 矩阵，并检查有实际提议与接受。
- Tensor d=2 与 Lorentz4 用 Pauli 基逐生成元验证酉等价，不计第二个结果。
- 共同 Wei 排除使用实际允许方向的简单实 Majorana 谱证书；推导针对指定 HS 表示，不是对全部可能 HS 的排除。
- 独立叶势破坏 Lorentz 正性：精确 -2800/81。
- 放开 tensor 共轭相关性：四模权重约 20.473942-1.254887i。
- 放开 symmetric-power 的相关Cartan限制：六模权重约 -7.14635870944e10；高精度参考、完整Fock与signed Pfaffian一致。
- Lorentz/tensor 真正零点单独记nan，不被算成sign1。

## 复现

按 [CONTEXT.md](../CONTEXT.md) 激活WSL oneAPI后：

```bash
build/python-tools/cmake/data/bin/cmake -S . -B build/context-check
build/python-tools/cmake/data/bin/cmake --build build/context-check --target sign_free_search -j 2
build/context-check/sign_free_search --model all --samples 128 --sweeps 16 --seed 9112026 --csv build/search_9112026.csv
build/context-check/sign_free_search --model all --samples 128 --sweeps 16 --seed 9112027 --csv build/search_9112027.csv
```

也可将 `--model all` 改成上表某一预设或 `controls`。默认 beta=.8、L_tau=3；独立场标准差轮换1.5、.5、1。
源代码：[search_models.h](../inc/catalog/search_models.h)、[catalog_search.cpp](../src/catalog_search.cpp)、[search.cpp](search.cpp)。
局域阵列的胞粒子数仍守恒，没有声称实现跨胞单费米子巡游。
这是验证器与通用dense walker，不据短扫描报告物理热化、自相关或性能提升。

- [种子9112026 CSV](search_results/round1_9112026.csv) / [日志](search_results/round1_9112026.log)
- [种子9112027 CSV](search_results/round1_9112027.csv) / [日志](search_results/round1_9112027.log)
- [CTest](search_results/ctest_round1.log) / [数据SHA-256](search_results/round1_data.sha256) / [源代码快照SHA-256](search_results/round1_sources.sha256)
- [计数与合并记录](search_registry.json) / [后续研究队列](SEARCH_QUEUE.md)

本轮是进展，不满足“超过20后停止并最终合并”的完成条件。后续构造必须保留同样的局域性、正号数值和重复性检查。

## 第二轮验算（2026-09-11）

新增 R004–R007 四项；四元数两个分支合并，充电／奇偶扇区合并。
九个 Gaussian 预设各使用两个种子 9112201/9112202，每种子 64 个独立场、8 次双向交替 sweep，
含初态合计 73 条；另重跑完整负／复相位／零控制。
投影验证种子 731/1907 各 32 samples；动态 Ising 种子 160507420/91012 各 64 samples。

| 记录分类 | 数量 |
| --- | ---: |
| 非零正权重（phase 检查通过） | 6074 |
| 零权重或硬限制排除 | 294 |
| 预期负号对照 | 24 |
| 一般复相位对照 | 2 |
| 合计 | 6394 |

最大权重相对误差 `6.48131e-13`，正权重最大 `|sign-1|=1.48332e-11`。
完整 [CTest 69/69 通过](search_results/ctest_round2.log)。
逐文件统计见 [round2_summary.json](search_results/round2_summary.json)，
原始 24 份 CSV／日志均为 [search_results/](search_results/) 中的 `round2_` 前缀。
[汇总脚本](summarize_search_round2.py) 只读 CSV、检查记录并生成统计，不修改搜索计数。

在已激活 oneAPI 的 WSL 中运行：

```bash
build/python-tools/cmake/data/bin/cmake -S . -B build/context-check
build/python-tools/cmake/data/bin/cmake --build build/context-check -j 2
build/context-check/sign_free_search --model quaternion_herm3 --samples 64 --sweeps 8 --seed 9112201 --csv build/quaternion_check.csv
build/context-check/sign_free_search --model lorentz_sym5 --samples 64 --sweeps 8 --seed 9112201 --csv build/lorentz_sym_check.csv
build/context-check/sign_free_projected --samples 32 --seed 731 --csv build/projected_check.csv
build/context-check/sign_free_dynamic_ising --samples 64 --seed 160507420 --csv build/ising_check.csv
```

另外六个四元数／三个 Lorentz-symmetric 预设完整名字见注册表，不能把它们算成九个新结果。
当前 `--model all` 会包括第二轮，故不再精确重现第一轮预设数量；历史文件不被覆盖。
Projected 和 DynamicIsing 是配置权重 API，不含生产 Monte Carlo 更新器；Gaussian 族使用实际 PfQMC sweep。
相对误差的分母：Gaussian 为 `max(1,abs(reference))`，两个专用验证器为 `1+abs(reference)`；
汇总取所记录的值。所有检查都保留 signed Pfaffian；正表达式仅作为独立 oracle。
[本轮数据校验和](search_results/round2_data.sha256)与[当时源码校验和](search_results/round2_sources.sha256)
用于保留证据快照；后续源码变更不会重写历史校验和。

## 第三轮验算（2026-09-11）

R008 Poisson 高阶顶点两种子 731/1907 各 512 samples，并枚举两胞全部 585 个 k<=3 顶点词；
R009 Albert 两种子 9120401/9120402 各 32 独立场配置/预设和 8 次真实双向 sweep。
Albert 覆盖单胞 27 模/27 channels 及邻胞 54 模/56 channels，beta=.45、两时间片。

共 3396 行：3386 正、4 零、6 预期负对照。最大记录相对误差 `3.85974e-13`，
最大正配置相位偏差 `4.19738e-16`。Poisson 使用 `abs(reference)` 作非零权重误差分母，
Albert 使用 `max(1,abs(reference))`，汇总未混淆这两种口径。
[完整 CTest 71/71 通过](search_results/ctest_round3.log)。
[逐文件统计](search_results/round3_summary.json)由[只读汇总脚本](summarize_search_round3.py)生成。
四份 CSV／日志为 `search_results/round3_{poisson,albert}_{seed}`。

```bash
build/context-check/sign_free_poisson --samples 512 --seed 731 --csv build/poisson_check.csv
build/context-check/sign_free_albert --model all --samples 32 --sweeps 8 --seed 9120401 --csv build/albert_check.csv
```

Poisson 的 oracle 是有限多项式与完整小 Fock，不把正的 Gaussian trace 单独当作完整权重。
Albert 用独立 Jordan quadratic representation 和保持三次形式检查，没有枚举 2^27 维 Fock 空间。
两者都用负控制检查 signed Pfaffian。群论证明和有限样本证据的范围见 MODELS.md R008/R009。
[第三轮数据 SHA-256](search_results/round3_data.sha256) /
[第三轮源码快照 SHA-256](search_results/round3_sources.sha256)。

## 第四轮验算（2026-09-11）

R010 MobileSpin 种子 160507421/91013，每尺寸每种子64随机样本；
R011 QuantumLink 种子230508917/230508918，开链/三点环/四点环各256样本。
量子链还核对非回文闭词、闭边带两个原始负 hopping 的词，以及三点环全部512态的物理转移相位证书。
三叉图和独立link flip明确作为负对照，未纳入sign-free域。

共1914条：1858正、46零、10预期负；最大相对误差 `1.99840e-15`，最大正相位误差 `4.48027e-15`。
[CTest73/73通过](search_results/ctest_round4.log)；[统计](search_results/round4_summary.json)由
[summarize_search_round4.py](summarize_search_round4.py)生成。数据和日志为 `search_results/round4_*.csv/.log`。
MobileSpin误差分母 `1+abs(reference)`，QuantumLink为 `max(1,abs(reference))`。
两个新模型都是完整CT配置权重API，不含生产更新器。

```bash
build/context-check/sign_free_mobile_spin --samples 64 --seed 160507421 --csv build/mobile_check.csv
build/context-check/sign_free_quantum_link --samples 256 --seed 230508917 --csv build/link_check.csv
```

[数据校验和](search_results/round4_data.sha256) / [当时源码校验和](search_results/round4_sources.sha256)。

## 第五轮验算（2026-09-11）

R012 Meron种子731/1907，每个64 samples；R013 RingExchange种子409112/409113，每个256 samples；
R014 ColorProjector种子240602926/240602927，每个256 samples；R015 FrustratedImpurity种子160805144/91015，每个64 samples/长度。
八份CSV共3454条：3146正、296零、12预期负；最大记录相对误差 `3.12312e-12`，
正权重最大相位偏差 `7.70770e-15`。[完整CTest77/77通过](search_results/ctest_round5.log)。
[逐文件统计](search_results/round5_summary.json)由[只读汇总脚本](summarize_search_round5.py)生成。
Meron使用非零reference的绝对值归一化误差，其余三者用max(1,abs(reference))。
四个接口均为CT配置权重验证器；Meron含实际signed路径重求和，其余直接计算完整闭词Pfaffian。
原始数据和日志为 `search_results/round5_*.csv/.log`。

```bash
build/context-check/sign_free_meron --samples 64 --seed 731 --csv build/meron_check.csv
build/context-check/sign_free_ring_exchange --samples 256 --seed 409112 --csv build/ring_check.csv
build/context-check/sign_free_color_projector --samples 256 --seed 240602926 --csv build/color_check.csv
build/context-check/sign_free_frustrated_impurity --samples 64 --seed 160805144 --csv build/impurity_check.csv
```

[数据校验和](search_results/round5_data.sha256) / [当时源码校验和](search_results/round5_sources.sha256)。

## 第六轮验算（2026-09-11）

新增R016–R021六项；严格排重将R007 Ising并入R010 XXZ，故累计是20而非21，尚不能停止。
Stabilizer种子731/1907各64 samples；J-Q 611343/611344各256；Bilayer 731/1907各64；
Zigzag 9707019/91017各64；t-J 13106496/13106497各128；Fredkin 160503842/160503843各256。
十二份CSV共15192条：14287正、869零、34预期负、2复相位对照。
最大记录相对误差 `5.18640e-14`，正相位最大偏差 `1.14353e-14`。
[完整CTest83/83通过](search_results/ctest_round6.log)，[统计](search_results/round6_summary.json)
由[只读汇总脚本](summarize_search_round6.py)生成。全部专用接口仍是配置权重API，无生产更新器。

```bash
build/context-check/sign_free_stabilizer --samples 64 --seed 731 --csv build/stabilizer_check.csv
build/context-check/sign_free_jq --samples 256 --seed 611343 --csv build/jq_check.csv
build/context-check/sign_free_bilayer --samples 64 --seed 731 --csv build/bilayer_check.csv
build/context-check/sign_free_zigzag --samples 64 --seed 9707019 --csv build/zigzag_check.csv
build/context-check/sign_free_tj_chain --samples 128 --seed 13106496 --csv build/tj_check.csv
build/context-check/sign_free_fredkin --samples 256 --seed 160503842 --csv build/fredkin_check.csv
```

原始数据/日志 `search_results/round6_*.csv/.log`；[数据校验和](search_results/round6_data.sha256) /
[当时源码校验和](search_results/round6_sources.sha256)。误差按各验证器CSV记录，汇总另独立检查
`abs(Pf-reference)/max(1,abs(reference))`及每个非零权重的实际phase，未把小负数放入误差容限掩盖。

## 第七轮与最终验算（2026-09-11）

R022 PXP种子731/1907各64 samples，两份CSV共1758条：1140正、612零、6预期负。
最大记录相对误差 `1.59807e-15`，正权重相位误差0；[完整CTest84/84通过](search_results/ctest_round7.log)。
[逐文件统计](search_results/round7_summary.json)由[只读脚本](summarize_search_round7.py)生成。
全部原始Fock空间被保留，未投影掉相邻激发簇；负对照包括真实条件符号挫折，单site负Omega另检查为正闭词。

```bash
build/context-check/sign_free_pxp --samples 64 --seed 731 --csv build/pxp_check.csv
build/python-tools/cmake/data/bin/ctest --test-dir build/context-check --output-on-failure
/home/han/miniconda3/bin/python sign_free_catalog/summarize_search_all.py
```

[数据校验和](search_results/round7_data.sha256) / [本轮源码校验和](search_results/round7_sources.sha256)。

## 最终汇总

严格合并Ising/XXZ后保留21个新构造；旧目录和参数/尺寸变体未计入。
七轮共有 **35306条记录**：33081条非零正权重、2123条零权重、96条预期负对照、6条复相位对照。
正权重最大 `|sign-1|=1.48332e-11`；所有记录的最大已报告误差 `3.12312e-12`，
分母按各轮说明，不能把这两个不同误差量混为一谈。
[机器可读总统计](search_results/all_rounds_summary.json)由[summarize_search_all.py](summarize_search_all.py)生成，
不改动旧数据或注册表。零点没有定义sign；负/复对照不是失败，也不计入正配置平均sign。

[MODELS.md](MODELS.md)含21类索引、完整H、允许条件、证明和排重边界。
Gaussian族包含真实PfQMC sweeps；其余专用模型目前提供Hamiltonian与配置权重API，未实现生产级更新器。
有限样本、有限beta与尺寸的通过是实现核验，不能取代文档中的正性证明。
[最终全部源码/文档快照](search_results/final_sources.sha256)与各轮历史数据快照分开保存；
历史源码校验和不会为后续修改而重写。
