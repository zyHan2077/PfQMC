# 对原始笔记的完整覆盖检查

范围是本次读取的 `/home/han/notes/sign_free` 下全部 **15 个 Markdown 文件**。
逐文件原始字节 SHA-256 在 [notes_manifest.txt](notes_manifest.txt)，文件未作修改。
这里只声称覆盖这些笔记中的具体例子和构造，不声称穷尽所有无符号 Hamiltonian。
数学表达式见 [MODELS.md](MODELS.md)，实现入口见 [README.md](README.md)。

记 `C` 为 `sign_free_classification/`，`L` 为 `sign_free_locality/`。
下表的路径均相对于上述原始目录；预设之外的无限系列用 C++ 参数化工厂表示。

| 原始文件 | 内容与逐项处理 |
| --- | --- |
| `README.md` | I、II Wei 排除；III 纯跃迁开链与奇配对环保留，interacting Kitaev chain 排除；IV 所有具体表示族保留；V 静态扇区保留；扩展 E 分别实现；局域组合不另当作基本机制 |
| `C/README.md` | 导航，无额外 Hamiltonian |
| `C/01_问题背景与Wei判据.md` | §§1.1–1.3 统一中心化迹和有符号 Spin 提升；§§1.4–1.5 Wei 两类是用户指定排除项 |
| `C/02_准一维充要条件及证明.md` | §2.1 一般共同占据态相位/Metzler 定理不是一份有限模型表；§2.2 开链 `open_hopping`；§2.3 奇配对环 `odd_pair_ring`、`odd_pair_ring_three`。开放的相邻 hopping+pairing 链归入排除的 Kitaev 家族。§2.4 是适用范围 |
| `C/03_Wei之外的具体无符号族.md` | §3.1 interacting Kitaev chain 排除；§3.2 `four_mode`；§3.3 `adjointOdd(d)`；§3.4 `spinTensor(jx,jy)`；§3.5 `creation_only_chain`、`creation_only_spin`（非厄米） |
| `C/04_群表示分类补充.md` | §§4.1–4.2 权格充要判据是抽象构造原则，落在随后具体系列中；§4.3 `exteriorSp(n,k)`，含 5、48、132 模和拒绝的 14 模负例；§4.4 `exteriorSL(2*k,k)`；§4.5 `sl2Orthogonal`，`2U2`、`2U3` 单独预设，`2U1` 不无符号，`4U1` 属于已知平方/Wei 例子；§4.6 紧群加倍由 Wei 按块覆盖且原文证明不完整，不作为额外非 Wei 例子 |
| `C/05_核对记录与结论边界.md` | 历史核对不是本轮运行证据；保留完整 Fock、固定通道、共同场及 Wei 排除证明的边界。新的数值原始记录单独存档 |
| `L/README.md` | 导航及更强二维巡游问题的限制，无额外模型 |
| `L/01_三轨道单胞与局域耦合定理.md` | §§1.1–1.4 `spin_cell_j1` 与零权重点；§§1.5–1.6 `spin_lattice_j1`；§1.7 `spin_lattice_mixed`，平方保留现场补偿；§1.8 各胞粒子数守恒；§1.9 `spinCell(2*j)`、`spinLattice(j,...)` 无限整数自旋系列 |
| `L/02_二维巡游费米子的局域耦合模型.md` | §§2.1–2.3 `fiveOrbital(lx,ly)`，两个 d 块同场并与 f 的 Sz 共享场；§2.4 一般局域组合可由 `Hamiltonian::channels` 的有限支撑实线性组合表达；§§2.5–2.6 是机制与排除范围说明 |
| `L/03_直接加跃迁的精确反例与受限不可能定理.md` | §§3.1–3.4 十脉冲负例实现为负号控制，并用于测试现有 PFQMC 追踪负号；§§3.5–3.7 是受限不可能定理和研究问题，不是遗漏的无符号 Hamiltonian |
| `L/04_文献搜索与下一步方向.md` | §§4.1–4.4 对应已列 Wei/局域单胞/五轨道/静态族；§4.5 fermion bag、meron、SSE 是求和表示或研究方向，没有给出可逐个实现的具体新 H；§4.6 `pairingLattice`；§4.7 `transformedHopping` 和 `PairHopping` 均落实；§4.8 为开放问题 |
| `L/05_验证记录与可复算证书.md` | 十脉冲有理数负权重重新计算，独立 Fock/Pfaffian 核对；不能把此历史记录当成本次模型已经通过的证据 |
| `L/06_静态守恒量与Falicov_Kimball路线.md` | §6.1 `StaticModel` 和静态 BdG 扩展；§6.2 Falicov–Kimball；§§6.3–6.5 明确只能采样静态扇区，不能任意时间片独立 HS |
| `L/07_复辅助场与充要条件的适用范围.md` | §§7.1–7.4 为假设与修正后的强连通块判据；§7.5 单向复跃迁 `directed_two_mode`、`directed_dag`，`directed_scc` 展示有非平凡强连通块时的构造；§7.6 成对虚势/共轭双块由 Wei 覆盖，排除；§7.7 自旋 1 复场负例实算 -9/2；§§7.8–7.9 为局域耦合与适用范围 |

“排除 Wei 的两类”在这里指不重复实现笔记单列的 Wei 基本族，以及明确由它们覆盖的
平方、紧群和共轭虚场例子。组合五轨道模型含有一个已知无符号的巡游块，但整个指定
共享场族是笔记特别讨论的组合例子，因此保留。外幂例子中某些严格排除 Wei 的证明在
笔记里仍未完成；数值正性不能补上那项证明，不能把它们一律标为“已证明严格非 Wei”。

零点不计为 sign=1；无符号模型也不允许任意拓宽辅助场。
非厄米生成元与世界线构造分别按原文含义实现，不能为满足字面“全部 Gaussian HS”
而编造一个未证明无符号的分解。
