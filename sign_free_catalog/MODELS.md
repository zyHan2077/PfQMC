# 非 Wei、非 interacting Kitaev chain 的模型清单

这里把 `/home/han/notes/sign_free` 的具体构造逐项落成可计算的 Hamiltonian。
来源逐文件、逐节对照见 [AUDIT.md](AUDIT.md)，运行与数值结果见
[README.md](README.md) 和 [RESULTS.md](RESULTS.md)。下述无符号结论的条件来自笔记；
本仓库的有限数值验算检查实现，不把抽样当作对所有参数的证明。

<a id="auxiliary-field-status"></a>

## 辅助场 PfQMC 确认状态（口径校正）

**此前“21个已验证构造”混合了不同正权重表示，不等于21个辅助场PfQMC无符号模型。**
本节按完整厄米H及其写出的物理费米子变量，检查是否给出常规Gaussian-HS分解、
完整辅助场配置的正性论证和signed Pfaffian数值核验。实际PfQMC sweep另行注明。

新增21类中，**确认6类：R001–R005、R009；其余15类未确认上述辅助场PfQMC正性**。
R007已经并入R010，不能重复计数。每节开头及下方索引都标明状态。

| 范围 | 当前确认内容 | 尚未确认/适用边界 |
| --- | --- | --- |
| 旧§1–7 | 26个直接Gaussian-HS预设；signed Pf和真实PfQMC sweep | 仅限各节允许的场、通道和参数 |
| 旧§8：static FK / BdG | 静态守恒背景权重与静态Metropolis | 不是逐虚时间片独立变化的HS场 |
| 旧§9：local_unitary_hopping | 经非Gaussian多体U逆变换后，采样基准H | 原物理费米子变量中的H'直接HS分解未确认 |
| 旧§10：二维pair hopping | occupation世界线路径 | 常规辅助场PfQMC未确认 |
| 旧§11：五个非厄米生成元 | 形式Gaussian-HS场配置迹 | 非厄米、无厄米Walker sweep，不算已确认厄米辅助场模型 |
| R001–R005、R009 | 完整Gaussian-HS权重、正性论证、实际双向PfQMC sweep | 已确认，限文档条件 |
| R006 | Gaussian-HS加守恒扇区求和 | 新增F(N,p)未作常规HS分解；Fourier项并非逐项正 |
| R008 | CT/Poisson Gaussian顶点词 | 常规HS分解和相应PfQMC sweep未确认；并非occupation世界线 |
| R007/R010、R011–R022 | CT、世界线、meron或SSE词，详见逐项标签 | 完整H的常规辅助场PfQMC未确认 |

“未确认”不表示已经证明不存在辅助场解法。相反，**把投影/世界线路径写成Gaussian算符乘积并用Pfaffian计算，
也不自动构成常规辅助场PfQMC验证**。专用验证器的枚举/随机闭词不是生产worldline/SSE Markov链采样。
静态背景、非Gaussian换基和非厄米Gaussian生成元分别具有上表所述证据，不能因加标签而抹去这些区别。
历史84项测试和全部相位统计仍有效，但不能全部归入“辅助场PfQMC测试”。

## 共同记号与实际采样的权重

令 $\delta n_i=c_i^\dagger c_i-1/2$，$D=\sum_i\delta n_i$。
对一个单粒子矩阵 $h$，定义中心化二次算符

$$Q(h)=\sum_{ij}c_i^\dagger h_{ij}c_j-\tfrac12\operatorname{tr}h.$$

若有配对，使用 $\gamma_{2i}=c_i+c_i^\dagger$、
$\gamma_{2i+1}=-i(c_i-c_i^\dagger)$，并记
$Q(A)=\gamma^T A\gamma/4$，$A^T=-A$。
`numberMajorana` 和 `pairMajorana` 把下文的 hopping、density、pairing
转换为这个约定。注意 $A$ 对应 **$Q$**，而 $H$ 中的动能是 **$-Q$**。

绝大多数构造统一实现为

$$H=-Q_0-\frac12\sum_C\lambda_C O_C^2+E_0,\qquad \lambda_C\ge0.$$

每个 $O_C$ 是下文明确允许的厄米二次算符。程序没有遗漏平方展开后的现场项。
例如 $-(O_r+\eta O_s)^2/2$ 包含 $-O_r^2/2$、$-O_s^2/2$ 以及
$-\eta O_rO_s$；只保留最后一项通常不满足这里的构造。
真实 Gaussian HS 恒等式为

$$e^{\Delta\tau\lambda O^2/2}
 =\int\frac{d\phi}{\sqrt{2\pi}}e^{-\phi^2/2}
   e^{\sqrt{\Delta\tau\lambda}\phi O}.$$

代码的每个时间片按时间顺序放置 $e^{\Delta\tau Q_0/2}$、各个 HS 因子、
$e^{\Delta\tau Q_0/2}$。不同 $O_C$ 不对易时仍有 Trotter 误差；该分解收敛到所写的 $H$，
不声称有限 $\Delta\tau$ 已精确等于 $e^{-\beta H}$。

$$W[\phi]=e^{-\beta E_0}\prod_{C,\ell}
 \frac{e^{-\phi_{C\ell}^2/2}}{\sqrt{2\pi}}\;
 Z_f[\phi],\qquad Z_f=\operatorname{Tr}_{\text{完整 Fock 空间}}\prod_\ell e^{Q_\ell}.$$

CSV 记录 $Z_f$，没有乘上前面的正标量；它与 $W$ 的 sign 相同。
当 $Z_f>0$ 时 sign 为 1；$Z_f=0$ 时 sign **没有定义**，不能填成 1。
所有结论都使用完整 Fock 迹，不能任意固定粒子数或奇偶扇区。

## 1. 纯跃迁开链，以及奇配对环

**纯跃迁：**

$$Q_0=\sum_i\mu_i\delta n_i+
 \sum_{i=1}^{N-1}t_i(e^{i\theta_i}c_i^\dagger c_{i+1}+\mathrm{h.c.}),
 \quad t_i\ge0,$$

$$H=-Q_0-\tfrac12\sum_C\lambda_C\left(\sum_i v_{Ci}\delta n_i\right)^2,
 \quad v_{Ci}\in\mathbb R.$$

边相位固定，现场实数场可以独立变化。图必须是若干条不相交的开链：
**一般树不行，三叉点也能产生交换负号。** 代码 `openHopping(N)` 的预设是
6 模，$t=.4$，$\theta_i=.19(i+1)$（这里 $i=0,\ldots,N-2$）；
采用现场平方和邻居密度差平方。`open_hopping` 无 pairing，保留它不重复实现被排除的 interacting Kitaev chain。
不相交链可以把这些矩阵作直和，允许密度平方共享场。

**奇配对环：** 在简单环的每条边上，固定选择一种通道

$$B_{ij}=e^{i\theta_{ij}}c_i^\dagger c_j+\mathrm{h.c.}
 \quad\text{或}\quad
 B_{ij}=e^{i\theta_{ij}}c_i^\dagger c_j^\dagger+\mathrm{h.c.},$$

$$Q_0=\sum_i\mu_i\delta n_i+\sum_{\langle ij\rangle}t_{ij}B_{ij},
 \qquad H=-Q_0-\tfrac12\sum_C\lambda_C(\sum_i v_{Ci}\delta n_i)^2.$$

要求配对边总数为奇数、$t_{ij}\ge0$、相位固定。
`oddPairRing(N,pairingEdges)` 实现任意合法奇数条数；预设
`odd_pair_ring`、`odd_pair_ring_three` 是五模环的 1 条、3 条配对边。
这不是同一边上任意独立控制 hopping 和 pairing 的结论。
这些图可以在一个固定占据态相位约定下保正；开链上的任意配对模式则属于本次排除的 Kitaev-chain 范围。

## 2. 四模本质配对模型

用编号 $1,2,3,4$，定义

$$O_1=n_1+n_2-1,\quad O_2=n_2+n_3-1,$$

$$B_1=c_2^\dagger c_3-c_4^\dagger(c_1-c_1^\dagger)+\mathrm{h.c.},$$

$$B_2=-c_1^\dagger c_2+e^{-i\pi/3}c_4^\dagger(c_3-c_3^\dagger)+\mathrm{h.c.}.$$

$$\boxed{H=-uO_1-vO_2-t_1B_1-t_2B_2
 -U\delta n_1\delta n_2-V\delta n_2\delta n_3},\qquad t_1,t_2,U,V\ge0.$$

因为 $\delta n_1\delta n_2=(O_1^2-1/2)/2$，代码的 $E_0=(U+V)/4$。
只对 $O_1,O_2$ 引入任意实 HS 场，不能对 $B_1,B_2$ 引入任意正负场。
`fourMode()` / `four_mode` 使用 $(u,v,t_1,t_2,U,V)=(.2,-.15,.4,.3,.7,1.1)$。

一个方便的独立校验是把每个脉冲对应到

$$X=\begin{pmatrix}
(2u+v)/3&t_1&0\\
t_1&(v-u)/3&t_2\\
0&t_2&(-u-2v)/3
\end{pmatrix},\quad g=\prod e^X,$$

$$Z_f=2[\operatorname{tr}g\operatorname{tr}g^{-1}-1]\ge16.$$

旧的 `sign_free_check` 继续独立检查这个字符公式；新模型另有 Fock/有符号 Pfaffian/实际 PFQMC 扫描核对。

## 3. 奇数阶伴随表示：四模构造的无限系列

取奇数 $d\ge3$。选 $d^2-1$ 个正交归一的无迹厄米矩阵 $T_a$，
$\operatorname{tr}(T_aT_b)=\delta_{ab}$。代码明确构造对角矩阵和实、虚对称组合，
无需外部表示论软件。给一个实对称三对角矩阵 $X$，定义

$$A(X)_{ab}=\operatorname{tr}\{T_a[X,T_b]\},\quad
 Q(X)=\tfrac14\sum_{ab}\gamma_a A(X)_{ab}\gamma_b.$$

这里 $A$ 是纯虚反对称矩阵，因此 $Q$ 厄米；物理复模数为 $N=(d^2-1)/2$。

$$H=-Q(X_0)-\tfrac12\sum_{a=1}^{d-1}\lambda_a Q(D_a)^2.$$

$X_0$ 的相邻非对角元必须非负；$D_a$ 是任意无迹实对角方向。
HS 改变对角势，不改变正跃迁的限制。预设 `adjoint_d3/d5/d7` 分别是 4、12、24 模，
$X_{0,i,i+1}=.3$、$\lambda_a=.2$。

若 $g=\prod e^{X_\ell}$ 的正本征值为 $\xi_i$，则

$$Z_f=2^{(d-1)/2}\prod_{i<j}(\xi_i+\xi_j)\ge2^N.$$

代码逐配置用小 $d\times d$ 矩阵这个公式交叉检查大 Majorana Pfaffian。

## 4. 整数自旋轨道单胞与二维张量积

这里的 $j$ 标记**轨道矩阵**，不是要求每个单胞只有一个粒子的量子自旋。
令 $f_m$ 有 $2j+1$ 个轨道，$m=j,j-1,\ldots,-j$，

$$(J_z)_{mm}=m,\qquad
 (J_x)_{m,m-1}=\tfrac12\sqrt{(j+m)(j-m+1)},$$

其余元由对称性或零补齐。记 $S^a=f^\dagger J_af$，$D=N_f-(2j+1)/2$。

$$H=-\mu D-aS^x-bS^z
 -\tfrac12(\lambda_D D^2+\lambda_x(S^x)^2+\lambda_z(S^z)^2),\quad j=1,2,\ldots.$$

$\mu,a,b$ 可任意实数，$\lambda\ge0$。
`spinCell(2*j)`，预设 `spin_cell_j1/j2` 是 3、5 模。
**只能改变这三个共同方向，不能加入任意独立轨道现场势。**

$j=1$ 时令 $g=\prod\exp[(a_\ell\sigma_x+b_\ell\sigma_z)/2]$、
$R=\exp\sum\mu_\ell>0$，便有

$$Z_f=R^{-3/2}(1+R)[(1-R)^2+R(\operatorname{tr}g)^2]\ge0.$$

该公式允许零点；也不要求 $Z_f\ge8$。完整 Fock 扇区间的相加是必要的。

二维张量积把轨道标记为 $(m_x,m_y)$，取

$$Q_0=Q\!\left(\mu I+(aJ_x+bJ_z)\otimes I+
 I\otimes(cJ_x+dJ_z)\right).$$

以这五个方向为 $O_C$，仍取 $H=-Q_0-\sum\lambda_C O_C^2/2$。
`spinTensor(jx,jy)` 要求两个自旋为正整数；`spin_tensor_3x3/5x5` 是 9、25 模。
跃迁图确为二维，但五个场在整个图上相关，平方相互作用通常不是短程相互作用。
不能改成每个二维格点独立随机势。

## 5. 局域单胞阵列、混合轨道交换和五轨道巡游模型

在任意有限度图上，每个胞放置上一节的整数自旋轨道，
$O_{r,a}\in\{D_r,S_r^x,S_r^z\}$。允许

$$H=-\sum_{r,a}b_{ra}O_{ra}
 -\tfrac12\sum_C\lambda_C\left(\sum_{(r,a)\in C}v_{C,ra}O_{ra}\right)^2.$$

$C$ 可为单胞、相邻两胞等有限支撑，$v$ 任意实数。
`spinLattice(j,lx,ly,mixed)` 的默认二维 2×2 阵列包括单胞平方、
$(S_r^a\pm S_s^a)^2$；`mixed=true` 包括 $(S_r^x\pm S_s^z)^2$。
预设 `spin_lattice_j1/j2/mixed` 为 12、20、12 模。
各胞 Gaussian 演化仍分块，逐配置的迹是各胞非负迹的乘积；共享场产生真实相互作用。
$N_r$ 守恒，因此这些 $f$ 费米子本身不跨胞巡游。

**五轨道例子**在每胞三个 $f$ 轨道之外，加上 $d_{r\uparrow},d_{r\downarrow}$：

$$H=H_f-t_d\sum_{\langle rs\rangle,\sigma}
(d_{r\sigma}^\dagger d_{s\sigma}+\mathrm{h.c.})
-\mu_d\sum_r C_r-\tfrac g2\sum_r(C_r+S_r^z)^2,$$

$$C_r=n^d_{r\uparrow}+n^d_{r\downarrow}-1,\quad t_d,g\ge0.$$

$H_f$ 是前面的三轨道阵列。`fiveOrbital(lx,ly)` / `five_orbital`
默认 2×2、20 模，$t_d=.35,\mu_d=.09,g=.25$。
同一个 HS 场必须同时作用于 $S_r^z$ 和两个 $d$ 自旋的现场势。
两个巡游块完全相同且为实数，因此它们贡献一个实行列式的平方；
乘上 $f$ 单胞迹仍非负。这是与已知无符号巡游块耦合的组合模型，保留在本次清单中。
不允许任意 $d^\dagger f$ 杂化或两个自旋使用独立场。

**四模配对单胞阵列**也完全一样：以第 2 节的 $H_r$ 为局部 Hamiltonian，加入

$$H=\sum_r H_r-\tfrac12\sum_{\langle rs\rangle,ab}
 g_{rs}^{ab}(O_{r,a}+\eta_{rs}^{ab}O_{s,b})^2,
 \quad a,b\in\{1,2\},\quad g\ge0.$$

`pairingLattice(lx,ly)` 的两个预设 `pairing_lattice/pairing_lattice_2d`
是 2×1、2×2 胞（8、16 模），展示 $O_{r,1}\pm O_{s,2}$ 的共享场耦合。
配对在每胞内部真实参与演化。

## 6. 外幂模型：把小矩阵的作用变成许多轨道的跃迁

这类构造可以不用抽象群语言来实现。给 $d\times d$ 矩阵 $X$，先枚举所有
含 $k$ 个不同编号的集合 $I$。给每个集合一个**独立物理轨道** $f_I$。
矩阵 $L_k(X)$ 的规则为：从集合 $J$ 中取出 $j$，放入尚未占据的 $i$，
矩阵元是 $X_{ij}$ 乘上重新排序的费米子符号；对角元是 $\sum_{i\in J}X_{ii}$。
`wedgeAction(X,k)` 直接按这个规则构造矩阵。
这只是定义轨道矩阵，随后仍取所有 $f_I$ 的完整 Fock 空间，**不是固定 $k$ 粒子投影**。

### SL 中间外幂

取 $d=2k,k\ge2$，$X,Y$ 为无迹实对称矩阵，

$$H=-Q(\alpha L_k(X)+\beta L_k(Y))
-\tfrac{\lambda_X}{2}Q(L_k(X))^2-\tfrac{\lambda_Y}{2}Q(L_k(Y))^2.$$

任意数量的这些方向均可。`exteriorSL(d,k)`；预设
`sl_exterior_4_2/6_3/8_4` 有 6、20、70 个物理复模。
对于所有相应实矩阵乘积 $g$，$\det(I+\Lambda^k g)\ge0$。
不能任意加入化学势 $\mu I$；非中间外幂也不享有这个结论。

### 辛群 primitive 外幂

从 $2n$ 个编号开始，使用实对称矩阵

$$X=\begin{pmatrix}A&B\\ B&-A\end{pmatrix},\qquad A=A^T,\ B=B^T.$$

先构造 $L_k(X)$，再删掉可缩并出一对 $(i,n+i)$ 的冗余分量：
定义 $C=\sum_i a_{n+i}a_i$ 在集合基上的矩阵，令 $P$ 的列为 $\ker C$ 的正交基，
取 $h(X)=P^\dagger L_k(X)P$。$a$ 仅是构造轨道矩阵的工具。
物理轨道数为

$$M=\binom{2n}{k}-\binom{2n}{k-2}.$$

Hamiltonian 是上述 SL 公式把 $L_k$ 换成 $h$。无符号的参数条件为

$$1\le k\le n,\qquad k\text{ 偶数，或 }
 k=2r+1\text{ 且 }\binom nr\text{ 为偶数}.$$

`exteriorSp(n,k)` 实现这个限制；`sp_exterior_2_2/4_3/5_5`
分别为 5、48、132 模。$(n,k)=(3,3)$ 的 14 模例子是负号对照，工厂会拒绝把它当作无符号模型。
代码还核对 $L_k(X)$ 确实保持 $\ker C$，不能随便截去一些轨道代替这个投影。
笔记没有完成这些表示全部排除 Wei 的证明；这里将其作为笔记中额外给出的正性系列保留，
不额外宣称“已证明全部严格在 Wei 之外”。

## 7. 整个 SL(2,R) 上的正交表示族

用 $U_m=\mathrm{Sym}^m\mathbb R^2$ 表示自旋 $m/2$ 的 $(m+1)$ 维轨道矩阵。
若 Majorana 表示写成 $\bigoplus_m n_m U_m$，笔记的条件是

$$n_m\text{ 全为非负偶数},\qquad
 \sum_{m\equiv1\pmod4}n_m/2\equiv0\pmod2.$$

实现中取每个 $m$ 的 $n_m/2$ 份物理轨道块，令 $X,Z$ 是这些块上的
$J_x,J_z$ 直和。所有块**共享同一对场**：

$$H=-aQ(X)-bQ(Z)-\lambda_x Q(X)^2/2-\lambda_z Q(Z)^2/2.$$

这里不允许任意化学势或每块独立的 HS 场。
`sl2Orthogonal(multiplicities)` 按 $m=0,1,\ldots$ 接收重数。
笔记列出的 $2U_2$ 是第 4 节 $j=1,\mu=0$，另给预设 `sl2_u2`；$2U_3$ 是
`sl2_u3` 的 4 模 Hamiltonian；`sl2_mixed` 组合 $2U_2+2U_3$，共 7 模。
若 $x=\operatorname{tr}g$，它们的迹分别为 $2x^2$、$(x+2)^2(x-1)^2$，以及两者乘积。
$4U_1$ 的实行列式平方已被 Wei 覆盖，本次不另立模型；$2U_1$ 不是无符号例子。

## 8. 静态守恒量：Falicov–Kimball 与静态 BdG 扩展

> **常规厄米辅助场 PfQMC：未确认/范围受限。** 未确认一般动态HS；已确认静态背景Pfaffian和静态标签Metropolis。FK及static BdG均适用，场不能逐时间片独立变化。

涵盖预设：`static_fk`, `static_bdg`。

令每点有巡游 $c_i$ 和完全不动的 $f_i$：

$$H=Q(K)+U\sum_i\delta n_i^c\delta n_i^f
 -\mu_f\sum_i\delta n_i^f.$$

$K$ 是任意厄米跃迁矩阵，可在二维图上有复相位/磁通。
$s_i=2n_i^f-1=\pm1$ 是**整个虚时间方向不变**的量。一个配置只需解

$$h_s=K+\tfrac U2\operatorname{diag}s,\quad
 W[s]=e^{\beta\mu_f\sum_i s_i/2}\prod_j2\cosh(\beta\epsilon_j(h_s)/2)>0.$$

`StaticModel` 默认 3×3 OBC，$t=.3,U=3,\mu_f=.2$，非零固定跃迁相位。
相位约定由 `generator(s)` 明确给出。实际权重使用单个厄米 Gaussian 的有符号 Pfaffian，
上式仅用作独立参考。Metropolis 每次翻转一个**静态** $s_i$。
可以在条件 $H_s$ 中加入厄米 $\Delta_{ij}c_i^\dagger c_j^\dagger+\mathrm{h.c.}$，
仍然是单个 $e^{-\beta H_s}$ 的正迹；`pairing=.17` 给出 `static_bdg` 检查。
不能把 $s_i$ 改成逐时间片独立场，再借用这个静态论证。

## 9. 有限深度局域非 Gaussian 酉变换

> **常规厄米辅助场 PfQMC：未确认/范围受限。** 仅在非Gaussian多体U逆变换回基准H后确认辅助场正性和sweep；原物理费米子变量中H'的直接Gaussian-HS分解未确认。

涵盖预设：`local_unitary_hopping`。

取第 1 节的开链 $H$，构造一层不相交的两点门

$$U=\prod_r e^{i\theta_r n_{2r}n_{2r+1}},\qquad H'=UHU^\dagger.$$

这是显式有限范围的相关跃迁 Hamiltonian。例如同一门内
$Uc_i^\dagger U^\dagger=c_i^\dagger e^{i\theta n_j}$。
`transformedHopping(N)` / `local_unitary_hopping` 取 $\theta=.37$。
`exactHamiltonian()` 产生物理 $H'$；采样时逆变换回 $H$ 的基底，逐配置迹不变。
可观测量也必须同样变换。这是等价表示，不是新正性机制。

## 10. 二维 pair hopping：世界线配置的 Pfaffian 实现

> **辅助场 PfQMC：未确认。** 只确认二维pair hopping的世界线配置；没有确认Gaussian-HS辅助场正性，也未实现生产世界线更新器。

涵盖预设：`pair_hopping`。

$$b_r^\dagger=f_{r\uparrow}^\dagger f_{r\downarrow}^\dagger,\qquad
 H=-t\sum_{\langle rs\rangle}(b_r^\dagger b_s+\mathrm{h.c.})+D(n),\quad t\ge0.$$

`PairHopping` 的默认 2×2 图、8 个复模采用

$$D(n)=V\sum_{\langle rs\rangle}(n_r-1)(n_s-1)-\mu\sum_r n_r,
 \quad (t,V,\mu)=(.4,.2,.1).$$

每点有空、双占据以及两个单占据阻塞态；异地点的偶费米子算符相互对易。
代码按原始四次费米子算符重建 Hamiltonian，并验证非对角 hopping 元为 $-t$。
采用离散转移矩阵

$$T_{\Delta\tau}=(I+\Delta\tau\,t\sum_{\langle rs\rangle}
 (b_r^\dagger b_s+\mathrm{h.c.}))e^{-\Delta\tau D},$$

它的矩阵元非负，且 $T_{\beta/L}^L\to e^{-\beta H}$。
每条闭合占据路径由正的 $\Delta\tau t$、$e^{-\Delta\tau D}$ 和偶算符跳跃构成。
为独立核对费米子符号，利用

$$b^\dagger=e^{\frac\pi2(b^\dagger-b)}P_{00},\qquad
 b=e^{-\frac\pi2(b^\dagger-b)}P_{11}.$$

纯占据投影 $P_v$ 是 Gaussian 态，$\operatorname{Tr}P_v=1$；其 Majorana 协方差
$G_{2i,2i+1}=i(1-2n_i)$。`pathPfaffian` 用整个时空反对称矩阵计算
投影与旋转的闭合乘积迹，而不是把 sign 写死。
`closedHopPfaffian` 是一次跳跃再返回的闭环检查，不能理解为非对角矩阵元。
验算覆盖所有合法两顶点闭环、二维绕圈、阻塞态，以及小系统所有三步路径之和。
这是世界线路径权重求值器和小系统完整性检查，**没有声称找到普通 Gaussian HS 下的无符号分解，也没有实现生产级世界线更新器**。

## 11. 笔记的非厄米扩展：如实实现生成元

> **常规厄米辅助场 PfQMC：未确认/范围受限。** 已确认形式上的非厄米Gaussian-HS生成元场配置迹；它们不是厄米Hamiltonian，厄米Walker拒绝它们，没有对应PfQMC sweep。

涵盖预设：`creation_only_chain`, `creation_only_spin`, `directed_two_mode`, `directed_dag`, `directed_scc`。

这些例子不是物理厄米 Hamiltonian。给它们强行补上 h.c. 会改变问题并破坏论证。
`makeGenerator(name)` 明确标记 `hermitian=false`，可以计算形式上的 $H=-Q_0-\sum\lambda O^2/2$
和有符号迹；厄米 `Walker` 会拒绝它们。

**仅创建配对：**

$$Q_\ell=Q(h_\ell)+\sum_{i<j}\Delta_{ij,\ell}c_i^\dagger c_j^\dagger.$$

$h_\ell$ 取第 1 节开链或第 4 节整数自旋单胞允许的场。
增加粒子数的严格上三角部分不改变任何完整乘积的迹。
`creation_only_chain/spin` 分别为 4、3 模，保留任意实数的允许 HS 场。

**单向复跃迁：** 两模例子允许 $c_1^\dagger c_2$ 和 $ic_1^\dagger c_2$
与任意实密度势一起出现；所有单粒子矩阵可在同一顺序下上三角化。
`directed_two_mode` 给出 $i t c_1^\dagger c_2$；
`directed_dag` 扩展到四点菱形有向无环图，可有分叉和无向回路。

**强连通块扩展：** `directed_scc` 有两个各含两点的实正双向 hopping 块，
块间只有任意复数单向 hopping。对角块自身无符号，块间单向项不影响迹。
这是笔记“每个强连通块有共同相位使其 Metzler”的一个显式实例，
不能推广成任意复矩阵都无符号。

成对虚现场势 $i\theta(n_1-n_2)$ 给出 $4\cos^2(\theta/2)$，属于共轭配对/Wei
覆盖的例子，按要求排除。任意复现场势则失败：单模 $i2\pi(n-1/2)$ 的迹为 $-2$，
程序把它保留为检测错误符号实现的反例。

## 12. 持续搜索登记：本轮新增构造（2026-09-11）

本轮计数范围经用户明确为：**本轮新找到、去除明显重复后的 Hamiltonian 构造**。
第 1–11 节的旧清单不计数，改变尺寸、参数、随机种子也不计数。
超过 20 个合格结果后停止扩展，统一去重、合并并整理每一类的详细证明。
当前登记由 [search_registry.json](search_registry.json) 管理；数值记录与复现命令见
[SEARCH_RESULTS.md](SEARCH_RESULTS.md)。历史候选、反例及最终重复性结论见
[SEARCH_QUEUE.md](SEARCH_QUEUE.md)，**未通过有符号数值检验的候选不计入合格数**。

本节的“新增”是相对于当前目录的构造新增，不等同于文献首创或新的基本正性机制。
前几个Gaussian构造使用表示／谱保正，后续构造也包括守恒扇区、世界线及重求和。
排除共同 Wei 结构针对下文指定的整个辅助场族；
不声称同一个物理 H 不存在另一种 Wei 分解。每类的详细证明和全清单的最终合并均已完成。

### 最终索引与统一计数

本轮历史宽口径保留 **21 个独立 Hamiltonian 构造**，其中 **仅6个确认常规辅助场PfQMC，15个未确认**。R007与R010是同一族的Ising/XXZ参数面，
保留两个实现和证明入口，合计一项；因此编号到R022，并不等于22项。
旧§1–11、尺寸、颜色/表示的已知等价分支、参数面、随机种子和同一H的不同展开都未额外计数。
这里的“新”仅指本目录中新加入的构造，不声称文献首创或21种互不相关的基本正性机制。

| 编号 | Hamiltonian的主要构造 | 有限支撑及正性依据 | 常规辅助场PfQMC |
| --- | --- | --- | --- |
| [R001](#r001) | 星形hopping及密度/boost负平方 | 固定轨道胞、邻胞；Lorentz锥 | **已确认** |
| [R002](#r002) | $K\otimes I+I\otimes\bar K$轨道负平方 | 胞与邻胞；共轭谱乘积 | **已确认** |
| [R003](#r003) | Jacobi hopping的对称幂轨道负平方 | 固定多轨道胞；全非负半群正谱 | **已确认** |
| [R004](#r004) | 四元数矩阵congruence的两个表示分支 | 胞与邻胞；共轭成对谱 | **已确认** |
| [R005](#r005) | Lorentz无迹二次多项式轨道 | 胞与邻胞；权值配对 | **已确认** |
| [R006](#r006) | 指定Gaussian族加$F(N_r,p_r)$ | 有限程守恒扇区能量；正扇区迹 | **未确认** |
| [R008](#r008) | $-\sum\kappa e^{aS_r+bS_s}$ | 局域高阶指数顶点；CT正标量与signed迹 | **未确认** |
| [R009](#r009) | Albert 27轨道Jordan乘法负平方 | 固定胞、邻胞；负实谱偶重 | **已确认** |
| [R010](#r010) | $H_{tV}$耦合动态Ising/移动XXZ pair-spin | 二分图；密度插入与闭合spin路径 | **未确认** |
| [R011](#r011) | $c_j^\dagger B_e^\dagger c_i+\mathrm{h.c.}$ | 一维链/环；动态link消去交换符号 | **未确认** |
| [R012](#r012) | $H_{tV}-JP_{ij}-h\sigma_i\delta n_iP_i$ | 二分图；meron分组或条件势路径 | **未确认** |
| [R013](#r013) | 双pair四点ring exchange | 方格plaquette；偶复合粒子路径 | **未确认** |
| [R014](#r014) | $-JP_{ij}-K\Pi_{ij}$颜色磁体 | 有限程图；Cartesian颜色正矩阵元 | **未确认** |
| [R015](#r015) | 三条AF链加边界受挫三角 | 一维链及impurity；成对闭线 | **未确认** |
| [R016](#r016) | $-\sum J_aK_a$，含toric code | 有界支撑；对易Pauli群无$-I$ | **未确认** |
| [R017](#r017) | $-\sum Q_a\prod_{(ij)\in M_a}P^s_{ij}$ | 局域matching；Marshall后的正投影 | **未确认** |
| [R018](#r018) | $J_\mu L_r^\mu L_s^\mu+K_\mu N_r^\mu N_s^\mu$ | 二分dimer底图；三条件下事件偶数 | **未确认** |
| [R019](#r019) | $J_0/J_1$双腿加单向交叉$J_3$ | 开zigzag ladder；显式KT相位 | **未确认** |
| [R020](#r020) | 投影hopping加AF交换的$t$–$J$ | 开链；有序粒子Marshall相位 | **未确认** |
| [R021](#r021) | $P^\uparrow P^s+P^sP^\downarrow$ Fredkin | 三点条件交换；非正非对角H | **未确认** |
| [R022](#r022) | $-\Omega P_{i-1}X_iP_{i+1}$ PXP | 三点条件单spin翻转；非正非对角H | **未确认** |

每个表项下方给出完整二次量子化定义、条件、证明、实现、数值边界与排重依据；证明针对其明确标注的表示。
正性针对指定的采样配置：R006先作扇区投影，R012可先作meron重求和；中间复数项或负resolved路径
不会被冒充成正样本。零权重不定义sign，另行记录。所有非零合格配置的signed Pfaffian相位数值为1。

Gaussian负平方族有实际PfQMC sweep；专用世界线/扇区接口当前是Hamiltonian与配置权重验证器，
没有生产级更新器。完整Fock编码和保留的阻塞态/复制在每节说明，不能混同于投影后的自旋热力学。


<a id="r001"></a>

### R001. Lorentz 星形轨道与局域密度相互作用

> **辅助场 PfQMC：已确认。** 已有完整Gaussian-HS配置正性论证、signed Pfaffian核验及真实双向sweep；只对本节给定通道和参数成立。

每胞有固定的 $M=m+1$ 个轨道，$m\ge3$。令

$$D_r=\sum_{a=0}^m\delta n_{ra},\qquad
 B_{ra}=f_{r0}^\dagger f_{ra}+f_{ra}^\dagger f_{r0},\quad a=1,\ldots,m.$$

一个直接可用的 Hamiltonian 是

$$\boxed{H=\sum_r\left[-\mu_rD_r-\sum_a t_{ra}B_{ra}
 -\frac{u_r}{2}D_r^2-\frac12\sum_a g_{ra}B_{ra}^2\right]
 -\frac12\sum_{\langle rs\rangle}\lambda_{rs}(D_r+\eta_{rs}B_{s1})^2.}$$

所有 $\mu,t,\eta$ 为实数；$u,g,\lambda\ge0$。
因为 $B_a^2=1/2-2\delta n_0\delta n_a$，单胞相互作用也可写成

$$\sum_a g_a\delta n_0\delta n_a-u\sum_{0\le i<j\le m}\delta n_i\delta n_j
 -\frac14\sum_a g_a-\frac{uM}{8}.$$

因此最小四轨道模型只是中心—叶之间的排斥、胞内均匀吸引，以及普通星形跃迁；
它不是一个需要抽象多体算符才能定义的模型。
固定 $m$ 时，胞内作用是 onsite，最后的共享平方只涉及相邻两胞，是严格有限范围四费米作用。
代码 `lorentzModel(m,lx,ly)` 给出单胞和二维阵列；预设 `lorentz4/5/8` 与 `lorentz4_lattice`。
阵列不含单费米子跨胞跃迁，$N_r$ 守恒，不把它误称为二维巡游费米子。

**逐配置正性的证明。** 实 Gaussian HS 只会使每片单粒子矩阵成为

$$h_\ell=\mu_\ell I_M+K(v_\ell),\qquad
 K(v)=\begin{pmatrix}0&v^T\\v&0_m\end{pmatrix},\quad v\in\mathbb R^m.$$

令 $J=\mathrm{diag}(1,-I_m)$。直接相乘得 $K^TJ+JK=0$，所以
$e^K$ 保持 Lorentz 二次型 $x_0^2-\sum_a x_a^2$。
它从单位阵连续出发，因而还保持未来光锥 $x_0>0$；任何有限乘积
$g=\prod e^{K(v_\ell)}$ 都保留这两个性质，且 $\det g=1$。

若 $g$ 有负实本征值 $\xi$，可选非零实本征向量 $x$。
保持二次型给出 $(\xi^2-1)x^TJx=0$。
当 $\xi\ne-1$ 时，$x$ 必须是零范向量。非零零范向量的 $x_0$ 不可能为零，
可把它定向为未来向量；然而 $gx=\xi x$ 会变成过去向量，矛盾。
所以唯一可能的负实本征值是 $-1$。
实矩阵的非实本征值成共轭对，而 $\det g=1$ 又迫使 $-1$ 的代数重数为偶数。
这个论证不要求可对角化。

写 $R=e^{\sum\mu_\ell}>0$。完整中心化 Fock 迹为

$$Z_r=R^{-M/2}\det(I+Rg).$$

正实本征值提供正因子；非实共轭对提供 $|1+R\xi|^2$；
$-1$ 只提供偶次的 $(1-R)$。因此 $Z_r\ge0$。
相邻平方的同一个 HS 场同时改变两个胞的 $\mu$ 或 $v$，但每片仍胞间分块，
$Z_f=\prod_rZ_r\ge0$，外面的 HS 密度严格为正。非零配置 sign=1，零点单列。

**排重与条件边界。** $m=2$ 正是既有整数自旋 $j=1$ 三轨道族，工厂不把它作为新结果。
$m=3$ 的 Lie 代数维数为 6，与已有四轨道 `sl2_u3` 的 3 维允许族不同。
$m=4$ 是 $\mathfrak{so}(1,4)$ 实形式，不是已有五轨道辛外幂的
$\mathfrak{so}(2,3)$，且本族允许任意共同密度 HS。
半填充／没有独立密度方向的某些子族已在 split-orthogonal/Wei 范围中，不能重复宣称新增。

对 $u>0,g_1>0$ 的完整族，$\pm D,\pm B_1$ 都允许。
若存在共同 Wei 结构，其锥不等式在 $A$ 和 $-A$ 上同时成立，就必须变成等式，
产生一个与这些 Majorana 矩阵对易且平方为 $-I$ 的反线性算符 $\Theta$。
因此每个实本征值的空间维数必须为偶数。
但 $2D+B_1$ 的单粒子谱为 $3,1,2,\ldots,2$，其 Majorana 矩阵有简单实本征值 $+3$，矛盾。
这同时排除 Wei 的两种共同结构，且实本征值重数不受固定正则换基影响。
若再有 $g_2>0$，不对易的 $\pm B_1,\pm B_2$ 也排除**任意固定共同 Metzler 基**：
一个算符及其负都 Metzler 就只能是对角矩阵，两个不对易算符不可能同时对角。
这里使用的共同结构条件见 [Wei 原文 §IV–V](https://arxiv.org/html/1712.09412v3)。

不能增加独立叶势 $P_1=n_1-1/2$：四轨道中取按时间顺序的脉冲
$[aB_1,2aP_1,-aB_1,-2aP_1]$、$a=\log3$，得到精确 $Z=-2800/81$。
代码用 signed Pfaffian 和完整 16 维 Fock 迹检查了这个负例。

<a id="r002"></a>

### R002. 共轭张量轨道：同一个复矩阵同时作用于两个指标

> **辅助场 PfQMC：已确认。** 已有完整Gaussian-HS配置正性论证、signed Pfaffian核验及真实双向sweep；只对本节给定通道和参数成立。

每胞用 $d^2$ 个物理轨道 $f_{ij}$，本轮新增范围为 $d\ge3$。
对任意无迹厄米 $d\times d$ 矩阵 $K$，定义

$$h(K)=K\otimes I_d+I_d\otimes\bar K,$$

$$O_r(K)=\sum_{ijk}\left[
 K_{ik}f_{r,ij}^\dagger f_{r,kj}
 +\bar K_{jk}f_{r,ij}^\dagger f_{r,ik}\right],\quad
 D_r=\sum_{ij}\delta n_{r,ij}.$$

$\bar K$ 是逐项复共轭。$K$ 厄米保证 $h(K)$ 厄米，且迹为零，无需隐含投影。
取

$$\boxed{H=\sum_r[-\mu_rD_r-O_r(K_0)-\tfrac{u_r}{2}D_r^2
 -\tfrac12\sum_a\lambda_{ra}O_r(K_a)^2]
 -\tfrac12\sum_{\langle rs\rangle}\nu_{rs}[D_r+\eta_{rs}O_s(K_1)]^2.}$$

$u,\lambda,\nu\ge0$，其余实系数任意；$K_a$ 可包括虚数非对角的厄米方向。
每个胞尺寸固定，以上仍是 onsite 和最近邻四费米相互作用。
`conjugateTensorModel(d,lx,ly)` 实现全部无迹厄米方向，含共同密度场。
预设 `conjugate_tensor3`、`conjugate_tensor3_lattice` 分别为 9、18 轨道。

**证明。** 在一个给定的实 HS 配置中，

$$B=R(g\otimes\bar g),\qquad g=\prod_\ell e^{K_\ell},\qquad
 R=e^{\sum_\ell\mu_\ell}>0.$$

因为 $\operatorname{tr}K_\ell=0$，中心化标量为 $R^{-d^2/2}$。
不论 $g$ 是否可对角化，张量积的本征值是 $\xi_i\bar\xi_j$；因而

$$\boxed{Z_r=R^{-d^2/2}
 \prod_i(1+R|\xi_i|^2)
 \prod_{i<j}|1+R\xi_i\bar\xi_j|^2\ge0.}$$

这里的模平方是独立证明／小矩阵参考公式；数值权重本身依然通过有符号 Pfaffian 算出，
不把根号正支当作 sign。邻胞共享平方保持块结构，正性逐胞相乘。

**排重。** 这是共轭**张量积**，不是两套独立费米子的共轭直和／行列式平方。
$d=2$ 与 R001 的四轨道 Lorentz 单胞严格重复：把两个指标看成 2×2 矩阵，用
$\{I,\sigma_x,\sigma_y,\sigma_z\}/\sqrt2$ 作正交基，
$h(\sigma_a/2)$ 就是中心与第 $a$ 个叶轨道之间的 $B_a$。
代码逐矩阵验证该酉变换，`conjugate_tensor2` 只作排重测试，**不增加结果数**。
这个 Kronecker–Lorentz 对应也有[原始数学物理文献](https://arxiv.org/abs/2110.15118)。
$d\ge3$ 的允许 Lie 代数维数为 $2(d^2-1)$，不同于现有两个独立实 SL(2) 的整数自旋张量族；
同为九轨道不意味着是同一允许场集合。

共同 Wei 的排除可使用与 R001 相同的简单谱矛盾：取
$K=\mathrm{diag}(1,-1,0,\ldots)$，则 $h(K)$ 的最大本征值 2 简单。
$3D+O(K)$ 的单粒子谱全正，最大值 5 简单；Majorana 的实本征值 $+5$ 因此也简单。
在独立 $\pm D,\pm O(K)$ 都允许时不可能存在共同 $\Theta^2=-I$。
至少两个不对易的双向 $O(K_a)$ 同时排除任意固定共同 Metzler 基。
这不排除同一个 H 的其他辅助场表示。

若把两个张量因子的场拆开独立变化，结论立即失败。
例如只在第一个因子按顺序作用 $e^{0.4\sigma_x},e^{0.4\sigma_y},e^{0.4\sigma_z}$，
第二个因子保持单位，四轨道权重约为 $20.473942-1.254887i$，不再 sign=1。
代码用完整 Fock 迹和 signed Pfaffian 交叉核对。

<a id="r003"></a>

### R003. 对称幂轨道：全正辅助半群的新表示实例

> **辅助场 PfQMC：已确认。** 已有完整Gaussian-HS配置正性论证、signed Pfaffian核验及真实双向sweep；只对本节给定通道和参数成立。

最小单胞有六个轨道 $a,b,c,d,e,f$，定义

$$O_1=2\delta n_a+\delta n_b+\delta n_c,\quad
 O_2=\delta n_b+2\delta n_d+\delta n_e,\quad
 O_3=\delta n_c+\delta n_e+2\delta n_f,$$

$$T_1=\sqrt2(a^\dagger b+b^\dagger d)+c^\dagger e+\mathrm{h.c.},\qquad
 T_2=b^\dagger c+\sqrt2(d^\dagger e+e^\dagger f)+\mathrm{h.c.}.$$

$$\boxed{H=-\sum_{r,i}u_{ri}O_{ri}-\sum_r(t_{r1}T_{r1}+t_{r2}T_{r2})
 -\frac12\sum_C\lambda_C\left(\sum_{(r,i)\in C}v_{C,ri}O_{ri}\right)^2.}$$

$t,\lambda\ge0$，$u,v$ 实数。$C$ 只包含一胞或相邻固定数目的胞。
这是受相关场约束的密度相互作用，固定六轨道每胞时严格局域。
不能把六个轨道的现场 HS 独立化。

**矩阵规则与推广。** 六轨道依次对应
$(200),(110),(101),(020),(011),(002)$。这些整数是**定义轨道的标签**，不是物理粒子数约束。
取辅助矩阵

$$X=\begin{pmatrix}u_1&t_1&0\\t_1&u_2&t_2\\0&t_2&u_3\end{pmatrix}.$$

一般对 $d\ge3,q\ge2$ 枚举 $\sum_i n_i=q$ 的标签，每个标签是一个物理费米子轨道，
共 $M=\binom{d+q-1}{q}$ 个。诱导单粒子矩阵的对角为 $\sum_i n_iX_{ii}$；
从 $n$ 把一个单位从 $j$ 移到 $i$ 的矩阵元为 $X_{ij}\sqrt{(n_i+1)n_j}$。
辅助 $X$ 仍只允许实三对角、相邻非对角元非负；HS 只改变其对角项。
`symmetricPowerModel(d,q,lx,ly)` 实现此规则；本轮验证 6、10 轨道单胞及 24 轨道二维阵列。

**证明。** 实 Jacobi 矩阵的指数是全非负矩阵：可通过正双对角／正对角分解的 Trotter
乘积得到，任意子式非负；乘法的 Cauchy–Binet 公式保持此性质。
对角 HS 指数的对角元为正，故辅助乘积 $g=\prod e^{X_\ell}$ 仍可逆且全非负，
其本征值 $\xi_i$ 均为正实数（不可约严格情形是经典全正矩阵正谱定理，
退化情形取极限，并由可逆性排除零）。全正矩阵理论及其群推广是已有结果，见
[Lusztig 的原始综述](https://arxiv.org/abs/0705.3842)。

对称幂表示的本征值为所有 $\prod_i\xi_i^{n_i}$，所以完整中心化迹为

$$\boxed{Z_r=\prod_{\sum n_i=q}
 2\cosh\left(\frac12\sum_i n_i\log\xi_i\right)\ge2^M>0.}$$

这里始终是 $M$ 个物理费米子的完整 Fock 迹；不能把辅助的 $q$ 当作固定物理粒子数。
有限支撑共享平方仅改变各胞对角实场，因此 $Z_f=\prod_rZ_r$。

**排重与限制。** $d=2$ 是已有整数自旋对称幂，$q=1$ 是原来的开链，不新增计数。
$d=3,q=2$ 不同于已有六轨道中间外幂：允许 Lie 代数和场限制均不同；
但它与已登记伴随构造使用同一个全正辅助半群思想，**不是新的基本正性机制**。
未来其他 Young 图的 Schur 表示应统一审查合并，不按每个标签增加结果数。

若三个 Cartan HS 均可正负变化，取 $X=\mathrm{diag}(1,3,9)$，物理单粒子谱为
$2,4,10,6,12,18$，Majorana 的正负实本征值全简单。
与前述共同 Wei 的反线性偶重要求矛盾。跳跃图的分叉还能产生两粒子交换负闭环，
排除固定占据态**相位**把全部 hopping 变正；这里没有声称排除了任意多体换基。

独立六轨道现场势确实会破坏正性。保持 hopping 为 $(T_1+T_2)/2$，按顺序添加六个对角势：

```text
[ 0, 1, 4, 0,-4, 2]
[-1, 1,-2,-3,-4, 2]
[ 3,-2,-2, 2, 1, 1]
[-2,-4,-1, 2, 0, 3]
[-1, 2,-4, 0,-2,-2]
[-1,-4, 3, 2, 3, 3]
```

完整中心化迹为约 $-7.14635870944\times10^{10}$。
高精度独立消元、64 维 Fock 指数和 signed Pfaffian 同号，说明共同场约束不是可省略的装饰。

<a id="r004"></a>

### R004：四元数矩阵轨道上的局域相互作用（两个分支合计一项）

> **辅助场 PfQMC：已确认。** 已有完整Gaussian-HS配置正性论证、signed Pfaffian核验及真实双向sweep；只对本节给定通道和参数成立。

这里的四元数只用来生成一组**实 hopping 系数**；物理粒子仍是普通复费米子。
一个固定单胞取 $M_-=n(2n-1)$ 或 $M_+=n(2n+1)$ 个轨道 $f_a$。
分别用四元数 Hermitian 或 anti-Hermitian 矩阵的实正交基 $T_a$ 标记它们，
内积为 $\langle X,Y\rangle=\operatorname{Re}\operatorname{tr}_{\mathbb H}(X^\dagger Y)$。
不熟悉四元数也可直接用代码中的四个二阶复矩阵

$$q_0=I_2,\quad q_1=\begin{pmatrix}i&0\\0&-i\end{pmatrix},\quad
q_2=\begin{pmatrix}0&1\\-1&0\end{pmatrix},\quad
q_3=\begin{pmatrix}0&i\\i&0\end{pmatrix}.$$

将 $q_\alpha/\sqrt2$ 放在第 $(i,j)$ 个块，上／下块互为 $\pm$ 共轭转置；
对角块在负分支用 $q_0$，正分支用 $q_{1,2,3}$。复矩阵内积需乘 $1/2$。
这些 $T_a$ 是轨道系数表，不是额外的量子算符。对任意无迹四元数 Hermitian $K$，定义

$$h(K)_{ab}=\langle T_a,KT_b+T_bK\rangle,\qquad
 O_r(K)=\sum_{ab}h(K)_{ab}f_{ra}^\dagger f_{rb},\qquad
 D_r=\sum_a(n_{ra}-\tfrac12).$$

$h(K)$ 实对称、无迹，故 $O$ Hermitian。允许的二次量子化 Hamiltonian 是

$$\boxed{H=-\sum_r[\mu_rD_r+O_r(K_{0r})]
-\frac12\sum_X\lambda_X\left[\sum_{r\in X}
 (d_{Xr}D_r+O_r(K_{Xr}))\right]^2,\quad\lambda_X\ge0.}$$

$X$ 只包含固定距离以内的有限个单胞。每个胞的 $n$ 固定，因而上述四费米子相互作用有 locality。
实系数 $d,\mu$ 任意；$K$ 可以包含全部四元数方向，不能仅保留实矩阵后宣称验证了完整族。
`quaternionModel` 的六个预设覆盖两分支、$n=2,3$、单胞与邻胞共享场。

**为什么非负。** 每个辅助时间片作用为 $X\mapsto e^KXe^K$，完整乘积为
$X\mapsto gXg^\dagger$。记其 $2n$ 阶复实现为 $G$，满足 $GJ=J\bar G$、$J^2=-I$。
于是复特征值成 $\lambda_i,\bar\lambda_i$ 对；实特征值也有偶代数重数。
这个标准四元数事实可见 [Loring, §§1–3](https://arxiv.org/html/1107.0500)。

映射 $Y=\chi(X)J$ 将负／正分支分别变成复反对称／复对称矩阵，
且 $Y\mapsto GYG^T$，所以物理传播矩阵的复化分别是 $\Lambda^2G$、$\operatorname{Sym}^2G$。
设密度场累计 $R=\exp(\sum\mu_\ell)>0$。完整中心化费米迹为

$$Z_-=R^{-M_-/2}\prod_i(1+R|\lambda_i|^2)
 \prod_{i<j}|1+R\lambda_i\lambda_j|^2
                 |1+R\lambda_i\bar\lambda_j|^2\ge0,$$
$$Z_+=R^{-M_+/2}\prod_i(1+R|\lambda_i|^2)
 \prod_{i<j}|1+R\lambda_i\lambda_j|^2
                 |1+R\lambda_i\bar\lambda_j|^2
 \prod_i|1+R\lambda_i^2|^2\ge0.$$

这只是 signed Fock 迹的解析表达，程序仍从 Pfaffian 算出原始符号后比较。
Jordan 块不改变特征值乘积公式。零因子允许出现：零权重的 sign 没有定义。
共享平方的 HS 仅使不同胞的 $K$ 相关，不破坏各胞公式。

**排重。** 两分支属于同一个四元数 congruence 构造，在本轮保守地合计一项。
$\mathrm{Herm}_{\mathbb H}(2)$ 的六轨道通过中心／差密度正交换基，完全等价于 R001 的
五叶 Lorentz 星；程序逐生成元检查这一点，该特例不新增计数。
计数依据是 $n\ge3$ Hermitian 分支及完整 anti-Hermitian 族，而不是六个预设的数量。
这不是 R002 的复矩阵共轭张量：允许辅助群和物理不可约表示不同。

对于完整允许场族，$K=\mathrm{diag}(1,-1)$（每项为二阶单位块）和 $3D$ 给出
$h(K)+3I$ 的最大实本征值 $5$，负／正分支重数分别为 $1,3$。
共同 Wei 的反线性平方为 $-I$ 要求这些实本征空间偶维，故此 HS 族不满足该共同条件。
不能由辅助 $G$ 有 Kramers 对称便把物理 $M_\pm$ 轨道重新算成 Wei 类。
这也不等于证明 Hamiltonian 不存在另一套 Wei 分解。

**数值保护。** 验证包括基正交、物理指数与直接 $gT_bg^\dagger$ 的矩阵比较、非对易生成元、
Pfaffian／辅助谱／普通 determinant 三方比较及真实 PfQMC 双向 sweep。
六轨道另做完整 64 维 Fock 比较。放宽为独立叶现场势时，四片给出
$Z=-11200/81$；该负号对照防止把 sign 强制设为 1。

<a id="r005"></a>

### R005：Lorentz 无迹对称矩阵轨道（多项式表示扩展合计一项）

> **辅助场 PfQMC：已确认。** 已有完整Gaussian-HS配置正性论证、signed Pfaffian核验及真实双向sweep；只对本节给定通道和参数成立。

取 $N\ge5$，令 $J=\mathrm{diag}(1,-1,\ldots,-1)$。
单胞轨道标记为实对称矩阵空间中满足 $\operatorname{tr}(JX)=0$ 的正交基 $T_a$，
共有 $M=N(N+1)/2-1$ 个普通费米子。只允许星形实矩阵
$K=\sum_{i=1}^{N-1}k_i(E_{0i}+E_{i0})$。定义

$$h(K)_{ab}=\operatorname{tr}[T_a(KT_b+T_bK)],\qquad
O_r(K)=\sum_{ab}h(K)_{ab}f^\dagger_{ra}f_{rb}.$$

Hamiltonian 与 R004 方框式相同，$K$ 改用上述星形限制，$D_r=\sum_a(n_{ra}-1/2)$。
这是明确的有限胞 hopping 加局域四费米子平方；`lorentzSymmetricModel` 给出
14、20 轨道和两胞 28 轨道实例，均只计这一个构造。

**证明。** 因 $KJ+JK=0$，映射 $X\mapsto KX+XK$ 保持 $\operatorname{tr}(JX)=0$，
并且在完整对称矩阵空间上湮灭 $J$。指数乘积是 $X\mapsto gXg^T$，
$g\in SO^+(1,N-1)$。R001 已证明 $g$ 的负实特征值有偶重数。
将其特征值写成正实单项及共轭成对项（负实的相等两项也作为一对）。
在 $\operatorname{Sym}^2g$ 的 $\lambda_i\lambda_j$ 特征值列表中，
交换所有共轭标签：非固定项成对，固定项是正实数，因此负实项仍必偶重。
去掉不变向量 $J$ 只除去特征值 $1$。所以

$$\boxed{Z=R^{-M/2}\frac{\prod_{i\le j}(1+R\lambda_i\lambda_j)}{1+R}\ge0,
\qquad R>0.}$$

同一标签配对论证适用于一般多项式 Schur 表示：权的重数在标签置换下相等。
这里将这个闭包以最简单的无迹二阶对称表示实现；未来不按 Young 图逐个增加计数。
$N=3$ 是旧整数自旋，$N=4$ 是 R002 的一个受限特例；实现从 $N=5$ 开始。
R003 使用 TN 辅助半群，其谱全正且 signed HS 仅在 Cartan；本构造允许所有正负 boost，
可产生复谱，允许场族并非 R003 的尺寸变化。

$N=5$ 时单个 boost 的物理谱为 $\{2^1,1^3,0^6,(-1)^3,(-2)^1\}$；
加 $3D$ 有简单实本征值，排除共同 Wei 偶重条件。两个 boost 不对易。
独立检查 $gT_bg^T$ 与物理指数乘积一致。
若把 $K=E_{00}$ 的作用强行正交投影回该 14 维空间，记为 $B$，另记合法 boost 为 $A$，
四片 $[aA,bB,-aA,-bB]$，$a=\log4,b=5\log(5/4)$，产生精确负迹

$$Z=-\frac{95232669932822733339285069632243}{4611686018427387904000000000000}
\simeq-20.65029352655228.$$

这是独立有理矩阵指数／消元结果；当前 Pfaffian 及 determinant 对照必须重现这个负号。
14 轨道未做指数大的完整 Fock 枚举，不能声称有这一 oracle。

<a id="r006"></a>

### R006：正超选择扇区上的局域充电／奇偶相互作用（合计一项）

> **辅助场 PfQMC：未确认。** 只确认基础Gaussian-HS场与守恒扇区的联合权重；新增F(N,p)经扇区精确处理，Fourier求和项可为复数，未确认完整H的常规HS分解。

若单胞粒子数 $N_r$ 与每个辅助 Gaussian 算符对易，可先按 $N_r$ 或奇偶
$p_r=(-1)^{N_r}$ 分块。**只有各被选扇区的迹都非负时**，才允许给它们乘任意正 Boltzmann 因子。
这带来不同于吸引负平方的局域相互作用，而不是把所有粒子数投影都误认为 sign-free。

充电实例以 R003 的六轨道胞为基础，加入

$$\boxed{H_{\rm charge}=H_{\rm TN}+\sum_r U_r(N_r-m_r)^2
 +\sum_{\langle rs\rangle}V_{rs}N_rN_s+\sum_r v_rN_r.}$$

$U,V,v$ 可以任意实，尤其允许排斥 $U>0$；固定单胞的谱有限。
也可以施加硬约束 $N_r=k_r$。若单胞传播矩阵的特征值全正，
固定 $k$ 粒子数的中心化迹为 $Z_{r,k}=e^{-\sum\operatorname{tr}h_\ell/2}e_k(B_r)>0$，
因为基本对称多项式 $e_k$ 是正实单项之和。于是

$$W(s,\{k_r\})=e^{-\beta F(\{k_r\})}\prod_r Z_{r,k_r}(s)\ge0.$$

这不是通过 imaginary HS 对排斥充电项声称每个场正；$F$ 在守恒扇区内被精确处理。

第二个实例从三个轨道的整数自旋 $j=1$ 胞出发，保持化学势为零，加入

$$\boxed{H_{\rm parity}=H_{j=1,\mu=0}
 +\sum_{\langle rs\rangle}J_{rs}p_rp_s+\sum_r v_rp_r,\quad
 p_r=\prod_{a=1}^3(1-2n_{ra}).}$$

$J,v$ 任意实，也可固定各胞奇偶。整数自旋传播谱包含 $1$，从而
$\operatorname{Tr}[(-1)^N\widehat B]=\det(I-B)=0$。
所以 $Z_{\rm even}=Z_{\rm odd}=Z/2\ge0$；邻胞 parity 乘积是明确的有限范围高阶费米子交互。
跨胞共享平方允许存在，但不能加入改变 $N_r$ 的胞间单粒子 hopping 后仍用同一证明。

**Pfaffian 实现。** `canonicalPfaffians` 用不少于 $M+1$ 个单位根计算完整有限多项式的系数：
$Z_k=L^{-1}\sum_j e^{-ik\theta_j}\operatorname{Tr}[e^{i\theta_jN}\widehat B]$。
每个求和项用 signed Pfaffian，保留中心化密度的相位；这些复数 Fourier 项本身不是采样权重。
所有项相加后的 sector weight 才要求 sign 为 1。`ProjectedModel` 实现有限程 $F$ 和硬限制。

**验证与边界。** 独立 Fock 粒子数块迹逐项比较；通过全局 fermion bit 操作重建实际 Hamiltonian，
并用 Gaussian-HS 数值积分及 2/4/8 片检查收敛，含非对易场。
$j=1$ 的固定 $N=1,2$ 负对照各为 $-1887/4096$，而奇偶和各为 $2209/4096$。
另有化学势破坏奇偶正性、零权重、Spin lift 负号、两种邻胞耦合符号和硬限制测试。
它是扇区权重 API／验证器，尚无生产级扇区 Monte Carlo 更新器。
充电与奇偶、不同约束及不同耦合符号共同算一个“正守恒扇区扩展”，不分开凑数。

<a id="r007"></a>

### R007：动态 Ising–fermion 的三轨道实现

> **辅助场 PfQMC：未确认。** 只确认CT-INT加自旋世界线路径的完整权重；用Pfaffian核验这些路径不等于验证Gaussian-HS辅助场。与R010合并。

**最终排重：本项是R010在 $J_{xy}=0,J_z=-J$ 的参数面，合并进R010，不再单独计数。**


每个二分格点放三个普通费米子 $c_i,a_i,b_i$，令
$\delta n_i^c=n_i^c-1/2$、$D_i=n_i^a+n_i^b-1$、$P_i=a_i^\dagger b_i^\dagger+b_i a_i$。

$$\boxed{H=-t\sum_{\langle ij\rangle}(c_i^\dagger c_j+\mathrm{h.c.})
 +V\sum_{\langle ij\rangle}\delta n_i^c\delta n_j^c
 -\frac J4\sum_{\langle ij\rangle}D_iD_j
 +\frac12\sum_i h_i\delta n_i^cP_i.}$$

$V\ge0$，$t,J,h_i$ 任意实；只连接相反子格，所有相互作用为 onsite 或 nearest-neighbor。
偶 $a,b$ 占据子空间中 $S_i^z=D_i/2,S_i^x=P_i/2$ 是 spin-1/2；
两个奇占据态则是不会翻转的 spin-zero 状态。**无需删除这些 Fock 态或施加物理投影**。

正性采用 [Huffman–Chandrasekharan §II](https://arxiv.org/html/1605.07420v2)
的 CT-INT 与自旋路径组合。子格旋转使交互变成 $h_i\sigma_i\delta n_i^cS_i^x$。
展开密度键与翻转交互，保留对角 Ising 演化。闭合自旋路径在每个点翻转偶数次，
因此 $(-1)^{l}\prod_i h_i^{l_i}$ 非负；对角传播是正指数，每次翻转的矩阵元是 $1/2$。
费米子因子带密度键展开的 $(-1)^m$ 以及 $\sigma_i\delta n_i^c$ 插入，
二分格点粒子—空穴配对给出非负迹。两者相乘得到非负的完整配置权重。
奇占据阻塞点不能翻转，错误／不闭合路径权重为零。
这是对文献模型的全费米子编码和 Pfaffian 实现，不作首创声明。

`DynamicIsing::configurationPfaffian` 不把上述正性公式当作计算符号的方法。
实际使用 $\delta n=(-i/2)e^{i\pi\delta n}$，再将指数拆成两个半角，
以 occupation Gaussian projector 及成对旋转实现 $S^x$ 路径。
完整 $3N$ 模时空 Pfaffian 乘上**全部**展开负号、$-i/2$、$h_i\sigma_i/2$、$V$ 和正对角传播。
Gaussian projector 不要求可逆 $B$，所以不能塞进普通连续 HS walker。

验证包含两点完整 64 维 Fock Hamiltonian／词迹、四点全部非零路径、随机事件时刻，
以及重新编号的断开二分图和三叉图。CSV 保留所有参数、边、事件、占据路径与原始相位。
零／负对照分别检查不闭合、阻塞、错误旋转和单独不满足正性条件的密度子迹。
该实现是可调用 CT-INT 配置权重 API，尚无生产级插入／删除更新器。
因为 $a,b$ 配对随虚时翻转，它不同于旧 FK/BdG 的静态守恒标签。
没有证明此 Hamiltonian 绝不可能存在另一种 Wei 分解；登记的是这套不同的显式正分解。

<a id="r008"></a>

### R008：局域 Gaussian 指数顶点产生的高阶相互作用

> **辅助场 PfQMC：未确认。** 已确认连续时间/Poisson展开的Gaussian顶点词权重非负；没有确认完整高阶H的常规HS分解，也没有运行对应辅助场PfQMC sweep。它不是occupation世界线验证。

每胞三个轨道 $a,b,c$，定义普通费米双线性

$$S_z=n_a-n_c,\quad S_x=(a^\dagger b+b^\dagger c+\mathrm{h.c.})/\sqrt2,
\quad D=N-3/2.$$

在完整八维单胞 Fock 空间上 $S_x^3=S_x,S_z^3=S_z$。对邻胞定义
$G_v=a_vS_{\alpha_v,r}+b_vS_{\beta_v,s}$，$\alpha_v,\beta_v\in\{x,z\}$，取

$$\boxed{H=-\sum_{v=(r,s)}\kappa_v e^{G_v}
 -\sum_r[\kappa_{r,+}e^{d_rD_r}+\kappa_{r,-}e^{-d_rD_r}],\quad
 \kappa_v,\kappa_{r,\pm}\ge0.}$$

所有 $a,b,d$ 实，边为最近邻；同时包括 $ZZ,XZ,ZX$ 及反号 $XZ$ 顶点。
指数是一个**有限的二次量子化多项式**，无需无穷多体项：

$$e^{aS}=1+\sinh(a)S+[\cosh(a)-1]S^2,$$
$$e^{aS_r+bT_s}=[1+\sinh(a)S_r+(\cosh(a)-1)S_r^2]
[1+\sinh(b)T_s+(\cosh(b)-1)T_s^2],$$
$$e^{dD}=\prod_{j=a,b,c}[\cosh(d/2)+(2n_j-1)\sinh(d/2)].$$

胞内项最高六费米子，邻胞项最高八费米子。特别是一个 $ZZ$ 顶点贡献
$n_{ra}n_{rc}n_{sa}n_{sc}$ 的系数 $-4\kappa(\cosh a-1)(\cosh b-1)$，一般非零。
程序用四个占据数的有限差分独立检查这个系数。
因此完整 Hamiltonian 不能重写成既有“二次动能加二次算符负平方”族：后者最高四费米子。
单独 onsite 指数与不同轴／参数不另计数。

**正性。** 用 $H_0=0$ 的连续时间展开，积分有序顶点时刻后，一个有序词的完整权重是

$$W(v_1,\ldots,v_k)=\frac{\beta^k}{k!}\prod_{j=1}^k\kappa_{v_j}
 \operatorname{Tr}(e^{G_{v_k}}\cdots e^{G_{v_1}}).$$

单胞因子均在 spin-1 的 Lorentz 群及正密度缩放中，故 R001 的同一谱配对给出完整迹非负。
多个顶点在共享胞上不对易，仍保持该群闭包。每项外部标量非负，所以配置非负。
归一化 Poisson 抽样若作为未来提议分布使用，要补对应常数；当前计算的就是上式未归一化展开测度。

局域 $\pm D$ 顶点的 Majorana 实本征值具有三重数，排除这套顶点族的共同 Wei 结构；
只用 $\mu=0$ 的 $S_x,S_z$ 不具有这个排除证据，不能把它们单独宣传为 Wei 之外的新模型。
反号 $XZ$ 顶点另给同一 hopping 矩阵元相反符号，排除固定占据态相位使每个顶点矩阵都非负。
这不是对任意多体换基或其他展开的绝对排除。

`PoissonModel` 将每个顶点生成元送入 signed Pfaffian，并保留 $\beta^k/k!$ 与耦合系数。
独立 oracle 使用上面的有限多项式，完整枚举两胞全部 $k\le3$ 的 585 个顶点词，
逐词比较并验证其和等于 $\beta^k\operatorname{Tr}(-H)^k/k!$。
另核对 Hamiltonian、多胞非对易性及 18 阶矩阵级数与 $e^{-\beta H}$。
负耦合对照和 Spin lift 负迹必须返回负号；真实零迹及零耦合单列。
提供两／三胞权重 API 和验证器，尚无生产级 CT 更新器。

<a id="r009"></a>

### R009：Albert 固定 27 轨道单胞

> **辅助场 PfQMC：已确认。** 已有完整Gaussian-HS配置正性论证、signed Pfaffian核验及真实双向sweep；只对本节给定通道和参数成立。

单胞取 27 个普通复费米子 $f_A$。轨道标签组织成三个实对角坐标及三组各八个实坐标：

$$X=\begin{pmatrix}x_1&u&v\\\bar u&x_2&w\\\bar v&\bar w&x_3\end{pmatrix},
\quad x_i\in\mathbb R,\quad u,v,w\in\mathbb O.$$

八元数只是一张生成实 hopping 系数的乘法表，物理场仍是复费米子。
实现用递归 Cayley–Dickson 规则 $(a,b)(c,d)=(ac-\bar d b,da+b\bar c)$，
从实数依次得到复数、四元数、八元数。所有矩阵乘法仅含两个因子；不假定八元数结合律。
取 $X$ 的实正交基 $T_A$（对角单位元，非对角元素及其共轭除以 $\sqrt2$），
定义 Jordan 乘积及实 hopping 系数

$$K\circ X=\tfrac12(KX+XK),\qquad
 h(K)_{AB}=\langle T_A,K\circ T_B\rangle,\qquad
 O_r(K)=\sum_{AB}h(K)_{AB}f_{rA}^\dagger f_{rB}.$$

$K$ 为无迹的上述矩阵，$h(K)$ 实对称、无迹。完整允许 26 个无迹方向，
Hamiltonian 是

$$\boxed{H=-\sum_r[\mu_rD_r+O_r(K_{0r})]
-\frac12\sum_X\lambda_X\left[\sum_{r\in X}(d_{Xr}D_r+O_r(K_{Xr}))\right]^2,
\quad D_r=\sum_{A=1}^{27}(n_{rA}-1/2),\quad\lambda_X\ge0.}$$

每胞固定 27 轨道，$X$ 只在 onsite／邻胞取支撑，所以它是有限程四费米子 Hamiltonian。
`makeAlbert` 给出单胞及邻胞阵列，使用全部 26 个方向和密度方向。

**正性所需的数学事实与证明。** 这部分比前面的矩阵谱配对更深；实现 Hamiltonian 无需先掌握群论。
关键结论是传播矩阵的任何负实本征值都成偶重数，故 $\det(I+Rg)\ge0$。
以下给出完整推导链，而不是用随机数值代替它。

Albert 空间有三次形式

$$\Delta(X)=x_1x_2x_3-x_1|w|^2-x_2|v|^2-x_3|u|^2
 +2\operatorname{Re}[(uw)\bar v].$$

其 adjugate 满足 $X^\#\circ X=\Delta(X)I$，导数
$d\Delta_X(Y)=\operatorname{Tr}(X^\#\circ Y)$。
trace 内积的关联性给出
$d\Delta_X(K\circ X)=\operatorname{Tr}K\,\Delta(X)$。
所以无迹 $K$ 的 $e^{h(K)}$ 保持 $\Delta$，连续地连接单位元。
该实群是 $E_{6(-26)}$，其最大紧子群是固定 $I$ 的 $F_4$；任何 Albert 元素可由 $F_4$ 对角化。
这些群及对角化事实见 [Chatterji §2](https://math.univ-cotedazur.fr/u/indira/papers/Dedicata1.pdf)。

对任意生成的 $g$ 作可交换分解 $g=ehu$：$e$ 椭圆、$h$ 双曲、$u$ 幺幂。
有限维表示中三者的谱分别在单位圆、正实轴及 $1$；
见 [Kostant, Proposition 2.1/2.3 与 Lemma 3.4](https://numdam.org/item/10.24033/asens.1254.pdf)。
经共轭可使 $h=e^{L_H}$，$H=\operatorname{diag}(H_1,H_2,H_3)$、$\sum H_i=0$，$L_H(X)=H\circ X$。

这里必须保证下一步的共轭不改变 $h$。其 centralizer $Z_G(h)$ 对 Cartan 对合稳定且 reductive；
$e$ 的幂的闭包是其中的紧子群，故可**在这个 centralizer 内部**共轭进入最大紧子群
$Z_G(h)\cap F_4$。这正是紧子群共轭定理在 centralizer 上的应用，不是对 $e,h$ 分别任意换基。
因为 $kL_Hk^{-1}=L_{kH}$ 且 $L_H(I)=H$，该紧子群恰为 $F_4$ 中固定 $H$ 的群。

$H_i$ 全不同、两个相同、全相同的稳定子分别是 $\operatorname{Spin}(8)$、
$\operatorname{Spin}(9)$、$F_4$，三个群都连通。
primitive idempotent／frame 稳定子的具体矩阵实现见
[Miyashita §1](https://www.math.okayama-u.ac.jp/mjou/mjou44/_43_Miyashita.pdf)。
因其连通，$e$ 在 $h$ 的每个实本征空间上的正交作用行列式都是 $+1$；
因此该空间中的 $-1$ 本征值必偶重。乘上同一正实 $h$ 本征值只把它变成负实 $-r$，
仍偶重；对易的幺幂因子不改变代数特征值重数。
结论适用于退化 $H$、非对角化 $g$，不局限于 generic 情形。

密度累计 $R=e^{\sum\mu_\ell}>0$，因为 $\operatorname{tr}h(K)=0$，完整中心化迹为

$$\boxed{Z_r=R^{-27/2}\det(I_{27}+Rg)\ge0.}$$

共享 HS 场使 $g_r$ 相关，但没有改变单胞所属群，故全权重仍是非负因子之积。

**去重、实现和证据边界。** 只用实 $3\times3$ 子代数时，27 轨道分解为
$\operatorname{Sym}^2\mathbb R^3\oplus7\Lambda^2\mathbb R^3$，可并入旧共轭张量与平方块，
不能作为新构造依据。这里允许全部 26 个方向，生成 exceptional 的不可约 27 维表示，
不是 R001/R002/R004 的尺寸或参数变化。$\pm D$ 有 27 重实 Majorana 本征空间，
已经排除本 HS 族的共同 Wei 偶重条件；方向也确实不对易。

程序的独立 oracle 用 Jordan 二次表示
$P(y)X=2y\circ(y\circ X)-(y\circ y)\circ X$、$y=\exp_\circ(K/2)$
计算 $e^{L_K}$，没有复用同一个 27 阶矩阵指数。
检查基、Jordan 恒等式、三次不变量、Pfaffian／determinant／Jordan 作用、实际左右 PfQMC sweep 与 Green。
放开为某个独立物理轨道势时，负控制为 $-23488102400/81$；独立两模 Fock 乘 $2^{25}$ 与完整 signed Pf 同号。
完整 27／54 模 Fock 空间未直接枚举，数值验证使用上述可控独立 oracle。
全参数正性来自群论推导，有限短虚时样本仅检验实现；不据此声称长虚时数值稳定性或文献首创。

补充中心化子步骤：Cartan involution满足 $\theta(h)=h^{-1}$，故 $Z_G(h)$ 是 $\theta$稳定的。
双曲元素共轭到正Cartan子群的步骤使用 [Kostant Proposition 2.4](https://numdam.org/item/10.24033/asens.1254.pdf)。
当 $H=\mathrm{diag}(a,a,-2a)$ 时，$L_H$权空间维数为10、16、1，Spin9作用分别为
$1\oplus9$、spin16和trivial；任何偶然权值碰撞合并这些实空间后，连通作用的行列式仍为正。

<a id="r010"></a>

### R010：有移动配对的 XXZ–fermion 局域模型

> **辅助场 PfQMC：未确认。** 只确认CT-INT加动态Ising/XXZ pair-spin世界线权重；包含R007的Ising参数面，未确认完整H的Gaussian-HS辅助场正性。

**统一计数包括R007的Ising面与下述移动XXZ面；共一项，保留各自参数域和实现。**


每点仍为 $c_i,a_i,b_i$，记 $B_i^\dagger=a_i^\dagger b_i^\dagger$、
$D_i=n_i^a+n_i^b-1$、$\delta n_i^c=n_i^c-1/2$。取

$$\boxed{H=-t\sum_{\langle ij\rangle}(c_i^\dagger c_j+\mathrm{h.c.})
 +V\sum_{\langle ij\rangle}\delta n_i^c\delta n_j^c
 +\frac{J_z}{4}\sum_{\langle ij\rangle}D_iD_j
 +\frac{J_{xy}}2\sum_{\langle ij\rangle}(B_i^\dagger B_j+\mathrm{h.c.})
 +\frac12\sum_i h_i\delta n_i^c(B_i^\dagger+B_i).}$$

二分图、$V,J_{xy}\ge0$，$t,J_z$ 实；同一可移动连通分量上的非零 $h_i$ 要同号。
所有项 onsite／nearest-neighbor，最高四费米子。$a,b$ 奇占据为阻塞的 spin-zero 状态，保留完整 Fock 空间。

**证明与新增内容。** 按 R007 的子格旋转，pair hopping 系数变成 $-J_{xy}/2$。
CT-INT 展开中的额外负号使其路径系数成为正 $J_{xy}/2$；$J_z$ 留在对角正传播中。
一条闭合路径可以先在 $i$ 产生 pair，再移动至 $j$，最后湮灭；
翻转次数只需在整个移动连通分量上为偶数，已不要求每个点分别偶数。
因此 $h_i$ 同号使 $\prod h_i$ 非负，费米子因子仍由二分粒子—空穴配对非负。
这将 [Huffman–Chandrasekharan §III](https://arxiv.org/html/1605.07420v2) 的 Heisenberg 例子
推广到独立实 $J_z$ 和非负 $J_{xy}$；其证明所需的对角／非对角符号不变。

$J_{xy}=0$ 时回到 R007（$J=-J_z$），该极限不新增计数；
Heisenberg、XY、XXZ 也只算这里一个移动配对构造。
新项 $B_i^\dagger B_j$ 改变了各胞的 pair 数和路径闭合条件，不是 Ising 的尺寸变化。
它也不是旧 `PairHopping` 的重命名，后者没有巡游 $c$ 及 onsite $\delta n^c(B^\dagger+B)$ 耦合。

`MobileSpin` 实现完整三轨道 signed Pfaffian、全部 CT 标量和动态路径。
每次 pair hop 用两个反向的 Gaussian pair 旋转与 occupation projector。
`signFreeParameters()` 单独检查定理条件；权重 API 仍允许不满足同号条件的真实负对照，不能强行给其 sign=1。
两点事件“$S^x_0@0.2$、pair hop@$0.5$、$S^x_1@0.9$”在 $\beta=1.2$、
$t=.8,V=.5,J_z=.7,J_{xy}=.9,h=(.6,.4)$ 时，全部自旋路径之和为
$0.0003586743649135766$。每点恰一次翻转，故不是原 Ising 路径。
把 $h$ 改成 $(.6,-.4)$，同一完整词为其负值。

验证含原未旋转基底的完整六模 Fock Hamiltonian／词迹、四点环、重新编号的二分图、阻塞／零控制。
只提供权重 API，没有生产级移动 CT 更新器。

<a id="r011"></a>

### R011：一维动态 quantum-link correlated hopping

> **辅助场 PfQMC：未确认。** 只确认matter/link占据世界线的CT权重；signed Pfaffian计算的是带occupation projector的路径。

一维开链或环上放 matter 费米子 $c_i$，每条有向键 $e=(i,j)$ 另放两个物理费米子 $a_e,b_e$。
令 $B_e^\dagger=a_e^\dagger b_e^\dagger$、$D_e=n_{a_e}+n_{b_e}-1$、$\delta n_i=n_i^c-1/2$。

$$\boxed{H=-t\sum_{e=(i,j)}(c_j^\dagger B_e^\dagger c_i+\mathrm{h.c.})
 +V\sum_{e=(i,j)}\delta n_i\delta n_j
 +\frac g2\sum_{e=(i,j)}D_e(n_j^c-n_i^c)+\mu\sum_i\delta n_i,\quad t\ge0.}$$

$V,g,\mu$ 可任意实；也可用任意有限范围实对角能量替换后三项。
向前移动一个 matter 粒子必须同时把键上 $00$ 变成 $11$；反向跃迁撤销该变化。
奇占据的键态阻塞跃迁而不被删除。全 Fock 空间不施加 Gauss-law 约束。
这是最高四费米子的局域 Hamiltonian；matter 数守恒，但总物理费米子数可改变 2。

它来自有限维 $U(1)$ quantum-link 的一维版本；文献模型与零绕行条件见
[Banerjee–Huffman, Eqs. (2)–(4) 及 Algorithm](https://arxiv.org/html/2305.08917v2)。
当前方向约定与该文交换键端点，$V=t,g=-t,\mu=0$ 对应其单 flavor 参数，差一个无关常数。
这里只验证一维允许的更一般对角参数，不将文献高维 meron 结论套到这些任意参数上。

**直接正性证明。** 将全部 $c$ 排在前，键上的 $a,b$ 邻接排列。
开链中最近邻 $c$ hopping 的 Jordan–Wigner 符号为正；pair 算符是偶算符，不额外产生交换符号。
在环上唯一闭边的符号为 $(-1)^{N_c-1}$。
给每个完整占据态 $s$ 乘一个固定实相位

$$U(s)=(-1)^{[N_c(s)-1]n_{a,\mathrm{boundary}}(s)}.$$

每次闭边 hopping 翻转该 $a$ 占据，同时 $N_c$ 不变，恰好消去闭边符号。
内部边既不改变这个 $a$ 也不改变 $N_c$。所以所有非对角 $UHU^\dagger$ 矩阵元均为 $-t$，
实对角项保持不变；全部闭合世界线权重非负。
等价地，闭合键状态要求每条边正反 hopping 次数相等，环上的 matter 绕行数为零。
这个相位是依赖多体占据的变换，不是普通 Gaussian 单粒子换基；不声称排除了 stoquastic 正性。

**去重与适用边界。** 该跃迁是四次项；不同于旧 open hopping／odd pairing ring 的二次动能。
相邻的 correlated hopping 项只共享一个 matter 轨道，而旧 pair hopping 的相邻项共享完整一对轨道，
所以并非固定轨道重标号、相位或逐轨道粒子—空穴的简单重复。
这里没有证明任意非线性多体幺正都不能联系这些模型。
高维 quantum-link 或多 flavor 后续若实现，应合并进本构造，不按维度／flavor 加计数。

`QuantumLink` 实现 Gaussian occupation projector 与
$\exp[\frac\pi2(c_j^\dagger c_i-c_i^\dagger c_j)\pm\frac\pi2(B_e^\dagger-B_e)]$
的有符号 Pfaffian 路径，乘每次 hopping 的 $t$ 及全部对角虚时传播。
独立 bit 算符 oracle 保留 fermion parity。两点开链、三点环还从完整 Fock creation/annihilation 矩阵重建 H，
逐条检查全部合法转移及上述相位变换；另有非回文闭词、闭边原始负 hopping 成对出现的词和完整二阶 CT 路径求和。

必要的负控制：三点环中两个 $c$ 粒子绕行交换后，键状态不闭合，完整合法权重为零；
若额外允许独立翻转键以恢复初态，同一扩展词的 signed Pfaffian 为 $-1$。
一维几何也必要：三叉树上 $1\to0\to3,2\to0\to1,3\to0\to2$ 可让所有键归位、
同时交换两个粒子，signed Pfaffian 为 $-1$。
因此 API 严格只接受链或环，不因零净键电流而宣称一般树／二维逐路径 sign-free。
这是路径权重验证器，无生产级更新器。

<a id="r012"></a>

### R012：动态 singlet-projector 与条件交错密度

> **辅助场 PfQMC：未确认。** 只确认meron重求和或旋转后的条件势/世界线表示；单条resolved路径可为负，c子系统的辅助场条件不等于完整H的HS分解。

每点 $c_i,a_i,b_i$，记 $B_i^\dagger=a_i^\dagger b_i^\dagger$、
$E_i=1-n_{a_i}-n_{b_i}+2n_{a_i}n_{b_i}$、$D_i=n_{a_i}+n_{b_i}-1$。
定义 $S_i^+=B_i^\dagger,S_i^z=D_i/2,S_i^x=(B_i^\dagger+B_i)/2$，以及

$$P_i=E_i/2+S_i^x,\qquad P_{ij}=E_iE_j/4-\boldsymbol S_i\cdot\boldsymbol S_j.$$

$P_i$ 是偶占据中 $|00\rangle+|11\rangle$ 的归一化 rank-one projector，
$P_{ij}$ 是两个 active pair-spin 的 singlet projector；两者都湮灭被触及的奇占据态。
取二分图及子格号 $\sigma_i=\pm1$：

$$\boxed{H=H_{tV}-J\sum_{\langle ij\rangle}P_{ij}
 -h\sum_i\sigma_i\delta n_i^c P_i,\quad V,J,h\ge0,}$$

其中 $H_{tV}$ 与 R007 相同，$t$ 任意实。所有项局域，最高八费米子。
$E_iE_j$ 的四个密度乘积在 H 中系数为 $-J$，不能丢掉它再宣称相同全 Fock Hamiltonian。
这使奇占据成为完整保留的阻塞态，而非隐含的 Hilbert-space 投影。

**两种正分解属于同一构造。** 原始 spin 基底下，对每个 singlet 顶点取水平连接、
每个 $P_i$ 顶点取两个端点，沿周期虚时连接。设触及部分产生 $N_c$ 条闭／开簇，
共有 $n_s$ 个 spin 顶点，未触及 $u$ 个点。若开簇端点位于相反子格，它是 meron，
两个方向的 spin 赋值抵消；其余每个簇贡献正因子2。因此 spin 迹是

$$Z_s=2^{N_c-n_s+2u}\,\mathbf1_{N_{\rm meron}=0}.$$

这里未触及点给四态，而非纯 spin 模型的两态。完整 CT-INT 权重
$W=V^{n_V}J^{n_J}h^{n_h}Z_sG_f\ge0$，其中 $G_f$ 是带密度键展开负号和
$\sigma_i\delta n_i^c$ 插入的 signed 费米迹。簇规则对应
[Huffman–Chandrasekharan §IV](https://arxiv.org/html/1605.07420v2)。
时刻已经指定，$W$ 是有序 simplex 上的密度，不额外乘 $\beta^k/k!$。

也可做现场 spin 旋转及子格 Marshall 变换，得到严格相同的 Hamiltonian

$$H'=H_{tV}+\frac J4\sum(D_iD_j-E_iE_j)
 -\frac J2\sum(B_i^\dagger B_j+\mathrm{h.c.})
 -h\sum_i\sigma_i\delta n_i^c n_{a_i}n_{b_i}.$$

此时 pair 世界线逐条非负，巡游 $c$ 看到系数非负的单向交错势。
$c$ 子系统的正性可用已有 Wei/粒子—空穴条件；不能据此声称这里发明了新基本正性机制，
也不能将“c 子系统满足 Wei”误等同于“完整含八费米子交互的 H 已是旧 Gaussian HS 模型”。
这两种展开只计一个完整动态 pair／conditional-density 构造。

**去重与数值。** 旧 PairHopping 无巡游 $c$；R010 无 $E_iE_j$ 八费米子项及条件投影势。
现场旋转不会消去这里的最高阶项。singlet rank-one 顶点的两费米子波函数不是 Slater，
也不能直接充当 R008 中未加常数的正 Gaussian 顶点；这不是对所有可能加常数／重分解的绝对排除。

`MeronModel` 既计算 cluster 分组权重，也在小系统中**实际求和完整三轨道 signed Pfaffian 路径**，
不把负 resolved 路径的 sign 置为 1。独立 oracle 为六模 Fock Hamiltonian／词迹，
另外逐矩阵核对上面的旋转式和八费米子系数。
一个合法 meron 词的单条路径为 $-0.00733930606247$，但分组后精确为零；
把 $P_i$ 错拆成 $S_i^x$ 会使完整词为 $-0.01467861212494$。
均匀四点环还存在额外的费米子零迹，用整数矩阵恒等式 $K^3=4K$ 及所有九个相关 trace
的严格零值确认，不能靠浮点阈值把未解析零点算成正样本。
此实现是完整配置权重／分组 API，没有生产级 cluster 更新器。

<a id="r013"></a>

### R013：四点复合 pair ring exchange

> **辅助场 PfQMC：未确认。** 只确认pair/ring交换的occupation世界线CT权重。

每个格点两个物理费米子，$b_r^\dagger=f_{ra}^\dagger f_{rb}^\dagger$，$N_r=n_{ra}+n_{rb}$。
按正方形顶点的循环顺序 $1,2,3,4$ 定义

$$\boxed{H=-J\sum_{\langle ij\rangle}(b_i^\dagger b_j+\mathrm{h.c.})
 -K\sum_\square(b_1^\dagger b_2b_3^\dagger b_4+\mathrm{h.c.})
 +V\sum_{\langle ij\rangle}(N_i-1)(N_j-1)-\mu\sum_i N_i,
\quad J,K\ge0.}$$

$V,\mu$ 实；可加其他有限程实对角项。ring 顶点同时交换两对粒子，最高八费米子。
偶占据 $00/11$ 是 hard-core pair，奇占据态完整保留并阻塞经过它的跃迁。
这是 [Melko–Sandvik 的 J–K 模型 Eqs. (1)–(3)](https://arxiv.org/html/cond-mat/0409112)
的显式全费米子实现；不冒称 quantum dimer 模型，因为这里没有额外 dimer-covering 约束。

**证明。** 不同胞的 $b,b^\dagger$ 都是偶费米算符，彼此对易。
合法 hopping/ring 在 pair 占据基中的非零矩阵元均为 $+1$，因此 CT 展开的非对角权重是
正 $J$ 或 $K$，对角演化是正实指数。任意闭合路径的完整权重非负。
有阻塞态时只删去相应可跃迁路径，没有引入负号。

$K=0$ 回到旧 pair hopping，不新增计数；$K>0$ 的单个 H 矩阵元可同时改变四胞，
任何仅由一／两胞 H 项相加的模型在该位置严格为零。
因此它是新的四胞交互构造；四点环的方向、格点尺寸及对角参数不各自计数。
R008 当前 nearest-neighbor 双胞指数顶点也不能作为单个 H 项产生这四胞变化；
没有声称排除了任意重新分组单胞或任意多体幺正的关系。

`RingExchange` 用每个发生变化的胞上的 $\pm\pi/2$ Gaussian pair rotation 及完整 occupation projector
计算 signed Pfaffian，再乘所有 $J,K$ 和对角虚时因子。
八模独立 Fock oracle 按写出的 fermion monomial 从右向左作用，保留每一步 CAR 符号。
全部物理转移、非回文“一个 ring 加两个 bond hop”、完整三顶点闭合路径之和及四／六胞随机路径均验证。
令 $K<0$ 而 $J>0$，上述三顶点闭环的完整权重为负，是不可被绝对值抹掉的对照。
零 $K$ 和奇占据阻塞单独记零。该接口允许计算负域权重，`signFreeParameters()` 才标记定理域。
当前提供 CT 路径权重 API，不含生产级 ring 更新器。

<a id="r014"></a>

### R014：SO(3) color-singlet 与 permutation 磁体

> **辅助场 PfQMC：未确认。** 只确认Cartesian颜色基中的投影/permutation世界线CT权重。

每点三个普通复费米子 $f_{ix},f_{iy},f_{iz}$，定义

$$S_i^a=-i\sum_{bc}\epsilon_{abc}f_{ib}^\dagger f_{ic},\quad
 N_i=\sum_a n_{ia},\quad A_i=N_i(3-N_i)/2,\quad Q_{ij}=\boldsymbol S_i\cdot\boldsymbol S_j.$$

完整八态是 $0\oplus1\oplus1\oplus0$：粒子数1和2各一个spin1，0和3为保留的阻塞态。
取任意有限程图（允许三角形），

$$P_{ij}=(Q_{ij}^2-A_iA_j)/3,\qquad
 \Pi_{ij}=Q_{ij}^2+Q_{ij}-A_iA_j,$$
$$\boxed{H=D(\{N_i\})-\sum_{\langle ij\rangle}(J_{ij}P_{ij}+K_{ij}\Pi_{ij}),\qquad J_{ij},K_{ij}\ge0.}$$

$D$ 是任意实有限程密度函数。API实现 $U\sum(N_i-3/2)^2+\mu\sum(N_i-3/2)
+V\sum(N_i-3/2)(N_j-3/2)+V_A\sum A_iA_j$。
这些是显式至多八费米子项；$A_iA_j$ 在完整Fock空间中不是常数。
特别地 $K=0,V_A=-J/3,D$的其余系数为零时，正好为负biquadratic磁体 $-(J/3)\sum Q_{ij}^2$。
原始模型见 [negative biquadratic §II](https://arxiv.org/html/1208.4133v1)
及 [SO(N) projector/permutation Eqs. (1)–(2)](https://arxiv.org/html/2406.02926)。

**证明。** 在两个active扇区分别取Cartesian颜色基
$|a\rangle_1=f_a^\dagger|0\rangle$、
$|a\rangle_2=\frac12\epsilon_{abc}f_b^\dagger f_c^\dagger|0\rangle$。
两者的spin矩阵都是 $(S^a)_{bc}=-i\epsilon_{abc}$。直接收缩两个epsilon，得到

$$P_{ij}=\frac13\sum_{a,b}|aa\rangle\langle bb|,\qquad
\Pi_{ij}=\sum_{a,b}|ab\rangle\langle ba|.$$

因此展开 $-H$ 的合法顶点矩阵元分别为 $J/3$ 和 $K$，对角传播为正实指数。
$N_i$ 沿世界线守恒，阻塞扇区仅删去触及的顶点；在每个保留扇区中任意闭合路径非负。
物理occupation基与颜色基只差固定现场相位（双占据101为负），每条跃迁的相位为
$u(s_{k+1})u(s_k)$，闭合词乘积严格为1。这不是把单步的原始负号删除。

**去重。** 所有颜色数、图、$J/K$分支统一一个构造。它包含真正的八费米子singlet/permutation项，
不能把它当作只有两态active空间的旧pair hopping或R010。
与R008直接正Gaussian顶点的区分可在未平移singlet上检查：两粒子波函数不是Slater，
对应反对称矩阵rank6；跨点partial transpose有本征值 $-1/3$。
这些检查排除该顶点本身属于正的product-Gaussian锥，不宣称排除任意加常数和重新分解。

`ColorProjector` 从occupation projector与带正确方向的Gaussian颜色旋转计算完整signed Pfaffian。
独立64维bond Fock矩阵核对两种多项式算符全部8192矩阵元；还验证完整三顶点词求和、
四种粒子/空穴组合的非回文三颜色闭环、三点三角形随机词和阻塞零权重。
$J<0$ 时奇数次singlet三颜色环为负。接口提供CT路径权重，不含生产更新器。

<a id="r015"></a>

### R015：三条反铁磁链的受挫三角 impurity

> **辅助场 PfQMC：未确认。** 只确认dimer基的CT闭合世界线；Gaussian态投影及转移是路径求值工具。

取三条等长开链 $a=1,2,3$，每条偶数个点 $i=0,\ldots,L-1$：

$$\boxed{H=\lambda\sum_{a=1}^3\sum_{i=0}^{L-2}\boldsymbol S_{a,i}\cdot\boldsymbol S_{a,i+1}
+\alpha[\boldsymbol S_{1,0}\cdot\boldsymbol S_{2,0}
 +(\boldsymbol S_{1,0}+\boldsymbol S_{2,0})\cdot\boldsymbol S_{3,0}],\quad\lambda,\alpha\ge0.}$$

它是局域四费米子Hamiltonian：每个spin用三个Majorana，
$S^a=-\frac i4\sum_{bc}\epsilon_{abc}\gamma^b\gamma^c$，$\{\gamma^a,\gamma^b\}=2\delta_{ab}$，
再代入 $\gamma_{2r}=f_r+f_r^\dagger,\gamma_{2r+1}=i(f_r^\dagger-f_r)$。
将两条对称臂的同位置spin配成三复模胞，第三臂相邻两spin配胞，得到 $9L/2$ 个普通复模。
此编码 $\boldsymbol S^2=3/4$ 恒成立；完整Fock空间是物理spin空间的 $2^{3L/2}$ 份复制，
无需单占据投影，也没有spin0阻塞态。两条对称臂的等耦合是证明条件。

**正性及边界。** 对配对臂采用现场dimer基
$|s\rangle=(|\uparrow\downarrow\rangle-|\downarrow\uparrow\rangle)/\sqrt2$、
$|t_0\rangle=(|\uparrow\downarrow\rangle+|\downarrow\uparrow\rangle)/\sqrt2$、
$|t_+\rangle=|\uparrow\uparrow\rangle,|t_-\rangle=-|\downarrow\downarrow\rangle$；第三臂采用 $S^z$ 基。
不能声称该基中每个局部 $-H$ 矩阵元都非负。正性来自完整闭合世界线：将两条对称臂折叠后，
将第三臂的up/down分别续接为 $t_+/t_0$ 线，$t_+$ 与 $t_-$ 视为正反向同类线，
另有 $t_0$ 和 $s$ 两类无向线。同类线不能自交；不同闭线在一维开链的虚时面上交叉偶数次。
局部符号还可分配为 $t_0t_0$ 与 $t_+t_-$ 的产生/湮灭各一个负号，$ss$为正；
每条闭线的产生与湮灭成对，因而交叉和转向的总负号都为偶数。
任意长度的论证及全部局部转移规则见
[Hann–Huffman–Chandrasekharan §§IV–V、Table 4](https://arxiv.org/html/1608.05144v2)。
实对角传播始终为正，各局部守恒parity复制具有相同闭词迹，故完整Fock求和也非负。
这里实现的是文献的spin Hamiltonian的精确费米编码；文中另一个有限 $U$ 三臂Hubbard模型
并没有被证明sign-free，不能把强耦合极限结论推广给它。

**实现与排重。** 该构造保留一个真实受挫三角形，并和任意长的动态spin链相连；
不同于旧无巡游介质的pair hopping、R007/R010的二分spin–density耦合或R014的八次颜色投影。
臂长、各parity复制和局部编码选择不另计数。
每个三复模胞的固定parity dimer态是纯Gaussian态；转移用Majorana双线性旋转实现，
保留实际复矩阵元的除数和每次 $-H$ 顶点，完整时空Pfaffian给出权重。
独立 $L=2$ 验证完整512维Fock H为64维spin谱的八份复制；$L=4$ 使用18模signed Pfaffian，
并核对局部矩阵词，不枚举 $2^{18}$ 维全矩阵。
普通 $S^z$ 基三角闭词和破坏两臂等耦合的词均有负号；小到 $10^{-33}$ 的合法权重仍保留并检查，
不以绝对截断把它归零。这是CT配置验证器，没有生产世界线更新器。

<a id="r016"></a>

### R016：局域 commuting-stabilizer 与 toric-code 费米子模型

> **辅助场 PfQMC：未确认。** 只确认对易Pauli词的SSE/CT展开权重；并非对高阶stabilizer相互作用作HS分解。

每条边 $e$ 放两个复费米子 $c_e,d_e$，定义

$$Z_e=1-2n_{c_e},\qquad
X_e=c_e^\dagger d_e+d_e^\dagger c_e+c_e^\dagger d_e^\dagger+d_ec_e.$$

它们在完整四态空间满足 $X^2=Z^2=1,XZ=-ZX$，即一个spin1/2的两份复制。
在方格环面上取star $s$和plaquette $p$：

$$\boxed{H=-\sum_sJ_s\prod_{e\ni s}X_e-\sum_pK_p\prod_{e\in\partial p}Z_e,\qquad J_s,K_p\ge0.}$$

各项有限支撑、最高八费米子，完整Fock空间没有投影。
原始toric code见 [Kitaev](https://arxiv.org/abs/quant-ph/9707021)。
统一的构造允许任意有界支撑Hermitian Pauli词 $K_a$，只要彼此对易且生成群不含 $-I$，
$H=-\sum_aJ_aK_a,J_a\ge0$。这包括带固定正负号的生成元；所有stabilizer格点/维度仅计一项。

**证明。** CT展开每个词的完整权重为

$$W(a_1,\ldots,a_k)=\frac{\beta^k}{k!}\prod_jJ_{a_j}
\operatorname{Tr}(K_{a_k}\cdots K_{a_1}).$$

对易且平方为1保证词为生成群元素。非平凡Pauli词的迹为0，恒等词的迹为完整Fock维数；
因不允许 $-I$，不会得到负迹。Toric star与plaquette共享0或2条边，因此对易；
只有各类全体生成元乘积为 $+I$，不存在负关系。
$E$条边的完整配分函数还可直接求得

$$Z=2^{2E}\prod_s\cosh(\beta J_s)\prod_p\cosh(\beta K_p)
[1+\prod_s\tanh(\beta J_s)][1+\prod_p\tanh(\beta K_p)].$$

一般生成群的条件用保留 $i$ 相位的GF(2)消元检查，不能只检查binary support。

**实现与去重。** 每个Pauli是 $-i\gamma_b\gamma_c=-i\exp(\pi\gamma_b\gamma_c/2)$，
拆成两个半角Gaussian以避单片零迹，所有 $-i$ 保留。实现按互异胞的偶算符因式分解，
每个胞内仍按原时序计算signed时空Pf；另用不因式分解的全局Pf检验Spin lift。
独立oracle直接作用完整occupation态，并用Pauli相位代数及小系统Fock谱核验。
$xzxz$给负、$zz$给正、单个$x$为零、$xyz$给纯虚；翻转一个star产生负恒等关系。
并未先将Pauli词化简为正数再冒称Pf验算。

Finite-depth CZ把独立X旋转成cluster-state不另计，属于旧§9；toric具有闭环多胞交互和全局关系。
这里不声称产生新基本正性机制，也不排除 $(I+K)/2$ 的其他正Gaussian分解。
本API计算完整展开词及其积分标量，无生产级stabilizer更新器。

<a id="r017"></a>

### R017：J-Q 多键 singlet-projector 模型

> **辅助场 PfQMC：未确认。** 只确认Marshall变换后的singlet-projector CT/SSE路径权重。

每个spin用两个普通复费米子 $a_i,b_i$：

$$S_i^z=\tfrac12-n_{a_i},\quad S_i^+=(b_i^\dagger+b_i)a_i,
\quad S_i^-=a_i^\dagger(b_i^\dagger+b_i),\quad
P_{ij}=\tfrac14-\boldsymbol S_i\cdot\boldsymbol S_j.$$

每个现场parity各含一份spin1/2，因此完整 $Z_f=2^{N_s}Z_{spin}$，无阻塞态或隐含单占据约束。
选二分图及局部互不相交的键集 $M_a$（matching），

$$\boxed{H=\sum_i h_iS_i^z-\sum_a Q_a\prod_{(ij)\in M_a}P_{ij},\qquad Q_a\ge0.}$$

常见J-Q2是单键J加方格上两条平行相邻键的乘积Q；Q3含三条相邻不相交键。
分别最高四、八、十二费米子；每个matching限定有限大小和直径。
原始模型见 [Sandvik Eq. (1)](https://arxiv.org/html/cond-mat/0611343)。
程序保留 $-J P$ 的常数，不能与 $J\boldsymbol S_i\cdot\boldsymbol S_j$ 的绝对能量混同。

**证明。** 在一个子格做 $S^\pm\mapsto-S^\pm$。
每个 $P_{ij}$ 在反平行子空间为 $\frac12\begin{pmatrix}1&1\\1&1\end{pmatrix}$，平行子空间为零。
不相交键乘积的矩阵元为这些非负数的乘积，CT顶点乘 $Q_a$ 仍非负。
任意 $h_i$ 的对角传播为正；完整闭词的子格相位首尾抵消。
对局部parity复制逐个应用同一论证即得全Fock正性。

**去重。** $Q=0$的普通pair/spin交换不新增。$Q>0$同时翻转四点，排除仅两胞color项；
它还使一条键的交换振幅依赖另一条键反平行还是平行，单纯pair hopping加XY ring不具备该矩阵元。
两种Q2方向的局部顶点不对易，因而不是R016 commuting-stabilizer特例。
所有matching长度、Q2/Q3、格点和颜色推广统一计一个多键投影构造。

`JQModel`在全部physical occupation态上使用Gaussian投影与两半角Majorana旋转，保留CT负号与系数。
四spin完整256维Fock H核对16份spin谱，逐项核对393216个顶点矩阵元和3712个二阶完整Fock词；
另有Q2/Q3随机词、零顶点、负Q奇词及正耦合三角交换的负号控制。无生产更新器。

<a id="r018"></a>

### R018：动态 singlet–triplet 各向异性受挫 bilayer

> **辅助场 PfQMC：未确认。** 只确认singlet-triplet/dimer基中的CT世界线权重。

每个二分底图格点$r$放两个spin，记 $\boldsymbol L_r=\boldsymbol S_{1r}+\boldsymbol S_{2r}$、
$\boldsymbol N_r=\boldsymbol S_{1r}-\boldsymbol S_{2r}$。每个spin用R015的三个Majorana，
故一胞三个普通复模，完整Fock为spin空间的 $2^{N_r}$ 份复制。取统一实耦合：

$$\boxed{H=\sum_{\langle rs\rangle}[J_zL_r^zL_s^z+J_\perp(L_r^xL_s^x+L_r^yL_s^y)
+K_zN_r^zN_s^z+K_\perp(N_r^xN_s^x+N_r^yN_s^y)]
+\sum_r[D_zS_{1r}^zS_{2r}^z+D_\perp(S_{1r}^xS_{2r}^x+S_{1r}^yS_{2r}^y)].}$$

所有项局域四费米子。正性域是 $K_z=0$ 或 $K_\perp=0$ 或 $J_\perp=0$；其余系数允许任意实号。
物理同层/交叉层耦合分别 $\mathcal J_\mu=J_\mu+K_\mu,\mathcal K_\mu=J_\mu-K_\mu$。

**证明。** 在每胞singlet/triplet基把非对角过程分为：triplet内磁量子数跃迁$m$、
$s\leftrightarrow t_0$交换$l$、$ss\leftrightarrow t_0t_0$的$p$、
$s$与$t_\pm$交换的$lm$、$ss\leftrightarrow t_+t_-$的$pm$。
闭词及底图二分性给出模2约束

$$N_p+N_{pm}=0,\quad N_m+N_{pm}+N_{lm}=0,\quad N_l+N_{lm}=0.$$

$K_\perp=0$删去$lm,pm$，$K_z=0$删去$l,p$；剩余每类次数都偶数。
$J_\perp=0$时$m$不存在，$pm$独占磁量子数pair产生/湮灭，所以它也出现偶数次，再用三式。
每类矩阵元的符号固定，故无论统一耦合的符号如何，闭词符号均为正。
现场项在该基对角，以实指数精确传播。论证不依赖配位数，适用任意二分有限程底图。
详见 [Alet–Damle–Pujari 的正权重证明](https://arxiv.org/html/1511.01586v2)。

**去重与边界。** 只有同时$K_z=K_\perp=0$才守恒每胞singlet；一般三种条件都保留动态$s$–$t$转换，
不能用静态sector特例替代。三种条件、各向异性和底图只计一个构造。
它含双层直/交叉交换，不是R015的三条开链和单个impurity三角；不是R010的巡游密度翻转。
具体而言R015配对双臂内部对应 $J_\perp=K_z=K_\perp=\lambda/2$，三个通道全部非零，
不在本节三个条件中；它依赖另外的一维闭合世界线证明。

`AnisotropicBilayer`保留原始dimer基signed顶点及每个parity Gaussian projector。
六模完整Fock谱、十二模独立逐Majorana occupation词、不同配位数图及三条件四种符号组合均核验。
三个条件均破坏的四胞方环词给 $-0.02997378689911806$，与独立物理Fock作用一致；
普通$S^z$基也有负闭词。所有复制和对角标量都保留，无生产更新器。

<a id="r019"></a>

### R019：Kennedy–Tasaki 可正化的 zigzag ladder

> **辅助场 PfQMC：未确认。** 只确认Bell基和非局域KT相位对应的CT世界线权重。

每个rung两个spin $\boldsymbol\sigma_i,\boldsymbol\tau_i$，两者都是spin1/2而非Pauli矩阵。
开边界上取

$$\boxed{H=\sum_{i=0}^{L-2}[J_0\boldsymbol\sigma_i\cdot\boldsymbol\sigma_{i+1}
+J_1\boldsymbol\tau_i\cdot\boldsymbol\tau_{i+1}+J_3\boldsymbol\tau_i\cdot\boldsymbol\sigma_{i+1}]
+J_2\sum_{i=0}^{L-1}\boldsymbol\sigma_i\cdot\boldsymbol\tau_i.}$$

$J_0,J_1,J_3\ge0,J_3\ge J_0+J_1$，$J_2$任意实。
每spin用三个Majorana，故每rung三复模，完整Fock有$2^L$份spin复制。
Hamiltonian自身始终是邻rung局域四费米子，非局域的是证明使用的基底相位。

**显式正化。** 每rung取Bell基 $v_1=\Phi_+,v_2=\Phi_-,v_3=\Psi_+,v_4=\Psi_-$。
记 $z(v)=1$当$v=v_1,v_2$，$x(v)=1$当$v=v_1,v_3$，其余为0。
给整条链状态乘

$$U(v_1',\ldots,v_L')=(-1)^{\sum_{i<j}z(v_i')x(v_j')}.$$

这是KT变换在Bell基中的对角写法，可线性扫描计算，不能误作有限深度现场相位。
将邻rung的16维矩阵直接乘此前后相位，所有非零非对角元素成为
$-a/4,-b/4,-c/4,-d/4$，其中

$$a=J_0+J_1+J_3,\quad b=-J_0+J_1+J_3,\quad
c=J_0-J_1+J_3,\quad d=-J_0-J_1+J_3.$$

给定条件使它们非正。$J_2$在Bell基对角，含完整两个端点上的rung项也不影响正性。
因此CT任意闭路径非负。推导和边界条件见 [Nakamura Eqs. (2)–(10)](https://arxiv.org/html/cond-mat/9707019v1)。
没有把周期边界或破坏spin-reversal的场自动加入定理域。

**去重。** 相邻rung只有一个定向交叉耦合，且允许两条腿不等强；R018同层对称且两条交叉同时存在，
写成 $L=\sigma+\tau,N=\sigma-\tau$，其mixed项系数为
$(J_0-J_1+J_3)L_i\cdot N_{i+1}/4+(J_0-J_1-J_3)N_i\cdot L_{i+1}/4$；
泛型$J_3>0$时不能同时为零，所以此族不只是R018改变格点大小或各向异性参数。R015只有一个边界受挫三角，
此处三角相互作用沿链延伸。KT长字符串也不是旧§9的有限深度受控相位例子。

`ZigzagLadder`直接计算未经取绝对值的Bell基Gaussian转移和所有CT系数。
两/三/四rung、全部parity、9模完整Fock H与spin谱复制、KT逐矩阵符号核验均通过；
包含不回走的三顶点闭环、条件边界的零、$J_2$负值以及$10^{-50}$非零权重。
普通$S^z$基和$J_3<J_0+J_1$均有真正负闭词。无生产更新器。

<a id="r020"></a>

### R020：完整Fock的 projected t–J 开链

> **辅助场 PfQMC：未确认。** 只确认squeezed-order相位下的occupation世界线CT权重；未确认投影hopping和交换项的Gaussian-HS分解。

每点普通费米子 $c_{i\uparrow},c_{i\downarrow}$，
$\tilde c_{i\sigma}=(1-n_{i\bar\sigma})c_{i\sigma}$，
$\boldsymbol S_i=\tfrac12c_i^\dagger\boldsymbol\sigma c_i$，$n_i=n_{i\uparrow}+n_{i\downarrow}$，
$d_i=n_{i\uparrow}n_{i\downarrow}$。取开链最近邻：

$$\boxed{H=-t\sum_{i,\sigma}(\tilde c_{i\sigma}^\dagger\tilde c_{i+1,\sigma}+\mathrm{h.c.})
+J\sum_i(\boldsymbol S_i\cdot\boldsymbol S_{i+1}-n_in_{i+1}/4)
+U\sum_i d_i+V\sum_i n_in_{i+1}-\mu\sum_i n_i,\quad t,J\ge0.}$$

$U,V,\mu$实。投影hopping展开至最高六费米子，每项只在一条键上。
完整Fock保留doublon：它是静态阻塞态，将链切成若干段，spin交换也不能移动它。
因此这不是只算无双占据扇区，也不是有限$U$普通Hubbard的近似等价。
在无doublon扇区，它恰为 [Wang–Ye 的t–J开链](https://arxiv.org/html/1310.6496v2)。

**简短完整证明。** 在每个doublon间的开链段，按空间顺序将单占据粒子编号$p=1,2,\ldots$，
给occupation态乘 $u(s)=(-1)^{\sum_{p:\downarrow}p}$。
相邻粒子与hole交换不会改变粒子顺序；按现场up/down交错排序的CAR基，合法hopping矩阵元本来就是$-t$，
且$u$不变。相邻反向spin交换的原矩阵元为$+J/2$，交换的粒子编号相邻，
$\sum_{\downarrow}p$改变奇数，故换基后为$-J/2$。
其余全为实对角项，传播指数正。因此每个完整静态doublon扇区的任意闭词非负，其总和亦非负。
相位依赖有序粒子而非固定格点，自动包含任意hole数；无需将spin-flip的原始CT负系数删除。

**去重与边界。** $J>0$产生两个spin flavor的局域交换，投影hopping由邻胞另一flavor的占据控制，
区别于旧单粒子open hopping及仅密度相互作用。R006给定Gaussian族的守恒扇区函数没有产生这类动态投影hopping。
所有密度、长度、spin/单hole等推广按一个projected t–J构造处理。
不允许由一维证明推广周期环：三点环两个同spin粒子绕行交换，或无hole时三次反铁磁spin交换，
都产生负的完整词。

`TJChain`用真实Gaussian hopping/spin rotation加occupation projector计算signed Pfaffian，
hopping乘$t$、spin交换乘原始$-J/2$，同时保留全部对角传播。
独立完整256维CAR Hamiltonian核对所有转移和squeezed相位；随机4/6点覆盖非回文及混合hopping/spin交换闭词。
阻塞和零J单独记零，两种环上负反例均保留。当前是CT路径API，没有生产更新器。

<a id="r021"></a>

### R021：Fredkin 条件交换链

> **辅助场 PfQMC：未确认。** 只确认条件交换的occupation世界线CT权重，包含resolved及合并顶点。

采用R017的完整双复模spin编码，$P_i^{\uparrow,\downarrow}=1/2\pm S_i^z$、
$P_{ij}^{s}=1/4-\boldsymbol S_i\cdot\boldsymbol S_j$。
开链上取

$$\boxed{H=\sum_{i=0}^{L-3}[g_i^L P_i^\uparrow P_{i+1,i+2}^{s}
+g_i^R P_{i,i+1}^{s}P_{i+2}^{\downarrow}]
+b_LP_0^\downarrow+b_RP_{L-1}^\uparrow+\sum_i h_iS_i^z,\qquad g_i^L,g_i^R\ge0.}$$

边界及$h_i$为实数。此H在完整Fock空间作用，每点两份spin复制；没有把迹投影到Dyck路径。
每个顶点支撑三个相邻spin，最高六费米子。
定义按 [Fredkin模型 §2 Eqs. (1)–(2)](https://arxiv.org/html/1605.03842) 的归一化projector，
不是忽略整体归一化的Pauli重写：$P_i^\uparrow P_{jk}^{s}=(1+\sigma_i^z)(1-\boldsymbol\sigma_j\cdot\boldsymbol\sigma_k)/8$。

**证明。** 原始$S^z$基中，singlet projector唯一非对角转移为
$|\uparrow\downarrow\rangle\leftrightarrow|\downarrow\uparrow\rangle$，矩阵元$-1/2$。
乘现场control的0/1及非负$g$后，H的非对角项仍非正。将全部对角项放入精确传播，
每个CT事件的权重为$g/2$或0，闭词乘积非负。
若多条Hamiltonian项产生同一跃迁，矩阵元为这些系数之和，仍非负。
每个parity复制使用同一论证，完整费米迹也正；Dyck零能基态仅是独立谱检查，不承担有限温正性证明。

**去重。** 条件单spin控制邻键交换使$Z_0X_1X_2$的系数为$-g/8$，
在同一spin–Majorana编码中是非零grade6项。
R017的不相交J-Q matching只有grades0、4、8、12（加场有2），线性CAR/Bogoliubov变换保留grade，
不能将此H归成现有J-Q族。它也不是R013的四点双pair运动。
等耦合时某些纯对角三密度项会抵消，因此不用这些已抵消项作为升阶证据。
所有有色/Fredkin条件变体统一计一个构造。

`FredkinModel`保留完整对角能量、resolved与合并跃迁两种CT权重；Gaussian旋转和occupation projector均走signed Pf。
核对262144个Fock/Pauli项矩阵元、320个完整二顶点词、spin谱复制及四spin Dyck零能态。
负对照将一个控制项耦合改为$-2$，六spin闭环$50\to42\to44\to52\to50$（数字为down位串）
在合并所有同跃迁顶点后系数依次为$-1/2,1,1/2,1/2$，仍为真实负配置；
完整resolved词的Pf之和等于合并词的负Pf，不能靠重新合并消掉此反例。
另有阻塞零与两个half-angle脉冲的符号检查。无生产更新器。

<a id="r022"></a>

### R022：PXP Rydberg 条件单spin翻转链

> **辅助场 PfQMC：未确认。** 只确认occupation世界线CT路径和SSE形式的kinetic词；Gaussian投影/旋转的signed Pfaffian不是常规HS场验证。

每点普通复费米子 $a_i,b_i$，定义

$$Z_i=1-2n_{a_i},\quad n_i=(1+Z_i)/2=1-n_{a_i},\quad P_i=1-n_i=n_{a_i},$$
$$X_i=a_i^\dagger b_i+b_i^\dagger a_i+a_i^\dagger b_i^\dagger+b_ia_i.$$

完整四态是两份spin1/2，$X_i$同时翻转$a_i,b_i$而保留现场parity。
开链上，缺失的边界control取$I$：

$$\boxed{H=-\sum_i\Omega_iP_{i-1}X_iP_{i+1}
+V\sum_i n_in_{i+1}+\Delta\sum_i n_i,\qquad\Omega_i\ge0,\quad V,\Delta\in\mathbb R.}$$

每项支撑至多三点，最高六费米子。这里$n_i$是被编码spin的激发数，不是$a,b$总费米数。
完整Fock空间和完整spin空间都保留；相邻激发造成的冻结区不会被隐含删去。
模型定义及full-space分解见 [Turner等 Eq. (5)、§II.2](https://arxiv.org/html/1806.10933v2)。
与原文正$X$号的差异由$\prod_iZ_i$现场相位消去；不是不同计数。

**证明。** occupation基中$X_i$的允许矩阵元均为$+1$，两个control都是对角0/1。
所以$H$的非对角元为$-\Omega_i$或0，对角传播为正实指数。
任意闭合CT路径权重是$\prod\Omega_i$乘这些正指数；每个现场parity复制相同，完整迹非负。
冻结激发簇对应完整保留的配置，亦满足同一论证。
若$\Omega_i$单独反号，每条闭路径上该点翻转偶数次，仍不能制造负环；实现明确检查这一点，
不把普通site-dependent耦合号当成必然有sign问题的反例。

**排重。** 本构造动态改变总激发数；Fredkin条件交换守恒它。
其controls会被邻近顶点改变，不能将它们当R006或FK的静态变量。
三个不同点的$Z_{i-1}X_iZ_{i+1}$系数为$-\Omega_i/4$，
既不是单/双点交换的各向异性参数，也不是旧§9仅改变矩阵元相位的有限深度twist。
相邻PXP顶点不对易，不能作为R016 commuting-stabilizer族的新名字。
所有facilitation/尺寸/激发扇区和符号约定只计这一种条件单spin翻转构造。

`PXPModel`实现occupation世界线路径，并另外直接用partial Gaussian density projector及
$X=-i\exp(\pi\gamma_1\gamma_2/2)$两个半角脉冲计算完整kinetic词；全部$-i$和
$\beta^k/k!$保留。含实对角项的路径接口使用指定时刻的密度，不重复乘该积分因子。
完整八模Fock谱是四spin谱的16份复制；枚举四点所有$k\le4$的341个有序kinetic词，
核对Fock乘积、非回文多点路径与实际零权重。

负对照给site0的振幅额外乘$(1-2n_2)$：激发位串
$0\to1\to5\to4\to0$的四次翻转产生负环。
它改变完整H的条件符号，实际路径Pf及完整kinetic词都为负；另保留原始Gaussian Spin-lift负号控制。
没有只验证一个不闭合跃迁的负幅度。当前无生产更新器。
