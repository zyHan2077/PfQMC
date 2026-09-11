# 持续搜索队列与排重证据

**辅助场口径校正：新增21类中，仅R001–R005、R009共6类确认常规Gaussian-HS辅助场PfQMC正性；其余15类未确认。**
“verified”、21类总数及历史全部sign统计涵盖不同表示，不能当作21类辅助场PfQMC结果。
逐例标签和旧目录中的静态/非Gaussian换基/非厄米例外见 [确认状态](MODELS.md#auxiliary-field-status)。

本文件保留候选推导及排重过程；每项开头的当前状态优先于其历史筛选笔记。
正式登记与计数见 [MODELS.md](MODELS.md) §12 和 [search_registry.json](search_registry.json)。
用户确认：只数本轮新找到、去除明显重复后的 Hamiltonian 构造；超过 20 后停止新增，统一合并推导。

## Q1. Quaternion-Hermitian 轨道

当前状态：已通过正式 signed Pfaffian 验证，与 Q2 合并登记 R004。

K 是 n×n quaternion-Hermitian 无迹矩阵，作用在 Herm_H(n) 的实轨道空间：
`h(K): X -> KX+XK`，轨道数 M=n(2n-1)。普通复费米子 H 为
`-mu D-Q(h(K0))-u D^2/2-sum lambda_a Q(h(Ka))^2/2`，允许固定胞的相邻共享平方。
所有 u/lambda>=0；任意实方向系数。不是把物理费米子本身定义成四元数。

普通复矩阵实现：每个 quaternion 用
`q0=I2, q1=diag(i,-i), q2=[[0,1],[-1,0]], q3=[[0,i],[i,0]]`。
Herm 基的每个对角放 q0，每个 i<j 的上块为 q_alpha/sqrt2，下块为其 adjoint。
内积是 ReTr(Ta^dag Tb)/2；物理 h_ab=ReTr[Ta^dag(KTb+TbK)]/2，为实对称。

复嵌入 G 的谱配成 lambda_i,bar(lambda_i)，本表示复化为 exterior^2 G。
任意 R>0 时，det(I+R rho(G)) 等于
`prod_i(1+R|lambda_i|^2) prod_i<j |1+R lambda_i lambda_j|^2 |1+R lambda_i bar(lambda_j)|^2`。
因此完整迹乘 R^(-M/2) 非负。Jordan 情形也适用。

n=2,M=6 与 Lorentz m=5 **重复**，只验证等价，不计新结果；从 n>=3 起登记。
此处是 SU*(2n) 的外二次幂实形式，不是原 SL(2n,R) 的中间外幂。
加入独立正负 D 和 generic 对角 K 后，物理谱 2kappa_i 的重数为1，
kappa_i+kappa_j 的重数为4；共同密度平移到正数后存在奇重实 Majorana 本征值，排除共同 Wei。
正负不对易 K 也排除共同 Metzler。

agent 独立 determinant/congruence 筛选 n=2,3,4 共600配置通过；还不是仓库 signed-Pf 验证。
数学来源：[Loring, Factorization of Matrices of Quaternions](https://arxiv.org/html/1107.0500)。

## Q2. Quaternion anti-Hermitian 轨道

当前状态：已通过正式 signed Pfaffian 验证，与 Q1 合并登记 R004，不另计一项。

相同的 K 改作用于 antiHerm_H(n)，轨道数 M=n(2n+1)，物理 h 同上仍实对称。
基的对角放 q1,q2,q3；i<j 上块 q_alpha/sqrt2，下块 **负** adjoint。
复化是 Sym^2 G，正性公式在 Q1 乘积基础上再乘 `prod_i |1+R lambda_i^2|^2`。
正负 D 与 generic K 给 2kappa_i 的重数3、其余重数4，仍与 Wei 偶重数矛盾。

n=1 是可对易密度退化例子，不计；n>=2 可以考虑。
n=2,M=10 可能等价于六维 Lorentz exterior^3 的一个 Hodge 分支；若以后引入该模型要合并，
不能反向因为维度不同就声称不重复。agent 又做了600配置 determinant/congruence 筛选，均正；待signed Pf。
Q1/Q2可视作同一底层quaternion构造的两个不等价表示分支，最终是否合并统一审查。

## Q3. 正超选择扇区与局域充电／奇偶相互作用

当前状态：已通过正式 signed Pfaffian、独立 sector Fock、Hamiltonian／Trotter 验证，登记 R006。

TN 正谱单胞的固定物理粒子数 k 迹为
`exp(-sum tr h/2) e_k(B)>0`，e_k 是基本对称多项式。
可以加入任意有限支撑真实函数 F({N_r})，例如
`U sum(N_r-m)^2 + V sum_<rs> N_r N_s`，U可正、V任意，或硬约束N_r=m。
整数 k_r 是静态扇区，但内部仍有时间变化、不对易的Cartan HS；不是单个静态自由H的FK公式。
要用 fugacity 多项式系数／完整Fock扇区迹验证，**不能把傅里叶投影的每个复相位项当成正配置**。

整数spin单胞在mu=0时还有正奇偶扇区：B恒有本征值1，所以det(I-B)=0，
偶／奇迹各是完整迹的一半。允许任意局域 `J sum p_r p_s + v sum p_r`，p_r=(-1)^N_r。
固定具体N可能负；共同mu非零也可能使一个奇偶扇区负。
这两条宜合为一个“正扇区相互作用扩展”，不按每个F或投影方式凑数。
背景：[Canonical AFQMC](https://arxiv.org/abs/2212.08654)。尚未接入仓库signed Pf。

## Q4. 局域 Gaussian 指数顶点／Poisson 相互作用

在相邻三轨道spin1单胞上取
`H=-sum_<rs> kappa_rs exp(a_rs Sx_r+b_rs Sz_s)`，kappa>=0，a/b实；
可同时加入x/z交换的另一组键，使不同项不对易。
CT展开每个顶点是原保正族的Gaussian指数，权重为正的kappa乘积乘signed Gaussian迹。
因为完整Fock中 S^3=S，顶点等于
`[1+sinh(a)Sx_r+(cosh(a)-1)(Sx_r)^2] [1+sinh(b)Sz_s+(cosh(b)-1)(Sz_s)^2]`。
是明确有限范围的六/八费米相互作用；单胞版本只是原负平方的重写，不另计。
已有指数二次算符designer Hamiltonian方法：
[Huffman–Chandrasekharan PRD96 114502](https://journals.aps.org/prd/abstract/10.1103/PhysRevD.96.114502)，
[论文第5章](https://arxiv.org/abs/1912.11776)。不是新基本机制，待signed Pf和物理H验证。

## Q5. 动态 Ising–fermion 三轨道模型

二分最近邻图，每点(c,a,b)三个复模，Dab=na+nb-1，P=a^dag b^dag+ba：

`H=-t sum(c_i^dag c_j+h.c.) + V sum dnci dncj -J/4 sum Dabi Dabj +1/2 sum_i h_i dnci P_i`。

V>=0，t/J/h实（h可异号），全部作用至多四费米且onsite/最近邻。
令 Sz=Dab/2,Sx=P/2。偶占据子空间是spin1/2，两个单占据态是静态spin0阻塞；
无需投影即可把完整Fock空间写成带自旋空缺的原模型扇区之和。
文献正性是CT-INT加自旋世界线，不是普通独立Gaussian HS：
[Huffman–Chandrasekharan §II Eq1–10](https://arxiv.org/html/1605.07420v2)。

按子格sigma_i作局部自旋旋转，使顶点为h_i sigma_i dnc_i Sx_i。
H0含c自由hop和-J SzSz，density顶点带CT-INT负号。
`dn=(-i/2) exp(i*pi*dn)`，指数拆为两次pi/2脉冲以避单片零迹；
Sx用pair creation/annihilation旋转乘纯占据投影，另乘1/2。
自旋对角传播给正数；每点闭合要求偶数次flip，从而h_i和子格号的总因子可正确处理。
若阻塞或路径不闭合，总权重为0，不能填sign1。

agent 的独立screen脚本已保存到 [research/dynamic_ising_screen.py](research/dynamic_ising_screen.py)，
原始日志 [research/dynamic_ising_screen.log](research/dynamic_ising_screen.log)。脚本不导入本库，
本轮暂作待移植证据。两点6mode全Fock词迹约4.88983042758992e-6；
方环4site正配置约4.649276296952922e-5；32随机配置正。
未合并成正式验算目标之前不计数。原审查脚本有硬编码子格和部分打印式核对，
移植时应改为一般二分图、断言和可复现逐配置记录，避免照搬为生产算法。

## 合并／排除记录

- conjugate_tensor2 与 Lorentz4：Pauli基酉等价，实际逐矩阵核对，不计两次。
- 对称幂 d=2 与整数spin、q=1与开链：重复。
- 一般 Schur 标签合计一个多项式表示族；d3 hook(2,1)=det×adjoint，
  无迹时其八复模迹是已有四复模配对伴随迹的平方，不因8与4不同就新增。
  d4 hook20可作为对R003的推广候选，不能每个Young图增加计数。
- 不等质量 Hubbard 的 split-orthogonal hopping HS 已由 Wei 覆盖：
  [Wang 等 Eq7,9](https://arxiv.org/html/1506.05349v2) 与 [Wei §V](https://arxiv.org/html/1712.09412v3)。
- 2021 meron键 A=-K-K^2=-2P_(K=1)，双自旋键是实二分Gaussian指数的正标量极限；
  加相应real spinflip Hubbard HS仍可纳入O(n,n)/Wei闭包，不能只看meron名称判成新类：
  [原式 Eq3–6](https://arxiv.org/html/2011.13208v2)。
- crossed horizontal/vertical open chains + onsite density coupling：原开链直和／共享场的实例，重复。
- 全部文献新颖性声明限于本目录比较；检索无命中不是首创证明，也不证明排除了所有未知HS。

## 第二轮状态补充

动态 Ising Q5 已登记 R007：完整三轨道 Pfaffian、Fock 及路径核验已通过。
Lorentz 无迹对称轨道登记 R005，N=5/6 和邻胞验证通过；其他 Schur 标签不另计。
Q1/Q2 之前关于待 signed Pf 或待合并的文字属于历史阶段，当前均以 R004 为准。
当前 7 项，详见 SEARCH_RESULTS.md 第二轮。Q4 Poisson 局域指数顶点在独立构建已通过，等待正式集成。
Albert 27 轨道的 exceptional Jordan 构造正在实现及核对群论证明；尚不计数。

## 第三轮状态补充

Q4 已登记 R008；各轴、反号顶点和密度顶点合并为一个高阶局域指数构造。
Albert 已登记 R009，完整 26 个 traceless Jordan 方向的 signed Pfaffian／Jordan oracle／真实 sweep 通过。
当前 9 项，以 registry 为准。第二轮状态段的 7 是历史快照。
SU(1,n) adjoint 的 ad_K 经第二张量因子 J 换基成为 R002 conjugate tensor 限制，
去掉不变 J 只是删除正 spectator；不计新增。
移动 XXZ pair 与 c-density 翻转耦合正在实现；Heisenberg/XY/XXZ 不按各向异性分开计数。
Singlet-projector/条件密度模型可用 meron 重求和，也可 spin 旋转加 Marshall 的 worldline 分解；
正在区分整个 Hamiltonian 是否新构造，不能只因某个 c 子系统满足 Wei 就判整个模型重复，
也不能只因换了展开名称就计新结果。

## 第四轮状态补充

MobileSpin 已登记 R010，XXZ/Heisenberg/XY合计一项；QuantumLink一维完整Fock correlated hopping登记R011。
当前11项；以上较早的7/9均为历史快照。meron E投影方案独立构建通过，待root正式集成；
其与旋转后的动态pair/条件staggered势为同一个Hamiltonian，不按求和表示计数。
受挫impurity的Majorana-spin编码、SO3 color-singlet投影磁体、四点pair环交换仍待核验。
R001跨胞hopping允许性有明确反例，见 research/lorentz_cross_cell_screen.py/.json；不新增计数。

## 第五轮状态补充

R012动态singlet/条件密度、R013四胞pair环交换、R014颜色singlet/permutation、R015受挫impurity已正式验证登记。
当前15项；既有计数段是历史快照。下一轮核验commuting stabilizer、J-Q多键投影、受挫bilayer、zigzag ladder。
Open spin1 bilinear-biquadratic的KT参数区使用R014相同P/Pi算符，严格按参数扩展合并，不再计数。

## 第六轮与严格合并

R016 stabilizer、R017 J-Q、R018 bilayer、R019 zigzag、R020 projected t-J、R021 Fredkin全部正式验证。
R007的H精确等于R010的Jxy=0,Jz=-J参数面；不同h符号域不改变这个恒等式，最终只计一项。
当前20项，目标仍进行中。PXP Rydberg条件单spin翻转待正式集成；Compass不作为新计数。

## 最终去重与停止新增（2026-09-11）

R022 PXP通过主构建、独立Fock及signed Pf两种子验算。最终21项，已停止新增。
历史“待核验/进行中”段落只记录当时状态，当前以registry及以下结论为准。

- R007精确并入R010的Ising参数面；任意h符号的Ising条件保留，两个API不删。
- R004四元数Hermitian/anti分支合一；R006 charge/parity/hard约束合一。
- R001低维与R002 tensor d=2、R003 d=2/q=1、R004 Herm2/Lorentz6只作等价检查。
- R012 meron和旋转条件势是同一H；R016所有stabilizer、R017 Q2/Q3、R018三个面分别只计一项。
- R006只扩展指定Gaussian族的逐因子守恒N/p及静态F(N,p)，不吞并任意interacting sector。
- R018 bilayer仅LL/NN；R019泛型zigzag有非零LN/NL，非参数面重复。
- R015配对双臂内部Jperp=Kz=Kperp=lambda/2全非零，不满足R018的三个条件；需一维闭线证明。
- R017无Majorana grade6，R021 Fredkin具有非零grade6；R022 PXP则为非守恒磁化的条件单spin翻转。
- Compass及open spin1 bilinear-biquadratic KT参数区不新增；前者交换轴变体歧义较强，后者归入R014算符族。
- 高阶交互、世界线正性、非Gaussian顶点均不证明不存在其他Wei分解；排除仅针对具体给定辅助场族。

没有将这些结果宣称为文献首创，也没有声称任意多体幺正等价性已被穷尽排除。
