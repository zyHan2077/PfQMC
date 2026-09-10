# 四复模厄米配对模型：验算 sign free

实现来源：`/home/han/notes/sign_free_classification/03_Wei之外的具体无符号族.md` §3.2。
本目录实现该节四复模模型及其实 Gaussian HS 分解，不包含 §3.4 的二维模型。
验算日期：2026-09-10。实际测试数据和结论见 [REPORT.md](REPORT.md)。

## 模型和约定

Fock 基用整数 `0,...,15` 编码，bit `j-1` 表示模式 `j` 的占据；
`c_j` 的费米子符号为其前方已占据模式数的奇偶性。

\[
O_1=n_1+n_2-1,\qquad O_2=n_2+n_3-1,
\]
\[
B_1=c_2^\dagger c_3-c_4^\dagger(c_1-c_1^\dagger)+\mathrm{h.c.},
\quad
B_2=-c_1^\dagger c_2+e^{-i\pi/3}c_4^\dagger(c_3-c_3^\dagger)+\mathrm{h.c.}
\]
\[
Q(u,v,t_1,t_2)=uO_1+vO_2+t_1B_1+t_2B_2,
\qquad u,v\in\mathbb R,\quad t_1,t_2\ge0.
\]

计算完整 Fock 迹 `Z = Tr(exp(Q_L)...exp(Q_1))`，不投影粒子数或宇称。
配对创建与湮灭必须同时保留，复相位及不同物理项之间的系数关联不能自行修改。

## 三条独立计算路径

1. 在 **16×16 Fock 空间**直接构造上述费米子算符，通过厄米谱分解求每片指数并相乘，最后取完整迹。
2. 在 **3×3 辅助空间**独立构造
   \[
   X=\begin{pmatrix}(2u+v)/3&t_1&0\\t_1&(v-u)/3&t_2\\0&t_2&(-u-2v)/3\end{pmatrix},
   \quad g=e^{X_L}\cdots e^{X_1},\quad Z=2(\operatorname{tr}g\operatorname{tr}g^{-1}-1).
   \]
   `g^{-1}` 通过逆指数的反序乘积生成，避免直接求逆三阶长乘积。
3. 以 `gamma_(2j)=c_j+c_j†`、`gamma_(2j+1)=-i(c_j-c_j†)` 构造 **8×8 Majorana 矩阵**。
   由 `Q=gamma^T A gamma/4` 得 `A_ab=-Tr(gamma_a gamma_b Q)/8`。
   程序逐个重构原始 Fock 算符，检查归一化。传给仓库的 Hamiltonian 为 `H=-A`。
   调用原有 `pfaffianForEta`、`pfaf` 递推完整带符号权重，并另用
   `pfaffianForSignOfEta`、`pfaffianForSignOfProduct` 跟踪相位。

第 3 条直接复用仓库 `src/skewMatUtils.cpp` 和 PFAPACK，并按
[PfQMC 论文式 (4)–(5)](https://arxiv.org/html/2408.10311v1)进行乘积合成。
注意仓库函数 `pfaffianForEta` 实际返回完整迹 `2^N eta`，这里 `N=4`。
因此每次合成需要除以 16，块 Pfaffian 的 `(-1)^N` 在本例为 +1。
不通过 `sqrt(det(I+B))` 决定符号，不取权重绝对值，也不将微小虚部提前抹掉。

该程序是独立辅助场配置的权重验证器，不是新增了带 Metropolis 更新和物理量测量的完整模型 walker。
16 维精确 Fock 计算提供独立参照；原有 PfQMC walker 的回归由仓库原测试负责。

## 无符号性的数学理由

定义 `E_i=(B_i+[O_i,B_i])/2`、`F_i=(B_i-[O_i,B_i])/2`，
`H_1=2O_1-O_2`、`H_2=-O_1+2O_2`。这些算符满足 `sl(3)` 的
Chevalley–Serre 关系。程序检查全部 Cartan、升降及 Serre 关系。
对角生成元 `uO_1+vO_2` 的权重为四个零，以及 `±u,±v,±(u+v)` 各两次，
所以 16 维 Fock 表示为两个八维伴随表示，得到上面的完整迹公式。
这里取实际表示的迹，没有平方根分支选择。

令 `J=diag(1,-1,1)`。所有允许的 `X_l` 和 `-J X_l J` 均为 Metzler 矩阵，
故指数逐项非负。保留在第 `i` 个辅助位置始终不移动的路径，令
`s_i=sum_l (X_l)_ii`，可得
\[
g_{ii}\ge e^{s_i},\qquad (g^{-1})_{ii}\ge e^{-s_i}.
\]
由于 `sum_i s_i=0`，有 `tr(g)>=3` 和 `tr(g^{-1})>=3`，从而
\[
\boxed{Z=2(\operatorname{tr}g\operatorname{tr}g^{-1}-1)\ge16.}
\]

这是任意有限片数和允许参数下的证明；随机抽样用于排查实现错误，不能替代全称证明。
“共同 Wei 判据之外”的排除论证仍按原笔记：允许任意正负对角方向，而一般对角方向的
伴随表示有简单实特征值，与共同 Wei 结构要求的偶重数矛盾。
该结论不意味着超出所有已知无符号机制，也不代表文献原创性已经核实。

## 真实相互作用与 Gaussian HS 场

取
\[
H=-Q_0-U\delta n_1\delta n_2-V\delta n_2\delta n_3,\quad
Q_0=Q(0.2,-0.15,0.4,0.3),\quad U,V\ge0.
\]
对 `dt=beta/L`，每片采用对称分解
\[
e^{dtQ_0/2}\,
e^{\sqrt{dtU}\phi_lO_1+\sqrt{dtV}\psi_lO_2}\,
e^{dtQ_0/2},\qquad \phi_l,\psi_l\overset{\mathrm{iid}}\sim N(0,1).
\]
Gaussian 平均之外的常数因子是 `exp[-beta(U+V)/4]`。
程序用 64 节点标准正态 Gauss–Hermite 求积，独立检查
`exp[-a/4] E exp(sqrt(a) phi O_i) = exp(a delta_n_i delta_n_(i+1))`，
覆盖 `a=0,0.7,8`。这里不是把随机场直接塞进一个错误的合并指数。

CSV 中的 `Z` 是不含 Gaussian 密度和该常数的完整费米子迹，满足 `Z>=16`。
相对于标准正态乘积测度的物理配置权重是
`positive_HS_factor * Z`，仍严格为正；相对于 Lebesgue 测度还需乘
`(2*pi)^(-L) exp[-sum_l(phi_l^2+psi_l^2)/2]`。
加入这些因子后不能继续声称物理权重本身下界是 16。
有限 `dt` 的 Trotter 误差与无符号性是不同问题；本次没有估计热力学可观测量或做连续时间外推。

## 随机抽样与对照

| CSV ensemble | 参数 |
|---|---|
| 0 | 单位元，4 片，精确 `Z=16` |
| 1 | 纯对角边界，4 片 |
| 2 | 仅第一条辅助边非零，4 片 |
| 3 | 每片 `u,v=s*N(0,1)`，`t1,t2=s*Uniform(0,1)`，`s=0.35/sqrt(L)` |
| 4 | 同上，`s=0.8/sqrt(L)` |
| 5 | Gaussian HS，`beta=1,U=0.7,V=1.1` |
| 6 | Gaussian HS，`beta=4,U=2,V=1.5` |

随机组的片数按样本序号轮换 `L=1,8,32,64`；HS 每个物理时间片实际有 3 个 Gaussian 算符。
`--samples N` 表示每个随机组 N 例，总计 `4*N+3`。
所有随机配置独立生成，没有按权重拒绝样本。
随机数为 `std::mt19937`，种子、样本编号、组别保存在 CSV；精确重现正态序列需使用相同 C++ 标准库。

另固定种子 731 构造允许族之外的实对称、无约束跃迁模型，要求发现负 Fock 迹且
Pfaffian 给出相同负号。删除配对的测试还要求迹确实变化，避免把只创建、不影响迹的三角构造混入本例。

断言包括：Fock 迹实部 `>=16-1e-9`，相对虚部 `<1e-9`，
三阶公式相对误差 `<1e-8`，Pfaffian 相对误差 `<1e-7`，相位偏差 `<1e-8`。
相对误差分母为 `max(16,abs(Z_Fock))`；非有限数直接失败。失败返回非零退出码。
程序用双精度和直接矩阵乘积/求逆，适合这里的小系统参数范围；没有实现原 walker 的长虚时间 UDT 稳定化。
更强场或更长虚时间下的数值失败需要先区分病态误差与真实负权重，不能静默改符号。

## 编译与复现

本机 WSL2 环境、依赖路径及恢复步骤统一记录在 [CONTEXT.md](../CONTEXT.md)。

依赖与原仓库相同：Intel oneAPI C++/Fortran/MKL/MPI、Eigen 3.4、GoogleTest、CMake。
在仓库根目录的 Bash 中执行（Eigen 与 GTest 路径按环境修改）：

```bash
source /opt/intel/oneapi/setvars.sh
export PATH="$PWD/build/python-tools/cmake/data/bin:$PATH"
make -C inc/pfapack/fortran mFC=ifx
make -C inc/pfapack/c_interface
cmake -S . -B build -DCMAKE_CXX_COMPILER=mpiicpx -DEIGEN3_INCLUDE_DIR="$PWD/build/eigen" \
  -DCMAKE_PREFIX_PATH=/home/han/miniconda3
cmake --build build -j 2
ctest --test-dir build --output-on-failure
build/sign_free_check --samples 10000 --seed 240810311 \
  --csv 'sign_free_verify/samples_240810311.csv'
build/sign_free_check --samples 1000 --seed 910 \
  --csv 'sign_free_verify/samples_910.csv'
```

本机 CMake 安装在 `build/python-tools/cmake/data/bin/`，可将以上 `cmake`、`ctest`
分别替换为对应绝对或相对路径。Eigen 源码位于 `build/eigen`，tag `3.4.0`。
构建工具下载方式（只需首次执行）：

```bash
/home/han/miniconda3/bin/python -m pip install --target build/python-tools cmake
git clone --depth 1 --branch 3.4.0 https://gitlab.com/libeigen/eigen.git build/eigen
```

CTest 的 `SignFree.FourMode` 每随机组取 24 例，并执行边界、代数、HS 与负权重对照。
详细随机样本不依赖 Python/NumPy/SciPy；新程序和原程序共用 C++/MKL/PFAPACK 构建链。
