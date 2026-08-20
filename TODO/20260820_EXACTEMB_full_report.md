# OpenMolcas EXACTEMB 精确电子密度静电嵌入：理论、实现、使用与数值验证

日期：2026-08-20

代码：OpenMolcasUP_ExactDenEmb，分支 `agent/exact-density-embedding`，基线提交 `65a3eeb11dc` 加本文所述的本地修正

验证平台：GV800，用户 `alpha`

验证体系：相距 3.5 Å 的平行乙烯二聚体，cc-pVDZ，$C_1$，每个中性片段采用 CAS(2e,2o)

## 摘要

`EXACTEMB` 在 OpenMolcas 内部由一个片段的最终 RASSCF 一阶约化密度矩阵构造它对另一片段产生的电子库仑势，并把该势作为一电子算符加入目标片段的 RASSCF 哈密顿量。

环境片段的核吸引势由 SEWARD 的标准 `XField` 功能提供，环境电子的库仑排斥势由 `EXACTEMB` 提供，因此目标片段感受到环境片段的完整经典静电势，而不需要把电子密度拟合为原子电荷或有限阶多极矩。

二电子积分通过 OpenMolcas 现有的 Cholesky/RI-CD 表示收缩，实际计算不依赖 PySCF；PySCF 只在小体系验证中用于直接四中心积分参考。

本实现已经在中性乙烯二聚体上完成无嵌入回归、清除功能、电子数与 AO 映射检查、双向库仑势构造、PySCF direct-J 对照、仅电子库仑势响应、完整静电嵌入以及六轮相互极化迭代，七项验证全部通过。

本文从量子化学原理、OpenMolcas 数据表示、程序实现、输入文件写法、数值验证和新体系复现要求六个方面给出独立完整的说明。

## 1. 方法的物理内容

### 1.1 要解决的问题

考虑由两个片段 $A$ 和 $B$ 组成的体系。

目标是分别优化片段 $A$ 和片段 $B$ 的多参考波函数，同时让每个片段的电子在另一个片段产生的静电场中弛豫。

以优化片段 $A$ 为例，片段 $B$ 的静电影响包含两个部分：片段 $B$ 的原子核对 $A$ 电子的吸引，以及片段 $B$ 的电子密度对 $A$ 电子的库仑排斥。

本方法保留片段 $B$ 的完整 AO 一阶密度矩阵，不把它压缩为点电荷、偶极矩或四极矩，因此能够保留有限距离下与电荷穿透有关的空间信息。

这里的“精确电子密度”是指在选定基组和选定 RASSCF 波函数内使用完整的一阶密度矩阵；二电子积分在实际实现中仍采用可控阈值的 Cholesky/RI-CD 表示。

### 1.2 术语与计算对象

本文把提供环境电子密度的片段称为“源片段”，把当前接受静电势并进行波函数优化的片段称为“目标片段”；从 B 向 A 构造势时，B 是源片段，A 是目标片段。

“片段计算文件组”是指由片段自身 GATEWAY、SEWARD、SCF 和 RASSCF 产生的几何、基组、积分、轨道与 RunFile；“二聚体积分文件组”是指完整 AB 几何上的 RunFile、ONEINT 和 Cholesky 文件。

`EXACTEMB` 必须在二聚体积分文件组处于当前状态时执行，因为跨片段库仑积分属于完整 AB 基组；最终目标 RASSCF 必须重新建立目标片段的积分文件，因为它只优化目标片段波函数。

RASSCF 是限制活性空间自洽场方法，CASSCF 是其活性空间没有 RAS1/RAS3 截断的特例；本文乙烯计算使用 CAS(2e,2o)，即两个活性电子分布在两个活性轨道中。

源 `D1ao` 是自旋求和的一阶 AO 密度，包含冻结、非活性和活性电子的总密度，而不是只有活性空间密度。

### 1.3 嵌入哈密顿量

片段 $A$ 的嵌入一电子哈密顿量为

$$
\hat h_A^{\mathrm{emb}}=\hat h_A+\hat V_{\mathrm{nuc}}^B+\hat J[P_B].
$$

$\hat h_A$ 是孤立片段 $A$ 的动能与自身核吸引算符，$\hat V_{\mathrm{nuc}}^B$ 是环境片段 $B$ 的原子核对 $A$ 电子产生的吸引势，$\hat J[P_B]$ 是由 $B$ 的一阶电子密度产生的库仑排斥势。

在 AO 基中，电子库仑势写为

$$
J_{\mu\nu}^{B}=\sum_{\kappa\lambda\in B}(\mu_A\nu_A\mid\kappa_B\lambda_B)P_{\lambda\kappa}^{B}.
$$

指标 $\mu,\nu$ 属于目标片段 $A$，指标 $\kappa,\lambda$ 属于源片段 $B$，$P^B$ 是从源片段最终 RASSCF RunFile 读取的一阶 AO 密度矩阵。

该算符是普通的一电子算符，因此必须在 RASSCF 构造核心哈密顿量、非活性 Fock 矩阵、活性空间一电子积分、轨道梯度和 CI 哈密顿量之前加入，而不能只在最终总能量上附加一个常数修正。

### 1.4 核与电子贡献必须各出现一次

环境原子核通过目标片段 SEWARD 输入中的 `XField` 点电荷提供。

对全电子计算，每个环境原子在其真实坐标处放置等于核电荷 $Z_I$ 的正点电荷，SEWARD 随后产生电子－环境核吸引积分，并把目标片段原子核与这些外部正电荷之间的排斥常数计入总能量。

环境电子通过 `EXACTEMB` 产生的 `Exact Emb Pot` 记录提供。

如果只加入 `XField`，则缺少环境电子排斥；如果只加入 `Exact Emb Pot`，则缺少环境核吸引；如果两种贡献中的任何一种被重复加入，能量与轨道响应都会错误。

本文验证中的“完整静电嵌入”始终表示 `XField` 与 `Exact Emb Pot` 各加入一次。

对于固定环境密度，目标片段总能量中的跨片段静电部分可写为

$$
E_{A\leftarrow B}^{\mathrm{elst}}=\operatorname{Tr}[P_A V_{\mathrm{nuc}}^B]+\operatorname{Tr}[P_A J[P_B]]+\sum_{I\in A}\sum_{J\in B}\frac{Z_IZ_J}{R_{IJ}}.
$$

这里的电子－电子交叉项不带 $1/2$，因为 $P_B$ 被视为给定的外部环境密度；这些片段嵌入能量用于波函数优化和同方法对照，不应直接当作完整二聚体总能量。

### 1.5 Cholesky/RI-CD 库仑收缩

OpenMolcas 将二电子排斥积分近似表示为

$$
(\mu\nu\mid\kappa\lambda)\approx\sum_Q L_{\mu\nu}^{Q}L_{\kappa\lambda}^{Q}.
$$

先将源密度与 Cholesky 向量收缩，

$$
d_Q=\sum_{\kappa\lambda\in B}L_{\kappa\lambda}^{Q}P_{\lambda\kappa}^{B},
$$

再构造目标库仑矩阵，

$$
J_{\mu\nu}^{B}=\sum_Q L_{\mu\nu}^{Q}d_Q.
$$

实现复用 `Fock_util_interface::CHORAS_DRV`，并设置 `ExFac=0`，因此只计算库仑项 $J$，不计算交换项 $K$。

Cholesky/RI-CD 阈值决定积分近似误差；乙烯的 `Medium Cholesky` 验证给出相对于直接四中心积分参考约 $1.8\times10^{-8}$ 的矩阵相对误差。

### 1.6 方法中不包含的物理项

当前版本不包含片段间交换、反对称化、Pauli 排斥投影、Huzinaga 投影、非加和动能、非加和交换相关泛函、色散、显式片段间电荷转移以及解析核梯度。

由于目标片段 RASSCF 只在目标片段自身 AO 空间中优化，环境片段 AO 不进入该片段的变分空间，所以片段间电荷转移不在这一步中描述。

该方法的作用是制备受到另一片段静电极化的片段波函数，而不是给出完整的二聚体相互作用能分解。

如果这些片段波函数随后用于 NOCI 或其他片段态耦合计算，最终矩阵元必须使用不含 `XField` 和 `Exact Emb Pot` 的真实二聚体哈密顿量；嵌入势只用于制备片段态，不能再次进入最终物理哈密顿量。

### 1.7 相互极化与 freeze-and-thaw

一次嵌入使用环境片段的固定密度，例如用孤立 $B$ 的密度优化 $A$。

相互极化则重复执行

$$
P_A^{(n+1)}=\mathcal F_A[P_B^{(n)}],\qquad P_B^{(n+1)}=\mathcal F_B[P_A^{(n)}],
$$

其中 $\mathcal F_A$ 表示在 $B$ 的固定静电势中完成一次片段 $A$ 的 RASSCF 优化。

本文使用 Jacobi 更新：同一轮的 $A$ 和 $B$ 都读取上一轮密度。

`EXACTEMB` 在每轮 RASSCF 之前调用一次，环境密度在该次 RASSCF 的全部宏迭代和微迭代中保持固定，因此外层迭代具有清晰的固定点含义。

## 2. OpenMolcas 中的数据表示和程序实现

### 2.1 计算步骤与数据传递

从源片段 $B$ 向目标片段 $A$ 构造势时，程序执行以下步骤：

1. 从当前二聚体 RunFile、ONEINT 和 Cholesky 文件读取二聚体 AO 维数及积分表示。
2. 从源片段 RunFile 读取最终 RASSCF `D1ao`、基函数数目、冻结轨道数、非活性轨道数、活性电子数和弛豫根。
3. 按用户给出的源片段首 AO 编号，把源密度嵌入二聚体 AO 下三角矩阵，其余元素置零。
4. 用二聚体重叠矩阵计算 $\operatorname{Tr}(PS)$，并与源 RunFile 的电子数比较。
5. 调用纯库仑 Cholesky 收缩构造完整二聚体 $J[P_B]$。
6. 提取目标片段 AO 对角块，并写入目标片段 RunFile 的 `Exact Emb Pot`。
7. 目标片段 RASSCF 的 `SGFCIN` 读取该记录，把它加入标准 `OneHam` 路径。

主要源码位置如下：

| 文件 | 功能 |
|---|---|
| `src/exactemb/rdinp_exactemb.F90` | 读取输入并验证短逻辑名、AO 起点和容差 |
| `src/exactemb/exactemb_data.F90` | 保存本次调用的输入参数 |
| `src/exactemb/exactemb.F90` | 读取密度、检查电子数、构造并写入库仑势 |
| `src/fock_util/cho_exactemb_coulomb.F90` | 调用 OpenMolcas Cholesky Fock 接口 |
| `src/rasscf/sgfcin.F90` | 在 RASSCF 中读取并加入嵌入算符 |
| `src/runfile_util/runfile_data.F90` | 注册 `Exact Emb *` RunFile 记录 |
| `src/Driver/exactemb.prgm.src` | 向 pymolcas 注册模块、积分文件和两个逻辑 RunFile |

### 2.2 密度与算符的 packed 约定

OpenMolcas 的对称 AO 矩阵通常以下三角 packed 数组保存，但密度和算符采用不同约定。

`D1ao` 是 folded packed 密度：对角元保持不变，非对角元保存为普通矩阵元的两倍。

`Exact Emb Pot` 和 ONEINT 重叠矩阵是 plain packed 算符：每个下三角元素就是普通矩阵元。

因此可以直接用 packed 点积计算期望值，

$$
D_{\mathrm{fold}}\cdot V_{\mathrm{LT}}=\operatorname{Tr}(PV).
$$

乙烯验证中，把 `D1ao` 按 folded 规则展开得到 $\operatorname{Tr}(PS)=16.0000000024$，按 plain 规则展开则得到 15.5875；`Exact Emb Pot` 按 plain 规则展开后与 PySCF direct-J 相符，而误按 folded 规则解释会产生约 44% 的偏差。

### 2.3 AO 块映射

当前版本要求 `Group=NoSymm`，即 $C_1$ 计算。

源片段和目标片段必须分别对应二聚体 AO 序中的连续块。

如果二聚体输入按“全部 $A$ 原子，然后全部 $B$ 原子”的顺序排列，并且片段与二聚体使用完全相同的原子坐标、原子顺序、基组定义和球谐约定，则通常有

$$
i_A=1,\qquad i_B=n_{\mathrm{AO}}^A+1.
$$

乙烯 cc-pVDZ 中每个片段有 48 个 AO，因此 $A$ 的首 AO 为 1，$B$ 的首 AO 为 49，二聚体共有 96 个 AO。

程序检查源块和目标块是否越界或重叠，并检查映射后的源密度电子数

$$
\left|\operatorname{Tr}(PS)-\left[2(n_{\mathrm{Fro}}+n_{\mathrm{Ish}})+n_{\mathrm{ActEl}}\right]\right|<\varepsilon.
$$

该不变量能有效发现错误的源 AO 起点、明显不一致的基组或源密度，但用户仍应从输入构造上保证目标块的几何和 AO 约定一致。

### 2.4 RunFile 记录

本实现注册了 11 个正式 RunFile 记录：

| 记录 | 类型 | 含义 |
|---|---|---|
| `Exact Emb Active` | iScalar | 非零时 RASSCF 加入嵌入势 |
| `Exact Emb Ver` | iScalar | 数据格式版本，当前为 1 |
| `Exact Emb Root` | iScalar | 源片段的 `Relax CASSCF root` |
| `Exact Emb SrcAO` | iScalar | 源片段在二聚体 AO 序中的首 AO |
| `Exact Emb TgtAO` | iScalar | 目标片段在二聚体 AO 序中的首 AO |
| `Exact Emb DimAO` | iScalar | 构造时二聚体 AO 总数 |
| `Exact Emb Nelec` | dScalar | 映射密度的 $\operatorname{Tr}(PS)$ |
| `Exact Emb JFrob` | dScalar | 目标 $J$ 矩阵的全矩阵 Frobenius 范数 |
| `Exact Emb JMax` | dScalar | 目标 $J$ 的最大绝对矩阵元 |
| `Exact Emb Pot` | dArray | plain packed 目标片段库仑势 |
| `Exact Emb Energy` | dScalar | 最近一次 SGFCIN 调用的 $D\cdot J$ |

### 2.5 RASSCF 中的算符加入位置

`SGFCIN` 首先从 ONEINT 读取标准 `OneHam`，随后检查 `Exact Emb Active`。

激活时，程序验证 `Exact Emb Pot` 的 packed 长度和 $C_1$ 条件，把该数组逐元素加入 `OneHam`，然后继续执行标准 RASSCF Fock、轨道优化和 CI 路径。

因此嵌入势影响波函数和密度，而不仅是输出能量。

程序禁止同时启用 `EXACTEMB` 和 `OFEMbedding`，以避免两个不同嵌入模型在同一 RASSCF 中未经定义地相加。

输出中的激活横幅只在当前 RASSCF 进程第一次调用 `SGFCIN` 时打印一次。

乙烯计算中横幅里的 `current <J[P_partner]>` 对应第一次调用时的非活性密度，不是最终收敛密度的期望值；`Exact Emb Energy` RunFile 记录则会在后续 SGFCIN 调用中继续更新。

### 2.6 CLEAR 的语义

`CLEAR` 将目标 RunFile 中的 `Exact Emb Active` 置为 0，把版本写为 1，并把 `Exact Emb Energy` 置零。

已有的 `Exact Emb Pot` 可以保留在 RunFile 中，但只要激活标志为 0，RASSCF 就不会读取或加入它。

这种设计允许安全停用嵌入势，而不要求修改 RunFile 的内部目录结构或删除数组记录。

### 2.7 输入一致性检查与停止条件

`EXACTEMB` 在以下情况下停止计算并给出诊断：

- SOURCE、TARGET 或 CLEAR 的逻辑名超过 8 个字符，或包含 `/`、反斜杠等路径分隔符。
- 二聚体、源片段或目标片段不是 $C_1$。
- 当前二聚体积分文件不是由 Cholesky/RI-CD 计算产生的。
- 源 RunFile 不含尺寸正确的最终 `D1ao` 或电子数元数据。
- 源或目标 AO 块越界，或两个 AO 块重叠。
- $\operatorname{Tr}(PS)$ 与 RunFile 电子数的差超过 `TOLERANCE`。
- RASSCF 发现 `Exact Emb Pot` 尺寸不符。
- RASSCF 同时启用 `OFEMbedding`。

## 3. 用户输入与文件管理

### 3.1 SOURCE 和 TARGET 是逻辑名，不是路径

OpenMolcas 的 `NameRun` 接口使用最多 8 个字符的逻辑 RunFile 名。

`src/Driver/exactemb.prgm.src` 已注册两个固定逻辑文件：

```text
EMBSRC -> $WorkDir/EMBSRC
EMBTGT -> $WorkDir/EMBTGT
```

因此 `&EXACTEMB` 输入中应使用 `EMBSRC` 和 `EMBTGT`，不能写 `A.RunFile`、`../A.RunFile` 或绝对路径。

在调用模块前，用 EMIL `>>COPY` 把实际源 RunFile 和目标 RunFile 复制到这两个工作目录短名；调用后，再把 `EMBTGT` 复制到需要保存的文件。

### 3.2 关键字

```text
&EXACTEMB
SOURCE
EMBSRC
<源片段在二聚体 AO 序中的首 AO>
TARGET
EMBTGT
<目标片段在二聚体 AO 序中的首 AO>
TOLERANCE
1.0d-6
PRINT
1
END
```

| 关键字 | 含义 |
|---|---|
| `SOURCE` | 下一行是源逻辑 RunFile 名，再下一行是源片段首 AO |
| `TARGET` | 下一行是目标逻辑 RunFile 名，再下一行是目标片段首 AO |
| `TOLERANCE` | 电子数不变量的绝对容差，默认 $10^{-6}$ |
| `PRINT` | 非负值打印构造摘要，2 或更高还打印完整目标矩阵 |
| `CLEAR` | 下一行给出待停用的目标逻辑 RunFile 名 |
| `END` | 输入结束标志，必须出现 |

### 3.3 最小的 B 到 A 构造片段

乙烯二聚体按 A 的 6 个原子在前、B 的 6 个原子在后排列，每个片段 48 AO，因此 B 到 A 的调用为

```text
>>COPY B_source.RunFile EMBSRC
>>COPY A_target.RunFile EMBTGT

&EXACTEMB
SOURCE
EMBSRC
49
TARGET
EMBTGT
1
TOLERANCE
1.0d-6
PRINT
1
END

>>COPY EMBTGT A_with_B_electronic_potential.RunFile
```

这个片段必须位于完整二聚体 GATEWAY/SEWARD Cholesky 计算之后，因为当前 RunFile、ONEINT、ChVec、ChRed、ChRst 和 ChMap 必须属于同一个二聚体几何与基组。

### 3.4 完整可运行的乙烯 B 到 A 输入

下面的输入从零开始生成源片段密度、目标片段核场参考、二聚体 Cholesky 积分、电子库仑势和最终嵌入 RASSCF，不依赖其他作业留下的 RunFile。

四个连续计算部分依次承担以下作用：第一部分得到源 B 的最终密度，第二部分建立含 B 核场的目标 A RunFile，第三部分在 AB Cholesky 积分上构造 $J[P_B]$ 并写入目标 RunFile，第四部分重建 A 的积分与 XField 后读取该势并完成最终 RASSCF。

输入文件名和 `$Project` 可以自行选择；最后两个 `$CurrDir` 副本用于在 OpenMolcas scratch 清理后保存结果。

```text
* Complete neutral-ethylene exact-density electrostatic embedding: B -> A

&GATEWAY
Title=Ethylene B source
Coord
6
ethylene B at z 3.5
C  -0.6695   0.0000   3.5000
C   0.6695   0.0000   3.5000
H  -1.2321   0.9230   3.5000
H  -1.2321  -0.9230   3.5000
H   1.2321   0.9230   3.5000
H   1.2321  -0.9230   3.5000
Basis=cc-pVDZ
Group=NoSymm

&SEWARD

&SCF
Title=Ethylene B RHF

&RASSCF
Title=Ethylene B CASSCF(2e,2o)
FileOrb=SCFORB
Charge=0
Spin=1
NActEl=2
Ras2=2
PRWF=1.0d-6
PRSD

>>COPY $Project.RunFile B0.RunFile

&GATEWAY
Title=Ethylene A with nuclei of B represented by XField
Coord
6
ethylene A at z 0
C  -0.6695   0.0000   0.0000
C   0.6695   0.0000   0.0000
H  -1.2321   0.9230   0.0000
H  -1.2321  -0.9230   0.0000
H   1.2321   0.9230   0.0000
H   1.2321  -0.9230   0.0000
Basis=cc-pVDZ
Group=NoSymm
XField
6 ANGSTROM 0 0 0 0
-0.6695  0.0000  3.5  6.0
 0.6695  0.0000  3.5  6.0
-1.2321  0.9230  3.5  1.0
-1.2321 -0.9230  3.5  1.0
 1.2321  0.9230  3.5  1.0
 1.2321 -0.9230  3.5  1.0

&SEWARD

&SCF
Title=Ethylene A RHF in the nuclear field of B

&RASSCF
Title=Ethylene A nuclei-only CASSCF(2e,2o)
FileOrb=SCFORB
Charge=0
Spin=1
NActEl=2
Ras2=2
PRWF=1.0d-6
PRSD

>>COPY $Project.RunFile A0.RunFile

&GATEWAY
Title=Ethylene AB Cholesky integrals
Coord
12
parallel ethylene dimer A then B at 3.5 Angstrom
C  -0.6695   0.0000   0.0000
C   0.6695   0.0000   0.0000
H  -1.2321   0.9230   0.0000
H  -1.2321  -0.9230   0.0000
H   1.2321   0.9230   0.0000
H   1.2321  -0.9230   0.0000
C  -0.6695   0.0000   3.5000
C   0.6695   0.0000   3.5000
H  -1.2321   0.9230   3.5000
H  -1.2321  -0.9230   3.5000
H   1.2321   0.9230   3.5000
H   1.2321  -0.9230   3.5000
Basis=cc-pVDZ
Group=NoSymm
RICD

&SEWARD
Medium Cholesky

>>COPY B0.RunFile EMBSRC
>>COPY A0.RunFile EMBTGT

&EXACTEMB
SOURCE
EMBSRC
49
TARGET
EMBTGT
1
TOLERANCE
1.0d-6
PRINT
1
END

>>COPY EMBTGT A1.RunFile

&GATEWAY
Title=Final embedded ethylene A
Coord
6
ethylene A at z 0
C  -0.6695   0.0000   0.0000
C   0.6695   0.0000   0.0000
H  -1.2321   0.9230   0.0000
H  -1.2321  -0.9230   0.0000
H   1.2321   0.9230   0.0000
H   1.2321  -0.9230   0.0000
Basis=cc-pVDZ
Group=NoSymm
XField
6 ANGSTROM 0 0 0 0
-0.6695  0.0000  3.5  6.0
 0.6695  0.0000  3.5  6.0
-1.2321  0.9230  3.5  1.0
-1.2321 -0.9230  3.5  1.0
 1.2321  0.9230  3.5  1.0
 1.2321 -0.9230  3.5  1.0

&SEWARD

&SCF
Title=Ethylene A RHF in the nuclear field of B

>>COPY A1.RunFile $Project.RunFile

&RASSCF
Title=Final exact-density embedded ethylene A CASSCF(2e,2o)
FileOrb=SCFORB
Charge=0
Spin=1
NActEl=2
Ras2=2
PRWF=1.0d-6
PRSD

>>COPY $Project.RunFile $CurrDir/A1_final.RunFile
>>COPY $Project.RasOrb.1 $CurrDir/A1_final.RasOrb
```

### 3.5 输出摘要的解释

一次成功的构造会打印源和目标逻辑名、源弛豫根、二聚体 AO 数、源和目标 AO 范围、期望电子数、$\operatorname{Tr}(PS)$、电子数偏差、$J$ 范数和最大矩阵元。

乙烯 B 到 A 的代表性输出为

```text
  Source RunFile        : EMBSRC
  Target RunFile        : EMBTGT
  Source relaxation root:                     1
  Dimer AO count        :                    96
  Source AO block       :                    49                   96
  Target AO block       :                     1                   48
 Expected electrons    :        16.0000000000
 Tr[P S]               :        16.0000000000
 |Delta N|             :   1.705303E-13
 packed ||J||_2        :   2.087182E+01
 full-matrix ||J||_F   :   2.505784E+01
 max |J_mn|            :   2.448000E+00
  Exchange contribution : DISABLED
  Partner nuclei        : supply with standard SEWARD XField
  Target RunFile record : Exact Emb Pot
```

`|Delta N|` 接近机器精度说明源密度与所选二聚体源 AO 块在重叠度量下相容。

`packed ||J||_2` 是 packed 数组自身的欧氏范数，`full-matrix ||J||_F` 是把 plain packed 算符展开成完整对称矩阵后的 Frobenius 范数；两者定义不同，不应混用。

### 3.6 CLEAR 输入

```text
>>COPY A1_final.RunFile EMBTGT

&EXACTEMB
CLEAR
EMBTGT
END

>>COPY EMBTGT $CurrDir/A1_cleared.RunFile
```

清除后 `Exact Emb Active=0`，后续 RASSCF 忽略已有的 `Exact Emb Pot`。

### 3.7 在 GV800 上运行

验证使用的关键环境为

```bash
export MOLCAS=/scratch/alpha/OpenMolcasUP_ExactDenEmb-65a3eeb11
export PATH=$MOLCAS:$MOLCAS/bin:/opt/mamba/envs/omc/bin:$PATH
export LD_LIBRARY_PATH=$MOLCAS/lib:/scratch/alpha/build-exactden-20260820_093910/wignernj-sys/lib:/opt/mamba/envs/omc/lib:$LD_LIBRARY_PATH
export MOLCAS_NPROCS=1
export MOLCAS_MEM=12000
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
```

验证计算使用 Slurm `molcas` 分区、1 个 MPI 进程和 12000 MiB 内存。

RunFile 的并行环境应与生成它时保持兼容；本文全部验证统一使用 1 个进程。

EMIL 裸文件名通常位于 OpenMolcas 工作目录并随 scratch 清理，必须保存的最终 RunFile、RasOrb、VecDet 和分析数据应复制到 `$CurrDir/`。

安装二进制依赖离线安装的 `wignernj-sys/lib`，该目录必须保留或通过 `LD_LIBRARY_PATH` 提供等价动态库。

## 4. 推广到其他分子体系

### 4.1 必需准备的数据

对新的二片段体系，至少需要准备以下内容：

1. 片段 $A$ 的完整几何、基组、电子态、活性空间和收敛的 RASSCF 输入。
2. 片段 $B$ 的对应输入。
3. 按连续片段顺序排列的完整二聚体几何。
4. 目标片段 SEWARD 输入中的环境核 `XField`。
5. 二聚体 `RICD` 和 SEWARD Cholesky 积分文件。

片段 RunFile 必须来自最终收敛的 RASSCF，因为 `EXACTEMB` 读取其中状态特定的最终 `D1ao`。

### 4.2 确定 AO 起点

先分别运行片段 A、片段 B 和二聚体的 GATEWAY/SEWARD，记录 $n_A$、$n_B$ 和 $n_{AB}$。

在二聚体原子顺序为 A 后接 B 且 AO 连续排列时，应满足

$$
n_{AB}=n_A+n_B,\qquad i_A=1,\qquad i_B=n_A+1.
$$

B 到 A 使用 `SOURCE EMBSRC, i_B` 和 `TARGET EMBTGT, i_A`；A 到 B 则交换源、目标及其起点。

首次新体系计算必须检查输出中的源/目标 AO 范围和 $\operatorname{Tr}(PS)$，不能只根据原子数猜测 AO 起点。

### 4.3 几何和基组一致性

片段在独立输入和二聚体输入中的绝对坐标、原子顺序、同位素或标签、基组、球谐/Cartesian 约定以及整体取向必须一致。

刚性平移的片段应在独立片段输入中使用其在二聚体中的实际平移后坐标；不要用位于原点的片段 RunFile 去代表二聚体中已经平移的源块。

混合基组可以使用，但必须保证片段和二聚体对应原子上的基组逐一相同。

当前实现只支持 $C_1$，因此三个 GATEWAY 输入都应使用 `Group=NoSymm`。

使用 ECP、不同核模型或非标准外场时，`XField` 中应采用何种有效核电荷需要另行核实；本文只验证了全电子 C/H 体系。

### 4.4 电子态和活性空间

源 RunFile 的电子数由

$$
N_e=2(n_{\mathrm{Fro}}+n_{\mathrm{Ish}})+n_{\mathrm{ActEl}}
$$

确定。

对于多根 RASSCF，应确认 `Relax CASSCF root` 与目标物理态一致，并在 freeze-and-thaw 各轮监测 root identity、自然占据数和活性轨道连续性。

如果环境极化导致活性轨道发生分支切换，只比较总能量可能不足以判断是否收敛到同一电子态。

### 4.5 数值参数与收敛检查

建议首次计算使用 `TOLERANCE=1.0d-6`，并要求实际电子数误差远小于该阈值。

Cholesky 阈值应按目标精度选择；乙烯中 `Medium Cholesky` 已足以使 $J$ 的相对误差达到约 $2\times10^{-8}$，但更大或更弥散的体系需要重新检查。

一次嵌入后至少检查 RASSCF 能量、自然占据数、目标态、轨道连续性以及嵌入前后的密度变化。

相互极化计算应同时监测相邻轮能量差和正交化 AO 度量下的密度差，例如

$$
\left\|S^{1/2}\Delta P S^{1/2}\right\|_F.
$$

该 Löwdin 度量把非正交 AO 基中的密度变化转换为正交空间中的 Frobenius 范数，因而比直接比较 packed AO 元素更适合跨轮收敛判断。

若无阻尼迭代振荡或发散，可在外层密度更新中加入线性混合，但混合方式和收敛阈值必须明确记录。

### 4.6 推荐的新体系复现顺序

1. 完成不含嵌入记录的普通 RASSCF，并保存能量、轨道和 RunFile 作为参考。
2. 分别生成源片段、目标片段和二聚体 Cholesky 积分文件。
3. 构造一个方向的 `Exact Emb Pot`，核对 AO 范围、电子数、$J$ 范数和 RunFile 记录。
4. 在小体系或代表性截断模型上，用直接四中心积分程序独立复核 $J$。
5. 先做仅电子库仑势计算，确认算符进入 RASSCF 并产生合理响应。
6. 加入环境核 `XField`，完成一次完整静电嵌入。
7. 构造反方向势，并进行相互极化迭代。
8. 如片段态用于 NOCI，另建不含任何嵌入势的二聚体积分与哈密顿量计算。

## 5. 乙烯二聚体数值验证

### 5.1 计算设置

片段 A 位于 $z=0$，片段 B 是沿 $z$ 方向平移 3.5 Å 的相同乙烯。

每个片段采用 cc-pVDZ、`Group=NoSymm`、RHF 初始轨道和 CAS(2e,2o) RASSCF，设置 `NActEl=2`、`Ras2=2`、`Spin=1` 和 `PRWF=1.0d-6`。

每个片段有 48 个基函数、14 个闭壳层电子和 2 个活性电子；二聚体有 96 个基函数。

二聚体通过 GATEWAY `RICD` 和 SEWARD `Medium Cholesky` 生成 Cholesky 积分。

### 5.2 验证 1–2：无嵌入回归与 CLEAR

普通中性乙烯 RASSCF 在没有 `Exact Emb Active` 记录时正常完成，root-1 总能量为

$$
E_{\mathrm{isolated}}=-78.06789017\ E_h,
$$

与既有中性乙烯参考在打印精度上完全一致，输出中没有嵌入激活横幅。

`CLEAR` 最终完整输入在 Slurm 作业 859 中正常进入 `exactemb.exe`，打印

```text
EXACTEMB: embedding flag cleared in target RunFile: EMBTGT
```

作业以 `RUN_RC=0` 和 `Happy landing` 结束，保存的清除后 RunFile 为 263192 B，SHA-256 为 `6b9803b308db8cdcf5c59589c3ab0996ff17c362134ce21fa37a574f15a21937`。

验证 1–2 通过。

### 5.3 片段与二聚体计算文件

Slurm 作业 860–862 分别生成 A、B 和 AB 计算文件，全部 `RUN_RC=0`。

A 和 B 都得到 48 个基函数、16 个电子和 $-78.06789017\ E_h$ 的 RASSCF root-1 能量，验证了刚性平移下的能量不变性。

AB 得到 96 个基函数以及非空的 RunFile、OneInt、ChVec1、ChRed、ChRst、ChMap 和 RICDLib 文件。

主要文件大小为：RunFile 358784 B、OneInt 2501024 B、ChVec1 24488736 B、ChRed 4157136 B、ChRst 31256 B、ChMap 3248 B、RICDLib 4692 B。

### 5.4 验证 3：双向库仑势和电子数不变量

Slurm 作业 863 构造 B 到 A 的势，源 AO 范围为 49–96，目标范围为 1–48。

测得

$$
N_{\mathrm{expected}}=16.0000000000,\qquad \operatorname{Tr}(PS)=16.0000000000,
$$

$$
|\Delta N|=1.705303\times10^{-13}.
$$

目标势的 packed 数组范数为 20.87182，全矩阵 Frobenius 范数为 25.05784，最大绝对矩阵元为 2.448000。

Slurm 作业 864 构造 A 到 B 的势，源 AO 范围为 1–48，目标范围为 49–96，电子数误差为 $1.776357\times10^{-13}$。

由于两个片段只是刚性平移，双向构造得到相同的三个 $J$ 统计量。

Slurm 作业 865 故意把源首 AO 设为 73，使源块变成 73–120 并超出二聚体 1–96 的范围；程序按设计打印明确诊断并停止计算，证明 AO 范围检查有效。

两个成功构造后的目标 RunFile 都为 273088 B，比未写势的片段 RunFile 增加 9896 B，其中 9408 B 正好对应 $48\times49/2$ 个双精度 packed 元素，其余为记录开销。

验证 3 通过。

### 5.5 验证 4：与 PySCF 直接四中心积分比较

为避免不同程序原生 AO 表示造成歧义，验证程序从 SEWARD 生成的 `guessorb.h5` 读取 `BASIS_FUNCTION_IDS`、`PRIMITIVE_IDS` 和 `PRIMITIVES`，在 PySCF 中重建与 OpenMolcas 完全相同的 AO 基组。

源 `D1ao` 和目标 `Exact Emb Pot` 由链接 `libmolcas` 的独立 Fortran 工具导出。

PySCF 使用直接四中心积分收缩

```python
get_jk(
    (molA, molA, molB, molB),
    dmB,
    scripts="ijkl,lk->ij",
    aosym="s4",
)
```

并与显式 `shls_slice` ERI contraction 互相核对，两种 PySCF 路径相差约 $6\times10^{-15}$。

除逐元素矩阵差外，验证还求解广义本征问题 $Jc=\lambda Sc$；其本征值在一致的非奇异 AO 表示变换下保持不变，可作为不依赖具体 AO 归一化和排列的补充检查。

最终比较为

```text
[check] Tr[P_B S_B] = 16.0000000024
[result] ||J_omc||_F = 25.05783926
[result] ||J_ref||_F = 25.05783942
[result] max|dJ| = 8.606e-08
[result] ||dJ||_F/||J_ref|| = 1.842e-08
[result] generalized eig max dev = 1.156e-06
[result] largest generalized eigenvalue = 3.114725 / 3.114725
```

该结果说明 OpenMolcas Cholesky 库仑势与 PySCF 直接四中心积分参考在约 $2\times10^{-8}$ 的矩阵相对误差内一致。

差值反映 Cholesky/RI-CD 截断、AO 重建和浮点数值误差的合成尺度；没有观察到可分辨的 AO 块映射错误。

验证 4 通过。

### 5.6 验证 5–6：RASSCF 对电子库仑势和完整静电势的响应

仅加入 $J[P_B]$ 而不加入环境核时，OpenMolcas 得到

$$
E_A[J_B]=-42.13929260\ E_h.
$$

同时加入环境核 `XField` 和 $J[P_B]$ 时，得到

$$
E_A[V_{\mathrm{nuc}}^B+J_B]=-41.66496895\ E_h.
$$

`XField` 单独作用时的参考能量为 $-78.82129253\ E_h$。

这些片段能量包含所选外部势及 OpenMolcas 的 XField 核－点电荷常数，不应直接解释为孤立片段能或完整二聚体相互作用能。

PySCF 参考使用相同 AO 基、相同源密度、直接积分 $J$、相同环境核势和相同 CASSCF 初始轨道。

OpenMolcas 的 XField 总能量还包含 A 原子核与 B 正点电荷之间的常数

$$
E_{\mathrm{nuc}}^{A-B\mathrm{\ charge}}=36.28788679\ E_h,
$$

比较时在 PySCF 一侧加入同一常数。

结果为

| 计算 | PySCF / $E_h$ | OpenMolcas / $E_h$ | 差值 / $E_h$ |
|---|---:|---:|---:|
| 仅环境核 XField | -78.82129253 | -78.82129253 | $-4.788\times10^{-10}$ |
| 仅环境电子 $J[P_B]$ | -42.13929286 | -42.13929260 | $-2.635\times10^{-7}$ |
| 完整静电嵌入 | -41.66496923 | -41.66496895 | $-2.803\times10^{-7}$ |

PySCF CASSCF 必须从与 OpenMolcas 相同的 XField-SCF 轨道开始；从嵌入 hcore 的另一组 SCF 轨道开始会收敛到不同的局部 CASSCF 解。

RASSCF 输出横幅第一次打印的 $D\cdot J$ 为 31.699010 和 33.382510 $E_h$，分别对应仅电子势和完整静电势计算的首次 SGFCIN 密度。

用最终密度直接计算得到约 35.425 和 36.470 $E_h$，这一区别来自打印时机，不表示嵌入算符缺失。

验证 5–6 通过。

### 5.7 验证 7：六轮相互极化

以孤立 A、B 密度作为第 0 轮，采用 Jacobi 更新完成六轮双向计算，Slurm 作业 870–875 全部 `RUN_RC=0`。

嵌入 RASSCF 能量为

| 轮次 $n$ | $E_A^{(n)}$ / $E_h$ | $E_B^{(n)}$ / $E_h$ | $\Delta E_A$ / $E_h$ | $\Delta E_B$ / $E_h$ |
|---:|---:|---:|---:|---:|
| 1 | -41.66496895 | -41.66496894 | — | — |
| 2 | -41.65985838 | -41.65985229 | $+5.111\times10^{-3}$ | $+5.117\times10^{-3}$ |
| 3 | -41.66066547 | -41.66065879 | $-8.071\times10^{-4}$ | $-8.065\times10^{-4}$ |
| 4 | -41.66059870 | -41.66059193 | $+6.68\times10^{-5}$ | $+6.69\times10^{-5}$ |
| 5 | -41.66060393 | -41.66059716 | $-5.23\times10^{-6}$ | $-5.23\times10^{-6}$ |
| 6 | -41.66060353 | -41.66059675 | $+4.0\times10^{-7}$ | $+4.1\times10^{-7}$ |

第 1 轮 A 能量与一次完整静电嵌入结果逐位相同，验证了两种输入路径的一致性。

正交化 AO 度量下的相邻轮密度差为

| 变化 | A | B |
|---|---:|---:|
| 孤立态 $\rightarrow$ 第 1 轮 | $1.38\times10^{-1}$ | $1.365\times10^{-1}$ |
| 第 1 轮 $\rightarrow$ 第 2 轮 | $1.239\times10^{-3}$ | $1.672\times10^{-3}$ |
| 第 2 轮 $\rightarrow$ 第 3 轮 | $9.120\times10^{-5}$ | $1.188\times10^{-4}$ |
| 第 3 轮 $\rightarrow$ 第 4 轮 | $7.035\times10^{-6}$ | $9.199\times10^{-6}$ |
| 第 4 轮 $\rightarrow$ 第 5 轮 | $5.475\times10^{-7}$ | $7.137\times10^{-7}$ |
| 第 5 轮 $\rightarrow$ 第 6 轮 | $4.235\times10^{-8}$ | $5.532\times10^{-8}$ |
| 第 1 轮 $\rightarrow$ 第 6 轮 | $1.158\times10^{-3}$ | $1.566\times10^{-3}$ |

能量与密度变化都呈稳定几何衰减，第 6 轮能量变化已低于 $1\ \mu E_h$，说明该中性乙烯测试在无阻尼条件下稳定收敛。

A、B 最终能量存在约 $6.8\times10^{-6}\ E_h$ 的小差异；两个方向由独立计算路径得到，现有计算没有通过改变 Cholesky 阈值或收敛参数进一步分解差异来源，因此本文只报告该观测值，不把它归因于特定物理或数值机制。

验证 7 通过。

## 6. 实现修正与可追溯性

### 6.1 编译与运行接口修正

从基线提交到最终验证版本共完成五项逻辑修正，其中两项是源码编译问题，三项是运行接口或程序集成问题：

| 类别 | 问题 | 修正 |
|---|---|---|
| 编译 | `Reset_ExactEmb_Input` 的 `PUBLIC` 声明位于 `contains` 之后 | 把声明移入模块 specification 部分 |
| 编译 | `rdinp_exactemb.F90` 使用 `wp` 但未导入 | 从 `Definitions` 导入 `wp` |
| 集成 | pymolcas 没有 `exactemb.prgm` | 新增 pymolcas 模块描述文件并注册积分文件 |
| 接口 | 把长物理路径传给 8 字符 `NameRun` | 改为 `EMBSRC`、`EMBTGT` 逻辑名并在解析阶段拒绝路径 |
| RunFile | 11 个 `Exact Emb *` 字段未注册 | 加入整数标量、实数标量和实数数组标签表 |

离线安装 `libwignernj` 是 GV800 无外网条件下的构建环境处理，不属于第三个源码编译错误。

### 6.2 早期开发作业的实际记录

早期 Slurm 作业 852–859 共 8 个，记录了从运行脚本到最终 CLEAR 成功的完整调试过程：

| Slurm 编号 | 结果 | 原因或结论 |
|---:|---|---|
| 852 | FAILED | 作业脚本未预建 OpenMolcas scratch 父目录 |
| 853 | COMPLETED | 无嵌入中性乙烯基线正常 |
| 854 | FAILED | pymolcas 缺少 `exactemb.prgm` |
| 855 | FAILED | 模块描述文件已生效，但测试输入缺少兼容的当前 RunFile |
| 856 | FAILED | 未注册的长相对 RunFile 名不能由 `NameRun` 解析 |
| 857 | FAILED | 绝对路径被旧接口截断为 `/scratch` |
| 858 | FAILED | `Exact Emb Active` 等新 RunFile 标签尚未注册 |
| 859 | COMPLETED | 短逻辑名与 11 个标签修正后，完整 CLEAR 输入通过 |

这些失败属于实现完善过程；最终代码和后续数值验证使用的是修正后的同一安装。

### 6.3 构建环境

GV800 最终安装位于

```text
/scratch/alpha/OpenMolcasUP_ExactDenEmb-65a3eeb11
```

构建使用 GCC/GFortran 13.3.0、Open MPI 4.1.6、MKL、HDF5、LibXC、GA 5.8.2 和离线 `libwignernj`，以 `make -j4` 完成。

`exactemb.exe` 的验证时 SHA-256 为

```text
b4c8244f9d5b7d49588453c7af6644dfc642250cd4ad35ba0234b56e0b6747d3
```

详细构建命令、四轮编译记录和后续增量重编译见 `TODO/20260820_GV800_exactemb_build_record.md`。

## 7. 数据、输入和输出位置

本文给出方法、完整代表性输入和全部关键数值，但不把数十兆字节的原始 OpenMolcas 输出逐字嵌入正文。

本地目录

```text
TODO/20260820_exactemb_ethylene_gv800/
```

保存了实际输入文件、运行脚本、RunFile 导出工具、PySCF 分析程序和逐步骤记录，包括：

- `S2_A_neutral.in`、`S2_B_neutral.in` 和 `S2_AB_ricd.in`。
- `S3_BtoA.in`、`S3_AtoB.in` 和故意错误的 `S3_guard_srcrange.in`。
- `S5_A_*.in` 所包含的仅电子库仑势输入和完整静电嵌入输入。
- 六轮相互极化模板 `S6_ft_template.in`。
- `dump_rf.F90`、`stage04_compare_final.py` 和 `stage05_verify.py`。
- `stage00_setup.md` 至 `stage06_freeze_thaw.md` 的原始分步记录。

GV800 原始输出、RunFile、轨道、Cholesky 文件、导出矩阵、NPY/NPZ 结果和 `.provenance` 环境记录保存在

```text
/scratch/alpha/exactemb-ethylene-20260820_103419/work/
```

构建日志保存在

```text
/scratch/alpha/build-exactden-20260820_093910/
```

## 8. 结论

`EXACTEMB` 已实现一种完全在 OpenMolcas 内执行的、基于完整片段一阶密度矩阵的静电嵌入方法。

环境核势由标准 `XField` 提供，环境电子势由二聚体 Cholesky 积分与源 RASSCF `D1ao` 收缩得到，目标势通过正式 RunFile 记录进入 RASSCF 的一电子哈密顿量。

乙烯二聚体的七项验证全部通过，其中最关键的矩阵级验证表明 OpenMolcas $J[P]$ 与 PySCF 直接四中心积分参考的相对差为 $1.842\times10^{-8}$，完整静电嵌入 CASSCF 能量相差 $2.803\times10^{-7}\ E_h$，六轮相互极化在能量和密度上均稳定几何收敛。

当前版本适用于 $C_1$、连续 AO 片段块、Cholesky/RI-CD 二聚体积分文件和无片段间交换/Pauli 投影的静电极化计算。

按照本文给出的完整输入、AO 映射规则、XField 构造、文件逻辑名和检查量，可以把同一方法移植到其他二片段体系；对新元素、ECP、弥散基组、更大活性空间和更强极化体系，应重新验证 AO 对应关系、电子态连续性、Cholesky 精度和相互极化收敛。
