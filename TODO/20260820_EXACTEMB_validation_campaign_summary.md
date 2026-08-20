# EXACTEMB 全局总结报告:GV800 编译与乙烯验证战役(2026-08-20)

分支:`agent/exact-density-embedding`(65a3eeb11 + 当日修复);执行:shuo/ZCode/Codex 三方接力;主机:GV800(alpha 账号);最终结论:**编译通过、验证 Tests 1–7 全部通过,纯 OpenMolcas 精确密度 M2 嵌入实现获得端到端独立验证**。

## 1. 一句话结论

`&EXACTEMB` 模块从源码编译、RunFile 接口、Coulomb 矩阵数值(对 PySCF 直接 ERI 1.8e-8)、J-only 与完整 one-shot M2 嵌入能量(2.6–2.8e-7)到互相 freeze-thaw 不动点(几何收敛,六轮达 4e-7 Eh / 密度 5e-8)全部通过验证;发现的全部问题(3 个编译错、2 个运行时接口错、若干输入/分析工具错)均已修复并留痕。

## 2. 时间线总览

| 时段 | 事件 | 记录 |
|---|---|---|
| 09:39–10:15 | 编译 4 轮:wignernj 离线方案 → exactemb_data.F90 public 位置 → rdinp_exactemb.F90 缺 wp → 成功安装 | `TODO/20260820_GV800_exactemb_build_record.md` |
| 10:32–11:09 | Stage 0/1:环境审计 + 无记录回归 + CLEAR 冒烟;修复 driver 注册(exactemb.prgm.src)与 RunFile 名字/标签接口(EMBSRC/EMBTGT + 11 个标签注册),job 859 通过 | `stage00–stage01_*.md` |
| 11:19–11:22 | Stage 2:A/B 片段 RunFile + AB RICD 上下文生成,全部不变量命中(零错) | `stage02_contexts.md` |
| 11:24–11:27 | Stage 3:B→A 与 A→B 构造(\|ΔN\|≈1.7e-13)+ AO 越界守卫硬中止验证(零错) | `stage03_exactemb_construction.md` |
| 11:30–12:10 | Stage 4:dump_rf 工具 + replica 基组桥 + direct-J 对比,元素级 1.842e-08 | `stage04_direct_j.md` |
| 12:10–12:40 | Stage 5:J-only 与 one-shot M2 三道门 2.6e-7–4.8e-10(一次输入级修复:$CurrDir 保存前缀) | `stage05_oneshot_m2.md` |
| 12:45–13:10 | Stage 6:六轮互相 freeze-thaw,几何收敛到亚 μEh(零错) | `stage06_freeze_thaw.md` |

## 3. 代码修改清单(全部尚未提交,HEAD 仍为 65a3eeb11)

| 文件 | 修改 | 作者 | 性质 |
|---|---|---|---|
| `src/exactemb/exactemb_data.F90` | public 语句移入模块说明部分(编译错修复) | 用户 | 编译错 |
| `src/exactemb/rdinp_exactemb.F90` | `use Definitions` 补 `wp`(编译错修复) | ZCode | 编译错 |
| `src/Driver/exactemb.prgm.src` | 新增,driver 注册(修 "Requested module not available") | Codex | 运行时接口 |
| `src/exactemb/exactemb.F90` | RunFile 短名契约(EMBSRC/EMBTGT,8 字符校验) | Codex | 运行时接口 |
| `src/runfile_util/runfile_data.F90` | 注册 11 个 `Exact Emb *` 标签(修临时字段拒绝) | Codex | 运行时接口 |
| `doc/openmolcas_..._20260820.md` | 实现文档同步(EMBSRC/EMBTGT 语义等) | Codex | 文档 |

验证产物工具(非源码):`dump_rf.F90`(RunFile/ONEINT 记录导出)、`stage04_compare_final.py`、`stage05_verify.py`、`S6_ft_template.in` 等均在 `TODO/20260820_exactemb_ethylene_gv800/`。

## 4. 各验证 Stage 结果汇总

| Stage / Test | 判据 | 关键数字 | 结果 |
|---|---|---|---|
| 1 无记录回归 | 与 main 行为一致 | 基线能量 -78.06789017 Eh 复现 | ✅ |
| 1 CLEAR | 清除 `Exact Emb Active` | job 859 `Happy landing` | ✅ |
| 2 上下文 | 48/48/96 BF、16e、镜像一致 | 能量逐位、J 统计 A↔B 完全相同 | ✅ |
| 3 构造+守卫 | Tr[PS]=Ne、错误输入硬中止 | \|ΔN\|=1.7e-13(容差 1e-6);守卫诊断精确命中 | ✅ |
| 4 direct-J | 对 PySCF 直接 ERI 元素级一致 | max\|ΔJ\|=8.6e-8,relF=1.842e-8,广义本征谱 1.2e-6 | ✅ |
| 5 J-only(Test 5) | 嵌入能量外部复算 | -42.13929260 vs -42.13929286,2.6e-7 | ✅ |
| 5 one-shot M2(Test 6) | 含 XField 记账对齐 | -41.66496895 vs -41.66496923(+Enn 常数),2.8e-7 | ✅ |
| 6 freeze-thaw(Test 7) | 不动点收敛 | 六轮,ΔE 4e-7,密度范数 4–6e-8,比率 ~0.073 | ✅ |

## 5. 建立的语义与约定知识(跨代码对齐必读)

1. **RunFile 存储约定**:`D1ao` 为 fold 打包(展开折半;Tr[PS]=16.0000000024 实证),`Exact Emb Pot` 为 plain 下三角(EXACTEMB 打印的 ‖J‖_F 仅在该解释下与参考一致);ONEINT 重叠为 plain。EXACTEMB 内部打包点积对 (D fold, S plain) 与 (D plain, S fold) 不敏感——本次靠 replica 基组外部判定消歧。
2. **XField 总能量构成** = 电子部分(h+V) + A核×B点电荷 常数(36.28788679 Eh 实证),无自能项;跨代码比较 M2 能量必须显式加这项。
3. **`<J[P_partner]>` 打印时机**:首次 SGFCIN 调用、仅 inactive(14e)密度,非终态期望(终态 35.4/36.5 vs 打印 31.7/33.4);算子正确性由能量级 10⁻⁷ 一致证明。建议在实现文档注明。
4. **OpenMolcas CC-PVDZ ↔ PySCF cc-pvdz 不是置换关系**(广收缩 vs 分割收缩):必须走精确交叉重叠变换或直接在 replica 基组内比较;`guessorb.h5`(SEWARD 自带 run2h5 元数据)+ `check_molcas_pyscf_ao_bridge.py` 是已验证通道。
5. **EMIL `>>COPY` 语义**:裸名落在 scratch workdir(随清理删除),跨任务保存必须 `$CurrDir/` 前缀;RunFile 短名 ≤8 字符且不得含路径分隔符。
6. **GV800 离线依赖**:ExternalProject 的 libwignernj clone 需网络,离线方案 = staging tarball 预装 + `CMAKE_PREFIX_PATH` 命中 `find_package`;二进制 RPATH 依赖 `wignernj-sys` 目录,不可清理。

## 6. 产物位置

- GV800 `/scratch/alpha/OpenMolcasUP_ExactDenEmb-65a3eeb11/`:安装产物(49 exe,含 exactemb.exe;含当日全部接口修复的增量重建)。
- GV800 `/scratch/alpha/build-exactden-20260820_093910/`:源码(已打 XMLLINT patch + 修复)、openmolcas-build-j4 增量构建树(.mod 文件供 dump_rf 用)、wignernj-sys、全部编译日志。
- GV800 `/scratch/alpha/exactemb-ethylene-20260820_103419/`:验证战役根(work/stage01–06 全部作业输入/输出/溯源、jobs 脚本、stage04 工具与 npz)。
- 本地 `TODO/`:本报告、编译记录、`20260820_exactemb_ethylene_gv800/`(README 索引 + stage00–06 + 全部脚本副本)。
- 本地 `~/develop/tmp/gv800_exactden_build_20260820_093910/`:编译日志副本与 BUILD_REPORT.md。

## 7. 未完成事项与建议

1. **尽快提交代码**:§3 清单中 6 个文件改动 + TODO/ 目录;GV800 安装对应修复后代码,提交后源码-产物可追溯关系才闭环。
2. **Test 8(生产目标)**:真实 pentacene Stage F 四路线对比(孤立 / LoProp 多极 / 精确冻结密度 / 互相精确密度极化),并将片段态接入不变裸哈密顿的 NOCI;本次战役的全部方法论(上下文模板、验证脚本、replica 桥)可直接复用。
3. 可选改进:`<J>` 诊断改为终态密度或注明打印时机;EXACTEMB 的 RunFile 短名契约与 XField 记账写入用户文档。
