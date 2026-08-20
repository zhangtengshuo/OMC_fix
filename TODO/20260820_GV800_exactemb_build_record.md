# GV800 编译全记录:agent/exact-density-embedding 分支(65a3eeb11)

日期:2026-08-20;执行:shuo(经 ZCode);GV800 账号:alpha;结果:**经 4 轮迭代编译成功**,安装于 `/scratch/alpha/OpenMolcasUP_ExactDenEmb-65a3eeb11`,新模块 `exactemb.exe` 正常产出。

## 1. 任务目标

把本仓库 `agent/exact-density-embedding` 分支(提交 65a3eeb11,即 "Add EXACTEMB execution and validation README")上传 GV800,按 2026-07-16 DownBas 实际构建配方编译,临时安装到 `/scratch/alpha` 下。要求:打包上传不带 `.git*` 隐藏文件;完整记录编译流程与所有输出;预计新功能代码可能编译失败。

## 2. 编译配方(复刻 GV800 上 2026-07-16 DownBas 实际构建)

| 项 | 值 |
|---|---|
| 工具链 | `/opt/mamba/envs/omc`(gcc/gfortran 13.3.0 + openmpi 4.1.6 + MKL + hdf5 2.1.0 shared + libxc + python 3.10.20) |
| GA | `/opt/ga-5.8.2-mambaomc`(libga.a,GA 5.8.2,blacs_openmpi_ilp64 链接) |
| 环境变量 | `GAROOT=/opt/ga-5.8.2-mambaomc`、`MKLROOT=/opt/mamba/envs/omc`、`LD_LIBRARY_PATH` 加 `$OMC/lib` |
| cmake 选项 | `-DCMAKE_BUILD_TYPE=Release -DLINALG=MKL -DHDF5=ON -DMPI=ON -DGA=ON -DBUILD_SHARED_LIBS=ON -DEXTERNAL_LIBXC=$OMC -DEXTRA_FFLAGS="-ffree-line-length-none -Wno-error" -DCMAKE_PREFIX_PATH="<wignernj-sys>;$GA;$OMC" -DCMAKE_INSTALL_PREFIX=/scratch/alpha/OpenMolcasUP_ExactDenEmb-65a3eeb11` |
| 编译器 | 全部指死 `/opt/mamba/envs/omc/bin/{mpifort,mpicc,mpicxx}`(含 `MPI_*_COMPILER` 三项) |
| Python | `-DPYTHON_EXECUTABLE=$OMC/bin/python3`(pymolcas) |
| XMLLINT patch | 分支不含,现场 `sed -i 's/if (XMLLINT)/if (FALSE)/'`,备份为 `CMakeLists.txt.orig-20260820_093910` |
| 并行度 | `make -j4`(用户指定 4 核) |

上传方式:本地 `git archive agent/exact-density-embedding --prefix=OpenMolcasUP_ExactDenEmb-65a3eeb11/ | gzip`(44M,天然不含 .git),解包时再加 `tar --exclude='.git*'` 双保险。

## 3. 完整编译过程时间线

- **09:39** 本地打包分支源码(44M);alpha 账号在 `/scratch/alpha` 下预建 `build-exactden-20260820_093910/`(工作区)与 `OpenMolcasUP_ExactDenEmb-65a3eeb11/`(安装目标)。注:alpha 无权 chown 给 shuo,故全程统一以 alpha 身份操作。
- **09:44 第 1 轮**(`driver.log`,`build_driver_exactden.sh`,make -j16,曾以 nohup 启动后按用户要求废弃):cmake 通过,make 推进至 2% 时卡死于 ExternalProject 从 GitHub clone libwignernj(GV800 无外网,详见报错 E1)。进程随后全部终止(仅杀本任务进程,未删任何文件)。
- **09:51 wignernj 离线解决**(`wignernj-{configure,build}.log`):用 7 月构建留下的离线包 `/scratch/shuo/noci-openmolcas-staging-20260715/libwignernj-0.6.0.tar.gz`(345K,world-readable)预装到 `build-exactden-.../wignernj-sys`;主 cmake 增加 `CMAKE_PREFIX_PATH` 使 `find_package(wignernj 0.6 CONFIG QUIET)`(CMakeLists.txt:1282)直接命中,configure 日志确认 `Using system libwignernj 0.6.0 (...)`,彻底跳过 clone。
- **09:53 第 2 轮**(`driver2-j4.log`,`build_driver2_j4.sh`,make -j4,前台 ssh 实时监控,新构建目录 `openmolcas-build-j4`):cmake 通过(Configuring done / Generating done),make 至 74% 时失败于 `exactemb_data.F90`(报错 E2)。安装未执行。
- **10:0x 第 3 轮**(`driver3-fix.log`):用户在本地修好 `exactemb_data.F90`(public 语句移入模块说明部分),检查确认后仅上传该单文件 + diff 溯源记录,增量重编;`exactemb_data.F90` 通过,但 74% 处 `rdinp_exactemb.F90` 接连报 5 个错(报错 E3,单一根因)。
- **10:1x 第 4 轮**(`driver4-wpfix.log`):ZCode 在本地工作区把 `rdinp_exactemb.F90` 第 13 行改为 `use Definitions, only: wp, iwp, u6`(纯机械补导入),上传单文件增量重编,**make 100% 通过,make install 成功**(10:15)。
- **10:15 验证**:`BUILD_STATUS=OK`;`.molcasversion = o65a3eeb11dcbde`(精确对应分支提交);`bin/` 49 个 exe 含 exactemb.exe、expbas.exe、guga.exe、gugadrt.exe、gugaci.exe;`ldd exactemb.exe` 零缺失。

## 4. 报错记录(全部)

### E1:libwignernj GitHub clone 失败(网络,第 1 轮)

原文:

```
fatal: 无法访问 'https://github.com/susilehtola/libwignernj.git/'：Failed to connect to github.com port 443 after 136206 ms: Couldn't connect to server
正克隆到 'wignernj'...
```

根因:OpenMolcas 未找到系统 libwignernj 时走 ExternalProject 从 `https://github.com/susilehtola/libwignernj.git` 现场 clone;GV800 无直连外网,AGENTS.md 记载的 6080 代理当时未监听(`ss -tln | grep 6080` 无结果),clone 反复超时重试,make 卡死。7 月 16 日 DownBas 构建时 clone 曾成功(当时网络条件不同),并在 staging 留下了离线包。

修复:离线预装 wignernj 0.6.0 到 `wignernj-sys`(用 `$OMC/bin/{gcc,gfortran}` 编),cmake 加 `-DCMAKE_PREFIX_PATH="<wignernj-sys>;$GA;$OMC"` 使 `find_package(wignernj 0.6 CONFIG)` 命中,输出 `Using system libwignernj 0.6.0`,完全绕开 clone。副作用见 §6 运行注意(RPATH 依赖 wignernj-sys 目录)。

### E2:`exactemb_data.F90:37` PUBLIC 语句位置非法(语法,第 2 轮)

原文:

```
/scratch/alpha/.../src/exactemb/exactemb_data.F90:37:7:
   37 | public :: Reset_ExactEmb_Input
      |       1
Error: PUBLIC statement at (1) is only allowed in the specification part of a module
make[2]: *** [CMakeFiles/exactemb/CMakeFiles/exactemb_obj.dir/build.make:78：.../exactemb_data.F90.o] 错误 1
```

根因:`public :: Reset_ExactEmb_Input` 写在了 `contains` 之后、`end subroutine` 之后的模块子程序区;Fortran 规定 `public` 语句只能出现在模块说明部分(`contains` 之前)。该文件是全库唯一编译错误,其余模块(caspt2 等)均正常。

修复(用户完成):把 `Reset_ExactEmb_Input` 并入模块说明部分的 public 声明(现第 23-24 行两条 `public ::`),删除 contains 后的非法行。

### E3:`rdinp_exactemb.F90` 缺 `wp` 导入(连锁 5 错,单一根因,第 3 轮)

原文(5 条,gfortran 上限截断):

```
rdinp_exactemb.F90:80:32:  if (ElectronTolerance <= 0.0_wp) then
Error: Missing kind-parameter at (1)
rdinp_exactemb.F90:83:3:   end if
Error: Expecting END SUBROUTINE statement at (1)
rdinp_exactemb.F90:123:14: real(kind=wp), intent(out) :: Value
Error: Symbol 'wp' at (1) has no IMPLICIT type
rdinp_exactemb.F90:121:32: subroutine ReadRealLine(LU,Value,What)
Error: Symbol 'value' at (1) has no IMPLICIT type
rdinp_exactemb.F90:47:60:  call ReadRealLine(LuSpool,ElectronTolerance,'TOLERANCE')
Error: Type mismatch in argument 'value' at (1); passed REAL(8) to UNKNOWN
```

根因分析:宿主子程序 `RdInp_ExactEmb` 第 13 行 `use Definitions, only: iwp, u6` 缺 `wp`。第 80 行的 `0.0_wp`(E3-1)和内部子程序 `ReadRealLine` 的 `real(kind=wp)`(E3-3)直接依赖 `wp`;E3-2 是 E3-1 失败后 if 块未打开导致 `end if` 孤立;E3-4、E3-5 是 E3-3 声明失败后哑元 `Value` 类型未知的连锁。文件里 `ReadTextLine/ReadIntLine/ReadRealLine` 是宿主的内部子程序(整个文件以 `end subroutine RdInp_ExactEmb` 收尾),`iwp`/`u6` 经 host association 已可用;`Abend`(src/system_util/abend.F90 独立例程)与 `RdNLst` 为 F77 外部例程,均无需 use 导入——故唯一缺口就是 `wp`。

修复(ZCode 完成):第 13 行改为 `use Definitions, only: wp, iwp, u6`。同轮已预检模块其余文件:`exactemb.F90`(wp/iwp/u6 等导入齐全)、`main.F90`(用 `DefInt`,确认是 Definitions 模块合法导出,definitions.F90:29/51),均无同类问题。

### 过程性插曲(非编译报错,如实记录)

第一次清理旧进程时 `pkill -f "build-exactden-..."` 的模式匹配到了承载该命令的 `bash -c` 自身命令行,ssh 会话自杀(exit 255);实际已顺带杀掉目标进程,复查确认无残留。教训:pkill/pgrep 匹配模式需用 `[b]uild_driver` 类字符类技巧避免自匹配。

## 5. 最终验证结果

- `BUILD_STATUS=OK`(GV800 `build-exactden-20260820_093910/STATUS-j4.txt`)
- `.molcasversion = o65a3eeb11dcbde`(对应提交 65a3eeb11)
- `bin/` 49 个 exe:含 **exactemb.exe**(1.5MB)、expbas.exe、guga.exe、gugadrt.exe、gugaci.exe
- `ldd exactemb.exe` 与 `ldd gateway.exe` 均零 "not found":libhdf5.so.320、libxcf03.so.15、libmkl_scalapack_ilp64.so.2 等解析到 `/opt/mamba/envs/omc/lib`,libwignernj 解析到 `wignernj-sys/lib`

## 6. 运行注意与复跑

- **wignernj-sys 不可清理**:已装二进制 RPATH 指向 `build-exactden-20260820_093910/wignernj-sys/lib/libwignernj_f03.so.0`,删该目录则运行时报库缺失。
- **运行前环境**:`export LD_LIBRARY_PATH=/opt/mamba/envs/omc/lib:$LD_LIBRARY_PATH`(MKL/hdf5/libxc 均在 omc env 内)。
- **复跑**:以 alpha 身份在 `/scratch/alpha/build-exactden-20260820_093910` 下执行 `bash build_driver2_j4.sh` 即可增量编译 + 安装;改源码后只需 scp 对应文件到 `src/OpenMolcasUP_ExactDenEmb-65a3eeb11/` 再跑。
- **本地未提交修复**(建议尽快提交到分支):`src/exactemb/exactemb_data.F90`(用户修,public 位置)+ `src/exactemb/rdinp_exactemb.F90`(ZCode 修,补 wp 导入)。

## 7. 产物与日志位置

- GV800:`/scratch/alpha/OpenMolcasUP_ExactDenEmb-65a3eeb11/`(安装产物);`/scratch/alpha/build-exactden-20260820_093910/`(源码树含 XMLLINT patch、openmolcas-build 与 openmolcas-build-j4 两个构建目录、wignernj-sys、全部 driver*.log / openmolcas-*.log / wignernj-*.log / env*.txt / STATUS*.txt 与三个驱动脚本,均未删除)。
- 本地:`~/develop/tmp/gv800_exactden_build_20260820_093910/`(全部日志与脚本副本 + BUILD_REPORT.md 过程报告);本文档为项目 TODO 下的正式记录。

## 8. 11:04–11:09 后续 RunFile 接口修复与重编译

- 修复 `SOURCE`、`TARGET`、`CLEAR` 把长物理文件名直接传给 8 字符 `NameRun` 的接口错误；解析器现先验证长度与路径分隔符，驱动注册 `$WorkDir/EMBSRC` 和 `$WorkDir/EMBTGT`。
- 全库动态 `NameRun` 调用审计确认：除 EXACTEMB 原有四处外，GENANO 和 AVERD 使用 6/7 字符 `RUNnnn`，SLAPAF 使用已注册的 `RUNOLD`/`RUNFILE`，未发现第二个长名字或物理路径调用。
- 11:04 首次后续重编译按本文 §2 原配方执行，CMake、`make -j4`、install、verify 全部返回 0；日志为实验目录 `logs/runfile-name-fix-*`。
- 完整输入作业 858 到达 CLEAR 后暴露相关遗漏：十一项 `Exact Emb *` 新记录均未加入 `RunFile_data` 的合法标签表，Release 构建在 `Exact Emb Active` 上以 temporary field 中止。
- 第二轮一次性把六项整数标量加入 `LabelsIS`、四项实数标量加入 `LabelsDS`、`Exact Emb Pot` 加入 `LabelsDA`，避免逐字段试错；11:08 再按同一配方重编译，四项状态均为 0，日志为 `logs/runfile-label-fix-*`。
- 完整输入作业 859 使用当前任务生成的 RunFile 和短名字 `EMBTGT`，以 `RUN_RC=0`、`Happy landing` 完成 CLEAR；输出没有 temporary field 警告，Stage 1 运行门通过。
