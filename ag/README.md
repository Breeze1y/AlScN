# Ag/AlN/AlScN/AlN/Pt：形貌优先基础 KMC

本目录是独立重建版本，不依赖 `ai4alscn/kmc`。参考该目录的 Zhao 2016 提炼、器件资料和静电推导；保留这里原有的 `papers/`。

GitHub副本包含本目录代码、配置与仿真输出；参考论文资料仍保留在本地。可从仓库根目录独立运行：

```powershell
python -m pip install -r ag/requirements.txt
python -m unittest discover -s ag -p test_physics.py -v
python ag/run.py --controls --output results/new_run
```

![60nm自然生长细丝](results/study/filament.png)

## 使用

从 `D:\Al4memorister` 运行（可复用已有环境）：

```powershell
.\ai4alscn\.venv\Scripts\python.exe ag/run.py --controls
.\ai4alscn\.venv\Scripts\python.exe -m unittest discover -s ag -p test_physics.py -v
```

独立环境安装 `ag/requirements.txt` 后也可运行。默认输出 `ag/results/baseline/`；目录已存在会停止，请用 `--output results/new_name` 指定新目录。`--seeds 7,17,27` 指定所有独立种子。配置在 `configs/baseline.json`，其余默认值见 `kmc_model.py` 的 `Config`，每条运行 JSON 会展开保存全部参数。

输出：`morphology.png/svg` 为末态对照，`evolution.png/svg` 为第一条轨迹的真实演化；每条 `npz` 保存逐帧占据状态、时刻、位点能和电势；JSON/CSV 保存停止原因、种子、事件数、计数与守恒。金色是中性 Ag，青色是 Ag⁺，虚线为材料界面。只有中性 Ag 四邻接贯通才算几何导电桥，青色离子不算。

## 文献和移植边界

优先来源：Zhao et al., *Atomic Monte-Carlo Simulation for CBRAM with Various Filament Geometries*, SISPAD 2016, pp.153–156，[会议原文](https://in4.iue.tuwien.ac.at/pdfs/sispad2016/SISPAD_2016_153-156.pdf)。图6明确比较 Ag/SiO₂/Pt 的模拟与 TEM 枝状细丝。原文提供氧化、离子输运、还原、金属迁移四类事件及相对速率影响形貌的依据；**本代码是面向氮化物的基础重建，不是原作者算法或图6参数的逐项复现。**

器件几何沿用 `ai4alscn/硕士课题完整执行方案-KMC主线-20261102.md`：Ag/AlN(20 nm)/AlScN(20 nm)/AlN(20 nm)/Pt。计算宽80 nm仅是局部截面窗口，不等于实验150 μm电极直径；横向绝缘边界。介电常数8.5/15来自既有静电工作假设，未作为本器件实测参数。

## 状态、事件和电场

二维有限体积格点，默认1 nm有效单元；三层各20行介质，额外两行作电极。每格为宿主、Ag⁺或Ag。粒子是有效占据单元，不是声称1 nm银单元只有一个原子。初始介质全空，无预置种子、定向模板或人为补桥。

1. 顶Ag库在第一层空格氧化供给Ag⁺；离子可返回库。入/出账本显式累计。
2. Ag⁺向四邻空格跳跃。AlN/AlScN体势垒分开，界面用两侧基准均值加对称跨界面屏障，离子能量为 `site + phi`（Ag⁺，eV和V数值对应）。
3. 仅有Pt或Pt连通Ag直接邻接的位置可还原，电子库电势0 V；裸Pt另加成核势垒。Ag可逆氧化回离子，配位键能稳定已有Ag。多个邻接金属构成同一个局部电子库，不重复计反应通道。
4. 表面Ag可中性迁移，按前后Ag–Ag配位能差计算聚集/等配位/脱离的不同率；不把Zhao的四个示例势垒硬套为氮化物实测值。此版宿主差异体现在输运和介电中，尚未单列Ag–宿主键能。

率采用 `ν exp[-(B + max(ΔG,0))/kT]`，正逆共享B；跳跃和中性迁移在冻结场下满足局部详细平衡。电极是开放库，阴极标准金属化能是可调工作假设。每次金属构型变化都重新识别连通团簇、求解电场；离子跳跃不更新电势（忽略离子空间电荷，假定背景屏蔽）。这不是已包含所有反应自由能和金属电荷的热力学完整模型，尤其不以局部速率比检查代替全体系闭合证明。

总率 `R=ΣΓ`；等待时间 `-ln(u1)/R`，独立均匀随机数u2选择事件。事件超过脉冲结束时间则不执行该事件。固定3 V、300 K，默认2 s或350000事件预算，到首次金属桥即停止。`event_budget` 表示轨迹截断，不算有效成桥或物理失败；`time_limit` 表示给定脉冲内未成桥。此时钟依赖未标定的有效势垒/前因子，不可拿来预测实验成丝秒数。

电势满足 `∇·(ε∇φ)=∇·P`，即 `D=−ε∇φ+P`。AlScN之外P=0；Ag源φ=V，Pt及Pt连通金属φ=0，孤立金属以高介电常数近似浮置包含物。面介电系数取调和平均，面极化取与两侧介电匹配的加权值；端电压和极化界面电荷在同一方程求解，因此外侧AlN自动出现补偿场。金属体内P=0，电极连通Ag对局部场有屏蔽和尖端增强。

默认 `P_eff=-0.004 C/m²` 是**屏蔽后的净有效极化**，既不是已测Pr，也不是从宏观开关电压唯一反演的场。z从Ag指向Pt；负P让AlScN退极化场朝+z，故称中层assist，对外层则可能减弱场。零偏、无细丝时，中层约+0.235 MV/cm，外层约−0.117 MV/cm，总压降为0。`--controls` 对比相同种子和无序场下的P、0、−P；比较的是整条随机过程，不能保证不同参数下事件一一配对，也不能假定assist必定缩短整体成桥时间。

## 当前范围

这是等温、首次forming的形貌模型。尚无电流连续网络、接触二极管、焦耳热、1 mA限流反馈、RESET/保持、极化翻转、氮空位或陷阱动力学。几何金属桥意味着模型中存在金属通路，尚未计算其电导，也不代表完整实验开关验证。停止于第一桥可以避免在尚无电流/限流闭环时继续非物理加厚。

本版所有有效动力学势垒（默认0.48/0.50 eV输运、0.68 eV源反应、0.52 eV还原/氧化基准、0.08 eV裸Pt成核附加、0.80 eV中性迁移）、0.06 eV键能、0.025 eV固定无序宽度都是探索参数。原文Ag/SiO₂参数不适用于直接预测AlN/AlScN。先看真实生长形貌和多种子稳定性，再受实测约束做输运/供给/成核比率、P、宽度/网格扫描。不可把一张枝状图作为Ag机制或真实材料势垒已证实的证据。

物理检查见 `test_physics.py`：三层极化解析解、零偏补偿与P反向、跳跃局部详细平衡、开放Ag库守恒、金属/离子连通区分。

## 已交付结果

汇总见 [九轨迹运行报告](results/study/运行报告.md)、[形貌对照](results/study/morphology_comparison.png)、[代表细丝](results/study/filament.png)、[分层场对照](results/study/built_in_field.png)、[自然演化](results/seed7/evolution.png)。三种子原始记录分别在 `results/seed7`、`results/seed17`、`results/seed27`；每批保存运行前源码快照及SHA256清单。

`summarize.py` 汇总上述三个固定目录，逐帧核对占据数与Ag源库账本，再输出 `results/study`；不挑选最佳种子，统一展示7、17、27。`results/pilot`是首条调通运行，汇总统计不重复计入。贯通终态NPZ内的电势是成桥事件前的场，见JSON的 `potential_state`，不能将其当作导通后的电势。

配位能中的接触Pt也计作有效金属邻接；本版尚未为Ag–Pt和Ag–Ag分别标定黏附能。中性迁移的独立基准较高，在当前室温短脉冲内可能几乎不发生，必须从事件计数检查其实际贡献。
