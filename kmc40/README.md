# KMC40 最终连续 100 圈

本版从初始化 RESET 检查点开始，连续计算 100 次 SET 和 RESET。断桥后，已开始退缩的 Ag 位点共享回流通量。每个阶段含 202 个 I–V 采样点，完整轨迹共 40,400 点。

| 查看内容 | 入口 |
|---|---|
| 逐圈实测与模拟 I–V、两组各自的 100 圈叠加 | [交互式 HTML](replay/cycles100/report/iv_comparison.html) |
| 最终连续 100 圈有符号 I–V | [独立 CSV](replay/reset68/continuous/iv/iv_all.csv) |
| 100 圈形貌和第 68、100 圈断桥快照 | [形貌 PDF](replay/cycles100/report/morphology_100cycles.pdf) |
| 结果解读与全周期读出 | [中文总结](replay/cycles100/report/summary.md) |
| 与原 KMC 第 68 圈的比较 | [局部放大与叠加](replay/reset68/report/reset68_comparison.html) |
| 各圈原始状态、检查点与运行说明 | [连续轨迹](replay/reset68/continuous/README.md) |

在本目录运行 `.\.venv\Scripts\python.exe -B run_100_cycles.py` 可从已有检查点续算。已完成的轨迹可以直接打开报告，或用 `$kmc-output` 从 `replay/reset68/continuous` 导出。

复算脚本仍从 `results/cycle_reset_fit40/mid_clock_set` 读取初始化检查点、参数和前五圈的参考输入，从 `replay/cycles100` 读取后续圈的原始参考事件、检查点和基线数据。最终模拟的 CSV 和逐圈状态位于 `replay/reset68/continuous`。
