# SET/RESET 机器学习数据集

创建日期：2026-11-04

## 文件

- cycle_features.csv：每个候选SET/RESET回环一行，只含扫描设置和SET/RESET波形派生输入特征。
- cycle_targets.csv：读态响应目标。缺少同偏压读数时，目标字段留空并由label_read_pair_available标记。
- cycle_waveform_points.csv：自包含的原始候选波形点表，保留原始带符号电流、R列、点序、数据块、源行号、限流设置和限流点标记。
- data_dictionary.csv：逐表逐字段的数据类型、单位和含义。
- split_assignments.csv：样本级分组和两种候选留出方案。
- split_group_manifest.csv：按完整条件目录分组的训练/验证/测试分配。
- common_waveform_protocols.csv：两批具有相同电压端点、步进和读偏压的协议交集；waveform_protocol_inventory.csv则列出全部协议。compliance单独保留，不能忽略其差异。
- source_file_manifest.csv：波形源文件的相对路径、SHA-256和点数。
- excluded_or_unpaired_records.csv：未配对、编号歧义和需复核的记录。
- model_input_feature_columns.json：模型输入列白名单和不应作为输入的标识字段。
- dataset_manifest.json：样本数、分组数和输出文件校验值。
- nature_data_source_verification.json、A_batch_source_verification.json、A_batch_duplicate_cleanup_log.csv：处理校验和A批次精确重复清理记录。

## 标签与输入

推荐的第一阶段监督目标为cycle_targets.csv中的target_log10_HRS_over_LRS（回归）或target_two_state_ge_20（分类）。只有SET后和RESET后读偏压相同时才提供目标；阈值20是现有方案中的暂定筛选线，不是普适物理判据。target_stable_calibration_provisional仅供敏感性分析，不应直接当成确定的器件真值。

模型输入从cycle_features.csv按model_input_feature_columns.json选取。不要把样本ID、批次/目录/运行分组、读态电阻、读阻比或目标标签加入X。原始长表与每行特征通过sample_id连接。CSV中的缺失数值为空；读取时应保留为空值，不要自动补零。

## 当前标签分布

共2381组有同偏压读态目标：2287组达到20倍阈值，94组未达到；另有1093组目标缺失。正负样本不平衡，训练分类器时需在训练折内处理类别权重或采样，并以按条件组划分的验证/测试数据评估。标签缺失不能解释为负类。

## 数据划分

directory_group_split_draft将整个批次/条件目录放入同一个训练、验证或测试集合，不能再按单个回环随机拆分。它是可复现的初始方案，不保证每个子集都覆盖所有标签或协议，模型训练前应检查分布并冻结最终划分。

cross_batch_holdout_candidate只标记两批间电压起止、步进和读偏压相同的条件；当前只有2种共享波形协议，涉及482组nature_data训练候选和1699组A批次外部验证候选。共享协议的主通道compliance并不完全匹配，须报告差异并检查限流删失。该划分只是候选方案，不等于已验证的跨器件泛化。

## 使用边界

本目录包含3474组候选回环和1,627,410个波形点；其中2381组有同读偏压标签，2287组读阻比达到20，988组带暂定稳定标定标记。两个批次分开标识。

SET/RESET配对依据目录、文件编号和时间戳推断，没有独立物理循环ID。达到B1500 compliance的电流点可能被限流。所有回环都保留在特征表和波形表中，质量筛选通过显式标志提供，不静默删除。读阻目标的量程和测量有效性仍需复核后才能作为高精度连续标签。
