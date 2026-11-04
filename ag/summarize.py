"""Collect the fixed three-seed/control experiment and render its evidence."""
from pathlib import Path
import json, csv
import numpy as np
import matplotlib.pyplot as plt
from run import panel
from kmc_model import Config

ROOT = Path(__file__).resolve().parent
def main():
    out = ROOT/'results/study'
    out.mkdir(parents=True, exist_ok=False)
    results = []
    fig, axes = plt.subplots(3, 3, figsize=(13.5, 13.8))
    for row, variant in enumerate(['assist', 'off', 'oppose']):
        for col, seed in enumerate([7, 17, 27]):
            prefix = ROOT/f'results/seed{seed}/{variant}_seed_{seed}'
            result = json.loads(prefix.with_suffix('.json').read_text())
            data = np.load(prefix.with_suffix('.npz'))
            states = data['snapshots']
            cfg = Config(**result['config'])
            # Independently validate every archived snapshot against its ledger.
            for state, record in zip(states, result['trace']):
                assert np.count_nonzero(state[1:-1]) == record['net_source']
            assert result['ag_balance_error'] == 0
            result['archive'] = str(prefix.relative_to(ROOT))
            results.append(result)
            panel(axes[row, col], states[-1], cfg,
                  f"{variant} | seed {seed} | {result['status']}\nt={result['time_s']:.3g} s, Ag={result['trace'][-1]['metal']}")
    fig.suptitle('Ag / AlN(20 nm) / AlScN(20 nm) / AlN(20 nm) / Pt\nNatural CT-KMC growth | gold: Ag, teal: Ag+ | coarse-grained, uncalibrated', fontsize=15)
    fig.tight_layout(rect=(0, 0, 1, .95))
    fig.savefig(out/'morphology_comparison.png', dpi=180)
    fig.savefig(out/'morphology_comparison.svg'); plt.close(fig)
    chosen = np.load(ROOT/'results/seed7/assist_seed_7.npz')
    fig, ax = plt.subplots(figsize=(6.4, 6))
    panel(ax, chosen['snapshots'][-1], Config(**results[0]['config']),
          'Natural branched Ag bridge | 60 nm film\n3 V, 300 K | seed 7 | effective 1 nm lattice')
    ax.text(2, 5, 'AlN', color='#596275')
    ax.text(2, 25, 'AlScN', color='#596275')
    ax.text(25, 44, 'AlN', color='#596275')
    fig.text(.5, .01, 'Gold: Ag | Teal: Ag+ | exploratory, uncalibrated', ha='center', fontsize=9)
    fig.tight_layout(rect=(0, .03, 1, 1))
    fig.savefig(out/'filament.png', dpi=200); fig.savefig(out/'filament.svg'); plt.close(fig)
    # Initial E uses layer-internal differences, avoiding interface mixing.
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.4))
    for variant in ['assist', 'off', 'oppose']:
        data = np.load(ROOT/f'results/seed7/{variant}_seed_7.npz')
        phi = data['initial_potential_V'][1:-1].mean(axis=1)
        z = np.arange(len(phi))+.5
        axes[0].plot(z, phi, label=variant)
        e = -np.diff(phi)*10  # V/nm -> MV/cm
        axes[1].plot((z[:-1]+z[1:])/2, e, label=variant)
    for ax in axes:
        ax.axvspan(20, 40, color='#d6c9eb', alpha=.35)
        ax.set_xlabel('Depth from Ag (nm)'); ax.legend(); ax.grid(alpha=.2)
    axes[0].set(ylabel='Initial potential (V)', title='Fixed 3 V electrode constraint')
    axes[1].set(ylabel='Initial E_z (MV/cm)', title='AlScN field and compensating AlN fields')
    fig.tight_layout(); fig.savefig(out/'built_in_field.png', dpi=180); fig.savefig(out/'built_in_field.svg'); plt.close(fig)
    (out/'summary.json').write_text(json.dumps(results, indent=2), encoding='utf-8')
    with (out/'summary.csv').open('w', encoding='utf-8-sig', newline='') as f:
        keys = ['variant', 'seed', 'status', 'time_s', 'events', 'metal', 'ions', 'ag_balance_error']
        writer = csv.DictWriter(f, keys); writer.writeheader()
        for r in results:
            writer.writerow({k: r['config']['seed'] if k=='seed' else r['trace'][-1].get(k, r.get(k)) for k in keys})
    lines = ['# 基础模型首轮运行记录', '', '完整60 nm目标厚度，80 nm局部截面，1 nm有效格点；3 V、300 K。九条轨迹均从空介质开始。', '',
             '|内建场组|种子|停止原因|模型时间/s|事件数|中性Ag有效单元|Ag守恒误差|', '|---|---:|---|---:|---:|---:|---:|']
    for r in results:
        lines.append(f"|{r['variant']}|{r['config']['seed']}|{r['status']}|{r['time_s']:.6g}|{r['events']}|{r['trace'][-1]['metal']}|{r['ag_balance_error']}|")
    lines += ['', '三组描述统计（只对已成桥轨迹取均值；事件预算截断不计入）：', '', '|组别|成桥轨迹/总数|模型成桥时间均值/s|', '|---|---:|---:|']
    for variant in ['assist', 'off', 'oppose']:
        group = [r for r in results if r['variant']==variant]
        times = [r['time_s'] for r in group if r['status']=='first_bridge']
        avg = f'{np.mean(times):.6g}' if times else 'NA'
        lines.append(f'|{variant}|{len(times)}/{len(group)}|{avg}|')
    lines += ['', '![真实末态对照](morphology_comparison.png)', '', '![初始分层场](built_in_field.png)', '',
              'assist/off/oppose分别为P_eff=−0.004/0/+0.004 C/m²。assist指中层场方向，不保证全器件生长更快。每组仅3种子，结果是探索性描述，不能估计可靠成丝概率或证明极化因果。模型时间是未标定参数下的CT-KMC输出，不是实验时间预测。', '',
              'oppose组有事件预算截断，已成桥样本均值存在选择偏差，不可用三组均值直接排列成丝快慢。', '',
              '每帧Ag占据数已独立核对源库净注入。电势文件在first_bridge终态记录的是最后一次成桥事件执行之前的场，JSON的potential_state已标明；桥接后的金属电流/压降未求解。', '',
              '形貌可用于检查枝状生长与局部竞争，不作为TEM定量复现。尚无电流、限流、热、RESET/保持或极化翻转，完整边界见../../README.md。']
    (out/'运行报告.md').write_text('\n'.join(lines)+'\n', encoding='utf-8')
    print((out/'summary.csv').read_text(encoding='utf-8-sig'))

if __name__ == '__main__': main()
