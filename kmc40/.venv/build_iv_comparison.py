from __future__ import annotations

import csv
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
TARGET = ROOT / 'replay/cycles100/report/iv_comparison.html'
PHASES = ('SET', 'RESET')
BRANCHES = ('outward', 'return')
SWEEPS = {'IV Sweep +': 'SET', 'IV Sweep -': 'RESET'}


def source_rows(path):
    with path.open(encoding='utf-8-sig', newline='') as stream:
        yield from csv.DictReader(stream)


def branch_for(point):
    return BRANCHES[0] if point < 101 else BRANCHES[1]


def split(points):
    assert len(points) == 202
    return dict(zip(BRANCHES, (points[:101], points[101:])))


reference = {phase: [] for phase in PHASES}
for row in source_rows(ROOT / 'inputs/experimental/cycle074.csv'):
    reference[row['phase']].append((float(row['voltage_V']), float(row['current_A'])))
grid = {phase: [p[0] for p in reference[phase]] for phase in PHASES}


def load_observed():
    curves = {}
    path = ROOT.parent / 'cycle_data/organized/all_iv_points.csv'
    for row in source_rows(path):
        if row['device_label'] != 'AlScN_40nm':
            continue
        assert row['source_file'] == 'D1-104.csv'
        n, phase = int(row['cycle']), SWEEPS[row['sweep']]
        section = curves.setdefault(n, {}).setdefault(phase, [])
        point = int(row['point_index']) - 1
        assert point == len(section) and float(row['voltage_V']) == grid[phase][point]
        section.append(float(row['current_I1_A']))
    assert sorted(curves) == list(range(1, 105))
    assert all(curves[74][phase] == [p[1] for p in reference[phase]] for phase in PHASES)
    return [{phase: split(curves[n][phase]) for phase in PHASES} for n in sorted(curves)]


def load_simulated():
    curves = {}
    for row in source_rows(ROOT / 'replay/reset68/continuous/iv/iv_all.csv'):
        n, phase = int(row['cycle']), row['phase']
        section = curves.setdefault(n, {}).setdefault(phase, [])
        point = int(row['point'])
        assert point == len(section) and float(row['voltage_V']) == grid[phase][point]
        assert row['branch'] == branch_for(point)
        section.append(float(row['current_A']))
    assert sorted(curves) == list(range(1, 101))
    return [{phase: split(curves[n][phase]) for phase in PHASES} for n in sorted(curves)]


voltage = {phase: split(grid[phase]) for phase in PHASES}
observed = load_observed()
simulated = load_simulated()
data = dict(voltage=voltage, observed=observed, simulated=simulated)
payload = json.dumps(data, separators=(',', ':'), allow_nan=False)
template = Path(__file__).with_name('iv_comparison_template.html').read_text(encoding='utf-8')
assert template.count('__IV_DATA__') == 1
TARGET.write_text(template.replace('__IV_DATA__', payload), encoding='utf-8')
print(f'{len(observed)} measured, {len(simulated)} simulated, 404 samples per cycle\n{TARGET}')
