"""Continue the exact reference-return RNG checkpoint to 1 s at zero bias."""
import os
os.environ['OPENBLAS_NUM_THREADS']='1'
os.environ['OMP_NUM_THREADS']='1'
import sys,json,hashlib,csv
from pathlib import Path
from dataclasses import asdict
import numpy as np
ROOT=Path(__file__).resolve().parents[2];sys.path.insert(0,str(ROOT))
from run_kmc import load_dynamic,save_dynamic,draw_morphology,draw_dynamic_results,write_json
from dynamic_kmc import stability
from morphology_metrics import metrics
OUT=Path(__file__).resolve().parent;SOURCE=ROOT/'results/morphology_reference_fast_return_20261112'
m,pos,records,states,labels=load_dynamic(SOURCE)
audit=json.loads((SOURCE/'manifest.json').read_text())
audit.update(continued_from=str(SOURCE),source_checkpoint_sha256=hashlib.sha256((SOURCE/'checkpoint.npz').read_bytes()).hexdigest(),
 driver_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
 extension='Fixed endpoint: total 1.000 s at zero; no parameter changes, no selection of an earlier stable-looking endpoint')
m.apply_command(0);pos.update(phase='zero_hold',deadline_s=None)
while pos['index']<50:
 deadline=m.time+.02
 m.advance(deadline,m.steps+2000000)
 assert abs(m.time-deadline)<1e-12
 m.record('extended_zero_hold_end');i=pos['index'];records.append(dict(phase='zero_hold',point_index=i,**m.iv[-1]))
 states.append(m.state.reshape(m.shape).copy());labels.append(f'zero_hold_{i:03d}');pos['index']+=1
 if pos['index']%5==0:
  print('zero_held_s',m.time-pos['zero_start_s'],'M',m.iv[-1]['metal'],'ions',m.iv[-1]['ions'],'connected',m.connected,flush=True)
  save_dynamic(OUT,m,pos,records,states,labels,'running_checkpoint',audit)
pos['phase']='settled'
zero=[r for r in records if r['phase']=='zero_hold' or (r['phase']=='waveform' and r['point_index']==65)]
start=zero[-6]['time_s'];continuous=[r for r in m.iv if start-1e-12<=r['time_s']<=m.time+1e-12]
g=np.array([r['network_conductance_S'] for r in continuous]);mcount=np.array([r['metal'] for r in continuous])
diagnostic=dict(all_connected=all(r['connected'] for r in continuous),
 network_conductance_relative_range=float(np.ptp(g)/np.median(g)),metal_relative_range=float(np.ptp(mcount)/np.median(mcount)),
 records=len(continuous),span_s=m.time-start)
summary=dict(status='reference_return_and_zero1s_complete',events=m.steps,time_s=m.time,zero_hold_s=m.time-pos['zero_start_s'],
 first_bridge=m.first_bridge,connected=m.connected,endpoint=m.iv[-1],counts=m.counts,net_source=m.net_source,
 ag_balance_error=int(np.count_nonzero(m.state[m.interior]))-m.net_source,
 reaction_balance_error=int(np.sum(m.state[m.interior]==2))-(m.counts['reduce']-m.counts['oxidize']),
 morphology=metrics(m.state.reshape(m.shape),m.c.spacing_nm,m.c.layer_cells),
 stability=stability(zero,m.topology_log),continuous_zero_bias_diagnostics=diagnostic,topology_log=m.topology_log,
 config=asdict(m.c),electronics=asdict(m.electronics),scope=audit['protocol_scope'])
save_dynamic(OUT,m,pos,records,states,labels,summary['status'],audit);write_json(OUT/'result.json',summary)
np.savez_compressed(OUT/'states.npz',state=m.state.reshape(m.shape),scan_states=states,labels=labels,first_bridge_state=m.first_bridge_state,
 site_energy_eV=m.site.reshape(m.shape),vacancy_mask=m.vacancy_mask.reshape(m.shape))
with (OUT/'forming_trace.csv').open('w',encoding='utf-8-sig',newline='') as f:
 w=csv.DictWriter(f,fieldnames=list(records[0]));w.writeheader();w.writerows(records)
draw_morphology(OUT,m);draw_dynamic_results(OUT,m,records,states,labels)
print(json.dumps({k:summary[k] for k in ('status','connected','events','zero_hold_s','stability','continuous_zero_bias_diagnostics')},indent=2),flush=True)
