import csv, json
import continuous_observer as m

P = m.RESULTS
P.mkdir(parents=True, exist_ok=True)
(P/'model_stress_protocol.json').write_text(json.dumps({'frozen_before_run':'2026-11-14','selection':'First four truths of each of the 12 unchanged confirmation physical configurations.','perturbation':'Loss multiplier fixed at 1.25 and conductance at 0.80, outside assumed [0.9,1.1]. Initial states, capacities, measurement errors remain unchanged.','purpose':'Explicit misspecified-physics sensitivity; not part of the stated bounded-uncertainty guarantee.'},indent=2))
rows=[]
for c in json.loads((m.DATA/'continuous_configs.json').read_text()):
 if c['split']!='confirmation':continue
 for j,q0 in enumerate(c['truths'][:4]):
  q={**q0,'loss':1.25,'conductance':.8};x,y=m.plant(c,q)
  for name in ['single_box','state_bank','state_fine_bank','parameter_bank']:
   lo,hi,cl,ch,n=m.observe(c,y,name);d=float(m.deadline(c,hi,cl,ch).min()) if len(hi) else -1
   rows.append({'config':c['id'],'truth':j,'method':name,'deadline_s':d,'cells':len(hi),**m.check(c,q,x,d)})
with (P/'model_stress.csv').open('w') as f:w=csv.DictWriter(f,fieldnames=rows[0]);w.writeheader();w.writerows(rows)
out={name:{'queries':48,'approvals':sum(r['approved'] for r in rows if r['method']==name),'empty_sets':sum(r['cells']==0 for r in rows if r['method']==name),'guard_crossings':sum(r['pre_cross']+r['post_cross'] for r in rows if r['method']==name)} for name in ['single_box','state_bank','state_fine_bank','parameter_bank']};(P/'model_stress_metrics.json').write_text(json.dumps(out,indent=2));print(out)
