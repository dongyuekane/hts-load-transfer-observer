import csv, json, random, time
import continuous_observer as m
P = m.RESULTS
P.mkdir(parents=True, exist_ok=True)
np = m.np
cfgs=[c for c in json.loads((m.DATA/'continuous_configs.json').read_text()) if c['split']=='confirmation'];rows=[]
for c in cfgs:
 for j,q in enumerate(c['truths']):
  x,y=m.plant(c,q)
  for dt in [2.,1.]:
   st=time.perf_counter();lo,hi,cl,ch,n=m.observe(c,y,'state_fine_bank',dt);ds=m.deadline(c,hi,cl,ch,dt) if len(hi) else np.array([-1]);d=float(ds.min());sec=time.perf_counter()-st
   retained=bool(np.any(np.all((x>=lo-1e-6)&(x<=hi+1e-6),axis=1)&np.all((np.array(q['cap'])>=cl-1e-9)&(np.array(q['cap'])<=ch+1e-9),axis=1))) if len(hi) else False
   rows.append({'config':c['id'],'truth':j,'method':'state_fine_bank','dt':dt,'deadline_s':d,'surviving_cells':len(hi),'truth_retained':int(retained),'cpu_seconds':sec,**m.check(c,q,x,d)})
 print(c['id'],flush=True)
with (P/'continuous_equal_budget_confirmation.csv').open('w') as f:w=csv.DictWriter(f,fieldnames=rows[0]);w.writeheader();w.writerows(rows)
out={}
for dt in [2.,1.]:
 rr=[r for r in rows if r['dt']==dt];out[str(dt)]={'mean_deadline_s':float(np.mean([max(0,r['deadline_s']) for r in rr])),'approvals':sum(r['approved'] for r in rr),'truth_not_retained':sum(1-r['truth_retained'] for r in rr),'pre_cross':sum(r['pre_cross'] for r in rr),'post_cross':sum(r['post_cross'] for r in rr)}
(P/'continuous_equal_budget_confirmation_metrics.json').write_text(json.dumps(out,indent=2));print(out,flush=True)
# Formal all-method decision-cost benchmark: one predeclared query from each configuration, three randomized paired repetitions.
inputs=[(c,*m.plant(c,c['truths'][0])) for c in cfgs];tr=[];rng=random.Random(2026111466)
for rep in range(3):
 for c,x,y in inputs:
  methods=['single_box','state_bank','state_fine_bank','parameter_bank'];rng.shuffle(methods)
  for name in methods:
   st=time.perf_counter();lo,hi,cl,ch,n=m.observe(c,y,name,2.);ds=m.deadline(c,hi,cl,ch,2.) if len(hi) else np.array([-1]);d=float(ds.min());sec=time.perf_counter()-st
   tr.append({'config':c['id'],'truth':0,'rep':rep,'method':name,'deadline_s':d,'decision_seconds':sec})
with (P/'continuous_isolated_timing.csv').open('w') as f:w=csv.DictWriter(f,fieldnames=tr[0]);w.writeheader();w.writerows(tr)
(P/'continuous_isolated_timing_metrics.json').write_text(json.dumps({name:{'median_ms':float(np.median([r['decision_seconds']*1000 for r in tr if r['method']==name])),'total_s':sum(r['decision_seconds'] for r in tr if r['method']==name)} for name in methods},indent=2))
