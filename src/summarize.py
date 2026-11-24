import csv, json, statistics
import continuous_observer as m

P = m.RESULTS
rows=list(csv.DictReader((P/'continuous_confirmation_dt2.csv').open()));fine=[r for r in csv.DictReader((P/'continuous_equal_budget_confirmation.csv').open()) if float(r['dt'])==2];rows+=fine
lookup={(r['config'],r['truth'],r['method']):r for r in rows};tim=json.loads((P/'continuous_isolated_timing_metrics.json').read_text());out={'confirmation_configurations':12,'continuous_truths_per_configuration':12,'total_queries':144,'mean_convention':'Abstention deadline -1 is mapped to zero only for mean/gain summaries; approvals are separately reported.','methods':{}}
for name in ['single_box','state_bank','state_fine_bank','parameter_bank']:
 rr=[r for r in rows if r['method']==name];d=[max(0,float(r['deadline_s'])) for r in rr];gain=[x-max(0,float(lookup[(r['config'],r['truth'],'state_bank')]['deadline_s'])) for r,x in zip(rr,d)]
 out['methods'][name]={'mean_deadline_s':statistics.mean(d),'approvals':sum(int(r['approved']) for r in rr),'truth_not_retained':sum(1-int(r['truth_retained']) for r in rr),'unsafe_replays':sum(int(r['pre_cross'])+int(r['post_cross']) for r in rr),'mean_gain_vs_six_state_s':statistics.mean(gain),'median_gain_vs_six_state_s':statistics.median(gain),'improved_queries_vs_six':sum(v>0 for v in gain),'worse_queries_vs_six':sum(v<0 for v in gain),'median_decision_ms':tim[name]['median_ms']}
rr1=list(csv.DictReader((P/'continuous_confirmation_dt1.csv').open()))+[r for r in csv.DictReader((P/'continuous_equal_budget_confirmation.csv').open()) if float(r['dt'])==1];lookup1={(r['config'],r['truth'],r['method']):r for r in rr1};out['step_refinement']={'classification_changes':sum(int(r['approved'])!=int(lookup1[(r['config'],r['truth'],r['method'])]['approved']) for r in rows),'max_absolute_deadline_change_s':max(abs(float(r['deadline_s'])-float(lookup1[(r['config'],r['truth'],r['method'])]['deadline_s'])) for r in rows),'truth_membership_failures':sum(1-int(r['truth_retained']) for r in rr1),'unsafe_replays':sum(int(r['pre_cross'])+int(r['post_cross']) for r in rr1)}
(P/'metrics.json').write_text(json.dumps(out,indent=2));print(json.dumps(out,indent=2))
# Descriptive worked case selected by largest six-to-36-state gain, not used to tune the method.
best=max(fine,key=lambda r: max(0,float(r['deadline_s']))-max(0,float(lookup[(r['config'],r['truth'],'state_bank')]['deadline_s'])))
c=next(c for c in json.loads((m.DATA/'continuous_configs.json').read_text()) if c['id']==best['config']);q=c['truths'][int(best['truth'])];x,y=m.plant(c,q);e={'selection':'Largest descriptive six-to-36-state gain in the unchanged confirmation set; no tuning.','config':c['id'],'load':c['load'],'cooling':c['cooling'],'delay':c['delay'],'truth':int(best['truth']),'actual_state':x.tolist(),'observations':y.tolist(),'methods':{}}
for name in out['methods']:
 lo,hi,cl,ch,n=m.observe(c,y,name);d=m.deadline(c,hi,cl,ch);e['methods'][name]={'cells':len(hi),'counts':n,'Tc_upper':float(hi[:,0].max()),'Tn_upper':float(hi[:,1].max()),'deadline':float(d.min())}
(P/'worked_observation_case.json').write_text(json.dumps(e,indent=2));print(e)
