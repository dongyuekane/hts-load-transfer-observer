"""Continuous uncertainty, observation-conditioned comparison-system bank.
Differential comparison is analytic; RK4 enclosure implementation is numerical.
"""
import argparse, os, json, csv, time, itertools, random
from pathlib import Path
for k in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS'):os.environ[k]='1'
import numpy as np
from scipy.integrate import solve_ivp
import thermal_model as b

P = Path(__file__).resolve().parents[1]
DATA = P / 'data'
RESULTS = P / 'results'
EPS=.1;GUARD=.15;MAXWAIT=1200.;DT=2.;LOSS=(.9,1.1);COND=(.9,1.1)

def field(x,I,u,clo,chi,upper=True,loss=LOSS,cond=COND):
 a,n=x[...,0],x[...,1];delta=a-n
 if upper:
  hc=np.where(delta>=0,cond[0],cond[1]);hn=np.where(delta>=0,cond[1],cond[0]);q=b.heat(I,a)*loss[1]
 else:
  hc=np.where(delta>=0,cond[1],cond[0]);hn=np.where(delta>=0,cond[0],cond[1]);q=b.heat(I,a)*loss[0]
 f=np.stack([q-hc*b.H*delta,hn*b.H*delta+b.Q0-u*b.G*(n-b.TS)],axis=-1)
 den=np.array([b.CC,b.CN])*np.where(f>=0,clo if upper else chi,chi if upper else clo)
 return f/den
def step(x,I,u,clo,chi,dt,upper=True,loss=LOSS,cond=COND):
 f=lambda z:field(z,I,u,clo,chi,upper,loss,cond)
 k1=f(x);k2=f(x+dt*k1/2);k3=f(x+dt*k2/2);k4=f(x+dt*k3)
 return x+dt*(k1+2*k2+2*k3+k4)/6
def holding(I,u,loss=LOSS,cond=COND):
 aa=np.arange(b.TS+.01,float(b.bound(I))-GUARD+1e-9,.01)
 q=b.heat(I,aa)*loss[1];nn=aa-q/(b.H*cond[0]);ok=(nn>=b.TS)&(b.H*cond[1]*(aa-nn)+b.Q0<=u*b.G*(nn-b.TS));aa=aa[ok];nn=nn[ok]
 if len(nn):nn=np.maximum.accumulate(nn[::-1])[::-1]
 return aa,nn
def deadline(cfg,z,clo,chi,dt=DT,loss=LOSS,cond=COND):
 I=b.electrical(cfg['load'])[0][4];Ib=b.electrical(cfg['load'])[0][2];aa,nn=holding(Ib,cfg['cooling'],loss,cond)
 x=z.copy();d=np.full(len(x),-1.);valid=x[:,0]<b.bound(I)-GUARD;prefix=valid.copy()
 if not len(aa):return d
 for tt in np.arange(0,MAXWAIT+cfg['delay']+dt/2,dt):
  if tt>=cfg['delay']-1e-7:
   idx=np.searchsorted(aa,x[:,0]);mem=(idx<len(aa))&(x[:,1]<=nn[np.minimum(idx,len(aa)-1)])
   prefix&=valid&mem;d[prefix]=max(0.,tt-cfg['delay'])
  if tt>=MAXWAIT+cfg['delay'] or not np.any(valid):break
  xn=step(x,I,cfg['cooling'],clo,chi,dt,True,loss,cond);valid&=np.maximum(x[:,0],xn[:,0])<b.bound(I)-GUARD;x=np.where(valid[:,None],xn,x)
 return d
def observe(cfg,y,method,dt=DT):
 nstate=1 if method=='single_box' else 36 if method=='state_fine_bank' else 6;ncc=3 if method=='parameter_bank' else 1;ncn=2 if method=='parameter_bank' else 1
 edges=[np.linspace(77.8,83.8,nstate+1),np.linspace(.8,1.2,ncc+1),np.linspace(.8,1.2,ncn+1)];cells=list(itertools.product(range(nstate),range(ncc),range(ncn)))
 lo=np.array([[edges[0][i],y[0]-EPS] for i,j,k in cells]);hi=np.array([[edges[0][i+1],y[0]+EPS] for i,j,k in cells]);clo=np.array([[edges[1][j],edges[2][k]] for i,j,k in cells]);chi=np.array([[edges[1][j+1],edges[2][k+1]] for i,j,k in cells]);alive=np.ones(len(cells),bool);I=b.electrical(cfg['load'])[0][4]
 history=[]
 for tt in np.arange(0,180,dt):
  lnew=step(lo,I,cfg['cooling'],clo,chi,dt,False);hnew=step(hi,I,cfg['cooling'],clo,chi,dt,True)
  # The ordinary numerical comparison is not outward-rounded interval arithmetic.
  lo=np.where(alive[:,None],lnew,lo);hi=np.where(alive[:,None],hnew,hi)
  for j,obs_t in enumerate([60,120,180],1):
   if abs(tt+dt-obs_t)<1e-8:
    lo[:,1]=np.maximum(lo[:,1],y[j]-EPS);hi[:,1]=np.minimum(hi[:,1],y[j]+EPS);alive&=np.all(lo<=hi+1e-10,axis=1)
    history.append(int(alive.sum()))
 return lo[alive],hi[alive],clo[alive],chi[alive],history
def plant(cfg,q):
 I=b.electrical(cfg['load'])[0][4]
 f=lambda t,x:b.rhs(x,I,cfg['cooling'],np.array(q['cap']),np.array(q['cap']),q['loss'],q['conductance'])
 sol=solve_ivp(f,[0,180],q['x0'],t_eval=[0,60,120,180],rtol=1e-10,atol=1e-11,max_step=2)
 return sol.y[:,-1],sol.y[1]+q['noise']
def check(cfg,q,x,wait):
 if wait<0:return {'approved':0,'pre_cross':0,'post_cross':0}
 I=b.electrical(cfg['load'])[0][4];Ib=b.electrical(cfg['load'])[0][2]
 def run(z,current,dur):
  def ev(t,v):return float(b.bound(current))-GUARD-v[0]
  ev.terminal=True;ev.direction=-1
  if z[0]>=b.bound(current)-GUARD:return z,True
  sol=solve_ivp(lambda t,v:b.rhs(v,current,cfg['cooling'],np.array(q['cap']),np.array(q['cap']),q['loss'],q['conductance']),[0,dur],z,events=ev,rtol=1e-10,atol=1e-11,max_step=10)
  return sol.y[:,-1],bool(len(sol.t_events[0]))
 z,pre=run(x,I,wait+cfg['delay']);z,post=run(z,Ib,3600) if not pre else (z,False)
 return {'approved':1,'pre_cross':int(pre),'post_cross':int(post)}
def run(split, dt=DT, *, config_limit=None, truth_limit=None, output=RESULTS):
 output.mkdir(parents=True, exist_ok=True)
 configs=[c for c in json.loads((DATA/'continuous_configs.json').read_text()) if c['split']==split][:config_limit];rows=[];evidence=[];rng=random.Random(2026111497)
 for cfg in configs:
  for j,q in enumerate(cfg['truths'][:truth_limit]):
   x,y=plant(cfg,q);methods=['single_box','state_bank','parameter_bank'];rng.shuffle(methods)
   for m in methods:
    tic=time.perf_counter();lo,hi,clo,chi,counts=observe(cfg,y,m,dt);d=deadline(cfg,hi,clo,chi,dt) if len(hi) else np.array([-1.]);value=float(d.min());elapsed=time.perf_counter()-tic
    truth_member=np.any(np.all((x>=lo-1e-6)&(x<=hi+1e-6),axis=1)&np.all((np.array(q['cap'])>=clo-1e-9)&(np.array(q['cap'])<=chi+1e-9),axis=1)) if len(hi) else False
    replay=check(cfg,q,x,value);r={'config':cfg['id'],'truth':j,'method':m,'dt':dt,'deadline_s':value,'surviving_cells':len(hi),'truth_retained':int(truth_member),'cpu_seconds':elapsed,'Tc_true':float(x[0]),'Tn_true':float(x[1]),**replay};rows.append(r)
    if split=='confirmation' and j==0:evidence.append({'config':cfg['id'],'method':m,'observations':y.tolist(),'actual_state':x.tolist(),'actual_parameters':q,'lower':lo.tolist(),'upper':hi.tolist(),'cap_lower':clo.tolist(),'cap_upper':chi.tolist(),'deadlines':d.tolist()})
  print(split,cfg['id'],'step',dt,flush=True)
 name=f'continuous_{split}_dt{dt:g}';
 with (output/(name+'.csv')).open('w') as f:w=csv.DictWriter(f,fieldnames=rows[0]);w.writeheader();w.writerows(rows)
 lookup={(r['config'],r['truth'],r['method']):r for r in rows};out={'configs':len(configs),'truths':len(rows)//3,'step_s':dt}
 for m in ['single_box','state_bank','parameter_bank']:
  rr=[r for r in rows if r['method']==m];g=np.array([max(0,r['deadline_s'])-max(0,lookup[(r['config'],r['truth'],'single_box')]['deadline_s']) for r in rr]);g2=np.array([max(0,r['deadline_s'])-max(0,lookup[(r['config'],r['truth'],'state_bank')]['deadline_s']) for r in rr]);out[m]={'mean_deadline_s':float(np.mean([max(0,r['deadline_s']) for r in rr])),'mean_gain_vs_single_s':float(np.mean(g)),'median_gain_vs_single_s':float(np.median(g)),'improved_vs_single':int(np.sum(g>0)),'mean_gain_vs_state_s':float(np.mean(g2)),'median_gain_vs_state_s':float(np.median(g2)),'improved_vs_state':int(np.sum(g2>0)),'approvals':sum(r['approved'] for r in rr),'truth_not_retained':sum(1-r['truth_retained'] for r in rr),'pre_cross':sum(r['pre_cross'] for r in rr),'post_cross':sum(r['post_cross'] for r in rr),'total_cpu_s':sum(r['cpu_seconds'] for r in rr),'median_cpu_ms':1000*float(np.median([r['cpu_seconds'] for r in rr]))}
 (output/(name+'_metrics.json')).write_text(json.dumps(out,indent=2));(output/(name+'_evidence.json')).write_text(json.dumps(evidence,indent=2));print(json.dumps(out,indent=2))
if __name__ == '__main__':
 parser = argparse.ArgumentParser(description=__doc__)
 parser.add_argument('split', choices=('development', 'confirmation'))
 parser.add_argument('step_s', nargs='?', default=DT, type=float, choices=(1., 2.))
 parser.add_argument('--smoke', action='store_true', help='Run one configuration and one synthetic query')
 args = parser.parse_args()
 if args.smoke:
  run(args.split, args.step_s, config_limit=1, truth_limit=1, output=RESULTS / 'smoke')
 else:
  run(args.split, args.step_s)
