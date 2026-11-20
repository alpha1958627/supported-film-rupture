"""Average five stopped hmin trajectories on common physical times, from full histories."""
from pathlib import Path
import argparse,json,subprocess,shlex,concurrent.futures
import numpy as np

p=argparse.ArgumentParser();p.add_argument('--workspace',required=True);args=p.parse_args()
root=Path(args.workspace);snapshot=json.loads((root/'bm_snapshot.json').read_text());selected=json.loads((root/'selected.json').read_text());campaign=json.loads((root/'campaign.json').read_text())
records=[r for r in snapshot['records'] if r['kind']=='production'];assert len(records)==25
heights={c['name']:c['initial'] for c in campaign['cases']};grids={}
for h in (.5,.6,.7,.8,.9):
 group=[r for r in records if heights[r['name']]==h];assert len(group)==5 and all(r['latest']['status']=='THRESHOLD' for r in group)
 end=max(r['latest']['time'] for r in group)
 grids[str(h)]=np.unique(np.r_[0,np.geomspace(1e-11,end,1600),[r['latest']['time'] for r in group]]).tolist()
byhost={}
for r in records:
 route=selected[r['name']]['winner'];d=r['latest']
 item=dict(name=r['name'],remote=route['remote'],height=str(heights[r['name']]),accepted=d['accepted'],time=d['time'],history_sha256=d['history_sha256'],checkpoint_sha256=d['checkpoint_sha256'])
 byhost.setdefault(route['host'],dict(route=route,items=[]))['items'].append(item)
code=r'''
from pathlib import Path
import json,hashlib,io
import numpy as np
items=ITEMS;grids=GRIDS;out=[]
for v in items:
 p=Path(v['remote'])/'run'
 assert hashlib.sha256((p/'checkpoint.npz').read_bytes()).hexdigest()==v['checkpoint_sha256']
 with (p/'history.csv').open('rb') as f:raw=b''.join(f.readline() for _ in range(v['accepted']+1))
 assert hashlib.sha256(raw).hexdigest()==v['history_sha256']
 a=np.loadtxt(io.BytesIO(raw),delimiter=',',skiprows=1,usecols=(0,2),ndmin=2)
 assert len(a)==v['accepted'] and a[-1,0]==v['time'] and np.all(np.diff(a[:,0])>0)
 # Stop at the accepted first passage; hold at the threshold thereafter.
 t=np.r_[0,a[:,0]];y=np.r_[float(v['height']),a[:,1]];y[-1]=.4
 assert np.all(y>=.4-1e-12)
 values=np.interp(grids[v['height']],t,y,left=float(v['height']),right=.4)
 out.append(dict(name=v['name'],height=v['height'],values=values.tolist()))
print(json.dumps(out))
'''
def fetch(group):
 route=group['route'];program=code.replace('ITEMS',repr(group['items'])).replace('GRIDS',repr(grids))
 command='PYTHONPATH='+shlex.quote(route.get('pythonpath',''))+' '+route['python']+' -'
 result=subprocess.run(['ssh','-o','BatchMode=yes',route['host'],'bash -lc '+shlex.quote(command)],input=program,capture_output=True,text=True,timeout=300)
 assert result.returncode==0,result.stderr[-2000:]
 return json.loads(result.stdout)
with concurrent.futures.ThreadPoolExecutor(max_workers=2) as pool:paths=[r for result in pool.map(fetch,byhost.values()) for r in result]
groups=[]
for h,t in grids.items():
 values=np.array([r['values'] for r in paths if r['height']==h]);assert values.shape==(5,len(t))
 mean=values.mean(axis=0);assert abs(mean[0]-float(h))<1e-12 and abs(mean[-1]-.4)<1e-12
 groups.append(dict(initial=float(h),n=5,time=t,mean=mean.tolist()))
result=dict(snapshot_utc=snapshot['utc'],definition='Arithmetic mean of five hmin paths at common physical times. Linear interpolation of full accepted-step histories; after first passage each path is held at .4. No smoothing; n=5 throughout. Last point is the latest event, not the mean first-passage time.',source_history_sha256={r['name']:r['latest']['history_sha256'] for r in records},groups=groups)
(root/'mean_curves.json').write_text(json.dumps(result,indent=2));print('Full-history five-seed means verified for all five thicknesses')
