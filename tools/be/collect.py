"""Read-only snapshots; plots stop at verified checkpoint accepted counts."""
from pathlib import Path
import argparse,json,subprocess,base64,hashlib,time

p=argparse.ArgumentParser();p.add_argument('--workspace',required=True);a=p.parse_args()
ROOT=Path(a.workspace).resolve();REPO=Path(__file__).resolve().parents[2]
campaign=json.loads((REPO/'campaign.json').read_text());selected=json.loads((ROOT/'selected.json').read_text()) if (ROOT/'selected.json').exists() else {}
routes=[v['winner'] for v in selected.values()]
script=r'''
import json,hashlib,subprocess,io
from pathlib import Path
import numpy as np
routes=ROUTES
ids=','.join(r['job'] for r in routes)
state=subprocess.run(['bash','-lc','squeue -h -j '+ids+' -o "%i|%T|%N|%M|%L"; sacct -nP -j '+ids+' -X -o JobID,State,NodeList,Elapsed,Timelimit,Start,End'],stdout=subprocess.PIPE,stderr=subprocess.PIPE,universal_newlines=True)
items=[]
for route in routes:
 item=dict(route=route);folder=Path(route['remote']);latest=folder/'run/latest.json'
 for name in ('validation.json','validation_failed.json','structural_pass.json','restart_gate.json'):
  if (folder/name).exists():item[name]=json.loads((folder/name).read_text())
 if latest.exists():
  for attempt in range(5):
   d=json.loads(latest.read_text());checkpoint=(folder/'run/checkpoint.npz').read_bytes()
   if hashlib.sha256(checkpoint).hexdigest()==d['checkpoint_sha256']:break
  else:raise RuntimeError('Checkpoint changed during snapshot')
  with (folder/'run/history.csv').open('rb') as f:raw=b''.join(f.readline() for _ in range(d['accepted']+1))
  data=np.loadtxt(io.BytesIO(raw),delimiter=',',skiprows=1,ndmin=2)
  if d['accepted']:
   assert len(data)==d['accepted'] and data[-1,0]==d['time']
   assert abs(data[-1,2]-d['hmin'])<1e-12 and np.isfinite(data).all()
   bins=np.minimum((data[:,0]/d['time']*1200).astype(int),1199);keep=set([0,len(data)-1])
   starts=np.r_[0,np.flatnonzero(np.diff(bins))+1];ends=np.r_[starts[1:],len(data)]
   for lo,hi in zip(starts,ends):keep.update([int(lo),int(hi-1),int(lo+np.argmin(data[lo:hi,2])),int(lo+np.argmax(data[lo:hi,2]))])
   item['points']=data[sorted(keep)][:,[0,2]].tolist()
   item['max_abs_mass']=float(abs(data[:,4]).max());item['max_full_residual']=float(data[:,5].max());item['min_hmin']=float(data[:,2].min())
   if item['max_abs_mass']>=1e-9 or item['max_full_residual']>=1e-8:item['numerical_gate_violation']=True
  else:item['points']=[]
  item['latest']=d;item['history_prefix_sha256']=hashlib.sha256(raw).hexdigest()
 items.append(item)
print(json.dumps(dict(items=items,scheduler=state.stdout,stderr=state.stderr)))
'''
groups={}
for route in routes:groups.setdefault(route['host'],[]).append(route)
result=dict(updated_utc=time.strftime('%Y-%m-%dT%H:%M:%SZ',time.gmtime()),campaign=campaign,groups=[],selected=selected)
for host,own in groups.items():
 code=script.replace('ROUTES',repr(own));command='echo '+base64.b64encode(code.encode()).decode()+' | base64 -d | '+own[0]['python']+' -'
 r=subprocess.run(['ssh','-o','BatchMode=yes',host,command],stdout=subprocess.PIPE,stderr=subprocess.PIPE,universal_newlines=True,encoding='utf-8',errors='replace',timeout=180)
 if r.returncode:result['groups'].append(dict(host=host,read_error=r.stderr));continue
 result['groups'].append(dict(host=host,**json.loads(r.stdout)))
(ROOT/'snapshot.json').write_text(json.dumps(result,separators=(',',':')),encoding='utf-8')
print(json.dumps(dict(updated_utc=result['updated_utc'],selected=len(selected),hosts=list(groups)),indent=2))
