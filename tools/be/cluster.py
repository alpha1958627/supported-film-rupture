"""SIPF-first submission, verified per-case race, validation-gated production."""
from pathlib import Path
import argparse,json,subprocess,shlex,tarfile,hashlib,time,traceback,base64

parser=argparse.ArgumentParser();parser.add_argument('--workspace',required=True);parser.add_argument('action',choices=['prepare','advance','watch']);args=parser.parse_args()
ROOT=Path(args.workspace).resolve();REPO=Path(__file__).resolve().parents[2]
campaign=json.loads((REPO/'campaign.json').read_text());sites=json.loads((ROOT/'sites.json').read_text())
def atomic(p,obj):
 temp=p.with_suffix('.tmp');temp.write_text(json.dumps(obj,indent=2),encoding='utf-8');temp.replace(p)
def run(command,timeout=180):
 r=subprocess.run(command,stdout=subprocess.PIPE,stderr=subprocess.PIPE,universal_newlines=True,encoding='utf-8',errors='replace',timeout=timeout)
 if r.returncode:raise RuntimeError((command[0],r.returncode,r.stderr[-1500:]))
 return r.stdout.strip()
def ssh(site,code):return run(['ssh','-o','BatchMode=yes','-o','ConnectTimeout=12',site['host'],'bash -lc '+shlex.quote(code)])
def query(routes):
 states={}
 for site in sites:
  own=[r for r in routes if r['server']==site['server']]
  if not own:continue
  ids=','.join(r['job'] for r in own)
  try:
   text=ssh(site,'squeue -h -j '+ids+' -o "%i|%T|%S|%N"; sacct -nP -j '+ids+' -X -o JobID,State,Start,NodeList')
   for line in text.splitlines():
    part=line.split('|')
    if len(part)>=4 and part[0] in [r['job'] for r in own]:states.setdefault(part[0],dict(state=part[1],start=part[2],node=part[3]))
  except Exception as e:print('QUERY_ERROR',site['server'],repr(e),flush=True)
 return states
def setup(site):
 receipts=ROOT/'setup';receipts.mkdir(exist_ok=True);receipt=receipts/(site['server']+'.json')
 digest=hashlib.sha256((ROOT/'runtime.tar.gz').read_bytes()).hexdigest()
 if receipt.exists() and json.loads(receipt.read_text()).get('runtime_sha256')==digest:return
 dest=site['base']+'/common';ssh(site,'mkdir -p '+shlex.quote(dest))
 run(['scp','-q',str(ROOT/'runtime.tar.gz'),site['host']+':'+dest+'/runtime.tar.gz'])
 command='cd '+shlex.quote(dest)+' && tar -xzf runtime.tar.gz'
 if site['server']!='SIPF':command+=' && export BASILISK='+shlex.quote(site['basilisk'])+' && bash build.sh'
 ssh(site,command);atomic(receipt,dict(site=site,ready=True,runtime_sha256=digest))
def submit(site,name,cfg,kind,validation=None):
 setup(site);remote=site['base']+'/'+('validation' if kind=='validation' else 'production')+'/'+name
 assert name.isalnum() or all(c.isalnum() or c=='_' for c in name)
 existing=ssh(site,'if test -f '+shlex.quote(remote+'/job.txt')+'; then cat '+shlex.quote(remote+'/job.txt')+'; fi')
 if existing:
  assert existing.isdigit() and json.loads(ssh(site,'cat '+shlex.quote(remote+'/config.json')))==cfg
  return dict(site,name=name,kind=kind,remote=remote,job=existing)
 ssh(site,'mkdir -p '+shlex.quote(remote)+' && cp '+shlex.quote(site['base']+'/common')+'/* '+shlex.quote(remote)+'/')
 local=ROOT/'inputs'/name;local.mkdir(parents=True,exist_ok=True)
 atomic(local/'config.json',cfg)
 names=['config.json']
 if validation is not None:atomic(local/'validation.json',validation);names.append('validation.json')
 for filename in names:run(['scp','-q',str(local/filename),site['host']+':'+remote+'/'+filename])
 py=site['python'];minutes='00:30:00' if kind=='validation' else '06:00:00'
 entry='validate.py' if kind=='validation' else 'launch.py'
 script='#!/bin/bash\n#SBATCH --job-name=BE_'+name+'\n#SBATCH --partition='+site['partition']+'\n#SBATCH --nodes=1\n#SBATCH --ntasks=1\n#SBATCH --cpus-per-task='+str(site['allocated_cpus'])+'\n#SBATCH --mem=4G\n#SBATCH --time='+minutes+'\n#SBATCH --requeue\n#SBATCH --signal=B:USR1@600\n#SBATCH --open-mode=append\n#SBATCH --output=slurm-%j.out\nset -euo pipefail\ncd '+remote+'\nsha256sum -c inputs.sha256 >/dev/null\nexport OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=4\ntouch runtime_gate.done\n'+py+' -u '+entry+' &\nchild=$!\ntrap \'kill -USR1 "$child" 2>/dev/null || true\' USR1 TERM\nset +e\nwait "$child"\nrc=$?\nif kill -0 "$child" 2>/dev/null; then wait "$child"; rc=$?; fi\nset -e\nif test "$rc" = 75; then scontrol requeue "$SLURM_JOB_ID"; exit 0; fi\nexit "$rc"\n'
 (local/'run.sbatch').write_bytes(script.encode());run(['scp','-q',str(local/'run.sbatch'),site['host']+':'+remote+'/run.sbatch'])
 command='cd '+shlex.quote(remote)+' && job=$(sbatch --parsable run.sbatch) && printf "%s" "${job%%;*}" > job.txt && cat job.txt'
 job=ssh(site,command);assert job.isdigit()
 route=dict(site,name=name,kind=kind,remote=remote,job=job);print('SUBMITTED',name,site['server'],job,flush=True);return route
def cleanup(route):
 folder=ROOT/'cancelled';folder.mkdir(exist_ok=True);receipt=folder/(route['name']+'_'+route['server']+'_'+route['job']+'.json')
 if receipt.exists():return
 try:
  info=ssh(route,'scontrol show job '+route['job']+' -o')
  assert 'WorkDir='+route['remote'] in info and 'UserId='+route['user']+'(' in info
 except RuntimeError:
  text=ssh(route,'sacct -nP -j '+route['job']+' -X -o JobID,User,WorkDir')
  assert any(v.split('|')[:3]==[route['job'],route['user'],route['remote']] for v in text.splitlines())
  assert ssh(route,'cat '+shlex.quote(route['remote']+'/job.txt'))==route['job']
 if ssh(route,'squeue -h -j '+route['job']):ssh(route,'scancel '+route['job'])
 for _ in range(30):
  if not ssh(route,'squeue -h -j '+route['job']):break
  time.sleep(2)
 else:raise RuntimeError('Cancellation not confirmed')
 archive=route['remote']+'.cancel_'+route['job']+'.tar.gz'
 digest=ssh(route,'tar -czf '+shlex.quote(archive)+' -C '+shlex.quote(route['remote'])+' . && sha256sum '+shlex.quote(archive)).split()[0]
 local=receipt.with_suffix('.tar.gz');run(['scp','-q',route['host']+':'+archive,str(local)])
 assert hashlib.sha256(local.read_bytes()).hexdigest()==digest
 code='from pathlib import Path;import shutil;p=Path('+repr(route['remote'])+');assert str(p.resolve())==str(p) and p.name=='+repr(route['name'])+' and p.parent.name in ("validation","production") and p.parent.parent.name=="wall_be_ensemble_20261119";shutil.rmtree(str(p));Path('+repr(archive)+').unlink()'
 ssh(route,route['python']+' -c '+shlex.quote(code));atomic(receipt,dict(route=route,cancelled=True,cleaned=True,archive_sha256=digest))
def cycle():
 routes=json.loads((ROOT/'candidates.json').read_text()) if (ROOT/'candidates.json').exists() else []
 selected=json.loads((ROOT/'selected.json').read_text()) if (ROOT/'selected.json').exists() else {}
 validations={}
 desired={}
 for height in (.5,.6,.7,.8,.9):
  cfg=dict(next(c for c in campaign['cases'] if c['initial']==height),dt_divisor=64)
  name='v2_h%02d'%round(height*10);desired[name]=(cfg,'validation',None)
 states=query(routes)
 for name in sorted({r['name'] for r in routes}):
  own=[r for r in routes if r['name']==name]
  if name not in selected:
   active=[r for r in own if states.get(r['job'],{}).get('state') in ('RUNNING','COMPLETED')]
   for winner in sorted(active,key=lambda r:states[r['job']].get('start','')):
    if ssh(winner,'test ! -f '+shlex.quote(winner['remote']+'/runtime_gate.done')+' || echo READY')=='READY':
     selected[name]=dict(winner=winner,selected_utc=time.strftime('%Y-%m-%dT%H:%M:%SZ',time.gmtime()),cleanup_complete=False);atomic(ROOT/'selected.json',selected);break
  if name in selected:
   choice=selected[name]
   if not choice['cleanup_complete']:
    for route in own:
     if route['job']!=choice['winner']['job']:cleanup(route)
    choice['cleanup_complete']=True;atomic(ROOT/'selected.json',selected)
   if choice['winner']['kind']=='validation':
    winner=choice['winner'];text=ssh(winner,'if test -f '+winner['remote']+'/validation.json; then cat '+winner['remote']+'/validation.json; elif test -f '+winner['remote']+'/validation_failed.json; then cat '+winner['remote']+'/validation_failed.json; fi')
    if text:
     validation=json.loads(text);validations[name]=validation
     (ROOT/'validation').mkdir(exist_ok=True);atomic(ROOT/'validation'/(name+'.json'),validation)
     if validation.get('passed'):
      for case in campaign['cases']:
       if case['initial']==validation['initial']:desired[case['name']]=(dict(case,dt_divisor=validation['selected_divisor']), 'production', validation)
 for name,(cfg,kind,validation) in desired.items():
  own=[r for r in routes if r['name']==name]
  if not own:
   route=submit(sites[0],name,cfg,kind,validation);routes.append(route);atomic(ROOT/'candidates.json',routes)
  elif name not in selected and not any(states.get(r['job'],{}).get('state') in ('FAILED','CANCELLED','TIMEOUT','OUT_OF_MEMORY') for r in own):
   if all(states.get(r['job'],{}).get('state')=='PENDING' for r in own):
    for site in sites[1:]:
     if not any(r['server']==site['server'] for r in own):
      route=submit(site,name,cfg,kind,validation);routes.append(route);atomic(ROOT/'candidates.json',routes)
 atomic(ROOT/'controller_status.json',dict(utc=time.strftime('%Y-%m-%dT%H:%M:%SZ',time.gmtime()),states=states,selected=selected,validations=validations,production_submitted=len({r['name'] for r in routes if r['kind']=='production'})))
 done=len(selected)==30 and all(v['cleanup_complete'] for v in selected.values())
 if done:atomic(ROOT/'controller_done.json',dict(all_25_production_selected=True,validation_passed_heights=5))
 return done
if args.action=='prepare':
 with tarfile.open(ROOT/'runtime.tar.gz','w:gz') as tar:
  names=[]
  for p in (REPO/'be').iterdir():
   if p.is_file():tar.add(p,arcname=p.name);names.append((p.name,hashlib.sha256(p.read_bytes()).hexdigest()))
  for p in [REPO/'campaign.json',ROOT/'wall_native.so']:
   tar.add(p,arcname=p.name);names.append((p.name,hashlib.sha256(p.read_bytes()).hexdigest()))
  # Config and height-specific validation are generated separately, fingerprinted in saved checkpoints.
  content='\n'.join(digest+'  '+name for name,digest in names)+'\n'
  import io
  info=tarfile.TarInfo('inputs.sha256');info.size=len(content.encode());tar.addfile(info,io.BytesIO(content.encode()))
else:
 while True:
  try:
   done=cycle()
   if done or args.action=='advance':break
  except Exception:print(traceback.format_exc(),flush=True)
  time.sleep(30)
