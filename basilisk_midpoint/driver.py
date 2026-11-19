"""Independent no-slip midpoint trajectory; every accepted step is recorded."""
from pathlib import Path
import os,json,time,signal,hashlib,csv,sys
import numpy as np
from basilisk_midpoint_candidate import Trajectory
from validate import structural,atomic

ROOT=Path.cwd();RUN=ROOT/'run';RUN.mkdir(exist_ok=True);cfg=json.loads((ROOT/'config.json').read_text());validation=json.loads((ROOT/'validation.json').read_text())
assert validation['passed'] and validation['initial']==cfg['initial'] and cfg['cap_native']==validation['selected_cap_native']
height=cfg['initial'];lib=ROOT/'compiled/basilisk_midpoint.so';checkpoint=RUN/'checkpoint.npz'
checks,film=structural(height);assert film.ranks==cfg['mpi_ranks']==4
sim=Trajectory.load(checkpoint) if checkpoint.exists() else Trajectory(height,cfg['grid'],cfg['seed'],dt_cap=cfg['cap_native'],threshold=cfg['stop_u']/height,library=lib,tiered=False)
rank=sim.film.rank;assert sim.film.slip==0. and sim.threshold==cfg['stop_u']/height
stop=[False]
for sig in (signal.SIGUSR1,signal.SIGTERM,signal.SIGINT):signal.signal(sig,lambda *_:stop.__setitem__(0,True))
started=lastsave=lastprint=time.monotonic();previous=RUN/'latest.json'
wall_before=json.loads(previous.read_text())['actual_wall_seconds'] if previous.exists() else 0.
history=RUN/'history.csv'
if rank==0:
 atomic(RUN/'allocation.json',dict(job=os.getenv('SLURM_JOB_ID'),node=os.uname().nodename,mpi_ranks=sim.film.ranks,cpus_per_task=1,rank0_affinity=sorted(os.sched_getaffinity(0))))
 atomic(ROOT/'runtime_structural_gate.json',checks)
 if history.exists():
  lines=history.read_bytes().splitlines(keepends=True);keep=sim.step_count+1
  if len(lines)>keep:
   (RUN/('uncheckpointed_history_%d.csv'%time.time_ns() if hasattr(time,'time_ns') else 'uncheckpointed_history.csv')).write_bytes(b''.join(lines[keep:]));history.write_bytes(b''.join(lines[:keep]))
  assert len(history.read_bytes().splitlines())==keep
 else:
  with history.open('w',newline='') as f:csv.writer(f).writerow(['tau','dt_tau','u_min','u_max','mass_relative','full_residual','accepted','rejected','actual_wall_seconds','native_time','native_hmin'])
 (ROOT/'runtime_gate.done').write_text('no-slip midpoint MPI runtime and operator gates passed\n')
sim.film.library.bm_sync()

def save():
 if rank==0:
  sim.save(checkpoint);status='THRESHOLD' if sim.status=='RUPTURE_BRACKET' else sim.status
  event=None
  if sim.event:
   event=dict(lower=height*sim.event['lower']/1e8,upper=height*sim.event['upper']/1e8,hmin_before=height*sim.event['hmin_before'],hmin_after=height*sim.event['hmin'],threshold=cfg['stop_u'])
  atomic(previous,dict(status=status,time=height*sim.time/1e8,hmin=float(sim.h.min())*height,hmax=float(sim.h.max())*height,accepted=sim.step_count,rejected=sim.ledger['rejections'],relative_mass=float(sim.h.mean()/sim.initial_mean-1),actual_wall_seconds=wall_before+time.monotonic()-started,updated_utc=time.strftime('%Y-%m-%dT%H:%M:%SZ',time.gmtime()),event=event,config=cfg,ledger=sim.ledger,checkpoint_sha256=hashlib.sha256(checkpoint.read_bytes()).hexdigest(),history_sha256=hashlib.sha256(history.read_bytes()).hexdigest(),source_sha256=hashlib.sha256((ROOT/'basilisk_midpoint.c').read_bytes()).hexdigest()))
 sim.film.library.bm_sync()

try:
 if not checkpoint.exists():save()
 while sim.status=='RUNNING':
  row=sim.advance();now=time.monotonic()
  if sim.status=='NUMERICAL_FAILURE':break
  assert sim.ledger['max_residual']<=1e-10 and sim.ledger['max_relative_volume_error']<1e-9 and sim.h.min()>0
  if rank==0:
   with history.open('a',newline='') as f:csv.writer(f).writerow([height*sim.time/1e8,height*row['dt']/1e8,float(sim.h.min())*height,float(sim.h.max())*height,float(sim.h.mean()/sim.initial_mean-1),row['residual'],sim.step_count,sim.ledger['rejections'],wall_before+now-started,sim.time,float(sim.h.min())])
   if now-lastprint>=60:print(json.dumps(dict(tau=height*sim.time/1e8,u_min=float(sim.h.min())*height,accepted=sim.step_count,rejected=sim.ledger['rejections'],actual_wall_seconds=wall_before+now-started)),flush=True)
  if now-lastprint>=60:lastprint=now
  if sim.film.library.bm_any(int(now-lastsave>=300)):save();lastsave=now
  if sim.film.library.bm_any(int(stop[0] or now-started>5.8*3600)):break
 save()
 if sim.status=='RUPTURE_BRACKET':raise SystemExit(0)
 if sim.status=='NUMERICAL_FAILURE':raise SystemExit(82)
 raise SystemExit(75)
except Exception as error:
 sim.status='NUMERICAL_FAILURE';save()
 if rank==0:atomic(RUN/'FAILED.json',dict(error=repr(error),tau=height*sim.time/1e8))
 raise
