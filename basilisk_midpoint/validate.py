"""No-slip MPI operator/restart gates and short same-path midpoint step screen."""
from pathlib import Path
import json,time,hashlib
import numpy as np
from midpoint_candidate import Film as ReferenceFilm
from basilisk_midpoint_candidate import Trajectory,Film

ROOT=Path.cwd();LIB=ROOT/'compiled/basilisk_midpoint.so'
def atomic(path,obj):
 temp=path.with_suffix('.tmp');temp.write_text(json.dumps(obj,indent=2));temp.replace(path)
def structural(height):
 film=Film(height,64,LIB);ref=ReferenceFilm(height,64)
 rng=np.random.default_rng(2026113001);h=1+.003*rng.normal(size=(64,64));v=rng.normal(size=(64,64));dw=.001*rng.normal(size=(2,64,64))
 errors=dict(force=float(np.max(abs(film.force(h)-ref.force(h)))),noise=float(np.max(abs(film.stochastic(h,dw)-ref.stochastic(h,dw)))))
 df,dn=film.derivatives(h,v,dw);rf,rn=ref.derivatives(h,v,dw)
 errors.update(dforce=float(np.max(abs(df-rf))),dnoise=float(np.max(abs(dn-rn))))
 assert max(errors.values())<1e-10,errors
 actual,diagnostic=film.step(h,.03,dw);expected,_=ref.step(h,.03,dw)
 errors['step']=float(np.max(abs(actual-expected)));assert errors['step']<1e-9 and diagnostic['residual']<=1e-10
 assert abs(actual.mean()/h.mean()-1)<1e-9 and film.library.bm_slip_value()==0.
 # With flat faces, removing slip changes noise by this exact amplitude ratio.
 flat=np.ones((64,64));old=ReferenceFilm(height,64);old.slip=.2/height
 ratio=float(np.max(abs(film.stochastic(flat,dw)-old.stochastic(flat,dw)/np.sqrt(1+.2/height))))
 assert ratio<1e-12
 return dict(passed=True,operator_errors=errors,flat_noise_scaling_error=ratio,slip_native=0.,mpi_ranks=film.ranks),film

def main():
 cfg=json.loads((ROOT/'config.json').read_text());height=cfg['initial'];policy=json.loads((ROOT/'campaign.json').read_text())['validation'];started=time.monotonic()
 checks,film=structural(height);rank=film.rank;assert film.ranks==4
 if rank==0:
  atomic(ROOT/'structural_pass.json',checks);(ROOT/'runtime_gate.done').write_text('source/operator gate passed; step screen pending\n')
 refs=[];refmins=[];replay=[];results={}
 for index in range(policy['seeds_per_height']):
  seed=2026113900+round(height*10)*10+index
  sim=Trajectory(height,64,seed,dt_cap=policy['reference_cap'],threshold=.01,library=LIB,tiered=False)
  frames=[];mins=[]
  for timepoint in np.linspace(policy['native_end']/policy['frames'],policy['native_end'],policy['frames']):
   boundary=sim.noise.ticks(timepoint)
   while sim.ticks<boundary:
    sim.advance(stop_ticks=boundary);assert sim.status=='RUNNING',sim.ledger
   frames.append(sim.h.copy());mins.append(float(sim.h.min())*height)
  refs.append(np.array(frames)*height);refmins.append(np.array(mins))
  checkpoint=ROOT/('replay_%d.npz'%index)
  if rank==0:sim.save(checkpoint)
  sim.film.library.bm_sync();sim.advance();expected=sim.h.copy()
  restored=Trajectory.load(checkpoint);restored.advance();error=float(np.max(abs(expected-restored.h)))*height
  assert error<1e-9;replay.append(error)
  if rank==0:print('REFERENCE',height,seed,sim.step_count,error,flush=True)
 for cap in policy['candidate_caps']:
  errors=[]
  for index in range(policy['seeds_per_height']):
   seed=2026113900+round(height*10)*10+index
   sim=Trajectory(height,64,seed,dt_cap=cap,threshold=.01,library=LIB,tiered=False)
   frames=[];mins=[]
   for timepoint in np.linspace(policy['native_end']/policy['frames'],policy['native_end'],policy['frames']):
    boundary=sim.noise.ticks(timepoint)
    while sim.ticks<boundary:
     sim.advance(stop_ticks=boundary);assert sim.status=='RUNNING',sim.ledger
    frames.append(sim.h.copy()*height);mins.append(float(sim.h.min())*height)
   fine=refs[index];rms=np.std(fine,axis=(1,2))
   field=float(np.max(np.sqrt(np.mean((np.array(frames)-fine)**2,axis=(1,2)))/rms));delta=float(np.max(abs(np.array(mins)-refmins[index])))
   errors.append(dict(seed=seed,max_field_rms_relative=field,max_abs_u_min_difference=delta,steps=sim.step_count,ledger=sim.ledger))
  passed=all(v['max_field_rms_relative']<=policy['max_field_rms_relative'] and v['max_abs_u_min_difference']<=policy['max_abs_u_min_difference'] for v in errors)
  results[str(cap)]=dict(passed=passed,seeds=errors)
  if rank==0:print('SCREEN',height,cap,passed,[(v['max_field_rms_relative'],v['max_abs_u_min_difference']) for v in errors],flush=True)
  if passed:
   if rank==0:atomic(ROOT/'validation.json',dict(passed=True,initial=height,selected_cap_native=cap,selected_cap_tau=cap*height/1e8,structural=checks,restart_max_errors=replay,screen=results,policy=policy,actual_seconds=time.monotonic()-started,source_sha256={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in ROOT.glob('*.py')},limitation='Two-seed short strong-error screen, not long-time or first-passage convergence. Validation only uses native stop .01 to avoid interrupting the common-time comparison. Production stop is u=.4.'))
   return
  sim.film.library.bm_sync()
 if rank==0:atomic(ROOT/'validation_failed.json',dict(passed=False,initial=height,screen=results,error='No candidate passed; no production allowed'))
 raise RuntimeError('No midpoint candidate passed')

if __name__=='__main__':
 try:main()
 except Exception as error:
  if int(__import__('os').environ.get('SLURM_PROCID','0'))==0 and not (ROOT/'validation_failed.json').exists():atomic(ROOT/'validation_failed.json',dict(passed=False,error=repr(error)))
  raise
