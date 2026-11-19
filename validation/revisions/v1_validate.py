"""Per-height, two-seed short strong-error screen; no FPT-convergence claim."""
from pathlib import Path
import json,sys,time,ctypes,math
import numpy as np
from worker import Simulation,atomic_json

ROOT=Path.cwd();cfg=json.loads((ROOT/'config.json').read_text());height=cfg['initial']
policy=json.loads((ROOT/'campaign.json').read_text())['validation']
ptr=ctypes.POINTER(ctypes.c_double)
def call(s,h,noise,dt):
 out=np.empty_like(h);d=np.zeros(8)
 code=s.lib.wall_step(s.ptr(np.ascontiguousarray(h)),s.ptr(np.ascontiguousarray(noise)),dt,s.ptr(out),s.ptr(d))
 assert code==0 and np.isfinite(out).all() and out.min()>0 and d[3]<1e-8 and abs(out.mean()/height-1)<1e-9,(code,d.tolist())
 return out,d

def structural():
 s=Simulation(dict(cfg,dt_divisor=64),ROOT/'structural')
 s.lib.wall_noise_probe.argtypes=[ctypes.c_double]*4+[ptr]
 z=np.zeros(2);s.lib.wall_noise_probe(height,1.,1/64,2.1e-9/8,s.ptr(z))
 assert np.allclose(z,[math.sqrt(2*cfg['phi'])*height**1.5,1.5*math.sqrt(2*cfg['phi']*height)],rtol=1e-13)
 s.lib.wall_set_phi(0.)
 n=s.n;x=np.arange(n)/n;h=height+1e-4*np.cos(2*np.pi*x[:,None])+np.zeros((n,n));zero=np.zeros((2,n,n))
 left,_=call(s,h,zero,1e-7);right,_=call(s,h.T.copy(),zero,1e-7)
 rotation=float(np.max(abs(left-right.T)));assert rotation<1e-9
 rates=[]
 for mode in (1,6):
  dt=1e-7;steps=128;wave=np.cos(2*np.pi*mode*x[:,None])+np.zeros((n,n))
  h=height+1e-4*wave
  for _ in range(steps):h,_=call(s,h,zero,dt)
  measured=math.log(float(np.mean((h-height)*wave)*2)/1e-4)/(steps*dt)
  q2=(2*n*math.sin(np.pi*mode/n))**2
  continuous_discrete=height**3*q2*(3*cfg['A']/height**4-q2)
  be=-math.log1p(-dt*continuous_discrete)/dt
  error=abs(measured/be-1);rates.append(dict(mode=mode,measured=measured,discrete_BE=be,relative_error=error))
  assert error<.02, rates[-1]
 s.lib.wall_set_phi(cfg['phi'])
 return dict(noise_value_derivative=z.tolist(),xy_rotation_error=rotation,normal_modes=rates)

started=time.monotonic()
try:
 checks=structural();atomic_json(ROOT/'structural_pass.json',checks)
 results={};reference_frames=[];reference_hmins=[];replay=[]
 for seed_index in range(policy['seeds_per_height']):
  seed=2026111990+int(round(height*10))*10+seed_index
  candidate_cfg=dict(cfg,seed=seed,dt_divisor=policy['reference_divisor'])
  s=Simulation(candidate_cfg,ROOT/('fine_%d'%seed));frames=[];hmins=[]
  target=policy['macro_intervals']*s.noise.width;framewidth=target//policy['frames'];boundary=framewidth
  while s.cursor<target:
   s.advance()
   if s.cursor==boundary:frames.append(s.h.copy());hmins.append(float(s.h.min()));boundary+=framewidth
  assert len(frames)==policy['frames'] and s.status=='RUNNING'
  s.save();s.advance();expected=s.h.copy();r=Simulation(candidate_cfg,s.out);r.advance()
  delta=float(np.max(abs(r.h-expected)));assert delta<1e-10;replay.append(delta)
  reference_frames.append(np.array(frames));reference_hmins.append(np.array(hmins))
  print('REFERENCE',height,seed,delta,flush=True)
 for divisor in policy['candidates']:
  errors=[]
  for seed_index in range(policy['seeds_per_height']):
   seed=2026111990+int(round(height*10))*10+seed_index
   s=Simulation(dict(cfg,seed=seed,dt_divisor=divisor),ROOT/('d%d_%d'%(divisor,seed)))
   frames=[];mins=[];target=policy['macro_intervals']*s.noise.width;framewidth=target//policy['frames'];boundary=framewidth
   while s.cursor<target:
    s.advance()
    if s.cursor==boundary:frames.append(s.h.copy());mins.append(float(s.h.min()));boundary+=framewidth
   assert len(frames)==policy['frames'] and s.status=='RUNNING'
   fine=reference_frames[seed_index];den=np.std(fine,axis=(1,2))
   field=float(np.max(np.sqrt(np.mean((np.array(frames)-fine)**2,axis=(1,2)))/den))
   herror=float(np.max(abs(np.array(mins)-reference_hmins[seed_index])))
   errors.append(dict(seed=seed,max_relative_field_rms_error=field,max_abs_hmin_difference=herror,steps=s.steps,rejected=s.rejected))
  passed=all(v['max_relative_field_rms_error']<=policy['max_relative_field_rms_error'] and v['max_abs_hmin_difference']<=policy['max_abs_hmin_difference'] for v in errors)
  results[str(divisor)]=dict(passed=passed,seeds=errors)
  print('SCREEN',height,divisor,passed,json.dumps(errors),flush=True)
  if passed:
   report=dict(passed=True,initial=height,selected_divisor=divisor,dt_cap=2.1e-9/divisor,structural=checks,restart_max_errors=replay,screen=results,policy=policy,elapsed_actual_seconds=time.monotonic()-started,limitation='Two-seed short-time strong-error screen, not long-time or first-passage convergence')
   atomic_json(ROOT/'validation.json',report);sys.exit(0)
 raise RuntimeError('No screened candidate passed; no production authorization for this height')
except Exception as error:
 atomic_json(ROOT/'validation_failed.json',dict(passed=False,initial=height,error=repr(error),elapsed_actual_seconds=time.monotonic()-started));raise
